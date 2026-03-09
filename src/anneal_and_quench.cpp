#include <algorithm>
#include <argparse.hpp>
#include <filesystem>
#include <iostream>
#include <ostream>


#include "MC.hpp"
#include "chain.hpp"
#include "lattice_IO.hpp"
#include "run_common.hpp"
#include "stats.hpp"
#include "pyrochlore_geometry.hpp"
#include "format_bits.hpp"

/*
This program performs classical Monte Carlo simulated annealing for a spin
model on the pyrochlore (diamond-based) lattice, using local Metropolis updates.
The simulation:

1. Constructs a cubic supercell of the diamond lattice.
2. Defines exchange couplings up to third neighbors.
3. Performs a burn-in at high temperature.
4. Anneals the system from a hot temperature to a target cold temperature.
5. Turns off the magnetic field, and re-equilibrates.
5. Samples static spin–spin correlation functions at the final temperature.
6. Optionally saves the final spin configuration.

The output consists of:

a. Static structure factor data (.ssf.h5)
b. Optionally, the final spin state (.spins.h5)

 */


using namespace std;
using namespace CMC;



int main (int argc, char *argv[]) {
    ///////////////////////////////////////////////////////////////////////////
    /// Setup and CLI
    ///////////////////////////////////////////////////////////////////////////
    /// Setup and CLI
    argparse::ArgumentParser prog("anneal");

    add_bookkeeping_args(prog);
    add_annealer_args(prog);
    prog.add_argument("--n_post_quench", "-q")
        .help("Number of sweeps to run after quench")
        .default_value(static_cast<size_t>(128))
        .scan<'i', size_t>();
    add_J1J2J3_args(prog);

    /// PHYSICAL
    prog.add_argument("--external_field", "-B")
        .help("Global magnetic field")
        .nargs(3)
        .default_value(std::vector<double>({0.,0.,0.}))
        .scan<'g', double>();

    try {
        prog.parse_args(argc, argv);
    } catch (const std::exception& err){
        cerr << err.what() << endl;
        cerr << prog;
        std::exit(1);
    }


    ///////////////////////////////////////////////////////////////////////////
    /// Input loading and validation
    

    /// Ensuring directories exist AHEAD of time (avoids heartbreak)
    std::string outdir_s = prog.get<std::string>("output_dir");
    filesystem::path outdir(outdir_s);
    if (! filesystem::exists(outdir) ){
        throw runtime_error("Cannot open outdir");
    }


    int L = prog.get<int>("L");
    auto supercell_spec = imat33_t::from_cols( {L,0,0}, {0,L,0},{0,0,L});

    auto seed = get_hex(prog, "--seed");

    const auto spec = DiamondC();
    CMC::Lattice lat(spec, supercell_spec); 

    auto J1 = prog.get<double>("--J1");
    auto J2 = prog.get<double>("--J2");
    auto J3 = prog.get<double>("--J3");
    vector3::vec3d global_field;
    { 
        auto B_tmp = prog.get<std::vector<double>>("-B");
        for (int i=0; i<3; i++){ global_field[i]= B_tmp[i]; }
        cout<<"B="<<global_field<<std::endl;
    }

    CMC::MC_runner mc(lat, seed, global_field);
    mc.define_coupling("J1", pyrochlore::nn1_dist, 
        mat33d::from_cols({J1, 0,0}, {0,J1, 0}, {0,0,J1})
        );
    mc.define_coupling("J2", pyrochlore::nn2_dist, J2*coupling::Heis);
    mc.define_coupling("J3a", pyrochlore::nn3a_dist, J3*coupling::Heis);
    mc.define_coupling("J3b", pyrochlore::nn3b_dist, J3*coupling::Heis);

    mc.settings.T_ref = prog.get<double>("--T_ref");
    mc.setup_lattice();

    const double T_hot = prog.is_used("--T_hot") ? 
        prog.get<double>("--T_hot") :
        sqrt(J1*J1 + J2* J2 + J3*J3)*10;
    const double T_cold = prog.get<double>("--T_cold");

    const size_t n_steps = prog.get<size_t>("--n_steps");
    const size_t n_sweep = prog.get<size_t>("--n_sweep");
    const size_t n_burn_in = prog.get<size_t>("--n_burn_in");
    const size_t n_sample = prog.get<size_t>("--n_sample");

    double T = T_hot;

    // Parameter specification complete. Set the name...

    std::string name = get_ofilename("anneal_quench", prog);

    printf("Burning in (%zu sweeps)...\n", n_burn_in);
    for (size_t i=0; i<n_burn_in; i++){
        mc.sweep_local_Metropolis(T_hot);
    }

    printf("Done. Begin anneal...\n");
    const double factor = pow(T_cold / T_hot, 1./n_steps);

    energy_manager e_manager(mc);

    // Cool down to low temperature
    for (size_t i=0; i <n_steps; i++){
        size_t accepted=0;
        e_manager.new_T(T);
        for (size_t n=0; n<n_sweep; n++){
            accepted += mc.sweep_local_Metropolis(T);
        }
        e_manager.sample();


        double E = e_manager.curr_E();
        printf("Iter %4zu T=%.3e E=%3e Acceptance rate: %.2f%%\n", 
                i, T, E, accepted*100.0/lat.links.size()/n_sweep);
        T *= factor;
    }

    //quench


    const size_t n_sweep_postQ = prog.get<size_t>("--n_post_quench");

    printf("Quenching field and re-equilibrating (%zu sweeps)\n",
            n_sweep_postQ);
    global_field = {0,0,0};

    for (size_t i=0; i<n_sweep_postQ; i++){
        mc.sweep_local_Metropolis(T_cold);
    }

    static_corr_3D ssf_manager(lat);

    ssf_manager.declare_observable("SdotS", 
            static_corr_3D::NEEDS_XX | static_corr_3D::NEEDS_YY | static_corr_3D::NEEDS_ZZ );
    ssf_manager.declare_observable("SzSz", static_corr_3D::NEEDS_ZZ);


    printf("Sampling at T=%lf (%zu sweeps)...\n", T, n_sample);
    for (size_t i=0; i<n_sample; i++){
        for (size_t n=0; n<n_sweep; n++){
            mc.sweep_local_Metropolis(T);
        }
        ssf_manager.sample();
    }

    auto file_path = outdir/( name + ".out.h5");

    hid_t file_id = H5Fcreate(file_path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        throw std::runtime_error("Failed to create HDF5 file: " + file_path.string());
    }
    ssf_manager.write_group(file_id, "/ssf");
    e_manager.write_group(file_id, "/energy");
    
    if (prog.get<bool>("--save_state")){
        auto f = outdir /( name + ".spins.h5" );
        printf("Saving spin state to %s\n", f.string().c_str());
        save_spin_state(lat, f);
    }



    return 0;
}
