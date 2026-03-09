#include "MC.hpp"
#include "argparse.hpp"
#include "pyrochlore_geometry.hpp"
#include "stats.hpp"
#include <iostream>
#include "format_bits.hpp"
#include "preset_cellspecs.hpp"


using mat33d=vector3::mat33<double>;
// short test program confirming that the
// spiral orders are picked up correctly.
//
// Initialises a spiral order by hand, computes the energy and exports a SSF

using namespace CMC;

using namespace std;
int main (int argc, char *argv[]) {
    /// Setup and CLI
    argparse::ArgumentParser prog(argv[0]);

    /// BOOK-KEEPING
    prog.add_argument("--output_dir", "-o")
        .help("Path to output")
        .required();


    /// PHYSICAL
    prog.add_argument("L")
        .help("Cubic side lenfth ")
        .scan<'i', int>();


    prog.add_argument("--J1")
        .help("Nearest-neighbour Heisenberg coupling strength")
        .default_value(-1.0)
        .scan<'g', double>();
    prog.add_argument("--J2")
        .help("Second-nearest-neighbour Heisenberg coupling strength")
        .default_value(0.)
        .scan<'g', double>();
    prog.add_argument("--J3")
        .help("Third-nearest-neighbour Heisenberg coupling strength")
        .default_value(0.)
        .scan<'g', double>();

    prog.add_argument("-Q")
        .help("Q of the spiral, specified as three integers (multiply the reciprocal lattice vectors)")
        .nargs(3)
        .scan<'i', int>()
        .default_value(std::vector<int>({0,0,1}));

    try {
        prog.parse_args(argc, argv);
    } catch (const std::exception& err){
        cerr << err.what() << endl;
        cerr << prog;
        std::exit(1);
    }


    /// Ensuring directories exist AHEAD of time (avoids heartbreak)
    std::string outdir_s = prog.get<std::string>("output_dir");
    filesystem::path outdir(outdir_s);
    if (! filesystem::exists(outdir) ){
        throw runtime_error("Cannot open outdir");
    }

    stringstream name;
    name << "ZZZ_spiraltest";

    auto L = prog.get<int>("L");

    const auto spec = DiamondC();
    CMC::Lattice lat(spec, imat33_t::from_cols({L,0,0},{0,L,0},{0,0,L})); 

    CMC::static_corr_3D ssf_manager(lat);

    cout<<"Primitive vectors: "<<lat.primitive_spec.latvecs <<"\n";
    cout<<"Cell vectors: "<<lat.cell_vectors <<"\n";
    cout<<"Index vectors: "<<lat.index_cell_vectors <<"\n";
    auto b = get_reciprocal_cell_vectors(lat.cell_vectors);
    cout<<"recip vectors: " << b <<std::endl;
    auto Q_int = prog.get<std::vector<int>>("Q");
    auto Q = b*vector3::vec3d(Q_int[0], Q_int[1], Q_int[2]);

    name << "?Q_2pi="<<(1./2/M_PI)*Q;

    // set the spin alignment by hand
    {
        for (const auto [link_id, l] : lat.links){
            l->S[0] = cos(dot(vector3::vec3d(l->position),Q));
            l->S[1] = sin(dot(vector3::vec3d(l->position),Q));
            l->S[2] = 0;
        }


    }


    ssf_manager.declare_observable("SdotS", 
            static_corr_3D::NEEDS_XX | static_corr_3D::NEEDS_YY | static_corr_3D::NEEDS_ZZ );
    ssf_manager.declare_observable("SzSz", static_corr_3D::NEEDS_ZZ);

    // one sample 
    ssf_manager.sample();
    
    ssf_manager.save(outdir/( name.str() + ".ssf.h5"));


    CMC::MC_runner mc(lat, 0);
    {
        double J1 = prog.get<double>("--J1");
        double J2 = prog.get<double>("--J2");
        double J3 = prog.get<double>("--J3");
        mc.define_coupling("J1", pyrochlore::nn1_dist, 
            mat33d::from_cols({J1, 0,0}, {0,J1, 0}, {0,0,J1})
            );
        mc.define_coupling("J2", pyrochlore::nn2_dist, 
            mat33d::from_cols({J2, 0,0}, {0,J2, 0}, {0,0,J2})
            );
        mc.define_coupling("J3a", pyrochlore::nn3a_dist, 
            mat33d::from_cols({J3, 0,0}, {0,J3, 0}, {0,0,J3})
            );
        mc.define_coupling("J3b", pyrochlore::nn3b_dist, 
            mat33d::from_cols({J3, 0,0}, {0,J3, 0}, {0,0,J3})
            );
        mc.setup_lattice();
    }




    
    printf(" E=%3e \n", mc.total_energy_per_unit_cell());

    auto f = outdir /( name.str() + ".spins.h5" );
    printf("Saving spin state to %s\n", f.string().c_str());
    save_spin_state(lat, f);
    
}
