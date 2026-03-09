#pragma once

#include "vec3.hpp"
#include <argparse/argparse.hpp>
#include "format_bits.hpp"


using mat33d=vector3::mat33<double>;
namespace coupling {
static const mat33d Heis = mat33d::from_cols(
        {1,0,0},{0,1,0},{0,0,1});

static const mat33d Isin = mat33d::from_cols(
        {0,0,0},{0,0,0},{0,0,1});
};


inline void add_bookkeeping_args(argparse::ArgumentParser& prog){

    prog.add_argument("--output_dir", "-o")
        .help("Path to output")
        .required();

    prog.add_argument("--seed", "-s")
        .required()
        .help("64-bit int to seed the RNG");
//        .store_into(seed_s);

    prog.add_argument("--save_state")
        .implicit_value(true)
        .default_value(false);

}

inline void add_annealer_args(argparse::ArgumentParser& prog){
    /// ANNEALING PROTOCOL
    prog.add_argument("--T_hot")
        .help("Temperature  to begin annealing from")
        .scan<'g',double>();
    prog.add_argument("--T_ref")
        .help("Reference temperature for proposal distribution")
        .default_value(1.0)
        .scan<'g',double>();
    prog.add_argument("--T_cold")
        .help("Final temperature")
        .required()
        .scan<'g',double>();
    prog.add_argument("--n_steps")
        .help("Number of annealing steps")
        .default_value(static_cast<size_t>(100))
        .scan<'i', size_t>();
    prog.add_argument("--n_sweep")
        .help("Number of sweeps to run per temperature step")
        .default_value(static_cast<size_t>(16))
        .scan<'i', size_t>();
    prog.add_argument("--n_burn_in")
        .help("Number of sweeps to run at T_hot before we start annealing")
        .default_value(static_cast<size_t>(16))
        .scan<'i', size_t>();
    prog.add_argument("--n_sample")
        .default_value(static_cast<size_t>(64))
        .help("Number of sweeps to run at T_cold while collecting statistics")
        .scan<'i', size_t>();
}

inline uint64_t get_hex(argparse::ArgumentParser& prog, const std::string& name){

    uint64_t seed; // ugly hack for loading seed as hex
    std::stringstream ss;
    std::string seed_s = prog.get<std::string>(name);
    ss << std::hex << seed_s;
    ss >> seed; 
    return seed;
}

inline void add_J1J2J3_args(argparse::ArgumentParser& prog){
    prog.add_argument("--J1")
        .help("Nearest-neighbour Heisenberg coupling strength")
        .required()
        .scan<'g', double>();
    prog.add_argument("--J2")
        .help("Second-nearest-neighbour Heisenberg coupling strength")
        .default_value(0.)
        .scan<'g', double>();
    prog.add_argument("--J3")
        .help("Third-nearest-neighbour Heisenberg coupling strength")
        .default_value(0.)
        .scan<'g', double>();

    prog.add_argument("L")
        .help("Linear dimension of cubic supercell (e.g. L=2 has 2^3 *4 = 32 primitive cells)")
        .scan<'i', int>(); // cubic unit cell hardcoded
}

inline std::string get_ofilename(const std::string& root, const argparse::ArgumentParser& prog) {

    std::stringstream name;
    name << root<<DELIM; // accumulates hashed options
                            
    auto B_tmp = prog.get<std::vector<double>>("-B");

    auto J1 = prog.get<double>("--J1");
    auto J2 = prog.get<double>("--J2");
    auto J3 = prog.get<double>("--J3");
    const double T_cold = prog.get<double>("--T_cold");

    name <<"L="<<prog.get<int>("L")<<DELIM<<
        "B="<<B_tmp[0]<<"_"<<B_tmp[1]<<"_"<<B_tmp[2]<<DELIM<<
        "J1="<<J1<<DELIM<<
        "J2="<<J2<<DELIM<<
        "J3="<<J3<<DELIM<<
        "seed="<<prog.get<std::string>("--seed")<<DELIM<<
        "T_c="<<T_cold<<DELIM;

    return name.str();

}

