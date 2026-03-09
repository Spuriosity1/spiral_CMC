#include "vec3.hpp"
#include "timeit.hpp"   
#include <random>

using namespace std;
using namespace vector3;

volatile double dump;



int main (int argc, char *argv[]) {
    
    int n_reps = atoi(argv[1]);
    std::mt19937 rng;
    std::normal_distribution<double> normdist;

    std::vector<vec3d> u, v;
    std::vector<mat33<double>> M;
    for (int i=0; i<n_reps; i++) {
        u.emplace_back(normdist(rng), normdist(rng), normdist(rng));
        v.emplace_back(normdist(rng), normdist(rng), normdist(rng));
        M.emplace_back(mat33<double>::from_cols(
                    {normdist(rng), normdist(rng), normdist(rng)},
                    {normdist(rng), normdist(rng), normdist(rng)},
                    {normdist(rng), normdist(rng), normdist(rng)}
                    ));

    }
    

    TIMEIT("Benchmark: vector-vector dot",
        for (int i=0; i<n_reps; i++) {
            dump = dot(u[i], v[i]);
        }
    )


    TIMEIT("Benchmark: vector-matrix-vector",
        for (int i=0; i<n_reps; i++) {
            dump = dot(u[i], M[i]* v[i]);
        }
    )

    
    return 0;
}
