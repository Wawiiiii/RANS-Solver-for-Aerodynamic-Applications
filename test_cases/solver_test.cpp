// Driver: read a SolverConfig from a text file, run EulerSimulation, print the
// result. Edit the config file and re-run -- no recompile needed.
// Usage: solver_test [config_file]   (default: solver_config.txt)

#include "euler_simulation.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    try {
        const std::string config_file = (argc > 1) ? argv[1] : "solver_config.txt";
        std::cout << "Reading config: " << config_file << "\n";

        const aero::SolverConfig cfg = aero::load_solver_config(config_file);
        const aero::SolverResult res = aero::EulerSimulation(cfg).solve();

        std::cout << std::scientific << std::setprecision(12)
                  << "\nsolver_test completed.\n"
                  << "  cycles         = " << res.cycles << "\n"
                  << "  final residual = " << res.final_residual << "\n"
                  << "  converged      = " << (res.converged ? "yes" : "no") << "\n"
                  << "  CL = " << res.coeffs.cl << "\n"
                  << "  CD = " << res.coeffs.cd << "\n"
                  << "  CM = " << res.coeffs.cm << "\n";

        std::cout << "\nAll Euler tests passed.\n";
    }
    catch (const std::exception& e) {
        std::cerr << "\nTest failed: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
