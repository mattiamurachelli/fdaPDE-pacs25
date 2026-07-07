#include <nlopt.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cmath>
#include <exception>
#include <random>
#include <string>

// ================================================================
// SAME 2D TEST CASE AS OUR ALM FILE
//
// f(x,y) = (x - 1)^2 + (y - 2)^2
// c(x,y) = x + y - 1 = 0
//
// Exact solution: x* = (0, 1)
// ================================================================

// g++ lagrangian_test1_nlopt.cpp -o lagrangian_test1_nlopt -lnlopt

double objective(
    const std::vector<double>& x,
    std::vector<double>& grad,
    void* data
) {
    if (!grad.empty()) {
        grad[0] = 2.0 * (x[0] - 1.0);
        grad[1] = 2.0 * (x[1] - 2.0);
    }

    return (x[0] - 1.0) * (x[0] - 1.0)
         + (x[1] - 2.0) * (x[1] - 2.0);
}

double equality_constraint(
    const std::vector<double>& x,
    std::vector<double>& grad,
    void* data
) {
    if (!grad.empty()) {
        grad[0] = 1.0;
        grad[1] = 1.0;
    }

    return x[0] + x[1] - 1.0;
}

int main() {
    // ============================================================
    // SIMULATION PARAMETERS
    // ============================================================

    constexpr int n = 2;

    const int num_simulations = 50;

    // Controls the size of the neighborhood around the optimum.
    // Starting points are generated inside the square:
    // [x*_0 - radius, x*_0 + radius] x [x*_1 - radius, x*_1 + radius]
    const double neighborhood_radius = 1.0;

    // Same seed as in the OUR ALM file.
    const unsigned int seed = 12345;

    // Exact solution
    std::vector<double> exact_solution = {0.0, 1.0};

    // Random generator for starting points.
    // Important: use the same generator, seed, distribution and order
    // as in the OUR ALM file.
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(
        -neighborhood_radius,
         neighborhood_radius
    );

    // ============================================================
    // CSV OUTPUT
    // ============================================================

    std::ofstream csv_file("lagrangian_test1_nlopt.csv");

    csv_file << "TestId,OptId,NumIterOrEval,CompTime,InitPoint,MinPoint,Fx,DistSol,Residual,Status"
             << std::endl;

    // ============================================================
    // RUN SIMULATIONS
    // ============================================================

    std::cout << "========================================================\n";
    std::cout << "NLOPT AUGLAG - N SIMULATIONS\n";
    std::cout << "Number of simulations      : " << num_simulations << "\n";
    std::cout << "Neighborhood radius        : " << neighborhood_radius << "\n";
    std::cout << "Seed                       : " << seed << "\n";
    std::cout << "Exact solution             : ("
              << exact_solution[0] << ", "
              << exact_solution[1] << ")\n";
    std::cout << "========================================================\n";

    for (int test_id = 0; test_id < num_simulations; ++test_id) {

        // ------------------------------------------------------------
        // Generate initial point around the optimum
        // ------------------------------------------------------------

        std::vector<double> x0 = exact_solution;

        x0[0] = exact_solution[0] + dis(gen);
        x0[1] = exact_solution[1] + dis(gen);

        // NLopt modifies x in-place, so keep x0 for CSV output.
        std::vector<double> x = x0;

        double f_min = 0.0;
        nlopt::result result = nlopt::FAILURE;
        int num_eval = -1;

        // ------------------------------------------------------------
        // Solve and measure computation time
        // ------------------------------------------------------------

        auto start = std::chrono::steady_clock::now();

        try {
            // Outer Augmented Lagrangian algorithm
            nlopt::opt opt(nlopt::AUGLAG, n);

            // Local optimizer used by AUGLAG to solve the auxiliary subproblems
            nlopt::opt local_opt(nlopt::LD_LBFGS, n);

            local_opt.set_xtol_rel(1e-8);
            local_opt.set_ftol_abs(1e-12);
            local_opt.set_maxeval(1000);

            opt.set_local_optimizer(local_opt);

            // Objective function
            opt.set_min_objective(objective, nullptr);

            // Equality constraint: x + y - 1 = 0
            opt.add_equality_constraint(
                equality_constraint,
                nullptr,
                1e-8
            );

            opt.set_xtol_rel(1e-8);
            opt.set_ftol_abs(1e-12);
            opt.set_maxeval(1000);

            result = opt.optimize(x, f_min);

            // Number of objective evaluations performed by NLopt
            num_eval = opt.get_numevals();
        }
        catch (const std::exception& e) {
            std::cerr << "NLopt failed at test "
                      << test_id + 1
                      << ": "
                      << e.what()
                      << std::endl;

            result = nlopt::FAILURE;
        }

        auto end = std::chrono::steady_clock::now();

        auto elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                end - start
            );

        long long elapsed_raw = elapsed.count();

        // ------------------------------------------------------------
        // Extract results
        // ------------------------------------------------------------

        double dist_solution =
            std::sqrt(
                (x[0] - exact_solution[0]) * (x[0] - exact_solution[0])
                +
                (x[1] - exact_solution[1]) * (x[1] - exact_solution[1])
            );

        double constraint_residual =
            std::abs(x[0] + x[1] - 1.0);

        // ------------------------------------------------------------
        // Write CSV line
        // ------------------------------------------------------------

        csv_file << test_id + 1 << ","
                 << "NLOPT_AUGLAG" << ","
                 << num_eval << ","
                 << elapsed_raw << ","
                 << x0[0] << ";" << x0[1] << ","
                 << x[0] << ";" << x[1] << ","
                 << f_min << ","
                 << dist_solution << ","
                 << constraint_residual << ","
                 << static_cast<int>(result)
                 << std::endl;

        // ------------------------------------------------------------
        // Print compact information to screen
        // ------------------------------------------------------------

        std::cout << "Test " << test_id + 1 << "/" << num_simulations << "\n";
        std::cout << "Initial point : (" << x0[0] << ", " << x0[1] << ")\n";
        std::cout << "Optimal point : (" << x[0] << ", " << x[1] << ")\n";
        std::cout << "f(x)          : " << f_min << "\n";
        std::cout << "Distance      : " << dist_solution << "\n";
        std::cout << "Residual      : " << constraint_residual << "\n";
        std::cout << "Num eval      : " << num_eval << "\n";
        std::cout << "Status        : " << static_cast<int>(result) << "\n";
        std::cout << "Time [ns]     : " << elapsed_raw << "\n";
        std::cout << "--------------------------------------------------------\n";
    }

    csv_file.close();

    std::cout << "========================================================\n";
    std::cout << "CSV written to lagrangian_test1_nlopt.csv\n";
    std::cout << "========================================================\n";

    return 0;
}