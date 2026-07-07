#include <nlopt.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cmath>
#include <exception>
#include <random>
#include <string>
#include <iomanip>

// ================================================================
// 4D TEST CASE WITH EQUALITY AND INEQUALITY CONSTRAINTS
//
// f(x) = (x0 - 1)^2 + (x1 - 2)^2 + (x2 + 1)^2 + (x3 - 3)^2
//      + sin(x0*x2) + cos(x1*x3)
//
// c1(x) = x0^2 + x1^2 + x2 - 3 = 0
// c2(x) = x0*x1 + x3^2 - 2 <= 0
// c3(x) = x1^2 + x2^2 + x0 - 4 <= 0
// ================================================================

// g++ lagrangian_test4_nlopt.cpp -o lagrangian_test4_nlopt -lnlopt

double objective(
    const std::vector<double>& x,
    std::vector<double>& grad,
    void* data
) {
    if (!grad.empty()) {
        grad[0] = 2.0 * (x[0] - 1.0)
                + x[2] * std::cos(x[0] * x[2]);

        grad[1] = 2.0 * (x[1] - 2.0)
                - x[3] * std::sin(x[1] * x[3]);

        grad[2] = 2.0 * (x[2] + 1.0)
                + x[0] * std::cos(x[0] * x[2]);

        grad[3] = 2.0 * (x[3] - 3.0)
                - x[1] * std::sin(x[1] * x[3]);
    }

    return (x[0] - 1.0) * (x[0] - 1.0)
         + (x[1] - 2.0) * (x[1] - 2.0)
         + (x[2] + 1.0) * (x[2] + 1.0)
         + (x[3] - 3.0) * (x[3] - 3.0)
         + std::sin(x[0] * x[2])
         + std::cos(x[1] * x[3]);
}

double equality_constraint_1(
    const std::vector<double>& x,
    std::vector<double>& grad,
    void* data
) {
    if (!grad.empty()) {
        grad[0] = 2.0 * x[0];
        grad[1] = 2.0 * x[1];
        grad[2] = 1.0;
        grad[3] = 0.0;
    }

    return x[0] * x[0] + x[1] * x[1] + x[2] - 3.0;
}

double inequality_constraint_2(
    const std::vector<double>& x,
    std::vector<double>& grad,
    void* data
) {
    if (!grad.empty()) {
        grad[0] = x[1];
        grad[1] = x[0];
        grad[2] = 0.0;
        grad[3] = 2.0 * x[3];
    }

    return x[0] * x[1] + x[3] * x[3] - 2.0;
}

double inequality_constraint_3(
    const std::vector<double>& x,
    std::vector<double>& grad,
    void* data
) {
    if (!grad.empty()) {
        grad[0] = 1.0;
        grad[1] = 2.0 * x[1];
        grad[2] = 2.0 * x[2];
        grad[3] = 0.0;
    }

    return x[1] * x[1] + x[2] * x[2] + x[0] - 4.0;
}

int main() {
    // ============================================================
    // SIMULATION PARAMETERS
    // ============================================================

    constexpr int n = 4;

    const int num_simulations = 50;

    // Must be equal to the OUR_ALM file.
    const double neighborhood_radius = 1.0;

    // Must be equal to the OUR_ALM file.
    const unsigned int seed = 12345;

    // Reference solution
    std::vector<double> reference_solution = {
        -0.07521872,
         1.90996816,
        -0.65363624,
         1.46412614
    };

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(
        -neighborhood_radius,
         neighborhood_radius
    );

    // ============================================================
    // CSV OUTPUT
    // ============================================================

    std::ofstream csv_file("lagrangian_test4_nlopt.csv");

    csv_file << "TestId,OptId,NumIterOrEval,CompTime,InitPoint,MinPoint,Fx,DistSol,Residual,Status"
             << std::endl;

    csv_file << std::setprecision(17);

    // ============================================================
    // RUN SIMULATIONS
    // ============================================================

    std::cout << "========================================================\n";
    std::cout << "NLOPT AUGLAG - 4D MIXED-CONSTRAINT TEST - N SIMULATIONS\n";
    std::cout << "Number of simulations      : " << num_simulations << "\n";
    std::cout << "Neighborhood radius        : " << neighborhood_radius << "\n";
    std::cout << "Seed                       : " << seed << "\n";
    std::cout << "Reference solution         : ";
    for (int j = 0; j < n; ++j) {
        std::cout << reference_solution[j] << " ";
    }
    std::cout << "\n";
    std::cout << "========================================================\n";

    for (int test_id = 0; test_id < num_simulations; ++test_id) {

        // ------------------------------------------------------------
        // Generate initial point around the reference solution
        // ------------------------------------------------------------

        std::vector<double> x0(n);

        for (int j = 0; j < n; ++j) {
            x0[j] = reference_solution[j] + dis(gen);
        }

        // NLopt modifies x in-place.
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

            // Local optimizer used by AUGLAG
            nlopt::opt local_opt(nlopt::LD_LBFGS, n);

            local_opt.set_xtol_rel(1e-8);
            local_opt.set_ftol_abs(1e-12);
            local_opt.set_maxeval(3000);

            opt.set_local_optimizer(local_opt);

            opt.set_min_objective(objective, nullptr);

            // Equality constraint: c1(x) = 0
            opt.add_equality_constraint(
                equality_constraint_1,
                nullptr,
                1e-8
            );

            // Inequality constraints: c_i(x) <= 0
            opt.add_inequality_constraint(
                inequality_constraint_2,
                nullptr,
                1e-8
            );

            opt.add_inequality_constraint(
                inequality_constraint_3,
                nullptr,
                1e-8
            );

            opt.set_xtol_rel(1e-8);
            opt.set_ftol_abs(1e-12);
            opt.set_maxeval(3000);

            result = opt.optimize(x, f_min);

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

        double dist_solution = 0.0;

        for (int j = 0; j < n; ++j) {
            dist_solution +=
                (x[j] - reference_solution[j])
              * (x[j] - reference_solution[j]);
        }

        dist_solution = std::sqrt(dist_solution);

        double c1 =
            x[0] * x[0]
          + x[1] * x[1]
          + x[2]
          - 3.0;

        double c2 =
            x[0] * x[1]
          + x[3] * x[3]
          - 2.0;

        double c3 =
            x[1] * x[1]
          + x[2] * x[2]
          + x[0]
          - 4.0;

        double constraint_residual =
            std::sqrt(
                c1*c1
              + std::max(0.0, c2) * std::max(0.0, c2)
              + std::max(0.0, c3) * std::max(0.0, c3)
            );

        // ------------------------------------------------------------
        // Write CSV line
        // ------------------------------------------------------------

        csv_file << test_id + 1 << ","
                 << "NLOPT_AUGLAG" << ","
                 << num_eval << ","
                 << elapsed_raw << ",";

        for (int j = 0; j < n; ++j) {
            if (j > 0) csv_file << ";";
            csv_file << x0[j];
        }

        csv_file << ",";

        for (int j = 0; j < n; ++j) {
            if (j > 0) csv_file << ";";
            csv_file << x[j];
        }

        csv_file << ","
                 << f_min << ","
                 << dist_solution << ","
                 << constraint_residual << ","
                 << static_cast<int>(result)
                 << std::endl;

        // ------------------------------------------------------------
        // Print compact information to screen
        // ------------------------------------------------------------

        std::cout << "Test " << test_id + 1 << "/" << num_simulations << "\n";

        std::cout << "Initial point : ";
        for (int j = 0; j < n; ++j) {
            std::cout << x0[j] << " ";
        }
        std::cout << "\n";

        std::cout << "Optimal point : ";
        for (int j = 0; j < n; ++j) {
            std::cout << x[j] << " ";
        }
        std::cout << "\n";

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
    std::cout << "CSV written to lagrangian_test4_nlopt.csv\n";
    std::cout << "========================================================\n";

    return 0;
}