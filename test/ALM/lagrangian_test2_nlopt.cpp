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
// 5D TEST CASE
//
// f(x) = (x0 - 1)^2 + (x1 - 2)^2 + (x2 - 3)^2
//      + (x3 - 2)^2 + (x4 - 1)^2
//
// c1(x) = x0 + x1 + x2 - 3 = 0
// c2(x) = x0^2 + x3 - 1 = 0
// c3(x) = x1*x4 - 2 = 0
// ================================================================

// g++ lagrangian_test2_nlopt.cpp -o lagrangian_test2_nlopt -lnlopt

double objective(
    const std::vector<double>& x,
    std::vector<double>& grad,
    void* data
) {
    if (!grad.empty()) {
        grad[0] = 2.0 * (x[0] - 1.0);
        grad[1] = 2.0 * (x[1] - 2.0);
        grad[2] = 2.0 * (x[2] - 3.0);
        grad[3] = 2.0 * (x[3] - 2.0);
        grad[4] = 2.0 * (x[4] - 1.0);
    }

    return (x[0] - 1.0)*(x[0] - 1.0)
         + (x[1] - 2.0)*(x[1] - 2.0)
         + (x[2] - 3.0)*(x[2] - 3.0)
         + (x[3] - 2.0)*(x[3] - 2.0)
         + (x[4] - 1.0)*(x[4] - 1.0);
}

double equality_constraint_1(
    const std::vector<double>& x,
    std::vector<double>& grad,
    void* data
) {
    if (!grad.empty()) {
        grad[0] = 1.0;
        grad[1] = 1.0;
        grad[2] = 1.0;
        grad[3] = 0.0;
        grad[4] = 0.0;
    }

    return x[0] + x[1] + x[2] - 3.0;
}

double equality_constraint_2(
    const std::vector<double>& x,
    std::vector<double>& grad,
    void* data
) {
    if (!grad.empty()) {
        grad[0] = 2.0 * x[0];
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 1.0;
        grad[4] = 0.0;
    }

    return x[0]*x[0] + x[3] - 1.0;
}

double equality_constraint_3(
    const std::vector<double>& x,
    std::vector<double>& grad,
    void* data
) {
    if (!grad.empty()) {
        grad[0] = 0.0;
        grad[1] = x[4];
        grad[2] = 0.0;
        grad[3] = 0.0;
        grad[4] = x[1];
    }

    return x[1]*x[4] - 2.0;
}

int main() {
    // ============================================================
    // SIMULATION PARAMETERS
    // ============================================================

    constexpr int n = 5;

    const int num_simulations = 50;

    // Must be equal to the OUR_ALM file.
    const double neighborhood_radius = 1.0;

    // Must be equal to the OUR_ALM file.
    const unsigned int seed = 12345;

    // Reference solution
    std::vector<double> reference_solution = {
        -0.08170078635350764,
         1.3278938903564708,
         1.7538068959970370,
         0.9933249815092187,
         1.5061444400976221
    };

    // Random generator for starting points.
    // Important: same seed, same distribution, same generation order
    // as in the OUR_ALM file.
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(
        -neighborhood_radius,
         neighborhood_radius
    );

    // ============================================================
    // CSV OUTPUT
    // ============================================================

    std::ofstream csv_file("lagrangian_test2_nlopt.csv");

    csv_file << "TestId,OptId,NumIterOrEval,CompTime,InitPoint,MinPoint,Fx,DistSol,Residual,Status"
             << std::endl;

    csv_file << std::setprecision(17);

    // ============================================================
    // RUN SIMULATIONS
    // ============================================================

    std::cout << "========================================================\n";
    std::cout << "NLOPT AUGLAG - 5D TEST CASE - N SIMULATIONS\n";
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
            local_opt.set_maxeval(2000);

            opt.set_local_optimizer(local_opt);

            opt.set_min_objective(objective, nullptr);

            opt.add_equality_constraint(
                equality_constraint_1,
                nullptr,
                1e-8
            );

            opt.add_equality_constraint(
                equality_constraint_2,
                nullptr,
                1e-8
            );

            opt.add_equality_constraint(
                equality_constraint_3,
                nullptr,
                1e-8
            );

            opt.set_xtol_rel(1e-8);
            opt.set_ftol_abs(1e-12);
            opt.set_maxeval(2000);

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

        double c1 = x[0] + x[1] + x[2] - 3.0;
        double c2 = x[0]*x[0] + x[3] - 1.0;
        double c3 = x[1]*x[4] - 2.0;

        double constraint_residual =
            std::sqrt(c1*c1 + c2*c2 + c3*c3);

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
    std::cout << "CSV written to lagrangian_test2_nlopt.csv\n";
    std::cout << "========================================================\n";

    return 0;
}