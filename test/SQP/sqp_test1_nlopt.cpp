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
// SIMPLE 2D TEST CASE
//
// f(x,y) = (x - 1)^2 + (y - 2)^2
// c(x,y) = x + y - 1 = 0
//
// Exact solution: x* = (0,1)
// ================================================================

// g++ sqp_test1_nlopt.cpp -o sqp_test1_nlopt -lnlopt

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

    // Must match the OUR SQP file
    const double neighborhood_radius = 1.0;

    // Must match the OUR SQP file
    const unsigned int seed = 12345;

    const std::vector<double> reference_solution = {
        0.0,
        1.0
    };

    std::mt19937 gen(seed);

    std::uniform_real_distribution<double> dis(
        -neighborhood_radius,
         neighborhood_radius
    );

    // ============================================================
    // CSV OUTPUT
    // ============================================================

    std::ofstream csv_file("sqp_test1_nlopt.csv");

    if (!csv_file.is_open()) {
        std::cerr << "Unable to open sqp_test1_nlopt.csv\n";
        return 1;
    }

    csv_file
        << "TestId,OptId,NumIterOrEval,CompTime,"
        << "InitPoint,MinPoint,Fx,DistSol,Residual,Status"
        << std::endl;

    csv_file << std::setprecision(17);

    // ============================================================
    // RUN SIMULATIONS
    // ============================================================

    std::cout << "========================================================\n";
    std::cout << "NLOPT SLSQP - 2D TEST CASE - N SIMULATIONS\n";
    std::cout << "Number of simulations : "
              << num_simulations << "\n";
    std::cout << "Neighborhood radius   : "
              << neighborhood_radius << "\n";
    std::cout << "Seed                  : "
              << seed << "\n";
    std::cout << "Reference solution    : ("
              << reference_solution[0] << ", "
              << reference_solution[1] << ")\n";
    std::cout << "========================================================\n";

    for (int test_id = 0;
         test_id < num_simulations;
         ++test_id) {

        // --------------------------------------------------------
        // Generate the same starting point as the OUR SQP file
        // --------------------------------------------------------

        std::vector<double> x0(n);

        for (int j = 0; j < n; ++j) {
            x0[j] =
                reference_solution[j] + dis(gen);
        }

        // NLopt modifies the vector in place
        std::vector<double> x = x0;

        double f_min = 0.0;
        int num_eval = -1;

        nlopt::result result = nlopt::FAILURE;

        // --------------------------------------------------------
        // Solve and measure computation time
        // --------------------------------------------------------

        const auto start =
            std::chrono::steady_clock::now();

        try {
            nlopt::opt opt(
                nlopt::LD_SLSQP,
                n
            );

            opt.set_min_objective(
                objective,
                nullptr
            );

            opt.add_equality_constraint(
                equality_constraint,
                nullptr,
                1e-8
            );

            opt.set_xtol_rel(1e-8);
            opt.set_xtol_abs(1e-10);
            opt.set_ftol_rel(1e-10);
            opt.set_ftol_abs(1e-12);
            opt.set_maxeval(1000);

            result = opt.optimize(
                x,
                f_min
            );

            num_eval = opt.get_numevals();
        }
        catch (const std::exception& e) {
            std::cerr
                << "NLopt failed at test "
                << test_id + 1
                << ": "
                << e.what()
                << std::endl;

            result = nlopt::FAILURE;
        }

        const auto end =
            std::chrono::steady_clock::now();

        const auto elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                end - start
            );

        const long long elapsed_raw =
            elapsed.count();

        // --------------------------------------------------------
        // Compute metrics
        // --------------------------------------------------------

        const double dist_solution =
            std::sqrt(
                (x[0] - reference_solution[0])
                    * (x[0] - reference_solution[0])
              + (x[1] - reference_solution[1])
                    * (x[1] - reference_solution[1])
            );

        const double constraint_residual =
            std::abs(
                x[0] + x[1] - 1.0
            );

        // --------------------------------------------------------
        // Write CSV row
        // --------------------------------------------------------

        csv_file
            << test_id + 1 << ","
            << "NLOPT_SLSQP" << ","
            << num_eval << ","
            << elapsed_raw << ","
            << x0[0] << ";" << x0[1] << ","
            << x[0] << ";" << x[1] << ","
            << f_min << ","
            << dist_solution << ","
            << constraint_residual << ","
            << static_cast<int>(result)
            << std::endl;

        // --------------------------------------------------------
        // Print results
        // --------------------------------------------------------

        std::cout
            << "Test "
            << test_id + 1
            << "/"
            << num_simulations
            << "\n";

        std::cout
            << "Initial point : ("
            << x0[0] << ", "
            << x0[1] << ")\n";

        std::cout
            << "Optimal point : ("
            << x[0] << ", "
            << x[1] << ")\n";

        std::cout
            << "f(x)          : "
            << f_min << "\n";

        std::cout
            << "Distance      : "
            << dist_solution << "\n";

        std::cout
            << "Residual      : "
            << constraint_residual << "\n";

        std::cout
            << "Num eval      : "
            << num_eval << "\n";

        std::cout
            << "Status        : "
            << static_cast<int>(result)
            << "\n";

        std::cout
            << "Time [ns]     : "
            << elapsed_raw << "\n";

        std::cout
            << "--------------------------------------------------------\n";
    }

    csv_file.close();

    std::cout << "========================================================\n";
    std::cout << "CSV written to sqp_test1_nlopt.csv\n";
    std::cout << "========================================================\n";

    return 0;
}