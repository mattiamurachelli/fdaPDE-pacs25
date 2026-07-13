// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <fdaPDE/core/fdaPDE/optimization.h>
#include <unsupported/Eigen/SparseExtra>
#include "../INCLUDE/obj_constr.h"

#include <cstdio>
#include <vector>
#include <fstream>
#include <chrono>
#include <cmath>
#include <string>
#include <random>
#include <iomanip>

int main() {
    // ============================================================
    // SIMPLE 2D TEST CASE
    //
    // f(x,y) = (x - 1)^2 + (y - 2)^2
    // c(x,y) = x + y - 1 = 0
    //
    // Exact solution: x* = (0,1)
    // ============================================================

    constexpr int n = 2;

    // ============================================================
    // SIMULATION PARAMETERS
    // ============================================================

    const int num_simulations = 50;

    // Starting points are generated inside the square
    // [x*_0 - radius, x*_0 + radius]
    // x
    // [x*_1 - radius, x*_1 + radius]
    const double neighborhood_radius = 1.0;

    // Use the same seed in the NLopt file
    const unsigned int seed = 12345;

    Eigen::Matrix<double, n, 1> reference_solution;
    reference_solution << 0.0, 1.0;

    std::mt19937 gen(seed);

    std::uniform_real_distribution<double> dis(
        -neighborhood_radius,
         neighborhood_radius
    );

    // ============================================================
    // OBJECTIVE FUNCTION
    // ============================================================

    ScalarField<n> objective;

    objective = [](Eigen::Matrix<double, n, 1> x) -> double {
        return (x[0] - 1.0) * (x[0] - 1.0)
             + (x[1] - 2.0) * (x[1] - 2.0);
    };

    // Exact objective gradient
    auto objective_gradient =
        [](const Eigen::Matrix<double, n, 1>& x)
        -> Eigen::Matrix<double, n, 1> {

        Eigen::Matrix<double, n, 1> g;

        g << 2.0 * (x[0] - 1.0),
             2.0 * (x[1] - 2.0);

        return g;
    };

    // ============================================================
    // EQUALITY CONSTRAINT
    // c(x) = x + y - 1 = 0
    // ============================================================

    ScalarField<n> constraint;

    constraint = [](Eigen::Matrix<double, n, 1> x) -> double {
        return x[0] + x[1] - 1.0;
    };

    struct ConstraintGradient {
        using vector_t = Eigen::Matrix<double, n, 1>;

        vector_t operator()(const vector_t& x) const {
            vector_t g;
            g << 1.0, 1.0;
            return g;
        }
    };

    constrFunction<n> constr_func(
        constraint,
        ConstraintGradient{},
        false
    );

    constrList<n> constraints = {constr_func};

    // ============================================================
    // CSV OUTPUT
    // ============================================================

    std::ofstream csv_file("SQP/sqp_test1.csv");

    if (!csv_file.is_open()) {
        std::cerr << "Unable to open SQP/sqp_test1.csv\n";
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

    printf("========================================================\n");
    printf("OUR SQP - 2D TEST CASE - N SIMULATIONS\n");
    printf("Number of simulations : %d\n", num_simulations);
    printf("Neighborhood radius   : %.4f\n", neighborhood_radius);
    printf("Seed                  : %u\n", seed);
    printf("Reference solution    : (%.6f, %.6f)\n",
           reference_solution[0],
           reference_solution[1]);
    printf("========================================================\n");

    for (int test_id = 0; test_id < num_simulations; ++test_id) {

        // --------------------------------------------------------
        // Generate starting point around the exact solution
        // --------------------------------------------------------

        Eigen::Matrix<double, n, 1> x0;

        for (int j = 0; j < n; ++j) {
            x0[j] = reference_solution[j] + dis(gen);
        }

        // New solver for every independent simulation
        fdapde::SQP<n> problem;

        // Give the solver the exact objective gradient
        objFunction<n> obj_func(
            objective,
            objective_gradient
        );

        // --------------------------------------------------------
        // Solve and measure computation time
        // --------------------------------------------------------

        const auto start =
            std::chrono::steady_clock::now();

        problem.solve(
            obj_func,
            constraints,
            x0
        );

        const auto end =
            std::chrono::steady_clock::now();

        const auto elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                end - start
            );

        const long long elapsed_raw = elapsed.count();

        // --------------------------------------------------------
        // Extract results
        // --------------------------------------------------------

        const auto& final_point =
            problem.optimum().back();

        const double final_value =
            problem.values().back();

        const double dist_solution =
            (final_point - reference_solution).norm();

        const double constraint_residual =
            std::abs(
                final_point[0]
              + final_point[1]
              - 1.0
            );

        const std::string status = "OK";

        // --------------------------------------------------------
        // Write CSV row
        // --------------------------------------------------------

        csv_file
            << test_id + 1 << ","
            << "OUR_SQP" << ","
            << problem.num_iter().size() << ","
            << elapsed_raw << ","
            << x0[0] << ";" << x0[1] << ","
            << final_point[0] << ";" << final_point[1] << ","
            << final_value << ","
            << dist_solution << ","
            << constraint_residual << ","
            << status
            << std::endl;

        // --------------------------------------------------------
        // Print results
        // --------------------------------------------------------

        printf("Test %d/%d\n",
               test_id + 1,
               num_simulations);

        printf("Initial point : (%.6f, %.6f)\n",
               x0[0],
               x0[1]);

        printf("Optimal point : (%.10f, %.10f)\n",
               final_point[0],
               final_point[1]);

        printf("f(x)          : %.10f\n",
               final_value);

        printf("Distance      : %.4e\n",
               dist_solution);

        printf("Residual      : %.4e\n",
               constraint_residual);

        printf("SQP steps     : %zu\n",
               problem.num_iter().size());

        printf("Time [ns]     : %lld\n",
               elapsed_raw);

        printf("--------------------------------------------------------\n");
    }

    csv_file.close();

    printf("========================================================\n");
    printf("CSV written to SQP/sqp_test1.csv\n");
    printf("========================================================\n");

    return 0;
}