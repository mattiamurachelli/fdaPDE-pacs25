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
#include <iostream>
#include <algorithm>

int main() {
    // ============================================================
    // 4D TEST CASE WITH EQUALITY AND INEQUALITY CONSTRAINTS
    //
    // f(x) = (x0 - 1)^2 + (x1 - 2)^2 + (x2 + 1)^2
    //      + (x3 - 3)^2 + sin(x0*x2) + cos(x1*x3)
    //
    // c1(x) = x0^2 + x1^2 + x2 - 3 = 0
    // c2(x) = x0*x1 + x3^2 - 2 <= 0
    // c3(x) = x1^2 + x2^2 + x0 - 4 <= 0
    // ============================================================

    constexpr int n = 4;

    // ============================================================
    // SIMULATION PARAMETERS
    // ============================================================

    const int num_simulations = 50;

    // Each component of the starting point is generated inside:
    // [x*_j - neighborhood_radius, x*_j + neighborhood_radius].
    const double neighborhood_radius = 1.0;

    // Must be identical in the corresponding NLopt file.
    const unsigned int seed = 12345;

    // Numerical reference solution.
    Eigen::Matrix<double, n, 1> reference_solution;

    reference_solution <<
        -0.0752187183004026,
         1.9099681631540490,
        -0.6536362398447993,
         1.4641261411596933;

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
        return
            (x[0] - 1.0) * (x[0] - 1.0)
          + (x[1] - 2.0) * (x[1] - 2.0)
          + (x[2] + 1.0) * (x[2] + 1.0)
          + (x[3] - 3.0) * (x[3] - 3.0)
          + std::sin(x[0] * x[2])
          + std::cos(x[1] * x[3]);
    };

    // Exact objective gradient.
    auto objective_gradient =
        [](const Eigen::Matrix<double, n, 1>& x)
        -> Eigen::Matrix<double, n, 1> {

        Eigen::Matrix<double, n, 1> g;

        g[0] =
            2.0 * (x[0] - 1.0)
          + x[2] * std::cos(x[0] * x[2]);

        g[1] =
            2.0 * (x[1] - 2.0)
          - x[3] * std::sin(x[1] * x[3]);

        g[2] =
            2.0 * (x[2] + 1.0)
          + x[0] * std::cos(x[0] * x[2]);

        g[3] =
            2.0 * (x[3] - 3.0)
          - x[1] * std::sin(x[1] * x[3]);

        return g;
    };

    // ============================================================
    // CONSTRAINT 1
    // c1(x) = x0^2 + x1^2 + x2 - 3 = 0
    // ============================================================

    ScalarField<n> constraint1;

    constraint1 = [](Eigen::Matrix<double, n, 1> x) -> double {
        return
            x[0] * x[0]
          + x[1] * x[1]
          + x[2]
          - 3.0;
    };

    struct Grad1 {
        using vector_t = Eigen::Matrix<double, n, 1>;

        vector_t operator()(const vector_t& x) const {
            vector_t g;

            g <<
                2.0 * x[0],
                2.0 * x[1],
                1.0,
                0.0;

            return g;
        }
    };

    constrFunction<n> constr_1(
        constraint1,
        Grad1{},
        false
    );

    // ============================================================
    // CONSTRAINT 2
    // c2(x) = x0*x1 + x3^2 - 2 <= 0
    // ============================================================

    ScalarField<n> constraint2;

    constraint2 = [](Eigen::Matrix<double, n, 1> x) -> double {
        return
            x[0] * x[1]
          + x[3] * x[3]
          - 2.0;
    };

    struct Grad2 {
        using vector_t = Eigen::Matrix<double, n, 1>;

        vector_t operator()(const vector_t& x) const {
            vector_t g;

            g <<
                x[1],
                x[0],
                0.0,
                2.0 * x[3];

            return g;
        }
    };

    constrFunction<n> constr_2(
        constraint2,
        Grad2{},
        true
    );

    // ============================================================
    // CONSTRAINT 3
    // c3(x) = x1^2 + x2^2 + x0 - 4 <= 0
    // ============================================================

    ScalarField<n> constraint3;

    constraint3 = [](Eigen::Matrix<double, n, 1> x) -> double {
        return
            x[1] * x[1]
          + x[2] * x[2]
          + x[0]
          - 4.0;
    };

    struct Grad3 {
        using vector_t = Eigen::Matrix<double, n, 1>;

        vector_t operator()(const vector_t& x) const {
            vector_t g;

            g <<
                1.0,
                2.0 * x[1],
                2.0 * x[2],
                0.0;

            return g;
        }
    };

    constrFunction<n> constr_3(
        constraint3,
        Grad3{},
        true
    );

    constrList<n> constraints = {
        constr_1,
        constr_2,
        constr_3
    };

    // ============================================================
    // CSV OUTPUT
    // ============================================================

    std::ofstream csv_file("SQP/sqp_test4.csv");

    if (!csv_file.is_open()) {
        std::cerr
            << "Unable to open SQP/sqp_test4.csv"
            << std::endl;

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
    printf("OUR SQP - 4D MIXED-CONSTRAINT TEST - N SIMULATIONS\n");
    printf("Number of simulations : %d\n", num_simulations);
    printf("Neighborhood radius   : %.4f\n", neighborhood_radius);
    printf("Seed                  : %u\n", seed);

    printf("Reference solution    : ");
    for (int j = 0; j < n; ++j) {
        printf("%.8f ", reference_solution[j]);
    }
    printf("\n");

    printf("========================================================\n");

    for (int test_id = 0;
         test_id < num_simulations;
         ++test_id) {

        // --------------------------------------------------------
        // Generate initial point around the reference solution
        // --------------------------------------------------------

        Eigen::Matrix<double, n, 1> x0;

        for (int j = 0; j < n; ++j) {
            x0[j] =
                reference_solution[j] + dis(gen);
        }

        // Fresh SQP object for each simulation.
        fdapde::SQP<n> problem;

        // Supply the analytical objective gradient.
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

        const long long elapsed_raw =
            elapsed.count();

        // --------------------------------------------------------
        // Extract results
        // --------------------------------------------------------

        const auto& final_point =
            problem.optimum().back();

        const double final_value =
            problem.values().back();

        const double dist_solution =
            (final_point - reference_solution).norm();

        const double c1 =
            final_point[0] * final_point[0]
          + final_point[1] * final_point[1]
          + final_point[2]
          - 3.0;

        const double c2 =
            final_point[0] * final_point[1]
          + final_point[3] * final_point[3]
          - 2.0;

        const double c3 =
            final_point[1] * final_point[1]
          + final_point[2] * final_point[2]
          + final_point[0]
          - 4.0;

        const double constraint_residual =
            std::sqrt(
                c1 * c1
              + std::max(0.0, c2)
                    * std::max(0.0, c2)
              + std::max(0.0, c3)
                    * std::max(0.0, c3)
            );

        const bool valid_result =
            final_point.allFinite()
            && std::isfinite(final_value);

        const std::string status =
            valid_result ? "OK" : "FAILURE";

        // --------------------------------------------------------
        // Write CSV row
        // --------------------------------------------------------

        csv_file
            << test_id + 1 << ","
            << "OUR_SQP" << ","
            << problem.num_iter().size() << ","
            << elapsed_raw << ",";

        for (int j = 0; j < n; ++j) {
            if (j > 0) {
                csv_file << ";";
            }

            csv_file << x0[j];
        }

        csv_file << ",";

        for (int j = 0; j < n; ++j) {
            if (j > 0) {
                csv_file << ";";
            }

            csv_file << final_point[j];
        }

        csv_file
            << ","
            << final_value << ","
            << dist_solution << ","
            << constraint_residual << ","
            << status
            << std::endl;

        // --------------------------------------------------------
        // Print compact information
        // --------------------------------------------------------

        printf(
            "Test %d/%d\n",
            test_id + 1,
            num_simulations
        );

        printf("Initial point : ");
        for (int j = 0; j < n; ++j) {
            printf("%.6f ", x0[j]);
        }
        printf("\n");

        printf("Optimal point : ");
        for (int j = 0; j < n; ++j) {
            printf("%.10f ", final_point[j]);
        }
        printf("\n");

        printf(
            "f(x)          : %.10f\n",
            final_value
        );

        printf(
            "Distance      : %.4e\n",
            dist_solution
        );

        printf(
            "Residual      : %.4e\n",
            constraint_residual
        );

        printf(
            "SQP steps     : %zu\n",
            problem.num_iter().size()
        );

        printf(
            "Time [ns]     : %lld\n",
            elapsed_raw
        );

        printf(
            "--------------------------------------------------------\n"
        );
    }

    csv_file.close();

    printf("========================================================\n");
    printf("CSV written to SQP/sqp_test4.csv\n");
    printf("========================================================\n");

    return 0;
}