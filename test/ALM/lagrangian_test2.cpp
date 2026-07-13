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
    // 5D TEST CASE
    //
    // f(x) = (x0 - 1)^2 + (x1 - 2)^2 + (x2 - 3)^2
    //      + (x3 - 2)^2 + (x4 - 1)^2
    //
    // c1(x) = x0 + x1 + x2 - 3 = 0
    // c2(x) = x0^2 + x3 - 1 = 0
    // c3(x) = x1*x4 - 2 = 0
    // ============================================================

    constexpr int n = 5;

    // ============================================================
    // SIMULATION PARAMETERS
    // ============================================================

    const int num_simulations = 50;

    // Controls the size of the neighborhood around the reference solution.
    // Starting points are generated inside:
    // x_ref[j] +/- neighborhood_radius
    const double neighborhood_radius = 1.0;

    // Fixed seed for reproducibility.
    // Use the same value in the NLopt file.
    const unsigned int seed = 12345;

    // Reference solution computed for this constrained problem.
    Eigen::Matrix<double, 5, 1> reference_solution;
    reference_solution << -0.08170078635350764,
                           1.3278938903564708,
                           1.7538068959970370,
                           0.9933249815092187,
                           1.5061444400976221;

    // Random generator for starting points
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(
        -neighborhood_radius,
         neighborhood_radius
    );

    // ============================================================
    // OBJECTIVE FUNCTION
    // ============================================================

    ScalarField<5> objective;
    objective = [](Eigen::Matrix<double, 5, 1> x) -> double {
        return (x[0] - 1)*(x[0] - 1)
             + (x[1] - 2)*(x[1] - 2)
             + (x[2] - 3)*(x[2] - 3)
             + (x[3] - 2)*(x[3] - 2)
             + (x[4] - 1)*(x[4] - 1);
    };

    // ============================================================
    // CONSTRAINT 1: x0 + x1 + x2 - 3 = 0
    // ============================================================

    ScalarField<5> constraint1;
    constraint1 = [](Eigen::Matrix<double, 5, 1> x) -> double {
        return x[0] + x[1] + x[2] - 3;
    };

    struct Grad1 {
        using vector_t = Eigen::Matrix<double, 5, 1>;

        vector_t operator()(const vector_t& x) const {
            vector_t g;
            g << 1.0, 1.0, 1.0, 0.0, 0.0;
            return g;
        }
    };

    constrFunction<5> constr_1(constraint1, Grad1{}, false);

    // ============================================================
    // CONSTRAINT 2: x0^2 + x3 - 1 = 0
    // ============================================================

    ScalarField<5> constraint2;
    constraint2 = [](Eigen::Matrix<double, 5, 1> x) -> double {
        return x[0]*x[0] + x[3] - 1;
    };

    struct Grad2 {
        using vector_t = Eigen::Matrix<double, 5, 1>;

        vector_t operator()(const vector_t& x) const {
            vector_t g;
            g << 2.0*x[0], 0.0, 0.0, 1.0, 0.0;
            return g;
        }
    };

    constrFunction<5> constr_2(constraint2, Grad2{}, false);

    // ============================================================
    // CONSTRAINT 3: x1*x4 - 2 = 0
    // ============================================================

    ScalarField<5> constraint3;
    constraint3 = [](Eigen::Matrix<double, 5, 1> x) -> double {
        return x[1]*x[4] - 2;
    };

    struct Grad3 {
        using vector_t = Eigen::Matrix<double, 5, 1>;

        vector_t operator()(const vector_t& x) const {
            vector_t g;
            g << 0.0, x[4], 0.0, 0.0, x[1];
            return g;
        }
    };

    constrFunction<5> constr_3(constraint3, Grad3{}, false);

    // Create a constrList object to hold the constraints
    constrList<5> constraints = {constr_1, constr_2, constr_3};

    // ============================================================
    // CSV OUTPUT
    // ============================================================

    std::ofstream csv_file("ALM/lagrangian_test2.csv");

    csv_file << "TestId,OptId,NumIterOrEval,CompTime,InitPoint,MinPoint,Fx,DistSol,Residual,Status"
             << std::endl;

    csv_file << std::setprecision(17);

    // ============================================================
    // RUN SIMULATIONS
    // ============================================================

    printf("========================================================\n");
    printf("OUR ALM - 5D TEST CASE - N SIMULATIONS\n");
    printf("Number of simulations      : %d\n", num_simulations);
    printf("Neighborhood radius        : %.4f\n", neighborhood_radius);
    printf("Seed                       : %u\n", seed);
    printf("Reference solution         : ");
    for (int j = 0; j < n; ++j) {
        printf("%.6f ", reference_solution[j]);
    }
    printf("\n");
    printf("========================================================\n");

    for (int test_id = 0; test_id < num_simulations; ++test_id) {

        // ------------------------------------------------------------
        // Generate initial point around the reference solution
        // ------------------------------------------------------------

        Eigen::Matrix<double, 5, 1> x0;

        for (int j = 0; j < n; ++j) {
            x0[j] = reference_solution[j] + dis(gen);
        }

        // ------------------------------------------------------------
        // Create optimizer and ALM solver
        // ------------------------------------------------------------

        fdapde::GradientDescent<5> optimizer;
        fdapde::Lagrangian<5, fdapde::GradientDescent<5>> problem(optimizer);

        objFunction<5> obj_func(objective);

        // ------------------------------------------------------------
        // Solve and measure computation time
        // ------------------------------------------------------------

        auto start = std::chrono::steady_clock::now();

        problem.solve(
            obj_func,
            constraints,
            x0
        );

        auto end = std::chrono::steady_clock::now();

        auto elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                end - start
            );

        long long elapsed_raw = elapsed.count();

        // ------------------------------------------------------------
        // Extract results
        // ------------------------------------------------------------

        const auto& final_point = problem.optimum().back();

        double final_value = problem.values().back();

        double dist_solution =
            (final_point - reference_solution).norm();

        double c1 = final_point[0] + final_point[1] + final_point[2] - 3.0;
        double c2 = final_point[0]*final_point[0] + final_point[3] - 1.0;
        double c3 = final_point[1]*final_point[4] - 2.0;

        double constraint_residual =
            std::sqrt(c1*c1 + c2*c2 + c3*c3);

        std::string status = "OK";

        // ------------------------------------------------------------
        // Write CSV line
        // ------------------------------------------------------------

        csv_file << test_id + 1 << ","
                 << "OUR_ALM" << ","
                 << problem.num_iter().size() << ","
                 << elapsed_raw << ",";

        for (int j = 0; j < n; ++j) {
            if (j > 0) csv_file << ";";
            csv_file << x0[j];
        }

        csv_file << ",";

        for (int j = 0; j < n; ++j) {
            if (j > 0) csv_file << ";";
            csv_file << final_point[j];
        }

        csv_file << ","
                 << final_value << ","
                 << dist_solution << ","
                 << constraint_residual << ","
                 << status
                 << std::endl;

        // ------------------------------------------------------------
        // Print compact information to screen
        // ------------------------------------------------------------

        printf("Test %d/%d\n", test_id + 1, num_simulations);

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

        printf("f(x)          : %.10f\n", final_value);
        printf("Distance      : %.4e\n", dist_solution);
        printf("Residual      : %.4e\n", constraint_residual);
        printf("Subproblems   : %zu\n", problem.num_iter().size());
        printf("Time [ns]     : %lld\n", elapsed_raw);
        printf("--------------------------------------------------------\n");
    }

    csv_file.close();

    printf("========================================================\n");
    printf("CSV written to ALM/lagrangian_test2.csv\n");
    printf("========================================================\n");

    return 0;
}