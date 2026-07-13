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

int main() {
    // ============================================================
    // SIMPLE 2D TEST CASE
    //
    // f(x,y) = (x - 1)^2 + (y - 2)^2
    // c(x,y) = x + y - 1 = 0
    //
    // Exact constrained solution: x* = (0,1)
    // ============================================================

    constexpr int n = 2;

    // ============================================================
    // SIMULATION PARAMETERS
    // ============================================================

    const int num_simulations = 50;

    // Controls the size of the neighborhood around the optimum.
    // Starting points are generated inside the square:
    // [x*_0 - radius, x*_0 + radius] x [x*_1 - radius, x*_1 + radius]
    const double neighborhood_radius = 1.0;

    // Fixed seed for reproducibility.
    // Use the same seed in the NLopt file to generate the same starting points.
    const unsigned int seed = 12345;

    // Exact solution
    Eigen::Matrix<double, 2, 1> exact_solution;
    exact_solution << 0.0, 1.0;

    // Random generator for starting points
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(
        -neighborhood_radius,
         neighborhood_radius
    );

    // ============================================================
    // OBJECTIVE FUNCTION
    // ============================================================

    ScalarField<2> objective;
    objective = [](Eigen::Matrix<double, 2, 1> x) -> double {
        return (x[0] - 1) * (x[0] - 1)
             + (x[1] - 2) * (x[1] - 2);
    };

    // ============================================================
    // CONSTRAINT
    // ============================================================

    // Constraint: x + y - 1 = 0
    ScalarField<2> constraint;
    constraint = [](Eigen::Matrix<double, 2, 1> x) -> double {
        return x[0] + x[1] - 1;
    };

    struct Grad {
        using vector_t = Eigen::Matrix<double, 2, 1>;

        vector_t operator()(const vector_t& x) const {
            vector_t g;
            g << 1.0, 1.0;
            return g;
        }
    };

    constrFunction<2> constr_func(constraint, Grad{}, false);

    // Create a constrList object to hold the constraints
    constrList<2> constraints = {constr_func};

    // ============================================================
    // CSV OUTPUT
    // ============================================================

    std::ofstream csv_file("ALM/lagrangian_test1.csv");

    csv_file << "TestId,OptId,NumIterOrEval,CompTime,InitPoint,MinPoint,Fx,DistSol,Residual,Status"
             << std::endl;

    // ============================================================
    // RUN SIMULATIONS
    // ============================================================

    printf("========================================================\n");
    printf("OUR ALM - N SIMULATIONS\n");
    printf("Number of simulations      : %d\n", num_simulations);
    printf("Neighborhood radius        : %.4f\n", neighborhood_radius);
    printf("Seed                       : %u\n", seed);
    printf("Exact solution             : (%.4f, %.4f)\n",
           exact_solution[0], exact_solution[1]);
    printf("========================================================\n");

    for (int test_id = 0; test_id < num_simulations; ++test_id) {

        // ------------------------------------------------------------
        // Generate initial point around the optimum
        // ------------------------------------------------------------

        Eigen::Matrix<double, 2, 1> x0;

        x0[0] = exact_solution[0] + dis(gen);
        x0[1] = exact_solution[1] + dis(gen);

        // ------------------------------------------------------------
        // Create optimizer and ALM solver
        // ------------------------------------------------------------

        fdapde::GradientDescent<2> optimizer;
        fdapde::Lagrangian<2, fdapde::GradientDescent<2>> problem(optimizer);

        // Important: create the objective wrapper inside the loop.
        // This avoids problems if the object is moved or modified during solve.
        objFunction<2> obj_func(objective);

        // ------------------------------------------------------------
        // Solve and measure computation time
        // ------------------------------------------------------------

        auto start = std::chrono::steady_clock::now();

        problem.solve(
            obj_func,
            constraints,
            x0,
            fdapde::BacktrackingLineSearch()
        );

        auto end = std::chrono::steady_clock::now();

        auto elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                end - start
            );

        long long elapsed_raw = elapsed.count();

        // std::chrono::steady_clock is monotonic and suitable for measuring intervals.
        // See cppreference documentation. 

        // ------------------------------------------------------------
        // Extract results
        // ------------------------------------------------------------

        const auto& final_point = problem.optimum().back();

        double final_value = problem.values().back();

        double dist_solution =
            (final_point - exact_solution).norm();

        double constraint_residual =
            std::abs(final_point[0] + final_point[1] - 1.0);

        std::string status = "OK";

        // ------------------------------------------------------------
        // Write CSV line
        // ------------------------------------------------------------

        csv_file << test_id + 1 << ","
                 << "OUR_ALM" << ","
                 << problem.num_iter().size() << ","
                 << elapsed_raw << ","
                 << x0[0] << ";" << x0[1] << ","
                 << final_point[0] << ";" << final_point[1] << ","
                 << final_value << ","
                 << dist_solution << ","
                 << constraint_residual << ","
                 << status
                 << std::endl;

        // ------------------------------------------------------------
        // Print compact information to screen
        // ------------------------------------------------------------

        printf("Test %d/%d\n", test_id + 1, num_simulations);
        printf("Initial point : (%.6f, %.6f)\n", x0[0], x0[1]);
        printf("Optimal point : (%.10f, %.10f)\n",
               final_point[0], final_point[1]);
        printf("f(x)          : %.10f\n", final_value);
        printf("Distance      : %.4e\n", dist_solution);
        printf("Residual      : %.4e\n", constraint_residual);
        printf("Subproblems   : %zu\n", problem.num_iter().size());
        printf("Time [ns]     : %lld\n", elapsed_raw);
        printf("--------------------------------------------------------\n");
    }

    csv_file.close();

    printf("========================================================\n");
    printf("CSV written to ALM/lagrangian_test1.csv\n");
    printf("========================================================\n");

    return 0;
}