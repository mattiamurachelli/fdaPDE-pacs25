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
#include "alm_obj_constr.h"
#include <cstdio>
#include <vector>

int main() {
    // SIMPLE 2D TEST CASE
    // f(x) = (x - 1)^2 + (y - 2)^2
    // c(x) = x + y - 1 = 0

    // Create the optimizer and the problem
    fdapde::GradientDescent<2> optimizer;
    fdapde::Lagrangian<2, fdapde::GradientDescent<2>> problem(optimizer);

    // Set up objective function and constraints using the Function wrapper
    // Objective function
    ScalarField<2> objective;
    objective = [](Eigen::Matrix<double, 2, 1> x) -> double {
        return (x[0] - 1) * (x[0] - 1) + (x[1] - 2) * (x[1] - 2);
    };
    objFunction<2> obj_func(objective);

    // Constraint: x + y - 1 = 0
    ScalarField<2> constraint;
    constraint = [](Eigen::Matrix<double, 2, 1> x) -> double {
        return x[0] + x[1]  - 1;
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
    //constrFunction<2> constr_func(constraint, false);

    // Create a constrList object to hold the constraints
    constrList<2> constraints = {constr_func};

    // Set up an initial point
    Eigen::Matrix<double, 2, 1> x0;
    x0 << -5.0, 1.0;

    // Solve the problem
    problem.solve(std::move(obj_func), constraints, x0, fdapde::BacktrackingLineSearch());

    // Print results
    printf("========================================================\n");
    printf("2D PROBLEM  : \n");
    // Initial point
    printf("Initial point : ");
    for(std::size_t i = 0; i < x0.size(); ++i) { printf("%.2f ", x0[i]); }
    printf("\n");
    // Number of subproblems and corresponding iterations
    printf("Number of subproblems solved : %d \n", problem.num_iter().size());
    printf("Number of iterations (for each subproblem): \n");
    std::vector<int> num_iter_ = problem.num_iter();
    for (std::size_t iter : num_iter_) { printf("%d ", iter); }
    printf("\n");
    // Values of f(x) at the optimal points
    printf("Values f(x) at optimal points : \n");
    std::vector<double> values = problem.values();
    for (double val : values) { printf("%.4f ", val); }
    printf("\n");
    // Optimal points
    printf("Optimal point : \n");
    const auto& opt_points = problem.optimum();
    for (const auto& point : opt_points) {
        printf("(");
        for (std::size_t i = 0; i < point.size(); ++i) { printf("%.4f, ", point[i]); }
        printf(") \t");
    }
    printf("\n");
    printf("========================================================\n");

    return 0;
}