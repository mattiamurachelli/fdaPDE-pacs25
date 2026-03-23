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
    // 4D TEST CASE WITH CONSTRAINT INEQUALITIES
    // f(x) = (x0 - 1)^2 + (x1 - 2)^2 + (x2 + 1)^2 + (x3 - 3)^2 + sin(x0x2) + cos(x1x3)
    // c1(x) = x0^2 + x1^2 + x2 - 3 = 0
    // c2(x) = x0x1 + x3^2 - 2 <= 0
    // c3(x) = x1^2 + x2^2 + x0 - 4 <= 0

    // Create the optimizer and the problem
    fdapde::GradientDescent<4> optimizer;
    fdapde::Lagrangian<4, fdapde::GradientDescent<4>> problem(optimizer);

    // Set up objective function and constraints using the Function wrapper
    // Objective function
    ScalarField<4> objective;
    objective = [](Eigen::Matrix<double, 4, 1> x) -> double {
        return (x[0] - 1) * (x[0] - 1) + (x[1] - 2) * (x[1] - 2) + (x[2] + 1) * (x[2] + 1) + 
               (x[3] - 3) * (x[3] - 3) + std::sin(x[0] * x[2]) + std::cos(x[1] * x[3]);
    };
    objFunction<4> obj_func(objective);

    // Constraint 1: x0^2 + x1^2 + x2 - 3 = 0
    ScalarField<4> constr_func_1;
    constr_func_1 = [](Eigen::Matrix<double, 4, 1> x) -> double {
        return x[0] * x[0] + x[1] * x[1] + x[2] - 3;
    };
    struct Grad1 {
        using vector_t = Eigen::Matrix<double, 4, 1>;

        vector_t operator()(const vector_t& x) const {
            vector_t g;
            g << 2.0*x[0], 2.0*x[1], 1.0, 0.0;
            return g;
        }
    };
    constrFunction<4> constr_1(constr_func_1, Grad1{}, false);
    //constrFunction<4> constr_1(constr_func_1, false);

    // Constraint 2: x0x1 + x3^2 - 2 <= 0
    ScalarField<4> constr_func_2;
    constr_func_2 = [](Eigen::Matrix<double, 4, 1> x) -> double {
        return  x[0] * x[1] + x[3] * x[3] - 2;
    };
    struct Grad2 {
        using vector_t = Eigen::Matrix<double, 4, 1>;

        vector_t operator()(const vector_t& x) const {
            vector_t g;
            g << x[1], x[0], 0.0, 2.0*x[3];
            return g;
        }
    };
    constrFunction<4> constr_2(constr_func_2, Grad2{}, true);
    //constrFunction<4> constr_2(constr_func_2, true);

    // Constraint 3: x1^2 + x2^2 + x0 - 4 <= 0
    ScalarField<4> constr_func_3;
    constr_func_3 = [](Eigen::Matrix<double, 4, 1> x) -> double {
        return x[1] * x[1] + x[2] * x[2] + x[0] - 4;
    };
    struct Grad3 {
        using vector_t = Eigen::Matrix<double, 4, 1>;

        vector_t operator()(const vector_t& x) const {
            vector_t g;
            g << 1.0, 2.0*x[1], 2.0*x[2], 0.0;
            return g;
        }
    };
    constrFunction<4> constr_3(constr_func_3, Grad3{}, true);
    //constrFunction<4> constr_3(constr_func_3, true);

    // Create a constrList object to hold the constraints
    constrList<4> constraints = {constr_1, constr_2, constr_3};

    // Set up an initial point
    Eigen::Matrix<double, 4, 1> x0;
    x0 << 0.0, 0.0, 0.0, 0.0;

    // Solve the problem
    problem.solve(std::move(obj_func), constraints, x0, fdapde::BacktrackingLineSearch());

    // Print results
    printf("========================================================\n");
    printf("4D PROBLEM  : \n");
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
    for (double val : values) { printf("%.6f ", val); }
    printf("\n");
    // Optimal points
    printf("Optimal point : \n");
    const auto& opt_points = problem.optimum();
    for (const auto& point : opt_points) {
        printf("(");
        for (std::size_t i = 0; i < point.size(); ++i) { printf("%.6f, ", point[i]); }
        printf(") \t");
    }
    printf("\n");
    printf("========================================================\n");

    return 0;
}