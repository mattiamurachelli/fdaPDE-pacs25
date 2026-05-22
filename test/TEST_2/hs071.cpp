// CODE TO TEST THE PERFORMANCE OF ALM AND SQP IMPLEMENTATION
// ON THE BENCHMARK TEST HS071

// Make necessary includes
#include <fdaPDE/core/fdaPDE/optimization.h>
#include <fdaPDE/fdapde.h>
#include <unsupported/Eigen/SparseExtra>
#include <filesystem>
#include <fstream>
#include <chrono>
#include "../INCLUDE/obj_constr.h"

// useful type aliases
using Eigen::VectorXd;
using fdapde::ScalarField;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

int main() {

    // Define the dimension of the problem and number of simulations to perform
    const int n = 4;
    int num_simulations = 10;

    // Define BOX constraints
    // n-dimensional cube with edge length 2*size centered in (3,3,...)
    double size = 2;
    std::vector<constrFunction<n>> constraints_vec;
    for (int i = 0; i < n; ++i) {

        // Constraint: x_i <= size + 3 <=> x_i - size - 3 <= 0 (x_i <= 5)
        ScalarField<n> c_upper;
        c_upper = [i, size](Eigen::Matrix<double, n, 1> x) -> double {
            return x[i] - size - 3;
        };

        constrFunction<n> f_upper(c_upper, true);
        constraints_vec.push_back(f_upper);

        // Constraint: x_i >= -size + 3 <=> -x_i - size + 3 <= 0 (x_i >= 1)
        ScalarField<n> c_lower;
        c_lower = [i, size](Eigen::Matrix<double, n, 1> x) -> double {
            return -x[i] - size + 3;
        };

        constrFunction<n> f_lower(c_lower, true);
        constraints_vec.push_back(f_lower);
    }
    // Define NON-BOX constraints
    // x0*x1*x2*x3 >= 25 <=>  - x0*x1*x2*x3 + 25 <= 0
    ScalarField<n> c_1;
    c_1 = [](Eigen::Matrix<double, n, 1> x) -> double {
        return -x[0]*x[1]*x[2]*x[3] + 25;
    };

    constrFunction<n> c_1_func(c_1, true);
    constraints_vec.push_back(c_1_func);
    // x0^2 * x1^2 * x2^2 * x3^2 = 40
    ScalarField<n> c_2;
    c_2 = [](Eigen::Matrix<double, n, 1> x) -> double {
        return x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 40;
    };

    constrFunction<n> c_2_func(c_2, false);
    constraints_vec.push_back(c_2_func);

    constrList<n> constraints(constraints_vec);

    // Objective function
    ScalarField<n> objective;
    objective = [](Eigen::Matrix<double, n, 1> x) -> double {
        return x[0]*x[3]*(x[0] + x[1] + x[2]) + x[2];
    };
    objFunction<n> obj_func(objective);
    // Set up output file stream
    std::ofstream file("TEST_2/simulations_hs071.csv");
    file << "TestId,OptId,NumIter,CompTime,MinPoint,Fx,DistSol,DistBetween" << std::endl;
    
    // Create minimizers
    // ALM
    fdapde::GradientDescent<n> optimizer;
    fdapde::Lagrangian<n, fdapde::GradientDescent<n>> solver_alm(optimizer, 1000, 1e-6);
    // SQP
    fdapde::SQP<n> solver_sqp;

    // Perform multiple simulations
    float dist_alm, dist_sqp, dist_between;
    Eigen::VectorXd solution(n);
    solution << 1.00000000, 4.74299963, 3.82114998, 1.37940829;

    for(int i=0; i < num_simulations; ++i) {
        // Random initial guess
        Eigen::VectorXd x0(n);
        x0 << 1.0, 5.0, 5.0, 1.0;

        // ALM
    fdapde::GradientDescent<n> optimizer;
    fdapde::Lagrangian<n, fdapde::GradientDescent<n>> solver_alm(optimizer, 1000, 1e-6);

        std::random_device rd;
        std::mt19937 gen(rd());
        // Generate in noise in (-1, 1)
        std::uniform_real_distribution<> dis(-1, 1);

        for(int i=0; i < x0.size(); ++i) {
            x0(i) = x0(i) + dis(gen);
        }

        // Now perform the minimization
        std::cout << "________________________________________________________________" << std::endl;
        std::cout << "TEST NUMBER : " << i+1 << std::endl;
        std::cout << "Initial guess = (";
        for(int i=0; i < x0.size(); ++i) {
            if(i == x0.size()-1) {
                std::cout << x0(i);
            } else {
                std::cout << x0(i) << ", ";
            }
        }
        std::cout << ")" << std::endl;
        std::cout << "________________________________________________________________" << std::endl;
        // Start timer (ALM)
        auto start_alm = std::chrono::steady_clock::now();
        solver_alm.solve(obj_func, constraints, x0, fdapde::BacktrackingLineSearch());
        // Stop timer (ALM)
        auto end_alm = std::chrono::steady_clock::now();
        // Convert timer to nanoseconds (ALM)
        auto elapsed_time_alm = std::chrono::duration_cast<std::chrono::nanoseconds>(end_alm - start_alm);
        long long elapsed_time_alm_raw = elapsed_time_alm.count();
        // Start timer (SQP)
        auto start_sqp = std::chrono::steady_clock::now();
        solver_sqp.solve(obj_func, constraints, x0);
        // Stop timer (SQP)
        auto end_sqp = std::chrono::steady_clock::now();
        // Convert timer to nanoseconds (SQP)
        auto elapsed_time_sqp = std::chrono::duration_cast<std::chrono::nanoseconds>(end_sqp - start_sqp);
        long long elapsed_time_sqp_raw = elapsed_time_sqp.count();
        // Compute required distances
        dist_alm = (solver_alm.optimum().back() - solution).norm();
        dist_sqp = (solver_sqp.optimum().back() - solution).norm();
        dist_between = (solver_alm.optimum().back() - solver_sqp.optimum().back()).norm();
        // Output
        // ALM Line
        std::string separator_alm = "";
        file << i+1 << "," << "ALM" << "," << solver_alm.num_iter().size() << "," << elapsed_time_alm_raw << ",";
        for(const auto value : solver_alm.optimum().back()) { file << separator_alm << value; separator_alm = ";";}
        file << "," << solver_alm.values().back() << "," << dist_alm << "," << dist_between << std::endl;
        // OUR line
        std::string separator_sqp = "";
        file << i+1 << "," << "SQP" << "," << solver_sqp.num_iter().size() << "," << elapsed_time_sqp_raw << ",";
        for(const auto value : solver_sqp.optimum().back()) { file << separator_sqp << value; separator_sqp = ";";}
        file << "," << solver_sqp.values().back() << "," << dist_sqp << "," << dist_between << std::endl;
    }
    
    return 0;
}