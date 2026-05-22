// CODE TO TEST THE PERFORMANCE OF OUR OPTIMIZERS IMPLEMENTATIONS
// AGAINST MIT'S LBFGSB IMPLEMENTATION
//
// METHODS TESTED:
// 1) MIT LBFGSB
// 2) OUR LBFGSB
// 3) OUR ALM
// 4) OUR SQP

// Make necessary includes
#include <LBFGSB.h>                                     // MIT version
#include <fdaPDE/core/fdaPDE/optimization.h>            // OUR version
#include <LBFGSpp/LineSearchBacktracking.h>
#include <fdaPDE/fdapde.h>
#include <unsupported/Eigen/SparseExtra>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <random>
#include <iostream>
#include "../INCLUDE/obj_constr.h"

// useful type aliases
using Eigen::VectorXd;
using fdapde::ScalarField;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

// Rosenbrock function for MIT optimizer
struct Rosenbrock {
    int n;

    double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
    {
        grad.setZero();

        double fx = 0.0;

        for (int i = 0; i < n - 1; i++) {
            double t1 = 1.0 - x[i];
            double t2 = x[i + 1] - x[i] * x[i];

            grad[i]   += -2.0 * t1 - 400.0 * x[i] * t2;
            grad[i+1] += 200.0 * t2;

            fx += t1 * t1 + 100.0 * t2 * t2;
        }

        return fx;
    }
};

int main() {

    // Define the dimension of the problem and number of simulations
    const int n = 5;
    int num_simulations = 10;

    // Define box constraints
    // n-dimensional cube with edge length 2*size centered in the origin
    double size = 10.0;

    // MIT VERSION BOUNDS
    VectorXd lb = VectorXd::Constant(n, -size);
    VectorXd ub = VectorXd::Constant(n,  size);

    // OUR VERSION CONSTRAINTS
    std::vector<constrFunction<n>> constraints_vec;

    for (int i = 0; i < n; ++i) {

        // Constraint: x_i <= size <=> x_i - size <= 0
        ScalarField<n> c_upper;

        c_upper = [i, size](Eigen::Matrix<double, n, 1> x) -> double {
            return x[i] - size;
        };

        constrFunction<n> f_upper(c_upper, true);
        constraints_vec.push_back(f_upper);

        // Constraint: x_i >= -size <=> -x_i - size <= 0
        ScalarField<n> c_lower;

        c_lower = [i, size](Eigen::Matrix<double, n, 1> x) -> double {
            return -x[i] - size;
        };

        constrFunction<n> f_lower(c_lower, true);
        constraints_vec.push_back(f_lower);
    }

    constrList<n> constraints(constraints_vec);

    // Define the function to minimize (n-dimensional Rosenbrock function)

    // MIT FUNCTION
    Rosenbrock fun_mit{n};

    // OUR FUNCTION
    ScalarField<n> obj_fun;

    obj_fun = [n](Eigen::Matrix<double, n, 1> x) -> double {

        double sum = 0.0;

        for (int i = 0; i < n - 1; ++i) {

            double t1 = x[i + 1] - x[i] * x[i];
            double t2 = 1 - x[i];

            sum += 100.0 * t1 * t1 + t2 * t2;
        }

        return sum;
    };

    auto obj_gradient = [n](Eigen::Matrix<double, n, 1> x)
        -> Eigen::Matrix<double, n, 1> {

        Eigen::Matrix<double, n, 1> grad =
            Eigen::Matrix<double, n, 1>::Zero();

        for (int i = 0; i < n - 1; i++) {

            double t1 = 1.0 - x[i];
            double t2 = x[i + 1] - x[i] * x[i];

            grad[i]   += -2.0 * t1 - 400.0 * x[i] * t2;
            grad[i+1] += 200.0 * t2;
        }

        return grad;
    };

    objFunction<n> obj(obj_fun, obj_gradient);

    // Define MIT minimizer
    LBFGSpp::LBFGSBParam<double> param;

    // Set of parameters to get convergence for this hard problem
    param.epsilon        = 1e-6;
    param.max_iterations = 1000;
    param.ftol           = 1e-4;
    param.wolfe          = 0.9;
    param.min_step       = 1e-20;
    param.max_step       = 1.0;

    LBFGSpp::LBFGSBSolver<double> solver_mit(param);

    // Define OUR LBFGSB minimizer
    fdapde::LBFGSBParams params;

    params.epsilon_  = 1e-6;
    params.max_iter_ = 1000;
    params.c1_       = 1e-4;
    params.c2_       = 0.9;
    params.alpha_min_= 1e-20;
    params.step_     = 1.0;

    fdapde::LBFGSB<n> solver_lbfgsb(params);

    // Define OUR ALM minimizer
    fdapde::GradientDescent<n> optimizer;
    fdapde::Lagrangian<n, fdapde::GradientDescent<n>> solver_alm(optimizer, 1000, 1e-6);

    // Define OUR SQP minimizer
    fdapde::SQP<n> solver_sqp;

    // Set up output file stream
    std::ofstream file("TEST_1/simulations_rosenbrock.csv");

    file << "TestId,OptId,NumIter,CompTime,MinPoint,Fx,DistSol,DistBetween" << std::endl;

    // Exact solution
    VectorXd solution = VectorXd::Constant(n, 1.0);

    // Perform multiple simulations
    for (int test_id = 0; test_id < num_simulations; ++test_id) {

        // Random initial guess
        Eigen::VectorXd x0(n);

        std::random_device rd;
        std::mt19937 gen(rd());

        std::uniform_real_distribution<> dis(-2 * size, 2 * size);

        for (int i = 0; i < x0.size(); ++i) {
            x0(i) = dis(gen);
        }

        // Output current test
        std::cout << "____________________________________________________________" << std::endl;

        std::cout << "TEST NUMBER : " << test_id + 1 << std::endl;

        std::cout << "Initial guess = (";

        for (int i = 0; i < x0.size(); ++i) {

            if (i == x0.size() - 1) {
                std::cout << x0(i);
            }
            else {
                std::cout << x0(i) << ", ";
            }
        }

        std::cout << ")" << std::endl;

        std::cout << "____________________________________________________________" << std::endl;

        // ============================================================
        // MIT LBFGSB
        // ============================================================

        double fx_mit;

        Eigen::VectorXd x_mit = x0;

        // Start timer
        auto start_mit = std::chrono::steady_clock::now();

        int niter_mit = solver_mit.minimize(fun_mit, x_mit, fx_mit, lb, ub);

        // Stop timer
        auto end_mit = std::chrono::steady_clock::now();

        // Convert timer to nanoseconds
        auto elapsed_time_mit = std::chrono::duration_cast<std::chrono::nanoseconds>(end_mit - start_mit);

        long long elapsed_time_mit_raw = elapsed_time_mit.count();

        // ============================================================
        // OUR LBFGSB
        // ============================================================

        auto start_lbfgsb = std::chrono::steady_clock::now();

        solver_lbfgsb.optimize(obj_fun, obj_gradient, x0, lb, ub, fdapde::MoreThuenteLineSearch());

        auto end_lbfgsb = std::chrono::steady_clock::now();

        auto elapsed_time_lbfgsb = std::chrono::duration_cast<std::chrono::nanoseconds>(end_lbfgsb - start_lbfgsb);

        long long elapsed_time_lbfgsb_raw = elapsed_time_lbfgsb.count();

        // ============================================================
        // OUR ALM
        // ============================================================

        auto start_alm = std::chrono::steady_clock::now();

        solver_alm.solve(obj, constraints, x0, fdapde::BacktrackingLineSearch());

        auto end_alm = std::chrono::steady_clock::now();

        auto elapsed_time_alm = std::chrono::duration_cast<std::chrono::nanoseconds>(end_alm - start_alm);

        long long elapsed_time_alm_raw = elapsed_time_alm.count();

        // ============================================================
        // OUR SQP
        // ============================================================

        auto start_sqp = std::chrono::steady_clock::now();

        solver_sqp.solve(obj, constraints, x0);

        auto end_sqp = std::chrono::steady_clock::now();

        auto elapsed_time_sqp = std::chrono::duration_cast<std::chrono::nanoseconds>(end_sqp - start_sqp);

        long long elapsed_time_sqp_raw = elapsed_time_sqp.count();

        // ============================================================
        // Compute distances
        // ============================================================

        // MIT
        double dist_mit = (x_mit - solution).norm();

        // OUR LBFGSB
        double dist_lbfgsb = (solver_lbfgsb.optimum() - solution).norm();

        double dist_between_lbfgsb = (solver_lbfgsb.optimum() - x_mit).norm();

        // OUR ALM
        double dist_alm =(solver_alm.optimum().back() - solution).norm();

        double dist_between_alm = (solver_alm.optimum().back() - x_mit).norm();

        // OUR SQP
        double dist_sqp = (solver_sqp.optimum().back() - solution).norm();

        double dist_between_sqp = (solver_sqp.optimum().back() - x_mit).norm();

        // ============================================================
        // OUTPUT CSV
        // ============================================================

        // ------------------------------------------------------------
        // MIT LINE
        // ------------------------------------------------------------

        std::string separator_mit = "";

        file << test_id + 1 << "," << "MIT-LBFGSB" << "," << niter_mit << "," << elapsed_time_mit_raw << ",";

        for (const auto value : x_mit) {
            file << separator_mit << value;
            separator_mit = ";";
        }

        file << "," << fx_mit << "," << dist_mit << "," << 0.0 << std::endl;

        // ------------------------------------------------------------
        // OUR LBFGSB LINE
        // ------------------------------------------------------------

        std::string separator_lbfgsb = "";

        file << test_id + 1 << "," << "OUR-LBFGSB" << "," << solver_lbfgsb.n_iter() << "," << elapsed_time_lbfgsb_raw << ",";

        for (const auto value : solver_lbfgsb.optimum()) {
            file << separator_lbfgsb << value;
            separator_lbfgsb = ";";
        }

        file << "," << solver_lbfgsb.value() << "," << dist_lbfgsb << "," << dist_between_lbfgsb << std::endl;

        // ------------------------------------------------------------
        // OUR ALM LINE
        // ------------------------------------------------------------

        std::string separator_alm = "";

        file << test_id + 1 << "," << "OUR-ALM" << "," << solver_alm.num_iter().size() << "," << elapsed_time_alm_raw << ",";

        for (const auto value : solver_alm.optimum().back()) {
            file << separator_alm << value;
            separator_alm = ";";
        }

        file << "," << solver_alm.values().back() << "," << dist_alm << "," << dist_between_alm << std::endl;
 
        // ------------------------------------------------------------
        // OUR SQP LINE
        // ------------------------------------------------------------

        std::string separator_sqp = "";

        file << test_id + 1 << "," << "OUR-SQP" << "," << solver_sqp.num_iter().size() << "," << elapsed_time_sqp_raw << ",";

        for (const auto value : solver_sqp.optimum().back()) {
            file << separator_sqp << value;
            separator_sqp = ";";
        }

        file << "," << solver_sqp.values().back() << "," << dist_sqp << "," << dist_between_sqp << std::endl;

        // ============================================================
        // EXPORT TRAJECTORIES AT LAST TEST
        // ============================================================

        if (test_id == num_simulations - 1) {

            // --------------------------------------------------------
            // LBFGSB trajectory
            // --------------------------------------------------------

            std::ofstream traj_lbfgsb("TEST_1/trajectory_lbfgsb_rosenbrock.csv");

            traj_lbfgsb << "Iter,Point,Fx" << std::endl;

            for (int j = 0; j < solver_lbfgsb.n_iter(); ++j) {

                std::string separator_traj = "";

                traj_lbfgsb << j + 1 << ",";

                for (const auto& value :
                     solver_lbfgsb.trajectory()[j]) {

                    traj_lbfgsb << separator_traj << value;
                    separator_traj = ";";
                }

                traj_lbfgsb << "," << solver_lbfgsb.values()[j] << std::endl;
            }

            // --------------------------------------------------------
            // ALM trajectory
            // --------------------------------------------------------

            std::ofstream traj_alm("TEST_1/trajectory_alm_rosenbrock.csv");

            traj_alm << "Iter,Point,Fx" << std::endl;

            // Add x0
            std::string separator_traj_alm = "";

            traj_alm << 0 << ",";

            for (int j = 0; j < x0.size(); ++j) {

                traj_alm << separator_traj_alm << x0[j];
                separator_traj_alm = ";";
            }

            traj_alm << "," << obj(x0) << std::endl;

            // Add trajectory
            for (int j = 0; j < solver_alm.num_iter().size(); ++j) {

                separator_traj_alm = "";

                traj_alm << j + 1 << ",";

                for (const auto& value :
                     solver_alm.optimum()[j]) {

                    traj_alm << separator_traj_alm << value;
                    separator_traj_alm = ";";
                }

                traj_alm << "," << solver_alm.values()[j] << std::endl;
            }

            // --------------------------------------------------------
            // SQP trajectory
            // --------------------------------------------------------

            std::ofstream traj_sqp("TEST_1/trajectory_sqp_rosenbrock.csv");

            traj_sqp << "Iter,Point,Fx" << std::endl;

            // Add x0
            std::string separator_traj_sqp = "";

            traj_sqp << 0 << ",";

            for (int j = 0; j < x0.size(); ++j) {

                traj_sqp << separator_traj_sqp << x0[j];
                separator_traj_sqp = ";";
            }

            traj_sqp << "," << obj(x0) << std::endl;

            // Add trajectory
            for (int j = 0; j < solver_sqp.num_iter().size(); ++j) {

                separator_traj_sqp = "";

                traj_sqp << j + 1 << ",";

                for (const auto& value :
                     solver_sqp.optimum()[j]) {

                    traj_sqp << separator_traj_sqp << value;
                    separator_traj_sqp = ";";
                }

                traj_sqp << "," << solver_sqp.values()[j] << std::endl;
            }
        }
    }

    return 0;
}