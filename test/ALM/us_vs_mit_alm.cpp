// CODE TO TEST THE PERFORMANCE OF OUR AUGMENTED LAGRANGIAN METHOD IMPLEMENTATION
// AGAINST MIT'S LBFGSB

// Make necessary includes
#include <LBFGSB.h>                                     // MIT version
#include <fdaPDE/core/fdaPDE/optimization.h>            // OUR version

#include <LBFGSpp/LineSearchBacktracking.h>
#include <fdaPDE/fdapde.h>
#include <unsupported/Eigen/SparseExtra>
#include <filesystem>
#include <fstream>
#include <chrono>

#include "alm_obj_constr.h"
// useful type aliases
using Eigen::VectorXd;
using fdapde::ScalarField;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

struct Rosenbrock {
    int n;
    double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
    {
        grad.setZero();
        double fx = 0.0;
        for (int i = 0; i < n - 1; i++) {
            double t1 = 1.0 - x[i];
            double t2 = x[i+1] - x[i]*x[i];
            grad[i]   += -2.0 * t1 - 400.0 * x[i] * t2;
            grad[i+1] += 200.0 * t2;
            fx += t1*t1 + 100.0*t2*t2;
        }
        return fx;
    }
};

int main() {

    // Define the dimension of the problem and number of simulations to perform
    const int n = 5;
    int num_simulations = 10;

    // Define box constraints
    // n-dimensional cube with edge length 2*size centered in the origin
    double size = 10;
    // MIT VERSION
    VectorXd lb = VectorXd::Constant(n, -size);
    VectorXd ub = VectorXd::Constant(n, size);
    // OUR VERSION
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

    // Define the function to minimize (n-dimensional Rosenbrock function in this case)
    // MIT FUNCTION
    Rosenbrock fun_mit{n};
    // OUR FUNCTION
    ScalarField<n> obj_fun;
    obj_fun = [n](Eigen::Matrix<double,n,1> x) -> double {
        double sum = 0.0;
        for (int i = 0; i < n - 1; ++i) {
            double t1 = x[i+1] - x[i]*x[i];
            double t2 = 1 - x[i];
            sum += 100.0 * t1 * t1 + t2 * t2;
        }
        return sum;
    };
    
    auto obj_gradient = [n](Eigen::Matrix<double,n,1> x) -> Eigen::Matrix<double,n,1> {
        Eigen::Matrix<double,n,1> grad = Eigen::Matrix<double,n,1>::Zero();
        for (int i = 0; i < n - 1; i++) {
            double t1 = 1.0 - x[i];
            double t2 = x[i+1] - x[i]*x[i];
            grad[i]   += -2.0 * t1 - 400.0 * x[i] * t2;
            grad[i+1] += 200.0 * t2;
        }
        return grad;
    };
    objFunction<n> obj(obj_fun, obj_gradient);
    // Set up output file stream
    std::ofstream file("ALM/simulations_rosenbrock.csv");
    file << "TestId,OptId,NumIter,CompTime,MinPoint,Fx,DistSol,DistBetween" << std::endl;
    
    // Perform multiple simulations
    float dist_mit, dist_our, dist_between;
    VectorXd solution = VectorXd::Constant(n, 1.0);

    for(int i=0; i < num_simulations; ++i) {
        // Random initial guess
        Eigen::VectorXd x0(n);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-2*size, 2*size);

        for(int i=0; i < x0.size(); ++i) {
            x0(i) = dis(gen);
        }

        // Define MIT minimizer
        LBFGSpp::LBFGSBParam<double> param;
        // Set of parameters to get convergence for this hard problem
        param.epsilon        = 1e-6;
        param.max_iterations = 1000;
        LBFGSpp::LBFGSBSolver<double> solver_mit(param);
        // Define OUR minimzer (using the same parameters)
        fdapde::GradientDescent<n> optimizer;
        fdapde::Lagrangian<n, fdapde::GradientDescent<n>> solver_our(optimizer, 1000, 1e-6);

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
        // MIT version
        double fx;                      // stores f(x) at the minimum
        Eigen::VectorXd x = x0;         // necessary since the solver overwrites the initial guess
        // Start timer (MIT)
        auto start_mit = std::chrono::high_resolution_clock::now();                                     
        int niter = solver_mit.minimize(fun_mit, x, fx, lb, ub);
        // Stop timer (MIT)
        auto end_mit = std::chrono::high_resolution_clock::now();
        // Convert timer to nanoseconds (MIT)
        auto elapsed_time_mit = std::chrono::duration_cast<std::chrono::nanoseconds>(end_mit - start_mit);
        long long elapsed_time_mit_raw = elapsed_time_mit.count();
        // Start timer (OUR)
        auto start_our = std::chrono::high_resolution_clock::now();
        solver_our.solve(obj, constraints, x0, fdapde::BacktrackingLineSearch());
        // Stop timer (OUR)
        auto end_our = std::chrono::high_resolution_clock::now();
        // Convert timer to nanoseconds (OUR)
        auto elapsed_time_our = std::chrono::duration_cast<std::chrono::nanoseconds>(end_our - start_our);
        long long elapsed_time_our_raw = elapsed_time_our.count();
        // Compute required distances
        dist_mit = (x - solution).norm();
        dist_our = (solver_our.optimum().back() - solution).norm();
        dist_between = (solver_our.optimum().back() - x).norm();
        // Output
        // MIT Line
        std::string separator_mit = "";
        file << i+1 << "," << "MIT" << "," << niter << "," << elapsed_time_mit_raw << ",";
        for(const auto value : x) { file << separator_mit << value; separator_mit = ";";}
        file << "," << fx << "," << dist_mit << "," << dist_between << std::endl;
        // OUR line
        std::string separator_our = "";
        file << i+1 << "," << "OUR" << "," << solver_our.num_iter().size() << "," << elapsed_time_our_raw << ",";
        for(const auto value : solver_our.optimum().back()) { file << separator_our << value; separator_our = ";";}
        file << "," << solver_our.values().back() << "," << dist_our << "," << dist_between << std::endl;
        // At the last iteration also export the trajectory of our optimizer
        if(i==num_simulations-1) {
            std::ofstream traj_("ALM/trajectory_rosenbrock.csv");
            traj_ << "Iter,Point,Fx" << std::endl;
            // Add x0
            std::string separator_traj = "";
            traj_ << 0 << ",";
            for(int j=0; j < x0.size(); ++j){ traj_ << separator_traj << x0[j]; separator_traj = ";";}
            traj_ << "," << obj(x0) << std::endl;
            // Add the rest of the trajectory
            for(int j=0; j < solver_our.num_iter().size(); ++j) {
                separator_traj = "";
                traj_ << j+1 << ",";
                for(const auto& value : solver_our.optimum()[j]) {traj_ << separator_traj << value; separator_traj = ";";}
                traj_ << "," << solver_our.values()[j] << std::endl;
            }
        }
    }
    
    return 0;
}