// CODE TO TEST THE PERFORMANCE OF OUR LBFGSB OPTIMIZER IMPLEMENTATION
// AGAINST MIT'S ONE (2D CASE)

// Make necessary includes
#include <LBFGSB.h>                                     // MIT version
#include <fdaPDE/core/fdaPDE/optimization.h>            // OUR version

#include <LBFGSpp/LineSearchBacktracking.h>
#include <fdaPDE/fdapde.h>
#include <unsupported/Eigen/SparseExtra>
#include <filesystem>
#include <fstream>
#include <chrono>

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

    // Define the dimension of the problem
    const int n = 10;

    // Define MIT minimizer
    LBFGSpp::LBFGSBParam<double> param;
    // Set of parameters to get convergence for this hard problem
    param.epsilon        = 1e-6;
    param.max_iterations = 10000;
    param.max_linesearch = 200;
    param.ftol           = 1e-3;
    param.wolfe          = 0.9;    
    param.min_step       = 1e-40;
    param.max_step       = 1.0;
    LBFGSpp::LBFGSBSolver<double> solver_mit(param);
    // Define OUR minimzer (using the same parameters)
    fdapde::LBFGSB<n> solver_our;

    // Define bounds vector (for constrained optimization)
    int min = -10, max = 10;
    VectorXd lb = VectorXd::Constant(n, min);
    VectorXd ub = VectorXd::Constant(n, max);

    // Define the function to minimize (n-dimensional Rosenbrock function in this case)
    ScalarField<n> fun;
    fun = [n](Eigen::Matrix<double,n,1> x) -> double {
        double sum = 0.0;
        for (int i = 0; i < n - 1; ++i) {
            double t1 = x[i+1] - x[i]*x[i];
            double t2 = 1 - x[i];
            sum += 100.0 * t1 * t1 + t2 * t2;
        }
        return sum;
    };
    Rosenbrock fun_mit{n};

    // Set up output file stream
    std::ofstream file("simulations_rosenbrock.csv");
    file << "TestId,OptId,NumIter,CompTime,MinPoint,Fx,DistSol,DistBetween" << std::endl;
    
    // Perform multiple simulations
    int num_simulations = 50;
    float dist_mit, dist_our, dist_between;
    VectorXd solution = VectorXd::Constant(n, 1.0);

    for(int i=0; i < num_simulations; ++i) {
        // Random initial guess
        Eigen::VectorXd x0(n);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(min, max);

        for(int i=0; i < x0.size(); ++i) {
            x0(i) = dis(gen);
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
        solver_our.optimize(fun, x0, fdapde::BacktrackingLineSearch());
        // Stop timer (OUR)
        auto end_our = std::chrono::high_resolution_clock::now();
        // Convert timer to nanoseconds (OUR)
        auto elapsed_time_our = std::chrono::duration_cast<std::chrono::nanoseconds>(end_our - start_our);
        long long elapsed_time_our_raw = elapsed_time_our.count();
        // Compute required distances
        dist_mit = (x - solution).norm();
        dist_our = (solver_our.optimum() - solution).norm();
        dist_between = (solver_our.optimum() - x).norm();
        // Output
        // MIT Line
        std::string separator_mit = "";
        file << i+1 << "," << "MIT" << "," << niter << "," << elapsed_time_mit_raw << ",";
        for(const auto value : x) { file << separator_mit << value; separator_mit = ";";}
        file << "," << fx << "," << dist_mit << "," << dist_between << std::endl;
        // OUR line
        std::string separator_our = "";
        file << i+1 << "," << "OUR" << "," << solver_our.n_iter() << "," << elapsed_time_our_raw << ",";
        for(const auto value : solver_our.optimum()) { file << separator_our << value; separator_our = ";";}
        file << "," << solver_our.value() << "," << dist_our << "," << dist_between << std::endl;
    }
    
    return 0;
}