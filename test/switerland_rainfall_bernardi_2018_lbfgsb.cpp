// TO BE runned including LBFGSpp !!!

#include <LBFGSB.h>
#include <LBFGSpp/LineSearchBacktracking.h>
#include <fdaPDE/fdapde.h>

#include <iostream>
#include <random>
#include <string>
#include <unsupported/Eigen/SparseExtra>
#include <vector>
using namespace fdapde;

int main() {
    std::cout << "\t --- Switzerland Rainfall ---" << std::endl;
    // geometry
    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
    constexpr int n_params = 2;

    auto K_ = [](double theta, double gamma) -> Eigen::Matrix<double, 2, 2> {
        Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Zero(2, 2);
        Q << std::cos(theta), std::sin(theta), -std::sin(theta), std::cos(theta);
        Eigen::Matrix<double, 2, 2> Sigma = Eigen::Matrix<double, 2, 2>::Zero(2, 2);
        Sigma << 1. / std::sqrt(gamma), 0., 0., std::sqrt(gamma);
        Eigen::Matrix<double, 2, 2> res = Q * Sigma * Q.inverse();
        return res;
    };

    std::string data_path = "../data/switzerland-rainfall/";
    Triangulation<2, 2> D(data_path + "points.csv", data_path + "cells.csv", data_path + "boundary.csv", true, true);

    // data
    GeoFrame data(D);
    auto& l = data.insert_scalar_layer<POINT>("layer", data_path + "locs.csv");
    l.load_csv<double>(data_path + "response.csv");

    std::cout << "\t --- GEOFRAME ---" << std::endl;
    std::cout << l << std::endl;

    double mu_old = 0.5;   // mu_{opt} = 1.;

    double rho = 0.5;
    vector_t rho_seq = vector_t::Zero(13, 1);
    rho_seq << 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99;

    int n_rho = rho_seq.rows();
    std::cout << "D.measure(): " << D.measure() << std::endl;
    std::cout << n_rho << std::endl;
    vector_t gcv_vec = vector_t::Zero(n_rho, 1);
    matrix_t mu_optim = matrix_t::Zero(n_rho, 2);

    int n_lambda = 50;
    matrix_t lambda_grid(n_lambda, 1);
    for (int i = 0; i < n_lambda; ++i) { lambda_grid(i, 0) = std::pow(10, -5.0 + 0.25 * i); }
    // std::cout << lambda_grid << std::endl;
    GridSearch<1> optimizer;

    vector_t x0 = vector_t::Zero(2, 1);
    x0 << std::numbers::pi / 2., 5.;
    // std::cout << "x0 " << x0 << std::endl;
    // std::cout << "K(x0): \n" << K_(x0(0, 0), x0(1, 0)) << std::endl;
    ScalarField<2> PC;
    PC = [&D, &data, &rho, &K_](vector_t x) -> double {
        double theta = x(0, 0);
        double gamma = x(1, 0);

        Eigen::Matrix<double, 2, 2> K = K_(theta, gamma);

        FeSpace Vh(D, P1<1>);
        TrialFunction f(Vh);
        TestFunction v(Vh);
        auto a = integral(D)(dot(K * grad(f), grad(v)));
        ZeroField<2> u;
        auto F = integral(D)(u * v);
        SRPDE m("V1 ~ f", data, fe_ls_elliptic(a, F));
        double lambda = rho / (D.measure() * (1. - rho));   // data[0].rows() / D.measure()
        //  AAA, nella versione corrente della libreria viene diviso per il numero di obs !!!
        // NON c'è bisogno di moltiplicare per n !!!
        m.fit(lambda);
        return ((m.fitted() - m.response()).squaredNorm());
    };
    
    struct function_to_optimize {
        function_to_optimize(const ScalarField<2>& Field) : Field_(Field) { }

        double operator()(const vector_t& x, vector_t& grad) {
            double Fx = Field_(x);
            grad = Field_.gradient()(x);
            return Fx;
        }

        const ScalarField<2>& Field_;
    };

    function_to_optimize F_(PC);

    // Set up optimizer
    LBFGSB<2> lbfgsb;

    // theta, gamma
    vector_t lower_bounds = vector_t::Zero(2, 1);
    lower_bounds << 0., 1.;
    vector_t upper_bounds = vector_t::Zero(2, 1);
    upper_bounds << std::numbers::pi, 1000.;

    for (int i = 0; i < rho_seq.size(); ++i) {
        // optimize SSE -----------------------------------
        std::cout << "\t --- rho " << rho_seq(i, 0) << " ---" << std::endl;
        rho = rho_seq(i, 0);
        vector_t x = vector_t::Zero(2, 1);
        if (i == 0)
            x = x0;
        else {
            // x = x0;
            x(0, 0) = mu_optim(i - 1, 0);   // Bernardi
            x(1, 0) = mu_optim(i - 1, 1);
        }

        lbfgsb.optimize(F_, x, lower_bounds, upper_bounds, BacktrackingLineSearch());
        std::cout << "num_iter : " << lbfgsb.n_iter() << std::endl;

        // see Bernardi
        if (x(0, 0) == 0.0) {
            x(0, 0) = std::numbers::pi;
            lbfgsb.optimize(F_, x, lower_bounds, upper_bounds, BacktrackingLineSearch());
            std::cout << "num_iter : " << lbfgsb.n_iter() << std::endl;
        } else if (x(0, 0) == std::numbers::pi) {
            std::cout << "sfigato 2" << std::endl;
            x(0, 0) = 0.0;
            lbfgsb.optimize(F_, x, lower_bounds, upper_bounds, BacktrackingLineSearch());
            std::cout << "num_iter : " << lbfgsb.n_iter() << std::endl;
        }

        std::cout << "iter: " << lbfgsb.n_iter() << std::endl;
        std::cout << "optimum: " << lbfgsb.optimum() << std::endl;

        mu_optim(i, 0) = x(0, 0);
        mu_optim(i, 1) = x(1, 0);

        auto K = K_(x(0, 0), x(1, 0));

        // evaluate outer criteria (based on GCV) ---------
        FeSpace Vh(D, P1<1>);
        TrialFunction f(Vh);
        TestFunction v(Vh);

        auto a = integral(D)(dot(K * grad(f), grad(v)));
        ZeroField<2> u;
        auto F = integral(D)(u * v);
        // modeling
        SRPDE m("V1 ~ f", data, fe_ls_elliptic(a, F));
        optimizer.optimize(m.gcv(100, 476813), lambda_grid);

        gcv_vec(i, 0) = optimizer.value();
        std::cout << "\t --- --- ---" << std::endl;
    }

    std::cout << "\ngcv_vec: \n" << gcv_vec << std::endl;

    Eigen::Index minRow;
    gcv_vec.col(0).minCoeff(&minRow);
    std::cout << "gcv min at " << minRow << "  value: " << gcv_vec(minRow, 0) << std::endl;

    std::cout << "Kappa: \n" << K_(mu_optim(minRow, 0), mu_optim(minRow, 1)) << std::endl;
    Eigen::saveMarket(K_(mu_optim(minRow, 0), mu_optim(minRow, 1)), data_path + "kappa.mtx");
    return 0;
}