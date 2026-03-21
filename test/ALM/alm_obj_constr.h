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

// THIS FILE INTRODUCES EXAMPLES OF DATA STRUCTURES THAT COMPLY WITH THE
// AUGMENTED LAGRANGIAN METHOD

#include <fdaPDE/core/fdaPDE/optimization.h>
#include <unsupported/Eigen/SparseExtra>
#include <cstdio>
#include <vector>

using fdapde::ScalarField;

// Objective Type
// Example of ObjectiveT object to be used within the Lagrangian class.
// We create a wrapper based on the ScalarField class that exposes a call operator and a gradient method
// Remark : Gradient can either be set manually or be computed with finite differences if not provided explicitly
template<int N>
struct objFunction {
    using vector_t = std::conditional_t<N == Eigen::Dynamic, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, N, 1>>;
    // Scalar Field
    ScalarField<N> objective_;
    // Constructor
    objFunction(const ScalarField<N>& objective) : objective_(objective) {}
    // Call operator
    double operator()(const vector_t& x) const {
        return objective_(x);
    }
    // Gradient method
    vector_t gradient(const vector_t& x) const {
        return objective_.gradient()(x);
    }
};
// Constraint Type (introduction)
// Example of constraint object to be used (inside a vector) within the Lagrangian class
// The main idea is the same as Objective type, but we also add a bool attribute to identify the type
// of constraint (true = inequality, false = equality)
template<int N>
struct constrFunction {
    using vector_t = std::conditional_t<N == Eigen::Dynamic, Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<double, N, 1>>;
    // Scalar Field
    ScalarField<N> constr_;
    // Boolean
    bool is_inequality_;
    // Constructor
    constrFunction(const ScalarField<N>& constr, const bool is_inequality) : constr_(constr), is_inequality_(is_inequality) {}
    // Call operator
    double operator()(const vector_t& x) const {
        return constr_(x);
    }
    // Gradient method
    vector_t gradient(const vector_t& x) const {
        return constr_.gradient()(x);
    }
};
// Constraint Type (final)
// In order to comply with the ALM interface an std::vector of constrFunctions needs to be passed as argument
// to the solve method. This struct exactly serves that purpouse
template<int N>
struct constrList {
    // Data
    std::vector<constrFunction<N>> constraints_;
    // Constructor (vector)
    constrList(const std::vector<constrFunction<N>>& constraints) : constraints_(constraints) {}
    // Constructor (initializer list)
    constrList(std::initializer_list<constrFunction<N>> list) : constraints_(list) {}
    // size() method
    std::size_t size() const {return constraints_.size();}
    // [] operator
    const constrFunction<N>& operator[](std::size_t i) const {return constraints_[i];}
};