// nn_preconditioner.hpp
#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <array>
#include <algorithm> 
#include <cassert>
#include "MLP_LU.hpp" //CODEJENN

namespace Eigen
{

    template <typename Scalar>
    class NNPreconditioner
    {
        static constexpr int M = 97;
        static constexpr int INPUT_DIM = 3256;   // number of stored non-zeros
        static constexpr int OUTPUT_DIM = M * M; // 9409

        using Matrix_ = Matrix<Scalar, M, M>;

    public:
        NNPreconditioner() = default;

        //eigen stuff
        template <typename MatType>
        NNPreconditioner &analyzePattern(const MatType &) { return *this; }

        //compute
        template <typename MatType>
        void compute(const MatType &A)
        {
            //assert input size
            Eigen::SparseMatrix<Scalar> A_ = A;
            A_.makeCompressed();
            // assert(A_.valueSize() == INPUT_DIM && "A does not match sparsity count");

            //pack all 3256 entries from A_.valuePtr()
            std::array<Scalar, INPUT_DIM> input_arr;
            std::copy_n(A_.valuePtr(), INPUT_DIM, input_arr.begin());

            //NN
            auto P = MLP_LU<Scalar>(input_arr);

            //copy into P_eig
            for (int i = 0, k = 0; i < M; ++i)
                for (int j = 0; j < M; ++j, ++k)
                    P_eig(i, j) = P[k];
        }

        //apply preconditioner
        template <typename Rhs>
        Rhs solve(const Rhs &b) const
        // { //INV_A
        //     return P_eig * b;
        // }
        { //LU
            Rhs y = P_eig.template triangularView<UnitLower>().solve(b);
            Rhs x = P_eig.template triangularView<Upper>().solve(y);
            return x;
        }

        //eigen stuff
        ComputationInfo info() const { return Success; }

    private:
        Matrix_ P_eig;
    };

}