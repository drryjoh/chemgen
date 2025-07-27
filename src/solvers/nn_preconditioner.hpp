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

        // using Matrix_ = Matrix<Scalar, M, M>;
        using Matrix_ = SparseMatrix<Scalar>;

    public:
        NNPreconditioner() = default;

        //eigen stuff
        template <typename MatType>
        NNPreconditioner &analyzePattern(const MatType &) { return *this; }

        //compute
        template <typename MatType>
        void compute(const MatType &A)
        {
            //permute input
            A.makeCompressed();
            // assert(A_.valueSize() == INPUT_DIM && "A does not match sparsity count");
            PermutationMatrix<n_variables, n_variables> perm; // column perm (input)
            perm.indices() = perm_indices;
            B = perm.transpose() * A * perm;

            //pack all 3256 entries from A_.valuePtr()
            std::array<Scalar, INPUT_DIM> input_arr;
            std::copy_n(B.valuePtr(), INPUT_DIM, input_arr.begin());

            //NN
            auto P = MLP_LU<Scalar>(input_arr);

            //dense matrix
            for (int i = 0, k = 0; i < M; ++i)
                for (int j = 0; j < M; ++j, ++k)
                    P_eig(i, j) = P[k];
            //sparse matrix
            P_eig.setZero();  
            for (int i = 0, k = 0; i < M; ++i)
                for (int j = 0; j < M; ++j, ++k)
                    if (P[k] != 0)  
                    {
                        P_eig.insert(i, j) = P[k];
                    }
            P_eig.makeCompressed();
        }

        //apply preconditioner
        template <typename Rhs>
        Rhs solve(const Rhs &b) const
        // { //INV_A
        //     return P_eig * b;
        // }
        // { //LU
        //     Rhs y = P_eig.template triangularView<UnitLower>().solve(b);
        //     Rhs x = P_eig.template triangularView<Upper>().solve(y);
        //     return x;
        // }
        { //LU with permutation
            x = perm.transpose() * b;
            P_eig.matrixLU().template triangularView<UnitLower>().solveInPlace(x);
            P_eig.matrixLU().template triangularView<Upper>().solveInPlace(x);
            x = perm * x;
        }

        //eigen stuff
        ComputationInfo info() const { return Success; }

    private:
        Matrix_ P_eig;
    };

}