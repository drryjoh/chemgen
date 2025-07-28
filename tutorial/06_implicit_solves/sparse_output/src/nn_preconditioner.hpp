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
        using Perm_ = Eigen::PermutationMatrix<M, M>;

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
            // PermutationMatrix<n_variables, n_variables> perm; // column perm (input)
            perm.indices() = perm_indices;
            B = perm.transpose() * A * perm;
            B.makeCompressed();
            // assert(B.valueSize() == INPUT_DIM && "A does not match sparsity count");

            //pack all 3256 entries from A_.valuePtr()
            std::array<Scalar, INPUT_DIM> input_arr;
            std::copy_n(B.valuePtr(), INPUT_DIM, input_arr.begin());

            //NN
            auto P = MLP_LU<Scalar>(input_arr);

            // //dense matrix
            // for (int i = 0, k = 0; i < M; ++i)
            //     for (int j = 0; j < M; ++j, ++k)
            //         P_eig(i, j) = P[k];
            //sparse matrix
            P_eig.resize(M, M);
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
            Rhs x = perm.transpose() * b;
            for (int i = 0; i < M; ++i) {
                Scalar sum = Scalar(0);
                for (typename Matrix_::InnerIterator it(P_eig, i); it; ++it) {
                    int j = it.col();
                    if (j < i) sum += it.value() * x[j];
                    else break;   
                }
                x[i] -= sum;
            }
            for (int i = M - 1; i >= 0; --i) {
                Scalar sum = Scalar(0), diag = Scalar(0);
                for (typename Matrix_::InnerIterator it(P_eig, i); it; ++it) {
                    int j = it.col();
                    if (j > i)      sum  += it.value() * x[j];
                    else if (j == i) diag = it.value();
                }
                x[i] = (x[i] - sum) / diag;
            }
            x = perm * x;
            return x;
        }

        //eigen stuff
        ComputationInfo info() const { return Success; }

    private:
        Matrix_ P_eig;
        Perm_ perm;
        std::array<int, M> perm_indices;
        Matrix_ B;
    };

}