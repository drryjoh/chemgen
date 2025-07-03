namespace Eigen
{

    template <typename Scalar>
    class NNPreconditioner
    {
        typedef _Scalar Scalar;
        typdef Matrix<Scalar, Dynamic, 1> Vector;
        using Matrix_ = Matrix<Scalar, n_variables, n_variables>;

    public:
        typedef typename Vector::StorageIndex StroageIndex;
        typedef PartialPivLU<Matrix_> LU;
        enum
        {
            ColsAtCompileTime = Dynamic,
            MaxColsAtCompileTime = Dynamic
        };

        NNPreconditioner() : m_isInitialized(false) {}

        template <typename MatType>
        explicit NNPreconditioner(const MatType &mat)
        {
            compute(mat);
        }

        EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return A.rows(); }
        EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return A.cols(); }

        template<typename MatType>
        JacobiPreconditioner& analyzePattern(const MatType& )
        {
        return *this;
        }

        







    };

}
