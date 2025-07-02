// This file provides user-defined preconditioners for Eigen.

namespace Eigen {

// Diagonal preconditioner - just used as an example
// This is just a copy of Eigen's DiagonalPreconditioner
template <typename _Scalar>
class JacobiPreconditioner
{
    typedef _Scalar Scalar;
    typedef Matrix<Scalar,Dynamic,1> Vector;
  public:
    typedef typename Vector::StorageIndex StorageIndex;
    enum {
      ColsAtCompileTime = Dynamic,
      MaxColsAtCompileTime = Dynamic
    };

    JacobiPreconditioner() : m_isInitialized(false) {}

    template<typename MatType>
    explicit JacobiPreconditioner(const MatType& mat) : m_invdiag(mat.cols())
    {
      compute(mat);
    }

    EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return m_invdiag.size(); }
    EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return m_invdiag.size(); }

    template<typename MatType>
    JacobiPreconditioner& analyzePattern(const MatType& )
    {
      return *this;
    }

    template<typename MatType>
    JacobiPreconditioner& factorize(const MatType& mat)
    {
      m_invdiag.resize(mat.cols());
      for(int j=0; j<mat.outerSize(); ++j)
      {
        typename MatType::InnerIterator it(mat,j);
        while(it && it.index()!=j) ++it;
        if(it && it.index()==j && it.value()!=Scalar(0))
          m_invdiag(j) = Scalar(1)/it.value();
        else
          m_invdiag(j) = Scalar(1);
      }
      m_isInitialized = true;
      return *this;
    }

    template<typename MatType>
    JacobiPreconditioner& compute(const MatType& mat)
    {
      return factorize(mat);
    }

    /** \internal */
    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
      // Elementwise multiplication
      x = m_invdiag.array() * b.array() ;
    }

    template<typename Rhs> inline const Solve<JacobiPreconditioner, Rhs>
    solve(const MatrixBase<Rhs>& b) const
    {
      eigen_assert(m_isInitialized && "JacobiPreconditioner is not initialized.");
      eigen_assert(m_invdiag.size()==b.rows()
                && "JacobiPreconditioner::solve(): invalid number of rows of the right hand side matrix b");
      return Solve<JacobiPreconditioner, Rhs>(*this, b.derived());
    }

    ComputationInfo info() { return Success; }

  protected:
    Vector m_invdiag;
    bool m_isInitialized;
};

// LU (exact) preconditioner - just used as an example
// Note: in practice, if LU decomposition is desired, should use SparseLU for sparse matrices
template <typename _Scalar>
class LUPreconditioner
{
    typedef _Scalar Scalar;
    typedef Matrix<Scalar,Dynamic,1> Vector;
    using Matrix_ = Matrix<Scalar, n_variables, n_variables>;
  public:
    typedef typename Vector::StorageIndex StorageIndex;
    typedef PartialPivLU<Matrix_> LU;
    enum {
      ColsAtCompileTime = Dynamic,
      MaxColsAtCompileTime = Dynamic
    };

    LUPreconditioner() : m_isInitialized(false) {}

    template<typename MatType>
    explicit LUPreconditioner(const MatType& mat)
    {
      compute(mat);
    }

    EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return lu.rows(); }
    EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return lu.cols(); }

    template<typename MatType>
    LUPreconditioner& analyzePattern(const MatType& )
    {
      return *this;
    }

    template<typename MatType>
    LUPreconditioner& factorize(const MatType& mat)
    {
      // Get LU decomposition
      lu = LU(mat);
      m_isInitialized = true;
      return *this;
    }

    template<typename MatType>
    LUPreconditioner& compute(const MatType& mat)
    {
      // analyzePattern(mat);
      factorize(mat);
      return *this;
    }

    /** \internal */
    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
      /* The decomposition PA = LU can be rewritten as A = P^{-1} L U.
       * So we proceed as follows:
       * Step 1: compute c = Pb.
       * Step 2: replace c by the solution x to Lx = c.
       * Step 3: replace c by the solution x to Ux = c.
       */

      // Step 1
      x = lu.permutationP() * b;

      // Step 2
      lu.matrixLU().template triangularView<UnitLower>().solveInPlace(x);

      // Step 3
      lu.matrixLU().template triangularView<Upper>().solveInPlace(x);

      /* Below: directly access lower and upper triangular parts and print out */
      // Note: triangularView creates a reference so can read/write
      // Matrix_ m_lu = lu.matrixLU(); // LU decomposition where L (excluding unit diagonal) is stored in lower triangular part and U is stored in upper triangular part
      // Matrix_ m_l = m_lu.triangularView<UnitLower>(); // lower triangular part of m_lu with unit diagonal
      // Matrix_ m_u = m_lu.triangularView<Upper>(); // upper triangular part of m_lu

      // std::cout << "m_lu = \n" << m_lu << std::endl;
      // std::cout << "m_l = \n" << m_l << std::endl;
      // std::cout << "m_u = \n" << m_u << std::endl;

      // std::exit(0);
    }

    template<typename Rhs> inline const Solve<LUPreconditioner, Rhs>
    solve(const MatrixBase<Rhs>& b) const
    {
      eigen_assert(m_isInitialized && "LUPreconditioner is not initialized.");
      eigen_assert(lu.cols()==b.rows()
                && "LUPreconditioner::solve(): invalid number of rows of the right hand side matrix b");
      return Solve<LUPreconditioner, Rhs>(*this, b.derived());
    }

    ComputationInfo info() { return Success; }

  protected:
    LU lu;
    bool m_isInitialized;
};

// Inverse (exact) preconditioner - just used as an example
// Should not actually be used in practice
template <typename _Scalar>
class InversePreconditioner
{
    typedef _Scalar Scalar;
    typedef Matrix<Scalar,Dynamic,1> Vector;
    using Matrix_ = Matrix<Scalar, n_variables, n_variables>;
  public:
    typedef typename Vector::StorageIndex StorageIndex;
    typedef PartialPivLU<Matrix_> LU;
    enum {
      ColsAtCompileTime = Dynamic,
      MaxColsAtCompileTime = Dynamic
    };

    InversePreconditioner() : m_isInitialized(false) {}

    template<typename MatType>
    explicit InversePreconditioner(const MatType& mat)
    {
      compute(mat);
    }

    EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return A.rows(); }
    EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return A.cols(); }

    template<typename MatType>
    InversePreconditioner& analyzePattern(const MatType& mat)
    {
      A = mat; // if mat is sparse, then converted to dense here
      return *this;
    }

    template<typename MatType>
    InversePreconditioner& factorize(const MatType& mat)
    {
      invA = A.inverse();
      m_isInitialized = true;
      return *this;
    }

    template<typename MatType>
    InversePreconditioner& compute(const MatType& mat)
    {
      analyzePattern(mat);
      factorize(mat);
      return *this;
    }

    /** \internal */
    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
      // Very roundabout way of doing A^{-1} * b for demonstration purposes
      // Create a copy of A^{-1} (not actually necessary - just for demonstration)
      // Get LU decomposition of A^{-1}, such that A^{-1} = P^{-1} * L * U = P^T * L * U (since P is orthogonal)
      // Do x = P^T * L * U * b

      Matrix_ invA_;
      invA_.setZero();

      // This is how you access and assign values to matrix entries
      // Note: since we're just copying values over, can just do invA_ = invA
      for (int i = 0; i < n_variables; ++i)
      {
        for (int j = 0; j < n_variables; ++j)
        {
          invA_(i,j) = invA(i,j);
        }
      }

      // Get LU decomposition
      LU lu_invA = LU(invA_);

      Matrix_ m_lu = lu_invA.matrixLU(); // LU decomposition of A^{-1} where L (excluding unit diagonal) is stored in lower triangular part and U is stored in upper triangular part
      Matrix_ m_l = m_lu.template triangularView<UnitLower>(); // lower triangular part of m_lu with unit diagonal
      Matrix_ m_u = m_lu.template triangularView<Upper>(); // upper triangular part of m_lu

      x = m_u * b;
      x = m_l * x;
      x = lu_invA.permutationP().transpose() * x;
    }

    template<typename Rhs> inline const Solve<InversePreconditioner, Rhs>
    solve(const MatrixBase<Rhs>& b) const
    {
      eigen_assert(m_isInitialized && "InversePreconditioner is not initialized.");
      eigen_assert(A.cols()==b.rows()
                && "InversePreconditioner::solve(): invalid number of rows of the right hand side matrix b");
      return Solve<InversePreconditioner, Rhs>(*this, b.derived());
    }

    ComputationInfo info() { return Success; }

  protected:
    bool m_isInitialized;
    Matrix_ A;
    Matrix_ invA;
};

} // end namespace Eigen