// This file is derived from Eigen, a lightweight C++ template library
// for linear algebra, to enable user-defined preconditioners.

namespace Eigen {

// Diagonal preconditioner - just used as an example
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
template <typename _Scalar>
class LUPreconditioner
{
    typedef _Scalar Scalar;
    typedef Matrix<Scalar,Dynamic,1> Vector;
  public:
    typedef typename Vector::StorageIndex StorageIndex;
    typedef PartialPivLU<MatrixXd> LU;
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

} // end namespace Eigen
