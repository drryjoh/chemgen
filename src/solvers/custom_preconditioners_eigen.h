// This file is derived from Eigen, a lightweight C++ template library
// for linear algebra, to enable user-defined preconditioners.

namespace Eigen {

// Diagonal preconditioner - just used as an example
template <typename _Scalar>
class CustomPreconditioner
{
    typedef _Scalar Scalar;
    typedef Matrix<Scalar,Dynamic,1> Vector;
  public:
    typedef typename Vector::StorageIndex StorageIndex;
    enum {
      ColsAtCompileTime = Dynamic,
      MaxColsAtCompileTime = Dynamic
    };

    CustomPreconditioner() : m_isInitialized(false) {}

    template<typename MatType>
    explicit CustomPreconditioner(const MatType& mat) : m_invdiag(mat.cols())
    {
      compute(mat);
    }

    EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT { return m_invdiag.size(); }
    EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT { return m_invdiag.size(); }

    template<typename MatType>
    CustomPreconditioner& analyzePattern(const MatType& )
    {
      return *this;
    }

    template<typename MatType>
    CustomPreconditioner& factorize(const MatType& mat)
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
    CustomPreconditioner& compute(const MatType& mat)
    {
      return factorize(mat);
    }

    /** \internal */
    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
      x = m_invdiag.array() * b.array() ;
    }

    template<typename Rhs> inline const Solve<CustomPreconditioner, Rhs>
    solve(const MatrixBase<Rhs>& b) const
    {
      eigen_assert(m_isInitialized && "CustomPreconditioner is not initialized.");
      eigen_assert(m_invdiag.size()==b.rows()
                && "CustomPreconditioner::solve(): invalid number of rows of the right hand side matrix b");
      return Solve<CustomPreconditioner, Rhs>(*this, b.derived());
    }

    ComputationInfo info() { return Success; }

  protected:
    Vector m_invdiag;
    bool m_isInitialized;
};

} // end namespace Eigen
