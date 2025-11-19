#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/IterativeSolvers>

namespace py = pybind11;

// Return object
struct GMRESResult {
    py::array_t<double> x;
    int iters;
    int info;     // 0 Success, 1 NoConvergence, 2 NumericalIssue
    double relres;
};

GMRESResult gmres_dense(py::array_t<double, py::array::c_style | py::array::forcecast> A_in,
                        py::array_t<double, py::array::c_style | py::array::forcecast> b_in,
                        double rtol = 1e-8,
                        int maxiter = 1000,
                        int restart = 50)
{
    auto A_buf = A_in.request();
    auto b_buf = b_in.request();
    if (A_buf.ndim != 2) throw std::runtime_error("A must be 2D");
    if (b_buf.ndim != 1) throw std::runtime_error("b must be 1D");

    const int n = static_cast<int>(A_buf.shape[0]);
    const int m = static_cast<int>(A_buf.shape[1]);
    if (n != m) throw std::runtime_error("A must be square");
    if (b_buf.shape[0] != n) throw std::runtime_error("b size mismatch");

    // Map RowMajor data from NumPy into Eigen
    using RowMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    auto* A_ptr = static_cast<double*>(A_buf.ptr);
    auto* b_ptr = static_cast<double*>(b_buf.ptr);
    Eigen::Map<RowMat> A(A_ptr, n, n);
    Eigen::Map<Eigen::VectorXd> b(b_ptr, n);

    Eigen::GMRES<RowMat> solver;
    solver.setTolerance(rtol);
    solver.setMaxIterations(maxiter);
    solver.set_restart(restart);
    solver.compute(A);

    Eigen::VectorXd x = solver.solve(b);

    // Pack output
    py::array_t<double> x_out(n);
    auto x_mut = x_out.mutable_unchecked<1>();
    for (int i = 0; i < n; ++i) x_mut(i) = x(i);

    double denom = std::max(1.0, b.norm());
    double relres = (A * x - b).norm() / denom;

    GMRESResult result;
    result.x = x_out;
    result.iters = static_cast<int>(solver.iterations());
    result.info = static_cast<int>(solver.info());
    result.relres = relres;
    return result;
}

PYBIND11_MODULE(eigen_gmres, m) {
    py::class_<GMRESResult>(m, "GMRESResult")
        .def_readonly("x", &GMRESResult::x)
        .def_readonly("iters", &GMRESResult::iters)
        .def_readonly("info", &GMRESResult::info)
        .def_readonly("relres", &GMRESResult::relres);

    m.def("gmres_dense", &gmres_dense,
          py::arg("A"), py::arg("b"),
          py::arg("rtol") = 1e-8,
          py::arg("maxiter") = 1000,
          py::arg("restart") = 50,
          "Eigen GMRES for dense A (RowMajor).");
}

