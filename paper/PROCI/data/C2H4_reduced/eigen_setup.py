from setuptools import setup, Extension
import pybind11, os

# Try common Eigen include paths if env var is not set
eigen_inc = os.environ.get("EIGEN3_INCLUDE_DIR")
if not eigen_inc:
    for p in ("/opt/homebrew/include/eigen3", "/usr/local/include/eigen3", "/usr/include/eigen3"):
        if os.path.isdir(p):
            eigen_inc = p
            break

ext_modules = [
    Extension(
        "eigen_gmres",
        sources=["eigen_gmres.cpp"],
        include_dirs=[pybind11.get_include()] + ([eigen_inc] if eigen_inc else []),
        language="c++",
        extra_compile_args=["-std=c++17", "-O3"],
    )
]

setup(
    name="eigen_gmres",
    version="0.1",
    ext_modules=ext_modules,
)

