"""The optional C++ extension rtt_fusion._diffraction (the per-ray loop of the grating block,
stage R12 in C++); everything else is in pyproject.toml. Optional: without a compiler the package
installs without it and the block keeps its Python loop."""

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

diffraction = Pybind11Extension(
    "rtt_fusion._diffraction",
    ["cpp/python/diffraction_module.cpp", "cpp/src/diffraction_table.cpp"],
    include_dirs=["cpp/include"],
    cxx_std=20,
)
diffraction.optional = True

setup(ext_modules=[diffraction], cmdclass={"build_ext": build_ext})
