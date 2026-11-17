// Python binding of the tabulated diffraction interaction (stage R12, C++) for the grating block:
// rtt_fusion._diffraction.prt_many – the per-ray loop of `graph.grating` in one call. Optional:
// without it (no compiler at install time) the block keeps its Python loop.
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstring>
#include <stdexcept>
#include <vector>

#include "rtt_fusion/diffraction_table.hpp"

namespace py = pybind11;
using cplx = std::complex<double>;
using Doubles = py::array_t<double, py::array::c_style | py::array::forcecast>;
using Complexes = py::array_t<cplx, py::array::c_style | py::array::forcecast>;
using Flags = py::array_t<bool, py::array::c_style | py::array::forcecast>;

namespace {

std::vector<double> vec(const Doubles& a) { return {a.data(), a.data() + a.size()}; }

rtt_fusion::DiffractionTable table(const std::vector<int>& orders, bool reflected, const Doubles& lam,
                                   const Doubles& sx, const Doubles& sy, double period, double n_in,
                                   double n_out, const Complexes& jones) {
  rtt_fusion::DiffractionTable t;
  t.orders = orders;
  t.reflected = reflected;
  t.wavelengths_um = vec(lam);
  t.sigma_x = vec(sx);
  t.sigma_y = vec(sy);
  t.period_um = period;
  t.n_in = n_in;
  t.n_out = n_out;
  const std::size_t cells = orders.size() * t.wavelengths_um.size() * t.sigma_x.size() * t.sigma_y.size();
  if (static_cast<std::size_t>(jones.size()) != 4 * cells)
    throw std::invalid_argument("jones must have the shape (orders, λ, σx, σy, 2, 2)");
  t.jones.resize(cells);
  std::memcpy(t.jones.data(), jones.data(), sizeof(cplx) * 4 * cells);
  t.validate();
  return t;
}

}  // namespace

PYBIND11_MODULE(_diffraction, m) {
  m.doc() = "Tabulated diffraction interaction in C++ (RTT-Fusion stage R12) for the grating block";
  m.def(
      "prt_many",
      [](const std::vector<int>& orders, bool reflected, const Doubles& lam_axis, const Doubles& sx,
         const Doubles& sy, double period, double n_in, double n_out, const Complexes& jones, int order,
         double lam, const Doubles& sigma, const Flags& alive) {
        const auto t = table(orders, reflected, lam_axis, sx, sy, period, n_in, n_out, jones);
        if (sigma.ndim() != 2 || sigma.shape(1) != 2) throw std::invalid_argument("sigma must be (n, 2)");
        const auto n = static_cast<std::size_t>(sigma.shape(0));
        if (static_cast<std::size_t>(alive.size()) != n) throw std::invalid_argument("alive must have n entries");
        Complexes prt({static_cast<py::ssize_t>(n), py::ssize_t{3}, py::ssize_t{3}});
        Flags ok(static_cast<py::ssize_t>(n));
        std::vector<unsigned char> in(n), flags(n);
        for (std::size_t i = 0; i < n; ++i) in[i] = alive.data()[i] ? 1 : 0;
        {
          py::gil_scoped_release free;
          rtt_fusion::prt_many(t, order, lam, n, sigma.data(), in.data(),
                               reinterpret_cast<rtt_fusion::Prt*>(prt.mutable_data()), flags.data());
        }
        for (std::size_t i = 0; i < n; ++i) ok.mutable_data()[i] = flags[i] != 0;
        return py::make_tuple(prt, ok);
      },
      py::arg("orders"), py::arg("reflected"), py::arg("wavelengths_um"), py::arg("sigma_x"), py::arg("sigma_y"),
      py::arg("period_um"), py::arg("n_in"), py::arg("n_out"), py::arg("jones"), py::arg("order"),
      py::arg("wavelength_um"), py::arg("sigma"), py::arg("alive"),
      "Local PRT matrices (n, 3, 3) of `order` for incidences sigma (n, 2) and where they are valid (n,)");
}
