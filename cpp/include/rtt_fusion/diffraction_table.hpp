// Tabulated diffraction interaction in C++ (stage R12): the rules of the ADR draft
// docs/adr-entwurf-beugungstabelle.md (approved E2), mirroring src/rtt_fusion/couplings/jones_table.py.
//
// Jones matrices J_m (power normalised, transverse Cartesian components of the local surface
// frame x = grating vector, y = grooves, z = normal) over uniform axes (λ, σ_x, σ_y), σ = k⊥/k₀;
// cubic Lagrange interpolation (4 points per axis, tensor product); Rayleigh anomalies detected
// over the interpolation stencil; PRT matrix P_m = L_m J_m Π (I − k kᵀ) + k_m kᵀ.
#pragma once

#include <array>
#include <complex>
#include <vector>

namespace rtt_fusion {

using Jones = std::array<std::complex<double>, 4>;  // row-major 2×2
using Prt = std::array<std::complex<double>, 9>;    // row-major 3×3

struct DiffractionTable {
  std::vector<int> orders;
  bool reflected = false;
  std::vector<double> wavelengths_um, sigma_x, sigma_y;  // uniform, ≥ 4 points each
  double period_um = 1.0;
  double n_in = 1.0, n_out = 1.0;
  // jones[((o·Nλ + a)·Nx + b)·Ny + c] for order index o
  std::vector<Jones> jones;

  void validate() const;  // throws std::invalid_argument (diagnostic surface.table_invalid)
  // J_m at (λ, σ_x, σ_y); NaN entries if a node of the stencil is not propagating.
  Jones interpolate(int order, double wavelength_um, double sigma_x, double sigma_y) const;
  // A Rayleigh anomaly (|σ + m′λ/Λ x̂| = n, |m′| ≤ max_order, both media) crosses the stencil.
  bool near_rayleigh_anomaly(double wavelength_um, double sigma_x, double sigma_y,
                             int max_order = 8) const;
};

// 3×3 PRT matrix of an order for unit wave vectors from (σ_in, σ_out) and the indices.
Prt prt_matrix(const Jones& j, std::array<double, 2> sigma_in, std::array<double, 2> sigma_out,
               double n_in, double n_out, bool reflected);

}  // namespace rtt_fusion
