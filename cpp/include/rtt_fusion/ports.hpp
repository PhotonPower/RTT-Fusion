// Port types of the Fusion (C++ prototype of `rtt-fusion`, Fusion ADR draft point 9).
//
// Mirrors src/rtt_fusion/ports.py and src/rtt_fusion/convert.py; conventions in
// docs/konventionen.md: mm, µm, rad; fields ∝ exp(i(k·r − ωt)); Im k_z ≥ 0 with the grazing rule;
// transverse fields (E_x, E_y), E_z from transversality. Acceptance: reference/*.result.json.
#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <optional>
#include <vector>

namespace rtt_fusion {

using cplx = std::complex<double>;
using Vec2c = std::array<cplx, 2>;
using Vec3 = std::array<double, 3>;
using Vec3c = std::array<cplx, 3>;

inline constexpr double kPi = 3.141592653589793;
inline constexpr double kGrazingTolerance = 1e-12;  // relative to |ñ²k₀²| (BWPM, Raytatouille)

// Port plane z = const with direction ±1, medium ñ and vacuum wavelength.
struct Port {
  double z_mm = 0.0;
  int direction = 1;
  cplx index{1.0, 0.0};
  double wavelength_um = 1.0;

  double k0() const { return 2.0 * kPi / (wavelength_um * 1e-3); }
  // k_z = sqrt(ñ²k₀² − k⊥²), Im ≥ 0 (Re ≥ 0 if real), grazing rule.
  cplx kz(double kx, double ky) const;
};

// Representation S: plane waves, phase 1 at the port origin; optional orders and dΩ weights.
struct PlaneWaves {
  Port port;
  std::vector<double> kx, ky;              // 1/mm
  std::vector<Vec2c> field;                // (E_x, E_y) per wave
  std::vector<std::array<int, 2>> orders;  // empty if none
  std::vector<double> solid_angle;         // empty if discrete waves

  std::size_t size() const { return kx.size(); }
  std::vector<cplx> ez() const;               // −direction·(k_x E_x + k_y E_y)/k_z, 0 if k_z = 0
  std::vector<double> power_density() const;  // Re(k_z)·|E|² (with E_z)
};

// Representation G: (E_x, E_y) on the centred grid x_m = (m − nx/2)·dx, row-major (ny, nx).
struct GriddedField {
  Port port;
  std::size_t nx = 0, ny = 0;
  double dx_mm = 0.0, dy_mm = 0.0;
  std::vector<Vec2c> field;  // size nx·ny, index j·nx + i
};

// G → S on the centred grid: Ẽ(k) = (1/N) Σ E(r) e^(−ik·r), numpy FFT order and rounding of k.
PlaneWaves to_plane_waves(const GriddedField& g);
// Σ Re(k_z)|Ẽ|² · N·dx·dy (Parseval).
double power(const GriddedField& g);
// (2π)²/k Σ |A|² dΩ for weighted spectra (A = field/dΩ), optionally within sin θ ≤ NA/n.
double spectrum_power(const PlaneWaves& s, std::optional<double> cone_na = std::nullopt);

// Gauss–Legendre in cos θ on (0, 1], uniform φ at (j + ½)·2π/n_phi; θ from +z or −z.
void hemisphere(int n_theta, int n_phi, bool forward, std::vector<Vec3>& directions,
                std::vector<double>& weights);

// Far field E∞ (E_s ≈ E∞ e^(ikr)/r, r in mm) of a body at `center` → plane waves at the port
// (Weyl; radial part of E∞ removed). Directions must point through the port.
PlaneWaves far_field_to_spectrum(const Port& port, const Vec3& center_mm,
                                 const std::vector<Vec3>& directions,
                                 const std::vector<double>& weights,
                                 const std::vector<Vec3c>& e_inf);

// Richards–Wolf beam (uniform pupil, pupil polarisation (e_x, e_y)) as weighted plane waves.
PlaneWaves richards_wolf(const Port& port, double na, const Vec3& focus_mm, const Vec2c& pupil,
                         int n_theta, int n_phi);

}  // namespace rtt_fusion
