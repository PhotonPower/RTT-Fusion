// Acceptance of the C++ port types against the reference files of the Python prototype (R7).
#include <cmath>
#include <cstdio>
#include <string>

#include "rtt_fusion/ports.hpp"
#include "rtt_fusion/results.hpp"

using namespace rtt_fusion;
namespace R = rtt_fusion::results;

static int failures = 0;

static void check(bool ok, const std::string& what, double got = 0, double limit = 0) {
  if (!ok) {
    ++failures;
    std::printf("FAIL %s (%.3e > %.3e)\n", what.c_str(), got, limit);
  } else {
    std::printf("ok   %s (%.3e)\n", what.c_str(), got);
  }
}

static double max_diff(const std::vector<cplx>& a, const std::vector<cplx>& b) {
  double m = 0;
  for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
  return a.size() == b.size() ? m : HUGE_VAL;
}

static std::vector<cplx> flat(const std::vector<Vec2c>& v) {
  std::vector<cplx> out;
  for (const auto& e : v) out.insert(out.end(), {e[0], e[1]});
  return out;
}

static std::vector<cplx> as_cplx(const std::vector<double>& v) { return {v.begin(), v.end()}; }

int main(int argc, char** argv) {
  const std::string dir = argc > 1 ? argv[1] : "reference";

  {  // G → S, E_z, power
    auto d = R::load(dir + "/ports_grid_to_spectrum.result.json");
    auto p = d.at("parameters");
    auto tol = d.at("tolerance");
    const double abs_tol = R::number(tol.at("abs"));
    GriddedField g;
    g.port = Port{p.at("z_mm"), p.at("direction"), {R::number(p.at("index")), 0.0}, p.at("wavelength_um")};
    g.nx = p.at("nx");
    g.ny = p.at("ny");
    g.dx_mm = p.at("dx_mm");
    g.dy_mm = p.at("dy_mm");
    auto f = R::complexes(d.at("inputs").at("field"));
    for (std::size_t q = 0; q < g.nx * g.ny; ++q) g.field.push_back({f[2 * q], f[2 * q + 1]});
    PlaneWaves s = to_plane_waves(g);
    auto e = d.at("expected");
    check(max_diff(as_cplx(s.kx), as_cplx(R::reals(e.at("kx_per_mm")))) < 1e-9, "G->S kx", max_diff(as_cplx(s.kx), as_cplx(R::reals(e.at("kx_per_mm")))), 1e-9);
    double dm = max_diff(flat(s.field), R::complexes(e.at("amplitude")));
    check(dm < abs_tol, "G->S amplitude", dm, abs_tol);
    dm = max_diff(s.ez(), R::complexes(e.at("ez")));
    check(dm < abs_tol, "G->S ez", dm, abs_tol);
    const double pw = R::scalar(e.at("power"));
    check(std::abs(power(g) - pw) < 1e-12 * std::abs(pw), "G power", std::abs(power(g) - pw) / pw, 1e-12);
  }
  {  // hemisphere quadrature
    auto d = R::load(dir + "/hemisphere_quadrature.result.json");
    std::vector<Vec3> dirs;
    std::vector<double> w;
    hemisphere(d.at("parameters").at("n_theta"), d.at("parameters").at("n_phi"), true, dirs, w);
    auto ed = R::reals(d.at("expected").at("forward_directions"));
    auto ew = R::reals(d.at("expected").at("forward_weights"));
    double m = 0;
    for (std::size_t j = 0; j < dirs.size(); ++j) {
      for (int c = 0; c < 3; ++c) m = std::max(m, std::abs(dirs[j][static_cast<std::size_t>(c)] - ed[3 * j + static_cast<std::size_t>(c)]));
      m = std::max(m, std::abs(w[j] - ew[j]));
    }
    check(m < 1e-13, "hemisphere nodes and weights", m, 1e-13);
  }
  {  // Richards-Wolf beam
    auto d = R::load(dir + "/focused_beam_spectrum.result.json");
    auto p = d.at("parameters");
    Port port{p.at("z_mm"), 1, {R::number(p.at("index")), 0.0}, p.at("wavelength_um")};
    auto fm = p.at("focus_mm");
    PlaneWaves s = richards_wolf(port, p.at("na"), {fm[0], fm[1], fm[2]}, {1.0, 0.0}, p.at("n_theta"), p.at("n_phi"));
    auto e = d.at("expected");
    const double abs_tol = R::number(d.at("tolerance").at("abs"));
    double dm = max_diff(flat(s.field), R::complexes(e.at("amplitude")));
    check(dm < abs_tol, "Richards-Wolf amplitude", dm, abs_tol);
    const double pw = spectrum_power(s), ref = R::scalar(d.at("reference").at("power"));
    check(std::abs(pw / ref - 1) < 1e-12, "Richards-Wolf power vs (2pi)^2/k pi NA^2", std::abs(pw / ref - 1), 1e-12);
  }
  {  // far field -> spectrum (dipole)
    auto d = R::load(dir + "/far_field_spectrum_dipole.result.json");
    auto p = d.at("parameters");
    Port port{p.at("z_mm"), p.at("direction"), {R::number(p.at("index")), 0.0}, p.at("wavelength_um")};
    auto c = p.at("center_mm");
    auto in = d.at("inputs");
    auto dr = R::reals(in.at("directions"));
    auto w = R::reals(in.at("weights"));
    auto ei = R::complexes(in.at("e_inf"));
    std::vector<Vec3> dirs;
    std::vector<Vec3c> e_inf;
    for (std::size_t j = 0; j < w.size(); ++j) {
      dirs.push_back({dr[3 * j], dr[3 * j + 1], dr[3 * j + 2]});
      e_inf.push_back({ei[3 * j], ei[3 * j + 1], ei[3 * j + 2]});
    }
    PlaneWaves s = far_field_to_spectrum(port, {c[0], c[1], c[2]}, dirs, w, e_inf);
    const double abs_tol = R::number(d.at("tolerance").at("abs"));
    double dm = max_diff(flat(s.field), R::complexes(d.at("expected").at("amplitude")));
    check(dm < abs_tol, "far field -> spectrum amplitude", dm, abs_tol);
    const double pw = spectrum_power(s), ref = R::scalar(d.at("reference").at("power"));
    check(std::abs(pw / ref - 1) < 1e-12, "dipole power k 4pi/3", std::abs(pw / ref - 1), 1e-12);
  }
  std::printf("%d failure(s)\n", failures);
  return failures == 0 ? 0 : 1;
}
