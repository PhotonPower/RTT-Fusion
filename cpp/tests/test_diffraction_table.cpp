// R12: tabulated diffraction interaction in C++ against reference/jones_table_lamellar (R6/R7):
// interpolation and anomaly flags from the stored table, PRT matrices from FourierGrating's direct
// Jones matrices, P k = k_m and energy for random polarisations, validation of the table.
#include <cmath>
#include <cstdio>
#include <string>

#include "rtt_fusion/diffraction_table.hpp"
#include "rtt_fusion/results.hpp"

using namespace rtt_fusion;
namespace R = rtt_fusion::results;
using cplx = std::complex<double>;

static int failures = 0;
static void check(bool ok, const char* what, double v, double lim) {
  failures += ok ? 0 : 1;
  std::printf("%s %s (%.3e, limit %.3e)\n", ok ? "ok  " : "FAIL", what, v, lim);
}

int main(int argc, char** argv) {
  const std::string dir = argc > 1 ? argv[1] : "reference";
  auto d = R::load(dir + "/jones_table_lamellar.result.json");
  auto p = d.at("parameters");
  auto in = d.at("inputs");
  DiffractionTable t;
  for (auto& o : p.at("orders")) t.orders.push_back(o.get<int>());
  t.reflected = p.at("side") == "reflected";
  t.wavelengths_um = R::reals(in.at("wavelengths_um"));
  t.sigma_x = R::reals(in.at("sigma_x"));
  t.sigma_y = R::reals(in.at("sigma_y"));
  t.period_um = p.at("period_um");
  t.n_in = p.at("n_in");
  t.n_out = p.at("n_out");
  auto jv = R::complexes(in.at("jones"));
  for (std::size_t q = 0; q < jv.size() / 4; ++q) t.jones.push_back({jv[4 * q], jv[4 * q + 1], jv[4 * q + 2], jv[4 * q + 3]});
  t.validate();
  const double abs_tol = R::number(d.at("tolerance").at("abs"));

  auto queries = R::reals(in.at("queries"));
  auto expected = R::complexes(d.at("expected").at("interpolated"));  // (query, order, 2, 2)
  auto flags = R::reals(d.at("expected").at("near_anomaly"));
  double m = 0;
  for (std::size_t q = 0; q < queries.size() / 3; ++q) {
    for (std::size_t o = 0; o < t.orders.size(); ++o) {
      Jones j = t.interpolate(t.orders[o], queries[3 * q], queries[3 * q + 1], queries[3 * q + 2]);
      for (std::size_t e = 0; e < 4; ++e) m = std::max(m, std::abs(j[e] - expected[(q * t.orders.size() + o) * 4 + e]));
    }
    const bool flag = t.near_rayleigh_anomaly(queries[3 * q], queries[3 * q + 1], queries[3 * q + 2]);
    check(flag == (flags[q] != 0.0), "anomaly flag as Python", flag, flags[q]);
  }
  check(m < abs_tol, "cubic interpolation vs Python", m, abs_tol);

  // PRT at query 0 from FourierGrating's direct Jones matrices (reference section)
  auto direct = R::complexes(d.at("reference").at("direct"));
  auto prt_ref = R::complexes(d.at("expected").at("prt_at_query_0"));
  const double lam = 0.6012, sx = 0.031, sy = 0.007;
  double mp = 0;
  const double kin[3] = {sx, sy, std::sqrt(1 - sx * sx - sy * sy)};
  double energy_max = 0;
  for (std::size_t o = 0; o < t.orders.size(); ++o) {
    Jones j{direct[o * 4], direct[o * 4 + 1], direct[o * 4 + 2], direct[o * 4 + 3]};
    Prt P = prt_matrix(j, {sx, sy}, {sx + t.orders[o] * lam / t.period_um, sy}, 1.0, 1.0, false);
    for (std::size_t e = 0; e < 9; ++e) mp = std::max(mp, std::abs(P[e] - prt_ref[o * 9 + e]));
    // P k_in = k_m
    const double so[2] = {sx + t.orders[o] * lam / t.period_um, sy};
    const double km[3] = {so[0], so[1], std::sqrt(1 - so[0] * so[0] - so[1] * so[1])};
    for (int r = 0; r < 3; ++r) {
      cplx v = 0;
      for (int c = 0; c < 3; ++c) v += P[static_cast<std::size_t>(3 * r + c)] * kin[c];
      energy_max = std::max(energy_max, std::abs(v - km[r]));
    }
  }
  check(mp < 1e-12, "PRT matrices vs Python", mp, 1e-12);
  check(energy_max < 1e-13, "P k_in = k_m", energy_max, 1e-13);

  // anomaly detection where order ±1 grazes (σ_x ≈ ±0.143 at λ = 0.6, as tests/test_jones_table.py)
  DiffractionTable a;
  a.orders = {0};
  for (int i = 0; i < 5; ++i) a.wavelengths_um.push_back(0.58 + 0.01 * i);
  for (int i = 0; i < 13; ++i) a.sigma_x.push_back(0.025 * i);
  for (int i = 0; i < 4; ++i) a.sigma_y.push_back(-0.06 + 0.04 * i);
  a.period_um = 0.7;
  const bool hit1 = a.near_rayleigh_anomaly(0.6, 0.140, 0.0), hit2 = a.near_rayleigh_anomaly(0.6, 0.150, 0.01);
  const bool miss1 = a.near_rayleigh_anomaly(0.6, 0.05, 0.0), miss2 = a.near_rayleigh_anomaly(0.6, 0.24, 0.0);
  check(hit1 && hit2 && !miss1 && !miss2, "anomaly cells as in Python (J5)", hit1 + hit2 - miss1 - miss2, 2);

  // table validation
  DiffractionTable bad = t;
  bad.sigma_x.pop_back();
  bool threw = false;
  try {
    bad.validate();
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  check(threw, "invalid table refused", threw, 1);
  std::printf("%d failure(s)\n", failures);
  return failures == 0 ? 0 : 1;
}
