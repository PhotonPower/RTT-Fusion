// Tabulated diffraction interaction, see cpp/include/rtt_fusion/diffraction_table.hpp.
#include "rtt_fusion/diffraction_table.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace rtt_fusion {

namespace {
using cplx = std::complex<double>;

// 4-point Lagrange window on a uniform axis, as `_cubic_weights` in Python.
void cubic(const std::vector<double>& axis, double x, std::array<std::size_t, 4>& idx,
           std::array<double, 4>& w) {
  const double h = axis[1] - axis[0];
  const auto last = static_cast<long long>(axis.size()) - 3;
  const long long i = std::clamp(static_cast<long long>(std::floor((x - axis[0]) / h)), 1LL, last);
  for (int j = 0; j < 4; ++j) idx[static_cast<std::size_t>(j)] = static_cast<std::size_t>(i - 1 + j);
  for (int j = 0; j < 4; ++j) {
    double p = 1.0;
    for (int k = 0; k < 4; ++k)
      if (k != j) p *= (x - axis[idx[static_cast<std::size_t>(k)]]) /
                       (axis[idx[static_cast<std::size_t>(j)]] - axis[idx[static_cast<std::size_t>(k)]]);
    w[static_cast<std::size_t>(j)] = p;
  }
}

bool uniform(const std::vector<double>& a) {
  if (a.size() < 4) return false;
  const double h = a[1] - a[0];
  if (!(h > 0)) return false;
  for (std::size_t i = 1; i < a.size(); ++i)
    if (std::abs((a[i] - a[i - 1]) - h) > 1e-9 * std::abs(h)) return false;
  return true;
}
}  // namespace

void DiffractionTable::validate() const {
  if (!uniform(wavelengths_um) || !uniform(sigma_x) || !uniform(sigma_y))
    throw std::invalid_argument("axes must be uniform with at least 4 points");
  if (std::set<int>(orders.begin(), orders.end()).size() != orders.size())
    throw std::invalid_argument("orders must be unique");
  if (jones.size() != orders.size() * wavelengths_um.size() * sigma_x.size() * sigma_y.size())
    throw std::invalid_argument("jones has the wrong size");
  if (!(period_um > 0)) throw std::invalid_argument("period must be positive");
}

Jones DiffractionTable::interpolate(int order, double lam, double sx, double sy) const {
  const auto it = std::find(orders.begin(), orders.end(), order);
  if (it == orders.end()) throw std::invalid_argument("order not in the table");
  const std::size_t o = static_cast<std::size_t>(it - orders.begin());
  std::array<std::size_t, 4> ia, ib, ic;
  std::array<double, 4> wa, wb, wc;
  cubic(wavelengths_um, lam, ia, wa);
  cubic(sigma_x, sx, ib, wb);
  cubic(sigma_y, sy, ic, wc);
  const std::size_t na = wavelengths_um.size(), nb = sigma_x.size(), nc = sigma_y.size();
  Jones out{};
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b)
      for (int c = 0; c < 4; ++c) {
        const double w = wa[static_cast<std::size_t>(a)] * wb[static_cast<std::size_t>(b)] *
                         wc[static_cast<std::size_t>(c)];
        const Jones& j = jones[((o * na + ia[static_cast<std::size_t>(a)]) * nb +
                                ib[static_cast<std::size_t>(b)]) * nc + ic[static_cast<std::size_t>(c)]];
        for (int e = 0; e < 4; ++e) out[static_cast<std::size_t>(e)] += w * j[static_cast<std::size_t>(e)];
      }
  return out;
}

bool DiffractionTable::near_rayleigh_anomaly(double lam, double sx, double sy, int max_order) const {
  std::array<std::size_t, 4> ia, ib, ic;
  std::array<double, 4> w;
  cubic(wavelengths_um, lam, ia, w);
  cubic(sigma_x, sx, ib, w);
  cubic(sigma_y, sy, ic, w);
  const double ls[2] = {wavelengths_um[ia[0]], wavelengths_um[ia[3]]};
  const double xs[2] = {sigma_x[ib[0]], sigma_x[ib[3]]};
  const double ys[2] = {sigma_y[ic[0]], sigma_y[ic[3]]};
  for (int m = -max_order; m <= max_order; ++m)
    for (double n : {n_in, n_out}) {
      double lo = HUGE_VAL, hi = -HUGE_VAL;
      for (double l : ls)
        for (double x : xs)
          for (double y : ys) {
            const double d = std::hypot(x + m * l / period_um, y) - n;
            lo = std::min(lo, d);
            hi = std::max(hi, d);
          }
      if (lo < 0 && 0 < hi) return true;
    }
  return false;
}

Prt prt_matrix(const Jones& j, std::array<double, 2> si, std::array<double, 2> so, double n_in,
               double n_out, bool reflected) {
  const double ci = std::sqrt(n_in * n_in - si[0] * si[0] - si[1] * si[1]);
  const double k_in[3] = {si[0] / n_in, si[1] / n_in, ci / n_in};
  const double cz = std::sqrt(n_out * n_out - so[0] * so[0] - so[1] * so[1]);
  const double sign = reflected ? -1.0 : 1.0;
  const double k_m[3] = {so[0] / n_out, so[1] / n_out, sign * cz / n_out};
  // L (3×2) · J (2×2) · Π (2×3) · (I − k kᵀ) + k_m k_inᵀ
  const double lift[3][2] = {{1, 0}, {0, 1}, {-sign * so[0] / cz, -sign * so[1] / cz}};
  cplx lj[3][2];
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 2; ++c) lj[r][c] = lift[r][0] * j[static_cast<std::size_t>(c)] + lift[r][1] * j[static_cast<std::size_t>(2 + c)];
  Prt p{};
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c) {
      cplx v = 0.0;
      for (int t = 0; t < 2; ++t) v += lj[r][t] * ((t == c ? 1.0 : 0.0) - k_in[t] * k_in[c]);
      p[static_cast<std::size_t>(3 * r + c)] = v + k_m[r] * k_in[c];
    }
  return p;
}

}  // namespace rtt_fusion
