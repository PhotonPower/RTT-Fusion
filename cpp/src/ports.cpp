// Port types and conversions, see cpp/include/rtt_fusion/ports.hpp (mirrors ports.py, convert.py).
#include "rtt_fusion/ports.hpp"

#include <cmath>
#include <stdexcept>

#include "pocketfft_hdronly.h"

namespace rtt_fusion {

cplx Port::kz(double kx, double ky) const {
  const double k = k0();
  const cplx n2k2 = index * index * k * k;  // as Python: ((ñ²)·k₀)·k₀
  cplx q = n2k2 - (kx * kx + ky * ky);
  if (std::abs(q) <= kGrazingTolerance * std::abs(n2k2)) q = 0.0;
  cplx r = std::sqrt(q);
  if (r.imag() < 0.0 || (r.imag() == 0.0 && r.real() < 0.0)) r = -r;
  return r;
}

std::vector<cplx> PlaneWaves::ez() const {
  std::vector<cplx> out(size());
  for (std::size_t j = 0; j < size(); ++j) {
    const cplx kz = port.kz(kx[j], ky[j]);
    out[j] = kz == cplx(0.0, 0.0)
                 ? cplx(0.0, 0.0)
                 : -static_cast<double>(port.direction) * (kx[j] * field[j][0] + ky[j] * field[j][1]) / kz;
  }
  return out;
}

std::vector<double> PlaneWaves::power_density() const {
  const std::vector<cplx> e_z = ez();
  std::vector<double> out(size());
  for (std::size_t j = 0; j < size(); ++j) {
    const double mag = std::norm(field[j][0]) + std::norm(field[j][1]) + std::norm(e_z[j]);
    out[j] = port.kz(kx[j], ky[j]).real() * mag;
  }
  return out;
}

namespace {
// 2π·numpy.fft.fftfreq(n, d) with numpy's rounding: 2π·(m·(1/(n·d)))
double fftfreq(std::size_t i, std::size_t n, double d) {
  const auto m = static_cast<long long>(i) - (i < (n + 1) / 2 ? 0LL : static_cast<long long>(n));
  const double val = 1.0 / (static_cast<double>(n) * d);
  return 2.0 * kPi * (static_cast<double>(m) * val);
}
}  // namespace

PlaneWaves to_plane_waves(const GriddedField& g) {
  const std::size_t n = g.nx * g.ny;
  if (g.field.size() != n) throw std::invalid_argument("field must have nx·ny entries");
  std::vector<cplx> data(2 * n);
  for (std::size_t q = 0; q < n; ++q) {
    data[2 * q] = g.field[q][0];
    data[2 * q + 1] = g.field[q][1];
  }
  const auto c = static_cast<std::ptrdiff_t>(sizeof(cplx));
  const pocketfft::shape_t shape{g.ny, g.nx, 2};
  const pocketfft::stride_t stride{static_cast<std::ptrdiff_t>(2 * g.nx) * c, 2 * c, c};
  pocketfft::c2c<double>(shape, stride, stride, {0, 1}, true, data.data(), data.data(), 1.0);
  PlaneWaves s;
  s.port = g.port;
  s.kx.resize(n);
  s.ky.resize(n);
  s.field.resize(n);
  const double inv_n = static_cast<double>(g.nx * g.ny);
  for (std::size_t j = 0; j < g.ny; ++j)
    for (std::size_t i = 0; i < g.nx; ++i) {
      const std::size_t q = j * g.nx + i;
      const double kx = fftfreq(i, g.nx, g.dx_mm), ky = fftfreq(j, g.ny, g.dy_mm);
      // spectrum / N, referred to the grid centre: divide by exp(−i(k_x·(nx/2)·dx + k_y·(ny/2)·dy))
      const cplx centre = std::exp(cplx(0.0, -(kx * static_cast<double>(g.nx / 2) * g.dx_mm +
                                                ky * static_cast<double>(g.ny / 2) * g.dy_mm)));
      s.kx[q] = kx;
      s.ky[q] = ky;
      s.field[q] = {data[2 * q] / inv_n / centre, data[2 * q + 1] / inv_n / centre};
    }
  return s;
}

double power(const GriddedField& g) {
  const PlaneWaves s = to_plane_waves(g);
  double sum = 0.0;
  for (double p : s.power_density()) sum += p;
  return sum * static_cast<double>(g.nx * g.ny) * g.dx_mm * g.dy_mm;
}

double spectrum_power(const PlaneWaves& s, std::optional<double> cone_na) {
  if (s.solid_angle.size() != s.size()) throw std::invalid_argument("plane waves without dΩ");
  const double n = s.port.index.real();
  const double k = n * s.port.k0();
  const std::vector<cplx> e_z = s.ez();
  double sum = 0.0;
  for (std::size_t j = 0; j < s.size(); ++j) {
    if (cone_na && std::hypot(s.kx[j], s.ky[j]) / k > *cone_na / n) continue;
    sum += (std::norm(s.field[j][0]) + std::norm(s.field[j][1]) + std::norm(e_z[j])) /
           s.solid_angle[j];
  }
  return (2.0 * kPi) * (2.0 * kPi) / k * sum;
}

namespace {
// Gauss–Legendre nodes (ascending) and weights on [−1, 1], Newton on P_n.
void leggauss(int n, std::vector<double>& x, std::vector<double>& w) {
  x.assign(static_cast<std::size_t>(n), 0.0);
  w.assign(static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; ++i) {
    double z = std::cos(kPi * (i + 0.75) / (n + 0.5));
    double dp = 0.0;
    for (int it = 0; it < 100; ++it) {
      double p0 = 1.0, p1 = z;
      for (int k = 2; k <= n; ++k) {
        const double p2 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      dp = n * (z * p1 - p0) / (z * z - 1.0);
      const double dz = p1 / dp;
      z -= dz;
      if (std::abs(dz) < 1e-16) break;
    }
    x[static_cast<std::size_t>(n - 1 - i)] = z;  // ascending
    w[static_cast<std::size_t>(n - 1 - i)] = 2.0 / ((1.0 - z * z) * dp * dp);
  }
}
}  // namespace

void hemisphere(int n_theta, int n_phi, bool forward, std::vector<Vec3>& directions,
                std::vector<double>& weights) {
  std::vector<double> x, w;
  leggauss(n_theta, x, w);
  directions.clear();
  weights.clear();
  const double sign = forward ? 1.0 : -1.0;
  for (int a = 0; a < n_theta; ++a) {
    const double mu = 0.5 * (x[static_cast<std::size_t>(a)] + 1.0);
    const double wmu = 0.5 * w[static_cast<std::size_t>(a)];
    const double st = std::sqrt(1.0 - mu * mu);
    for (int b = 0; b < n_phi; ++b) {
      const double phi = 2.0 * kPi * (b + 0.5) / n_phi;
      directions.push_back({st * std::cos(phi), st * std::sin(phi), sign * mu});
      weights.push_back(wmu * (2.0 * kPi / n_phi));
    }
  }
}

PlaneWaves far_field_to_spectrum(const Port& port, const Vec3& center_mm,
                                 const std::vector<Vec3>& directions,
                                 const std::vector<double>& weights,
                                 const std::vector<Vec3c>& e_inf) {
  const std::size_t n = directions.size();
  if (weights.size() != n || e_inf.size() != n) throw std::invalid_argument("size mismatch");
  const double k = port.index.real() * port.k0();
  const Vec3 p0{-center_mm[0], -center_mm[1], port.z_mm - center_mm[2]};
  PlaneWaves s;
  s.port = port;
  for (std::size_t j = 0; j < n; ++j) {
    const Vec3& d = directions[j];
    if ((d[2] > 0 ? 1 : -1) != port.direction)
      throw std::invalid_argument("directions must point through the port");
    cplx radial = 0.0;
    for (int c = 0; c < 3; ++c) radial += e_inf[j][static_cast<std::size_t>(c)] * d[static_cast<std::size_t>(c)];
    const cplx factor = cplx(0.0, k / (2.0 * kPi)) * weights[j] *
                        std::exp(cplx(0.0, k * (d[0] * p0[0] + d[1] * p0[1] + d[2] * p0[2])));
    s.kx.push_back(d[0] * k);
    s.ky.push_back(d[1] * k);
    s.field.push_back({factor * (e_inf[j][0] - radial * d[0]), factor * (e_inf[j][1] - radial * d[1])});
    s.solid_angle.push_back(weights[j]);
  }
  return s;
}

PlaneWaves richards_wolf(const Port& port, double na, const Vec3& focus_mm, const Vec2c& pupil,
                         int n_theta, int n_phi) {
  const double n = port.index.real();
  if (!(na > 0.0 && na < n)) throw std::invalid_argument("0 < NA < n required");
  const double cmax = std::sqrt(1.0 - (na / n) * (na / n));
  std::vector<double> x, w;
  leggauss(n_theta, x, w);
  const double k = n * port.k0();
  const Vec3 shift{focus_mm[0], focus_mm[1], focus_mm[2] - port.z_mm};
  PlaneWaves s;
  s.port = port;
  for (int a = 0; a < n_theta; ++a) {
    const double ct = cmax + (1.0 - cmax) * 0.5 * (x[static_cast<std::size_t>(a)] + 1.0);
    const double wmu = (1.0 - cmax) * 0.5 * w[static_cast<std::size_t>(a)];
    const double st = std::sqrt(1.0 - ct * ct);
    for (int b = 0; b < n_phi; ++b) {
      const double phi = 2.0 * kPi * (b + 0.5) / n_phi;
      const double cp = std::cos(phi), sp = std::sin(phi);
      const double weight = wmu * (2.0 * kPi / n_phi);
      const cplx e_phi = -sp * pupil[0] + cp * pupil[1];
      const cplx e_rho = cp * pupil[0] + sp * pupil[1];
      const double sq = std::sqrt(ct);
      const cplx ax = sq * (e_phi * (-sp) + e_rho * (ct * cp));
      const cplx ay = sq * (e_phi * cp + e_rho * (ct * sp));
      const Vec3 d{st * cp, st * sp, ct};
      const cplx ph = std::exp(cplx(0.0, -k * (d[0] * shift[0] + d[1] * shift[1] + d[2] * shift[2])));
      s.kx.push_back(d[0] * k);
      s.ky.push_back(d[1] * k);
      s.field.push_back({ax * weight * ph, ay * weight * ph});
      s.solid_angle.push_back(weight);
    }
  }
  return s;
}

}  // namespace rtt_fusion
