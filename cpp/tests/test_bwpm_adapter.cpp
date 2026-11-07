// R11: the BWPM adapter in C++ through the registry against BWPM's own reference file
// (bwpm3d_tilted: tilted glass interface and silver, absorbing mask, 2 iterations) and energy
// conservation through rtt-fusion's power for a lossless stack.
#include <cmath>
#include <cstdio>
#include <string>

#include "bwpm_adapter.hpp"
#include "rtt_fusion/results.hpp"

using namespace rtt_fusion;
namespace R = rtt_fusion::results;

int main(int argc, char** argv) {
  const std::string bwpm_reference = argc > 1 ? argv[1] : "reference";
  int failures = 0;
  adapters::register_bwpm_adapter();
  {
    auto d = R::load(bwpm_reference + "/bwpm3d_tilted.result.json");
    auto p = d.at("parameters");
    auto in = d.at("inputs");
    auto adapter_base = AdapterRegistry::instance().make("BWPM");
    auto* adapter = dynamic_cast<adapters::BwpmAdapter*>(adapter_base.get());
    auto shape = R::shape(in.at("labels"));
    std::vector<std::int32_t> labels;
    for (double v : R::reals(in.at("labels"))) labels.push_back(static_cast<std::int32_t>(v));
    adapter->set_labels(labels, shape[0]);
    adapter->set_mask(R::reals(in.at("mask")));
    nlohmann::json idx = nlohmann::json::array();
    for (auto c : R::complexes(in.at("indices"))) idx.push_back({c.real(), c.imag()});
    const double ib = R::scalar(in.at("index_before"));
    nlohmann::json s = {{"nx", p.at("nx")}, {"ny", p.at("ny")}, {"dx_mm", p.at("dx_mm")},
                        {"dy_mm", p.at("dy_mm")}, {"dz_mm", p.at("dz_mm")}, {"indices", idx},
                        {"index_before", {ib, 0.0}}, {"z0_mm", 0.0}, {"iterations", 2}};
    auto handle = adapter->prepare(s, {p.at("wavelength_um").get<double>()});
    GriddedField g{Port{0.0, 1, {ib, 0.0}, p.at("wavelength_um")}, p.at("nx"), p.at("ny"),
                   p.at("dx_mm"), p.at("dy_mm"), {}};
    auto f = R::complexes(in.at("field"));
    for (std::size_t q = 0; q < g.nx * g.ny; ++q) g.field.push_back({f[2 * q], f[2 * q + 1]});
    PortFields out = adapter->apply(*handle, g);
    auto tol = d.at("tolerance");
    const double abs_tol = R::number(tol.at("abs"));
    for (const char* side : {"front", "back"}) {
      const auto& got = std::get<GriddedField>(out.at(side)).field;
      auto ref = R::complexes(d.at("expected").at(std::string(side) == "front" ? "reflected" : "transmitted"));
      double m = 0;
      for (std::size_t q = 0; q < got.size(); ++q)
        m = std::max({m, std::abs(got[q][0] - ref[2 * q]), std::abs(got[q][1] - ref[2 * q + 1])});
      const bool ok = m < abs_tol;
      failures += ok ? 0 : 1;
      std::printf("%s  BWPM adapter %s vs bwpm3d_tilted (%.3e, limit %.3e)\n", ok ? "ok  " : "FAIL", side, m, abs_tol);
    }
    std::printf("     describe: %s\n", adapter->describe().at("name").get<std::string>().c_str());
  }
  {  // energy: lossless slab, plane wave on the grid, many iterations
    adapters::BwpmAdapter a;
    const std::size_t n = 32, nz = 24;
    std::vector<std::int32_t> labels(nz * n * n, 0);
    for (std::size_t z = 4; z < 18; ++z)
      for (std::size_t q = 0; q < n * n; ++q) labels[z * n * n + q] = 1;
    a.set_labels(labels, nz);
    const double lam = 1.0, step = lam * 1e-3 / 16;
    nlohmann::json s = {{"nx", n}, {"ny", n}, {"dx_mm", step}, {"dy_mm", step}, {"dz_mm", step},
                        {"indices", {{1.0, 0.0}, {1.5, 0.0}}}, {"index_before", {1.0, 0.0}},
                        {"z0_mm", 0.0}, {"iterations", 60}};
    auto h = a.prepare(s, {lam});
    GriddedField g{Port{0.0, 1, {1.0, 0.0}, lam}, n, n, step, step, {}};
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t i = 0; i < n; ++i) {
        const double x = (static_cast<double>(i) - n / 2) * step;
        const cplx ph = std::exp(cplx(0.0, 2 * kPi * 1 / (n * step) * x));
        g.field.push_back({0.4 * ph, cplx(0.0, 0.7) * ph});
      }
    PortFields out = a.apply(*h, g);
    const double balance = (power(std::get<GriddedField>(out.at("front"))) +
                            power(std::get<GriddedField>(out.at("back")))) / power(g);
    const bool ok = std::abs(balance - 1.0) < 1e-12;
    failures += ok ? 0 : 1;
    std::printf("%s  energy through rtt-fusion power (%.3e)\n", ok ? "ok  " : "FAIL", std::abs(balance - 1.0));
  }
  std::printf("%d failure(s)\n", failures);
  return failures == 0 ? 0 : 1;
}
