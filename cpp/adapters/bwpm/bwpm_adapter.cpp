// BWPM adapter, see bwpm_adapter.hpp.
#include "bwpm_adapter.hpp"

#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>

#include "bwpm/solver.hpp"

namespace rtt_fusion::adapters {

namespace {
struct BwpmHandle final : Handle {
  nlohmann::json structure;
  std::vector<std::int32_t> labels;
  std::size_t nz = 0;
  std::vector<double> mask;
  std::vector<cplx> indices;
  cplx index_before;
};

cplx complex_of(const nlohmann::json& v) { return {v.at(0).get<double>(), v.at(1).get<double>()}; }
}  // namespace

void BwpmAdapter::set_labels(std::vector<std::int32_t> labels, std::size_t nz) {
  labels_ = std::move(labels);
  nz_ = nz;
}

void BwpmAdapter::set_mask(std::vector<double> mask) { mask_ = std::move(mask); }

nlohmann::json BwpmAdapter::describe() const {
  return {{"name", "BWPM"},
          {"adapter", "rtt_fusion::adapters::BwpmAdapter (C++)"},
          {"representation", "G"},
          {"ports", {{{"name", "front"}, {"in", "+z"}, {"out", "-z"}}, {{"name", "back"}, {"out", "+z"}}}},
          {"linear", true},
          {"reciprocal", true},
          {"deterministic", {{"cpu", "bitwise"}}},
          {"conventions_of_solver", {"none (BWPM uses the conventions of Raytatouille)"}},
          {"validity", {"reflections quantitative only at z-normal interfaces (BWPM E5)",
                        "trapped waves do not converge; limit iterations"}}};
}

std::unique_ptr<Handle> BwpmAdapter::prepare(const nlohmann::json& s,
                                             const std::vector<double>& wavelengths_um) {
  const std::size_t nx = s.at("nx"), ny = s.at("ny");
  if (labels_.size() != nz_ * nx * ny || nz_ == 0) throw std::invalid_argument("labels (nz·ny·nx) missing");
  if (!mask_.empty() && mask_.size() != nx * ny) throw std::invalid_argument("mask must have ny·nx entries");
  if (wavelengths_um.empty()) throw std::invalid_argument("at least one wavelength is required");
  auto h = std::make_unique<BwpmHandle>();
  h->structure = s;
  h->labels = labels_;
  h->nz = nz_;
  h->mask = mask_;
  for (const auto& v : s.at("indices")) h->indices.push_back(complex_of(v));
  for (auto l : labels_)
    if (l < 0 || static_cast<std::size_t>(l) >= h->indices.size()) throw std::invalid_argument("label out of range");
  h->index_before = complex_of(s.at("index_before"));
  h->wavelengths_um = wavelengths_um;
  std::ostringstream key;
  key << "bwpm-cpp " << s.dump();
  std::size_t hash = std::hash<std::string>{}(key.str());
  for (auto l : labels_) hash = hash * 1099511628211ULL ^ static_cast<std::size_t>(l);
  for (double m : mask_) hash = hash * 1099511628211ULL ^ std::hash<double>{}(m);
  for (double w : wavelengths_um) hash = hash * 1099511628211ULL ^ std::hash<double>{}(w);
  h->content_hash = std::to_string(hash);
  return h;
}

PortFields BwpmAdapter::apply(Handle& handle, const PortField& incident) {
  auto& h = dynamic_cast<BwpmHandle&>(handle);
  const auto* g = std::get_if<GriddedField>(&incident);
  if (g == nullptr) throw std::invalid_argument("BWPM takes a gridded field (G)");
  const auto& s = h.structure;
  const std::size_t nx = s.at("nx"), ny = s.at("ny");
  const double dx = s.at("dx_mm"), dy = s.at("dy_mm"), dz = s.at("dz_mm"), z0 = s.value("z0_mm", 0.0);
  if (g->nx != nx || g->ny != ny || g->dx_mm != dx || g->dy_mm != dy)
    throw std::invalid_argument("incident grid differs from the structure grid");
  if (g->port.direction != 1 || g->port.z_mm != z0 || g->port.index != h.index_before)
    throw std::invalid_argument("incident field must be on the front port with direction +1");
  bool known = false;
  for (double w : h.wavelengths_um) known = known || w == g->port.wavelength_um;
  if (!known) throw std::invalid_argument("wavelength was not prepared");
  const double k0 = g->port.k0();
  bwpm::Solver solver({ny, nx, dy, dx, dz, k0}, h.indices, h.labels, h.nz, h.index_before, h.mask,
                      s.value("threads", 1));
  std::vector<cplx> in(2 * nx * ny);
  for (std::size_t q = 0; q < nx * ny; ++q) {
    in[2 * q] = g->field[q][0];
    in[2 * q + 1] = g->field[q][1];
  }
  const bwpm::RunResult r = solver.run(in, s.value("iterations", 3), {});
  auto to_grid = [&](const std::vector<cplx>& v, Port port) {
    GriddedField out{port, nx, ny, dx, dy, {}};
    for (std::size_t q = 0; q < nx * ny; ++q) out.field.push_back({v[2 * q], v[2 * q + 1]});
    return out;
  };
  const std::int32_t last = h.labels[(h.nz - 1) * nx * ny + (ny / 2) * nx + nx / 2];
  PortFields out;
  out["front"] = to_grid(r.reflected, Port{z0, -1, h.index_before, g->port.wavelength_um});
  out["back"] = to_grid(r.transmitted, Port{z0 + static_cast<double>(h.nz) * dz, 1,
                                            h.indices[static_cast<std::size_t>(last)], g->port.wavelength_um});
  return out;
}

nlohmann::json BwpmAdapter::estimate(const Handle& handle) const {
  const auto& h = dynamic_cast<const BwpmHandle&>(handle);
  const double points = static_cast<double>(h.nz) * h.structure.at("nx").get<double>() *
                        h.structure.at("ny").get<double>();
  const double passes = 2.0 * h.structure.value("iterations", 3);
  return {{"seconds", 1.5e-7 * points * passes}, {"memory_bytes", 4.0 * points + 300.0 * points / h.nz}};
}

void register_bwpm_adapter() {
  AdapterRegistry::instance().add("BWPM", [] { return std::make_unique<BwpmAdapter>(); });
}

}  // namespace rtt_fusion::adapters
