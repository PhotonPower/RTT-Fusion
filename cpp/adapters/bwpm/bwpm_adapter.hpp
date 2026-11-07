// BWPM adapter in C++ (stage R11, Fusion ADR draft point 9, way 2): built outside Raytatouille
// against the public headers of rtt-fusion and BWPM's C++ core (`bwpm::Solver`, BWPM B4c).
// Mirrors src/rtt_fusion/adapters/bwpm.py → bwpm.adapter (representation G).
#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "rtt_fusion/adapter.hpp"

namespace rtt_fusion::adapters {

class BwpmAdapter final : public SolverAdapter {
 public:
  // Structure JSON: nx, ny, dx_mm, dy_mm, dz_mm, indices [[re, im], …], index_before [re, im],
  // z0_mm, iterations, threads. Arrays: labels (nz·ny·nx, row-major (z, y, x)) and an optional
  // absorbing mask (ny·nx) through the setters before `prepare`.
  void set_labels(std::vector<std::int32_t> labels, std::size_t nz);
  void set_mask(std::vector<double> mask);

  nlohmann::json describe() const override;
  std::unique_ptr<Handle> prepare(const nlohmann::json& structure,
                                  const std::vector<double>& wavelengths_um) override;
  PortFields apply(Handle& handle, const PortField& incident) override;
  nlohmann::json estimate(const Handle& handle) const override;

 private:
  std::vector<std::int32_t> labels_;
  std::size_t nz_ = 0;
  std::vector<double> mask_;
};

// Registers "BWPM" in the AdapterRegistry.
void register_bwpm_adapter();

}  // namespace rtt_fusion::adapters
