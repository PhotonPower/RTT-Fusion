// Adapter interface of the Fusion in C++ (Fusion ADR draft point 5 and 9, way 2: adapters in
// C++ outside Raytatouille, built against these headers). Mirrors docs/schnittstelle.md.
#pragma once

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <variant>
#include <vector>

#include <nlohmann/json.hpp>

#include "rtt_fusion/ports.hpp"

namespace rtt_fusion {

// A field at a port in one of the representations (G or S).
using PortField = std::variant<GriddedField, PlaneWaves>;

// Result of `apply`: reflected field at the front port ("front", direction −1) and transmitted
// field at the back port ("back", direction +1); adapters may add further ports.
using PortFields = std::map<std::string, PortField>;

// Opaque prepared structure; adapters derive from it.
struct Handle {
  virtual ~Handle() = default;
  std::string content_hash;
  std::vector<double> wavelengths_um;
};

class SolverAdapter {
 public:
  virtual ~SolverAdapter() = default;
  // name, version, representation ("G"/"S"), ports, linear, reciprocal, deterministic,
  // conventions_of_solver, validity – as data (JSON).
  virtual nlohmann::json describe() const = 0;
  // Checks the structure (solver-specific JSON plus solver-specific arrays passed by the caller
  // through the adapter's own constructor or setters) and returns a handle; computes nothing.
  virtual std::unique_ptr<Handle> prepare(const nlohmann::json& structure,
                                          const std::vector<double>& wavelengths_um) = 0;
  // Field incident on the front port (direction +1) → fields at the ports.
  virtual PortFields apply(Handle& handle, const PortField& incident) = 0;
  // {"seconds": …, "memory_bytes": …}
  virtual nlohmann::json estimate(const Handle& handle) const = 0;
};

// Registry of adapters by name (the Fusion core looks adapters up here).
class AdapterRegistry {
 public:
  using Factory = std::function<std::unique_ptr<SolverAdapter>()>;
  static AdapterRegistry& instance() {
    static AdapterRegistry r;
    return r;
  }
  void add(const std::string& name, Factory f) { factories_[name] = std::move(f); }
  std::unique_ptr<SolverAdapter> make(const std::string& name) const {
    auto it = factories_.find(name);
    return it == factories_.end() ? nullptr : it->second();
  }
  std::vector<std::string> names() const {
    std::vector<std::string> out;
    for (const auto& [k, v] : factories_) out.push_back(k);
    return out;
  }

 private:
  std::map<std::string, Factory> factories_;
};

}  // namespace rtt_fusion
