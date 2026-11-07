// Reader of reference files in Raytatouille's result format (ADR 0023), for the C++ tests.
// Arrays: {"dtype", "shape", "values"}, complex values as [re, im]; non-finite as strings.
#pragma once

#include <complex>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace rtt_fusion::results {

using json = nlohmann::json;

inline json load(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("cannot open " + path);
  json j = json::parse(f);
  if (j.at("format") != "raytatouille-result") throw std::runtime_error("not a result file");
  return j.at("data");
}

inline double number(const json& v) {
  if (v.is_string()) {
    const std::string s = v.get<std::string>();
    if (s == "nan") return std::nan("");
    if (s == "inf") return HUGE_VAL;
    if (s == "-inf") return -HUGE_VAL;
    throw std::runtime_error("unknown number " + s);
  }
  if (v.is_boolean()) return v.get<bool>() ? 1.0 : 0.0;
  if (v.is_object() && v.contains("float")) return number(v.at("float"));
  return v.get<double>();
}

inline std::vector<std::size_t> shape(const json& a) { return a.at("shape").get<std::vector<std::size_t>>(); }

inline std::vector<double> reals(const json& a) {
  std::vector<double> out;
  for (const auto& v : a.at("values")) out.push_back(number(v));
  return out;
}

inline std::vector<std::complex<double>> complexes(const json& a) {
  std::vector<std::complex<double>> out;
  for (const auto& v : a.at("values")) out.emplace_back(number(v[0]), number(v[1]));
  return out;
}

inline double scalar(const json& v) {
  if (v.is_object() && v.contains("complex")) return number(v.at("complex")[0]);
  return number(v);
}

}  // namespace rtt_fusion::results
