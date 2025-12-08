#ifndef UNITS_H
#define UNITS_H

#include "units_data.h"
#include <optional>
#include <string_view>

namespace combaero::units {

struct UnitInfo {
    std::string_view input;
    std::string_view output;
};

// Query units by function name (linear search, ~100 entries)
inline std::optional<UnitInfo> get_units(std::string_view name) {
    for (const auto& e : function_units) {
        if (e.name == name) {
            return UnitInfo{e.input, e.output};
        }
    }
    return std::nullopt;
}

// Convenience: get input units only (empty string if not found)
inline std::string_view input_units(std::string_view name) {
    if (auto u = get_units(name)) return u->input;
    return "";
}

// Convenience: get output units only (empty string if not found)
inline std::string_view output_units(std::string_view name) {
    if (auto u = get_units(name)) return u->output;
    return "";
}

// Check if a function is registered
inline bool has_units(std::string_view name) {
    return get_units(name).has_value();
}

// Iterate all entries (for validation, documentation generation, etc.)
inline constexpr const Entry* begin() { return function_units; }
inline constexpr const Entry* end() { return function_units + function_count; }

} // namespace combaero::units

#endif // UNITS_H
