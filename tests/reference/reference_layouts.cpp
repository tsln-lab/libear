#include <string>
#include "ear/bs2051.hpp"
#include "ear/exceptions.hpp"
#include "objects_reference.hpp"

namespace ear {
  namespace reference {

    Layout layoutForSpec(const std::string& spec) {
      const std::string marker = ":SC=";
      auto pos = spec.find(marker);
      if (pos == std::string::npos) return getLayout(spec);

      Layout layout = getLayout(spec.substr(0, pos));
      double az = std::stod(spec.substr(pos + marker.size()));
      for (Channel& channel : layout.channels()) {
        if (channel.name() == "M+SC") channel.polarPosition({az, 0.0, 1.0});
        if (channel.name() == "M-SC") channel.polarPosition({-az, 0.0, 1.0});
      }
      return layout;
    }

  }  // namespace reference
}  // namespace ear
