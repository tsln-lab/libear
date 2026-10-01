#pragma once
#include <Eigen/Core>
#include <utility>
#include "ear/layout.hpp"
#include "ear/screen.hpp"

namespace ear {

  /** @brief Internal screen representation for scaling polar coordinates.
   *
   * This stores the azimuths of the right and left edges, and the elevations
   * of the top and bottom edges. See `PolarEdges` in the reference
   * implementation (`ear.core.screen_common`).
   */
  struct PolarEdges {
    double leftAzimuth;
    double rightAzimuth;
    double bottomElevation;
    double topElevation;

    /// @throws invalid_argument for screens that extend past -y or +-z
    static PolarEdges fromScreen(const Screen& screen);
  };

  /** @brief Modify az and el so that vertical panning in allocentric
   * coordinates produces vertical source positions in the given layout.
   */
  std::pair<double, double> compensatePosition(double azimuth,
                                               double elevation,
                                               const Layout& layout);

}  // namespace ear
