#pragma once
#include <Eigen/Core>
#include <boost/optional.hpp>
#include <utility>
#include "ear/layout.hpp"
#include "ear/screen.hpp"
#include "screen_common.hpp"

namespace ear {

  /** @brief Modifies the azimuth and elevation of a position independently,
   * preserving the distance:
   *
   * - azimuth is interpolated between the screen edges and +-180
   * - elevation is interpolated between the screen edges and +-90
   *
   * See Rec. ITU-R BS.2127 section 7.3.3 and `PolarScreenScaler` in the
   * reference implementation (`ear.core.screen_scale`).
   */
  class PolarScreenScaler {
   public:
    PolarScreenScaler(const Screen& referenceScreen,
                      const Screen& reproductionScreen);

    std::pair<double, double> scaleAzEl(double azimuth, double elevation) const;
    Eigen::Vector3d scalePosition(const Eigen::Vector3d& position) const;

   private:
    PolarEdges _refScreenEdges;
    PolarEdges _repScreenEdges;
  };

  /// Apply screen scaling to a position if screenRef is set.
  class ScreenScaleHandler {
   public:
    ScreenScaleHandler(boost::optional<Screen> reproductionScreen,
                       const Layout& layout);

    /** @param position source position (polar space if !cartesian,
     *   allocentric space if cartesian)
     * @param screenRef block format screenRef flag
     * @param referenceScreen audioProgrammeReferenceScreen
     * @param cartesian block format cartesian flag
     */
    Eigen::Vector3d handle(const Eigen::Vector3d& position, bool screenRef,
                           const Screen& referenceScreen,
                           bool cartesian) const;

   private:
    boost::optional<Screen> _reproductionScreen;
    Layout _layout;
  };

}  // namespace ear
