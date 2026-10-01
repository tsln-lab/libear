#pragma once
#include <Eigen/Core>
#include <boost/optional.hpp>
#include <tuple>
#include <utility>
#include "ear/layout.hpp"
#include "ear/metadata.hpp"
#include "ear/screen.hpp"
#include "screen_common.hpp"

namespace ear {

  /** @brief Lock positions to the edges of the reproduction screen.
   *
   * See Rec. ITU-R BS.2127 section 7.3.4 and `ScreenEdgeLockHandler` in the
   * reference implementation (`ear.core.screen_edge_lock`).
   */
  class ScreenEdgeLockHandler {
   public:
    ScreenEdgeLockHandler(boost::optional<Screen> reproductionScreen,
                          const Layout& layout);

    std::pair<double, double> handleAzimuthElevation(
        double azimuth, double elevation,
        const ScreenEdgeLock& screenEdgeLock) const;

    /** @param pos source position (polar space if !cartesian, allocentric
     *   space if cartesian)
     */
    std::tuple<double, double, double> handleVector(
        const Eigen::Vector3d& pos, const ScreenEdgeLock& screenEdgeLock,
        bool cartesian = false) const;

    Eigen::Vector3d handle(const Eigen::Vector3d& pos,
                           const ScreenEdgeLock& screenEdgeLock,
                           bool cartesian = false) const;

   private:
    bool shouldModifyPosition(const ScreenEdgeLock& screenEdgeLock) const;
    std::pair<double, double> lockToScreenEdge(
        double azimuth, double elevation,
        const ScreenEdgeLock& screenEdgeLock) const;

    boost::optional<PolarEdges> _repScreenEdges;
    Layout _layout;
  };

}  // namespace ear
