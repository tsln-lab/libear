#pragma once
#include <Eigen/Core>
#include <vector>

namespace ear {

  /// Source positions and weights produced by object divergence.
  struct DivergedPositions {
    /// weight of each position; these sum to 1
    std::vector<double> gains;
    /// Cartesian source positions
    std::vector<Eigen::Vector3d> positions;
  };

  /** @brief Implement polar object divergence by duplicating and modifying
   * source directions.
   *
   * See Rec. ITU-R BS.2127 section 7.3.7 and `diverge` in the reference
   * implementation (`ear.core.objectbased.gain_calc`).
   *
   * @param position Cartesian source position
   * @param divergence objectDivergence value, 0 to 1
   * @param azimuthRange azimuthRange of the objectDivergence in degrees
   * @return weights and positions for the left, centre and right sources
   *   (or just the original position if divergence is zero)
   */
  DivergedPositions divergePolar(const Eigen::Vector3d& position,
                                 double divergence, double azimuthRange);

}  // namespace ear
