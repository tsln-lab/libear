#pragma once
#include <Eigen/Core>

namespace ear {

  /** @brief Gains for an object with extent in allocentric (Cartesian)
   * coordinates.
   *
   * See Rec. ITU-R BS.2127 section 7.3.10 and `get_gains` in the reference
   * implementation (`ear.core.objectbased.allo_extent`).
   *
   * @param channelPositions allocentric loudspeaker positions, one row per
   *   channel; this must not contain excluded channels
   * @param position allocentric source position
   * @param sizeX block format width
   * @param sizeY block format height
   * @param sizeZ block format depth
   * @return normalised gains, one per channel
   */
  Eigen::VectorXd alloExtentGains(const Eigen::MatrixXd& channelPositions,
                                  const Eigen::Vector3d& position,
                                  double sizeX, double sizeY, double sizeZ);

}  // namespace ear
