#pragma once
#include <Eigen/Core>
#include "ear/layout.hpp"
#include "ear/metadata.hpp"

namespace ear {

  /** @brief Implementation of channel locking as a position transformation.
   *
   * Egocentric (polar) variant: the real, normalised loudspeaker positions
   * are used, and the distance calculation is unweighted. See Rec. ITU-R
   * BS.2127 section 7.3.6 and `EgoChannelLockHandler` in the reference
   * implementation (`ear.core.objectbased.gain_calc`).
   */
  class ChannelLockHandler {
   public:
    /// @param layout layout to lock to; this should not contain LFE channels
    explicit ChannelLockHandler(const Layout& layout);

    /** @brief Apply channel lock to a position.
     *
     * @param position Cartesian source position
     * @param channelLock channel lock information
     * @return the position of the selected loudspeaker, or \p position if
     *   channel lock is not enabled or no loudspeaker is close enough
     */
    Eigen::Vector3d handle(const Eigen::Vector3d& position,
                           const ChannelLock& channelLock) const;

   private:
    /// normalised loudspeaker positions, one row per channel
    Eigen::MatrixXd _positions;
    /// priority of each channel, used to break ties; lower is better
    Eigen::VectorXi _priority;
  };

}  // namespace ear
