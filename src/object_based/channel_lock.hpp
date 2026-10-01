#pragma once
#include <Eigen/Core>
#include <vector>
#include "ear/layout.hpp"
#include "ear/metadata.hpp"

namespace ear {

  /** @brief Implementation of channel locking as a position transformation.
   *
   * See Rec. ITU-R BS.2127 section 7.3.6 and `ChannelLockHandlerBase` and
   * its subclasses in the reference implementation
   * (`ear.core.objectbased.gain_calc`).
   */
  class ChannelLockHandler {
   public:
    enum class Mode {
      /// real normalised loudspeaker positions are used, and the distance
      /// calculation is unweighted
      Egocentric,
      /// allocentric loudspeaker positions are used, and the distance
      /// calculation is weighted
      Allocentric,
    };

    /** @param layout layout to lock to; this should not contain LFE channels
     * @param mode egocentric (polar) or allocentric (Cartesian)
     * @param positions loudspeaker positions to use, one row per channel; if
     *   not given, they are derived from the layout (which must be known to
     *   allocentric::positionsForLayout in allocentric mode)
     */
    ChannelLockHandler(const Layout& layout, Mode mode = Mode::Egocentric);
    ChannelLockHandler(const Layout& layout, Mode mode,
                       Eigen::MatrixXd positions);

    /** @brief Apply channel lock to a position.
     *
     * @param position Cartesian source position
     * @param channelLock channel lock information
     * @param excluded if given, channels with excluded[i] set are not
     *   considered
     * @return the position of the selected loudspeaker, or \p position if
     *   channel lock is not enabled or no loudspeaker is close enough
     */
    Eigen::Vector3d handle(const Eigen::Vector3d& position,
                           const ChannelLock& channelLock,
                           const std::vector<bool>* excluded = nullptr) const;

   private:
    void initPriority(const Layout& layout);

    Mode _mode;
    /// loudspeaker positions, one row per channel
    Eigen::MatrixXd _positions;
    /// priority of each channel, used to break ties; lower is better
    Eigen::VectorXi _priority;
  };

}  // namespace ear
