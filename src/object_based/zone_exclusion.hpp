#pragma once
#include <Eigen/Core>
#include <vector>
#include "ear/layout.hpp"
#include "ear/metadata.hpp"

namespace ear {

  /** @brief Calculate downmix coefficients to route output away from a given
   * set of loudspeakers.
   *
   * For each channel, this stores a list of groups of other channels, sorted
   * in priority order. The downmix matrix generated is such that the energy
   * from each excluded channel is distributed between the non-excluded
   * channels in the highest priority group containing at least one
   * non-excluded channel.
   *
   * See Rec. ITU-R BS.2127 section 7.3.12 and `ZoneExclusionDownmix` in the
   * reference implementation (`ear.core.objectbased.zone`).
   */
  class ZoneExclusionDownmix {
   public:
    /// @param layout layout without LFE channels
    explicit ZoneExclusionDownmix(const Layout& layout);

    /** @brief Calculate a downmix matrix for a given set of excluded
     * channels.
     *
     * @param excluded if excluded[i], channel i is excluded
     * @return Downmix matrix M; M(i, j) is the coefficient from channel i to
     *   channel j
     */
    Eigen::MatrixXd downmixForExcluded(const std::vector<bool>& excluded) const;

   private:
    int _numChannels;
    /// for each channel, groups of channels in priority order
    std::vector<std::vector<std::vector<int>>> _channelGroups;
  };

  /// Apply zone exclusion to loudspeaker gains; see `ZoneExclusionHandler` in
  /// the reference implementation.
  class ZoneExclusionHandler {
   public:
    /// @param layout layout without LFE channels
    explicit ZoneExclusionHandler(const Layout& layout);

    /// Channels excluded by the given zones
    std::vector<bool> getExcluded(const ZoneExclusion& zoneExclusion) const;

    /// Gains with the energy of excluded channels redistributed
    Eigen::VectorXd handle(const Eigen::VectorXd& gains,
                           const ZoneExclusion& zoneExclusion) const;

   private:
    int _numChannels;
    Eigen::MatrixXd _positions;
    Eigen::VectorXd _azimuths;
    Eigen::VectorXd _elevations;
    ZoneExclusionDownmix _zed;
  };

}  // namespace ear
