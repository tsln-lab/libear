#pragma once
#include <Eigen/Core>
#include <boost/optional.hpp>
#include <vector>
#include "ear/layout.hpp"

namespace ear {
  /// Allocentric (Cartesian) loudspeaker positions, see Rec. ITU-R BS.2127
  /// section 7.3.9 and `ear.core.allocentric` in the reference
  /// implementation.
  namespace allocentric {

    /** @brief Allocentric position for a polar screen loudspeaker (M+SC or
     * M-SC); corresponds to `_screen_spk_position_to_cart`.
     */
    Eigen::Vector3d screenSpeakerPositionToCart(const PolarPosition& position);

    /** @brief Allocentric positions for the channels of a layout, one row per
     * channel; corresponds to `positions_for_layout`.
     *
     * @throws invalid_argument if the layout name or a channel name is not in
     * the table of allocentric positions
     */
    Eigen::MatrixXd positionsForLayout(const Layout& layout);

    /** @brief positionsForLayout, or none if the layout is not known; this
     * allows polar rendering of custom layouts.
     */
    boost::optional<Eigen::MatrixXd> positionsForLayoutIfKnown(
        const Layout& layout);

    /** @brief Extend a set of excluded channels so that the remaining layout
     * works with the allocentric panner; corresponds to `get_excluded`.
     */
    std::vector<bool> getExcluded(const Eigen::MatrixXd& channelPositions,
                                  std::vector<bool> isExcluded);

  }  // namespace allocentric
}  // namespace ear
