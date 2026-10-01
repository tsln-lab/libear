#pragma once
#include <boost/optional.hpp>
#include <memory>
#include <vector>
#include "../common/point_source_panner.hpp"
#include "../common/screen_edge_lock.hpp"
#include "../common/screen_scale.hpp"
#include "channel_lock.hpp"
#include "ear/common_types.hpp"
#include "ear/helpers/output_gains.hpp"
#include "ear/layout.hpp"
#include "ear/metadata.hpp"
#include "ear/warnings.hpp"
#include "polar_extent.hpp"
#include "zone_exclusion.hpp"

namespace ear {

  class GainCalculatorObjectsImpl {
   public:
    GainCalculatorObjectsImpl(const Layout& layout);
    void calculate(const ObjectsTypeMetadata& metadata, OutputGains& direct,
                   OutputGains& diffuse,
                   const WarningCB& warning_cb = default_warning_cb);
    template <typename T>
    void calculate(const ObjectsTypeMetadata& metadata, std::vector<T>& direct,
                   std::vector<T>& diffuse,
                   const WarningCB& warning_cb = default_warning_cb) {
      OutputGainsT<T> direct_wrap(direct);
      OutputGainsT<T> diffuse_wrap(diffuse);

      calculate(metadata, direct_wrap, diffuse_wrap, warning_cb);
    }

   private:
    /// gains for a Cartesian source with extent, excluding some channels;
    /// the panner for the non-excluded channels is cached
    void alloExtentPan(const Eigen::Vector3d& position, double width,
                       double height, double depth,
                       const std::vector<bool>& excluded,
                       Eigen::Ref<Eigen::VectorXd> out);

    Layout _layout;
    std::shared_ptr<PointSourcePanner> _pointSourcePanner;
    ScreenEdgeLockHandler _screenEdgeLockHandler;
    ScreenScaleHandler _screenScaleHandler;
    ChannelLockHandler _egoChannelLockHandler;
    PolarExtent _polarExtentPanner;
    ZoneExclusionHandler _zoneExclusionHandler;
    Eigen::Array<bool, Eigen::Dynamic, 1> _isLfe;

    /// allocentric state; only available for layouts with known allocentric
    /// loudspeaker positions
    boost::optional<Eigen::MatrixXd> _alloChannelPositions;
    boost::optional<ChannelLockHandler> _alloChannelLockHandler;
    std::vector<bool> _alloPannerExcluded;
    std::unique_ptr<AllocentricPanner> _alloPanner;

    Eigen::VectorXd _pvTmp;
    Eigen::VectorXd _pvTmpDiverged;
  };

}  // namespace ear
