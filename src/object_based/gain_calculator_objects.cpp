#include "gain_calculator_objects.hpp"

#include <cmath>
#include <functional>
#include "../common/allocentric.hpp"
#include "../common/geom.hpp"
#include "../common/helpers/eigen_helpers.hpp"
#include "allo_extent.hpp"
#include "divergence.hpp"
#include "ear/exceptions.hpp"

namespace ear {

  namespace {
    /// Get a Cartesian position vector given the ADM position; corresponds
    /// to `coord_trans` in the reference implementation
    struct CoordTrans : public boost::static_visitor<Eigen::Vector3d> {
      Eigen::Vector3d operator()(const PolarPosition& position) const {
        return toCartesianVector3d(position);
      }
      Eigen::Vector3d operator()(const CartesianPosition& position) const {
        return toCartesianVector3d(position).cwiseMax(-1.0).cwiseMin(1.0);
      }
    };

    /// divergence parameters from the objectDivergence variant, for the
    /// coordinate system of the block format; a mismatched type uses the
    /// default range, with a warning, as in the reference implementation
    struct DivergenceParams : public boost::static_visitor<void> {
      DivergenceParams(bool cartesian, const WarningCB& warning_cb)
          : cartesian(cartesian), warning_cb(warning_cb) {}

      void operator()(const PolarObjectDivergence& divergence) {
        value = divergence.divergence;
        if (cartesian) {
          range = 0.0;
          if (value != 0.0)
            warning_cb({Warning::Code::DIVERGENCE_AZIMUTHRANGE_IGNORED,
                        "azimuthRange specified for blockFormat in Cartesian "
                        "mode; using Cartesian divergence"});
        } else {
          range = divergence.azimuthRange;
        }
      }
      void operator()(const CartesianObjectDivergence& divergence) {
        value = divergence.divergence;
        if (cartesian) {
          range = divergence.positionRange;
        } else {
          range = 45.0;
          if (value != 0.0)
            warning_cb({Warning::Code::DIVERGENCE_POSITIONRANGE_IGNORED,
                        "positionRange specified for blockFormat in polar "
                        "mode; using polar divergence"});
        }
      }

      bool cartesian;
      const WarningCB& warning_cb;
      double value = 0.0;
      double range = 0.0;
    };
  }  // namespace

  GainCalculatorObjectsImpl::GainCalculatorObjectsImpl(const Layout& layout)
      : _layout(layout),
        _pointSourcePanner(configurePolarPanner(_layout.withoutLfe())),
        _screenEdgeLockHandler(_layout.screen(), _layout.withoutLfe()),
        _screenScaleHandler(_layout.screen(), _layout.withoutLfe()),
        _egoChannelLockHandler(_layout.withoutLfe(),
                               ChannelLockHandler::Mode::Egocentric),
        _polarExtentPanner(_pointSourcePanner),
        _zoneExclusionHandler(_layout.withoutLfe()),
        _isLfe(copy_vector<decltype(_isLfe)>(layout.isLfe())),
        _pvTmp(_pointSourcePanner->numberOfOutputChannels()),
        _pvTmpDiverged(_pointSourcePanner->numberOfOutputChannels()) {
    // Cartesian rendering needs allocentric loudspeaker positions, which are
    // only defined for the BS.2051 layouts; polar rendering of other layouts
    // is still possible
    _alloChannelPositions =
        allocentric::positionsForLayoutIfKnown(_layout.withoutLfe());
    if (_alloChannelPositions)
      _alloChannelLockHandler.emplace(_layout.withoutLfe(),
                                      ChannelLockHandler::Mode::Allocentric,
                                      *_alloChannelPositions);
  }

  void GainCalculatorObjectsImpl::alloExtentPan(
      const Eigen::Vector3d& position, double width, double height,
      double depth, const std::vector<bool>& excluded,
      Eigen::Ref<Eigen::VectorXd> out) {
    const Eigen::MatrixXd& allPositions = *_alloChannelPositions;

    Eigen::Index numIncluded = 0;
    for (bool ex : excluded)
      if (!ex) numIncluded++;

    Eigen::MatrixXd positions(numIncluded, 3);
    Eigen::Index j = 0;
    for (Eigen::Index i = 0; i < allPositions.rows(); i++)
      if (!excluded[i]) positions.row(j++) = allPositions.row(i);

    Eigen::VectorXd gains;
    if (width == 0.0 && height == 0.0 && depth == 0.0) {
      if (!_alloPanner || _alloPannerExcluded != excluded) {
        _alloPanner.reset(new AllocentricPanner(positions));
        _alloPannerExcluded = excluded;
      }
      gains = _alloPanner->handle(position).get();
    } else {
      gains = alloExtentGains(positions, position, width, height, depth);
    }

    out.setZero();
    j = 0;
    for (Eigen::Index i = 0; i < allPositions.rows(); i++)
      if (!excluded[i]) out(i) = gains(j++);
  }

  void GainCalculatorObjectsImpl::calculate(const ObjectsTypeMetadata& metadata,
                                            OutputGains& direct,
                                            OutputGains& diffuse,
                                            const WarningCB& warning_cb) {
    Eigen::Vector3d position =
        boost::apply_visitor(CoordTrans(), metadata.position);

    position = _screenScaleHandler.handle(position, metadata.screenRef,
                                          metadata.referenceScreen,
                                          metadata.cartesian);

    position = _screenEdgeLockHandler.handle(position, metadata.screenEdgeLock,
                                             metadata.cartesian);

    std::vector<bool> excluded;
    if (metadata.cartesian) {
      if (!_alloChannelPositions)
        throw invalid_argument(
            "Cartesian rendering is not possible for layout '" +
            _layout.name() + "': no allocentric loudspeaker positions");

      excluded = allocentric::getExcluded(
          *_alloChannelPositions,
          _zoneExclusionHandler.getExcluded(metadata.zoneExclusion));

      position = _alloChannelLockHandler->handle(position, metadata.channelLock,
                                                 &excluded);
    } else {
      position = _egoChannelLockHandler.handle(position, metadata.channelLock);
    }

    DivergenceParams divergenceParams(metadata.cartesian, warning_cb);
    boost::apply_visitor(divergenceParams, metadata.objectDivergence);

    DivergedPositions diverged =
        metadata.cartesian
            ? divergeCartesian(position, divergenceParams.value,
                               divergenceParams.range)
            : divergePolar(position, divergenceParams.value,
                           divergenceParams.range);

    // pan each diverged source separately and combine with power
    // normalisation, as in the reference implementation
    _pvTmp.setZero();
    for (size_t i = 0; i < diverged.positions.size(); i++) {
      if (metadata.cartesian)
        alloExtentPan(diverged.positions[i], metadata.width, metadata.height,
                      metadata.depth, excluded, _pvTmpDiverged);
      else
        _polarExtentPanner.handle(diverged.positions[i], metadata.width,
                                  metadata.height, metadata.depth,
                                  _pvTmpDiverged);
      _pvTmp.array() += diverged.gains[i] * _pvTmpDiverged.array().square();
    }
    _pvTmp = _pvTmp.array().sqrt();

    if (!metadata.cartesian)
      _pvTmp = _zoneExclusionHandler.handle(_pvTmp, metadata.zoneExclusion);

    // equivalent to np.nan_to_num in the reference implementation
    _pvTmp = _pvTmp.unaryExpr(
        [](double x) { return std::isnan(x) ? 0.0 : x; });

    _pvTmp *= metadata.gain;

    Eigen::VectorXd pv_full = Eigen::VectorXd::Zero(_isLfe.size());
    mask_write(pv_full, !_isLfe, _pvTmp);

    // apply diffuse split
    direct.write_vector(pv_full * std::sqrt(1.0 - metadata.diffuse));
    diffuse.write_vector(pv_full * std::sqrt(metadata.diffuse));
  }

}  // namespace ear
