#include "gain_calculator_objects.hpp"

#include <cmath>
#include <functional>
#include "../common/geom.hpp"
#include "../common/helpers/eigen_helpers.hpp"
#include "divergence.hpp"
#include "ear/exceptions.hpp"

namespace ear {

  // throws an exception if the given component is not implemented
  struct throw_if_not_implemented : public boost::static_visitor<void> {
    void operator()(const CartesianPosition&) const {
      throw not_implemented("cartesian");
    }
    template <typename T>
    void operator()(const T&) const {}
  };

  namespace {
    /// polar divergence parameters from the objectDivergence variant; a
    /// Cartesian divergence in polar mode uses the default azimuthRange,
    /// with a warning, as in the reference implementation
    struct PolarDivergenceParams : public boost::static_visitor<void> {
      PolarDivergenceParams(const WarningCB& warning_cb)
          : warning_cb(warning_cb) {}

      void operator()(const PolarObjectDivergence& divergence) {
        value = divergence.divergence;
        azimuthRange = divergence.azimuthRange;
      }
      void operator()(const CartesianObjectDivergence& divergence) {
        value = divergence.divergence;
        azimuthRange = 45.0;
        if (value != 0.0)
          warning_cb({Warning::Code::DIVERGENCE_POSITIONRANGE_IGNORED,
                      "positionRange specified for blockFormat in polar mode; "
                      "using polar divergence"});
      }

      const WarningCB& warning_cb;
      double value = 0.0;
      double azimuthRange = 45.0;
    };
  }  // namespace

  GainCalculatorObjectsImpl::GainCalculatorObjectsImpl(const Layout& layout)
      : _layout(layout),
        _pointSourcePanner(configurePolarPanner(_layout.withoutLfe())),
        _channelLockHandler(_layout.withoutLfe()),
        _polarExtentPanner(_pointSourcePanner),
        _isLfe(copy_vector<decltype(_isLfe)>(layout.isLfe())),
        _pvTmp(_pointSourcePanner->numberOfOutputChannels()),
        _pvTmpDiverged(_pointSourcePanner->numberOfOutputChannels()){};

  void GainCalculatorObjectsImpl::calculate(const ObjectsTypeMetadata& metadata,
                                            OutputGains& direct,
                                            OutputGains& diffuse,
                                            const WarningCB& warning_cb) {
    if (metadata.cartesian) throw not_implemented("cartesian");
    boost::apply_visitor(throw_if_not_implemented(), metadata.position);
    if (metadata.zoneExclusion.zones.size())
      throw not_implemented("zoneExclusion");
    if (metadata.screenRef) throw not_implemented("screenRef");

    Eigen::Vector3d position =
        toCartesianVector3d(boost::get<PolarPosition>(metadata.position));

    position = _channelLockHandler.handle(position, metadata.channelLock);

    PolarDivergenceParams divergenceParams(warning_cb);
    boost::apply_visitor(divergenceParams, metadata.objectDivergence);

    if (divergenceParams.value == 0.0) {
      _polarExtentPanner.handle(position, metadata.width, metadata.height,
                                metadata.depth, _pvTmp);
    } else {
      // pan each diverged source separately and combine with power
      // normalisation, as in the reference implementation
      DivergedPositions diverged = divergePolar(
          position, divergenceParams.value, divergenceParams.azimuthRange);
      _pvTmp.setZero();
      for (size_t i = 0; i < diverged.positions.size(); i++) {
        _polarExtentPanner.handle(diverged.positions[i], metadata.width,
                                  metadata.height, metadata.depth,
                                  _pvTmpDiverged);
        _pvTmp.array() += diverged.gains[i] * _pvTmpDiverged.array().square();
      }
      _pvTmp = _pvTmp.array().sqrt();
    }

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
