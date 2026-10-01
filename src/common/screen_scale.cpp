#include "screen_scale.hpp"
#include "ear/conversion.hpp"
#include "geom.hpp"
#include "helpers/eigen_helpers.hpp"

namespace ear {

  PolarScreenScaler::PolarScreenScaler(const Screen& referenceScreen,
                                       const Screen& reproductionScreen)
      : _refScreenEdges(PolarEdges::fromScreen(referenceScreen)),
        _repScreenEdges(PolarEdges::fromScreen(reproductionScreen)) {}

  std::pair<double, double> PolarScreenScaler::scaleAzEl(double az,
                                                         double el) const {
    Eigen::Vector4d azX(-180.0, _refScreenEdges.rightAzimuth,
                        _refScreenEdges.leftAzimuth, 180.0);
    Eigen::Vector4d azY(-180.0, _repScreenEdges.rightAzimuth,
                        _repScreenEdges.leftAzimuth, 180.0);
    Eigen::Vector4d elX(-90.0, _refScreenEdges.bottomElevation,
                        _refScreenEdges.topElevation, 90.0);
    Eigen::Vector4d elY(-90.0, _repScreenEdges.bottomElevation,
                        _repScreenEdges.topElevation, 90.0);

    return {interp(az, azX, azY), interp(el, elX, elY)};
  }

  Eigen::Vector3d PolarScreenScaler::scalePosition(
      const Eigen::Vector3d& position) const {
    double az = azimuth(position), el = elevation(position),
           dist = position.norm();
    std::pair<double, double> newAzEl = scaleAzEl(az, el);
    return cart(newAzEl.first, newAzEl.second, dist);
  }

  ScreenScaleHandler::ScreenScaleHandler(
      boost::optional<Screen> reproductionScreen, const Layout& layout)
      : _reproductionScreen(reproductionScreen), _layout(layout) {}

  Eigen::Vector3d ScreenScaleHandler::handle(const Eigen::Vector3d& position,
                                             bool screenRef,
                                             const Screen& referenceScreen,
                                             bool cartesian) const {
    if (!(screenRef && _reproductionScreen)) return position;

    PolarScreenScaler scaler(referenceScreen, *_reproductionScreen);

    if (cartesian) {
      PolarPosition polar = conversion::pointCartToPolar(
          CartesianPosition(position(0), position(1), position(2)));
      std::pair<double, double> scaled =
          scaler.scaleAzEl(polar.azimuth, polar.elevation);
      std::pair<double, double> comp =
          compensatePosition(scaled.first, scaled.second, _layout);
      CartesianPosition c = conversion::pointPolarToCart(
          PolarPosition(comp.first, comp.second, polar.distance));
      return {c.X, c.Y, c.Z};
    } else {
      return scaler.scalePosition(position);
    }
  }

}  // namespace ear
