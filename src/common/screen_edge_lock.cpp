#include "screen_edge_lock.hpp"
#include "ear/conversion.hpp"
#include "geom.hpp"

namespace ear {

  ScreenEdgeLockHandler::ScreenEdgeLockHandler(
      boost::optional<Screen> reproductionScreen, const Layout& layout)
      : _layout(layout) {
    if (reproductionScreen)
      _repScreenEdges = PolarEdges::fromScreen(*reproductionScreen);
  }

  bool ScreenEdgeLockHandler::shouldModifyPosition(
      const ScreenEdgeLock& screenEdgeLock) const {
    return _repScreenEdges &&
           (screenEdgeLock.horizontal || screenEdgeLock.vertical);
  }

  std::pair<double, double> ScreenEdgeLockHandler::lockToScreenEdge(
      double az, double el, const ScreenEdgeLock& screenEdgeLock) const {
    if (screenEdgeLock.horizontal) {
      if (*screenEdgeLock.horizontal == "left")
        az = _repScreenEdges->leftAzimuth;
      if (*screenEdgeLock.horizontal == "right")
        az = _repScreenEdges->rightAzimuth;
    }
    if (screenEdgeLock.vertical) {
      if (*screenEdgeLock.vertical == "top") el = _repScreenEdges->topElevation;
      if (*screenEdgeLock.vertical == "bottom")
        el = _repScreenEdges->bottomElevation;
    }
    return {az, el};
  }

  std::pair<double, double> ScreenEdgeLockHandler::handleAzimuthElevation(
      double azimuth, double elevation,
      const ScreenEdgeLock& screenEdgeLock) const {
    if (shouldModifyPosition(screenEdgeLock))
      return lockToScreenEdge(azimuth, elevation, screenEdgeLock);
    return {azimuth, elevation};
  }

  Eigen::Vector3d ScreenEdgeLockHandler::handle(
      const Eigen::Vector3d& pos, const ScreenEdgeLock& screenEdgeLock,
      bool cartesian) const {
    if (!shouldModifyPosition(screenEdgeLock)) return pos;

    if (cartesian) {
      PolarPosition polar = conversion::pointCartToPolar(
          CartesianPosition(pos(0), pos(1), pos(2)));
      std::pair<double, double> locked =
          lockToScreenEdge(polar.azimuth, polar.elevation, screenEdgeLock);
      std::pair<double, double> comp =
          compensatePosition(locked.first, locked.second, _layout);
      CartesianPosition c = conversion::pointPolarToCart(
          PolarPosition(comp.first, comp.second, polar.distance));
      return {c.X, c.Y, c.Z};
    } else {
      double az = azimuth(pos), el = elevation(pos), dist = pos.norm();
      std::pair<double, double> locked = lockToScreenEdge(az, el, screenEdgeLock);
      return cart(locked.first, locked.second, dist);
    }
  }

  std::tuple<double, double, double> ScreenEdgeLockHandler::handleVector(
      const Eigen::Vector3d& pos, const ScreenEdgeLock& screenEdgeLock,
      bool cartesian) const {
    Eigen::Vector3d ret = handle(pos, screenEdgeLock, cartesian);
    return std::make_tuple(ret(0), ret(1), ret(2));
  }

}  // namespace ear
