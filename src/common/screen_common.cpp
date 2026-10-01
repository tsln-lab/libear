#include "screen_common.hpp"
#include <algorithm>
#include <cmath>
#include "ear/exceptions.hpp"
#include "geom.hpp"
#include "helpers/eigen_helpers.hpp"

namespace ear {

  namespace {
    struct ScreenGeometry {
      Eigen::Vector3d centre;
      Eigen::Vector3d xVec;
      Eigen::Vector3d zVec;
    };

    struct ScreenGeometryVisitor : public boost::static_visitor<ScreenGeometry> {
      ScreenGeometry operator()(const PolarScreen& screen) const {
        Eigen::Vector3d centre = toCartesianVector3d(screen.centrePosition);
        double width = screen.centrePosition.distance *
                       std::tan(radians(screen.widthAzimuth / 2.0));
        double height = width / screen.aspectRatio;

        Eigen::Matrix3d axes = localCoordinateSystem(
            screen.centrePosition.azimuth, screen.centrePosition.elevation);

        return {centre, axes.row(0).transpose() * width,
                axes.row(2).transpose() * height};
      }

      ScreenGeometry operator()(const CartesianScreen& screen) const {
        Eigen::Vector3d centre = toCartesianVector3d(screen.centrePosition);
        double width = screen.widthX / 2.0;
        double height = width / screen.aspectRatio;

        return {centre, Eigen::Vector3d(width, 0.0, 0.0),
                Eigen::Vector3d(0.0, 0.0, height)};
      }
    };
  }  // namespace

  PolarEdges PolarEdges::fromScreen(const Screen& screen) {
    // Determine the Cartesian position, angle and size of the screen, then
    // use that to determine the edge azimuths and elevations. The screen
    // surface is given by
    //
    //     centre + x * x_vec + y * y_vec
    //
    // for x and y in the range [-1, 1], where x=1, y=1 is the top right
    // corner, and x=1, y=-1 is the bottom right corner.
    ScreenGeometry g = boost::apply_visitor(ScreenGeometryVisitor(), screen);

    double leftAzimuth = azimuth(g.centre - g.xVec);
    double rightAzimuth = azimuth(g.centre + g.xVec);
    if (rightAzimuth > leftAzimuth)
      throw invalid_argument(
          "invalid screen specification: screen must not extend past -y");

    if ((azimuth(g.centre - g.zVec) - azimuth(g.centre + g.zVec)) > 1e-3)
      throw invalid_argument(
          "invalid screen specification: screen must not extend past +z or "
          "-z");

    return {leftAzimuth, rightAzimuth, elevation(g.centre - g.zVec),
            elevation(g.centre + g.zVec)};
  }

  std::pair<double, double> compensatePosition(double az, double el,
                                               const Layout& layout) {
    auto names = layout.channelNames();
    if (std::find(names.begin(), names.end(), "U+045") != names.end()) {
      Eigen::Vector3d elX(0.0, 30.0, 90.0);
      Eigen::Vector3d elY(30.0, 30.0 * (30.0 / 45.0), 30.0);
      double rightAz = interp(el, elX, elY);

      Eigen::Vector4d azX(-180.0, -30.0, 30.0, 180.0);
      Eigen::Vector4d azY(-180.0, -rightAz, rightAz, 180.0);
      double newAz = interp(az, azX, azY);

      return {newAz, el};
    }
    return {az, el};
  }

}  // namespace ear
