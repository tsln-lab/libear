#include "divergence.hpp"
#include "../common/geom.hpp"

namespace ear {

  namespace {
    // Find gains g_l, g_c, g_r for the left, centre and right objects for
    // divergence value x such that:
    // - g_l + g_r + g_c = 1 for all x
    // - g_l = g_r = 0 and g_c = 1 for x = 0
    // - g_l = g_r = g_c = 1/3 for x = 0.5
    // - g_l = g_r = 0.5 and g_c = 0 for x = 1
    std::vector<double> divergenceGains(double divergence) {
      double g_lr = divergence / (divergence + 1.0);
      double g_c = (1.0 - divergence) / (divergence + 1.0);
      return {g_lr, g_c, g_lr};
    }
  }  // namespace

  DivergedPositions divergePolar(const Eigen::Vector3d& position,
                                 double divergence, double azimuthRange) {
    if (divergence == 0.0) return {{1.0}, {position}};

    double dist = position.norm();
    Eigen::Vector3d p_l = cart(azimuthRange, 0.0, dist);
    Eigen::Vector3d p_r = cart(-azimuthRange, 0.0, dist);

    Eigen::Matrix3d M =
        localCoordinateSystem(azimuth(position), elevation(position))
            .transpose();

    return {divergenceGains(divergence), {M * p_l, position, M * p_r}};
  }

  DivergedPositions divergeCartesian(const Eigen::Vector3d& position,
                                     double divergence, double positionRange) {
    if (divergence == 0.0) return {{1.0}, {position}};

    Eigen::Vector3d offset(positionRange, 0.0, 0.0);
    auto clip = [](const Eigen::Vector3d& p) -> Eigen::Vector3d {
      return p.cwiseMax(-1.0).cwiseMin(1.0);
    };

    return {divergenceGains(divergence),
            {clip(position + offset), clip(position), clip(position - offset)}};
  }

}  // namespace ear
