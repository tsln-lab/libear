#include "channel_lock.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <tuple>
#include "../common/geom.hpp"

namespace ear {

  ChannelLockHandler::ChannelLockHandler(const Layout& layout) {
    const auto& channels = layout.channels();
    const Eigen::Index n = static_cast<Eigen::Index>(channels.size());

    _positions.resize(n, 3);
    std::vector<double> azimuths(n), elevations(n);
    for (Eigen::Index i = 0; i < n; i++) {
      PolarPosition position = channels[i].polarPosition();
      azimuths[i] = position.azimuth;
      elevations[i] = position.elevation;
      _positions.row(i) = toNormalisedVector3d(position).transpose();
    }

    // define a priority for channels, used to select a single channel when
    // multiple channels are the same distance from the position. Channels
    // with the lowest absolute elevation have the highest priority, with
    // ties broken by elevation, absolute azimuth then azimuth.
    std::vector<Eigen::Index> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](Eigen::Index a, Eigen::Index b) {
                       return std::make_tuple(std::abs(elevations[a]),
                                              elevations[a],
                                              std::abs(azimuths[a]),
                                              azimuths[a]) <
                              std::make_tuple(std::abs(elevations[b]),
                                              elevations[b],
                                              std::abs(azimuths[b]),
                                              azimuths[b]);
                     });

    _priority.resize(n);
    for (Eigen::Index i = 0; i < n; i++) {
      _priority(order[i]) = static_cast<int>(i);
    }
  }

  Eigen::Vector3d ChannelLockHandler::handle(
      const Eigen::Vector3d& position, const ChannelLock& channelLock) const {
    if (!channelLock.flag) return position;

    const double tol = 1e-5;
    const Eigen::Index n = _positions.rows();

    Eigen::VectorXd distances =
        (_positions.rowwise() - position.transpose()).rowwise().norm();

    // find possible channels closer than the maxDistance if given
    std::vector<Eigen::Index> possible;
    for (Eigen::Index i = 0; i < n; i++) {
      if (!channelLock.maxDistance ||
          distances(i) < *channelLock.maxDistance + tol)
        possible.push_back(i);
    }

    // if there are no possible channels, don't channel lock
    if (possible.empty()) return position;

    // find the minimum distance, and of the channels with the same distance,
    // return the position of the channel with the lowest priority
    double minDistance = distances(possible[0]);
    for (Eigen::Index i : possible) minDistance = std::min(minDistance, distances(i));

    Eigen::Index closest = -1;
    for (Eigen::Index i : possible) {
      if (distances(i) < minDistance + tol &&
          (closest < 0 || _priority(i) < _priority(closest)))
        closest = i;
    }

    return _positions.row(closest).transpose();
  }

}  // namespace ear
