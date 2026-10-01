#include "channel_lock.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <tuple>
#include "../common/allocentric.hpp"
#include "../common/geom.hpp"

namespace ear {

  ChannelLockHandler::ChannelLockHandler(const Layout& layout, Mode mode)
      : _mode(mode) {
    const auto& channels = layout.channels();
    const Eigen::Index n = static_cast<Eigen::Index>(channels.size());

    if (mode == Mode::Egocentric) {
      _positions.resize(n, 3);
      for (Eigen::Index i = 0; i < n; i++)
        _positions.row(i) =
            toNormalisedVector3d(channels[i].polarPosition()).transpose();
    } else {
      _positions = allocentric::positionsForLayout(layout);
    }

    initPriority(layout);
  }

  ChannelLockHandler::ChannelLockHandler(const Layout& layout, Mode mode,
                                         Eigen::MatrixXd positions)
      : _mode(mode), _positions(std::move(positions)) {
    initPriority(layout);
  }

  void ChannelLockHandler::initPriority(const Layout& layout) {
    const auto& channels = layout.channels();
    const Eigen::Index n = static_cast<Eigen::Index>(channels.size());

    std::vector<double> azimuths(n), elevations(n);
    for (Eigen::Index i = 0; i < n; i++) {
      PolarPosition position = channels[i].polarPosition();
      azimuths[i] = position.azimuth;
      elevations[i] = position.elevation;
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
      const Eigen::Vector3d& position, const ChannelLock& channelLock,
      const std::vector<bool>* excluded) const {
    if (!channelLock.flag) return position;

    const double tol = 1e-5;
    const Eigen::Index n = _positions.rows();

    // unweighted distances, used to find channels within maxDistance
    Eigen::VectorXd distances =
        (_positions.rowwise() - position.transpose()).rowwise().norm();

    // weighted distances, used to find the closest loudspeaker
    Eigen::VectorXd distancesWeighted;
    if (_mode == Mode::Egocentric) {
      distancesWeighted = distances;
    } else {
      Eigen::RowVector3d w(1.0 / 16.0, 4.0, 32.0);
      distancesWeighted =
          ((_positions.rowwise() - position.transpose()).array().square().rowwise() * w.array())
              .rowwise()
              .sum()
              .sqrt();
    }

    // find possible channels: not excluded, and closer than the maxDistance
    // if given
    std::vector<Eigen::Index> possible;
    for (Eigen::Index i = 0; i < n; i++) {
      if (excluded && (*excluded)[i]) continue;
      if (!channelLock.maxDistance ||
          distances(i) < *channelLock.maxDistance + tol)
        possible.push_back(i);
    }

    // if there are no possible channels, don't channel lock
    if (possible.empty()) return position;

    // find the minimum weighted distance, and of the channels with the same
    // distance, return the position of the channel with the lowest priority
    double minDistance = distancesWeighted(possible[0]);
    for (Eigen::Index i : possible)
      minDistance = std::min(minDistance, distancesWeighted(i));

    Eigen::Index closest = -1;
    for (Eigen::Index i : possible) {
      if (distancesWeighted(i) < minDistance + tol &&
          (closest < 0 || _priority(i) < _priority(closest)))
        closest = i;
    }

    return _positions.row(closest).transpose();
  }

}  // namespace ear
