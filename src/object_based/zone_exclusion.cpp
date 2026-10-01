#include "zone_exclusion.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include "../common/geom.hpp"
#include "ear/helpers/assert.hpp"

namespace ear {

  namespace {
    const double epsilon = 1e-6;

    /// sign function with some tolerance around zero
    int sign(double x) {
      if (x > epsilon) return 1;
      if (x < -epsilon) return -1;
      return 0;
    }

    /// layer number of a channel
    int layer(const Channel& channel) {
      double elevation = channel.polarPositionNominal().elevation;
      if (elevation < -10.0) return 0;
      if (elevation < 10.0) return 1;
      if (elevation < 75.0) return 2;
      return 3;
    }

    using Key = std::array<double, 4>;

    /// calculate a key for this channel; channels with a lower key
    /// (lexicographically) have a higher priority. If two channels have the
    /// same key, then the energy may be split between them
    Key calcKey(const Channel& from, const Channel& to) {
      // priority when moving between layers; prefer to move up before down
      static const int layerPrio[4][4] = {
          {0, 1, 2, 3},  // B
          {3, 0, 1, 2},  // M
          {3, 2, 0, 1},  // U
          {3, 2, 1, 0},  // T
      };

      Eigen::Vector3d fromPos = toCartesianVector3d(from.polarPositionNominal());
      Eigen::Vector3d toPos = toCartesianVector3d(to.polarPositionNominal());

      // prefer channels on the same layer (see layerPrio above)
      double layerPriority = layerPrio[layer(from)][layer(to)];

      // prefer to keep sources behind/in front of the listener; this results
      // in less extreme front/back movement when one side (left/right) is
      // excluded
      double frontBackChange = std::abs(sign(fromPos(1)) - sign(toPos(1)));

      // prefer closer speakers
      double cartDist = (fromPos - toPos).norm();

      // break ties by the front/back distance; this eliminates splitting that
      // is not symmetrical around +x or +y
      double frontBackDist = std::abs(fromPos(1) - toPos(1));

      return {layerPriority, frontBackChange, cartDist, frontBackDist};
    }

    bool keysEqual(const Key& a, const Key& b) {
      for (size_t i = 0; i < a.size(); i++)
        if (!(std::abs(a[i] - b[i]) < epsilon)) return false;
      return true;
    }
  }  // namespace

  ZoneExclusionDownmix::ZoneExclusionDownmix(const Layout& layout)
      : _numChannels(static_cast<int>(layout.channels().size())) {
    for (const Channel& channel : layout.channels())
      ear_assert(!channel.isLfe(), "lfe channel passed to zone exclusion panner");

    const auto& channels = layout.channels();
    for (int i = 0; i < _numChannels; i++) {
      // groups of channels with (approximately) the same key, in the order
      // they were first encountered
      std::vector<std::pair<Key, std::vector<int>>> groups;
      for (int j = 0; j < _numChannels; j++) {
        Key key = calcKey(channels[i], channels[j]);
        bool found = false;
        for (auto& group : groups) {
          if (keysEqual(key, group.first)) {
            group.second.push_back(j);
            found = true;
            break;
          }
        }
        if (!found) groups.emplace_back(key, std::vector<int>{j});
      }

      std::stable_sort(groups.begin(), groups.end(),
                       [](const std::pair<Key, std::vector<int>>& a,
                          const std::pair<Key, std::vector<int>>& b) {
                         return a.first < b.first;
                       });

      ear_assert(groups[0].second == std::vector<int>{i},
                 "channel should always be mapped to itself if possible");

      std::vector<std::vector<int>> channelGroups;
      for (auto& group : groups) channelGroups.push_back(std::move(group.second));
      _channelGroups.push_back(std::move(channelGroups));
    }
  }

  Eigen::MatrixXd ZoneExclusionDownmix::downmixForExcluded(
      const std::vector<bool>& excluded) const {
    ear_assert(static_cast<int>(excluded.size()) == _numChannels,
               "wrong number of channels in excluded");

    bool all = true, none = true;
    for (bool ex : excluded) {
      all = all && ex;
      none = none && !ex;
    }
    if (all || none) return Eigen::MatrixXd::Identity(_numChannels, _numChannels);

    Eigen::MatrixXd downmix = Eigen::MatrixXd::Zero(_numChannels, _numChannels);

    for (int i = 0; i < _numChannels; i++) {
      bool done = false;
      for (const auto& group : _channelGroups[i]) {
        std::vector<int> notExcluded;
        for (int j : group)
          if (!excluded[j]) notExcluded.push_back(j);
        if (!notExcluded.empty()) {
          for (int j : notExcluded) downmix(i, j) = 1.0 / notExcluded.size();
          done = true;
          break;
        }
      }
      ear_assert(done, "no non-excluded channels found");
    }

    return downmix;
  }

  ZoneExclusionHandler::ZoneExclusionHandler(const Layout& layout)
      : _numChannels(static_cast<int>(layout.channels().size())),
        _positions(layout.channels().size(), 3),
        _azimuths(layout.channels().size()),
        _elevations(layout.channels().size()),
        _zed(layout) {
    const auto& channels = layout.channels();
    for (int i = 0; i < _numChannels; i++) {
      PolarPosition nominal = channels[i].polarPositionNominal();
      _positions.row(i) = toCartesianVector3d(nominal).transpose();
      _azimuths(i) = nominal.azimuth;
      _elevations(i) = nominal.elevation;
    }
  }

  namespace {
    struct ExcludedVisitor : public boost::static_visitor<void> {
      ExcludedVisitor(const Eigen::MatrixXd& positions,
                      const Eigen::VectorXd& azimuths,
                      const Eigen::VectorXd& elevations,
                      std::vector<bool>& excluded)
          : positions(positions),
            azimuths(azimuths),
            elevations(elevations),
            excluded(excluded) {}

      void operator()(const CartesianExclusionZone& zone) const {
        for (Eigen::Index i = 0; i < positions.rows(); i++) {
          if (positions(i, 0) - epsilon < zone.maxX &&
              positions(i, 1) - epsilon < zone.maxY &&
              positions(i, 2) - epsilon < zone.maxZ &&
              positions(i, 0) + epsilon > zone.minX &&
              positions(i, 1) + epsilon > zone.minY &&
              positions(i, 2) + epsilon > zone.minZ)
            excluded[i] = true;
        }
      }

      void operator()(const PolarExclusionZone& zone) const {
        for (Eigen::Index i = 0; i < positions.rows(); i++) {
          if (elevations(i) - epsilon < zone.maxElevation &&
              elevations(i) + epsilon > zone.minElevation &&
              // speakers at the poles have indeterminate elevation and
              // should match any range
              (std::abs(elevations(i)) > 90.0 - epsilon ||
               insideAngleRange(azimuths(i), zone.minAzimuth, zone.maxAzimuth,
                                epsilon)))
            excluded[i] = true;
        }
      }

      const Eigen::MatrixXd& positions;
      const Eigen::VectorXd& azimuths;
      const Eigen::VectorXd& elevations;
      std::vector<bool>& excluded;
    };
  }  // namespace

  std::vector<bool> ZoneExclusionHandler::getExcluded(
      const ZoneExclusion& zoneExclusion) const {
    std::vector<bool> excluded(_numChannels, false);
    ExcludedVisitor visitor(_positions, _azimuths, _elevations, excluded);
    for (const auto& zone : zoneExclusion.zones)
      boost::apply_visitor(visitor, zone);
    return excluded;
  }

  Eigen::VectorXd ZoneExclusionHandler::handle(
      const Eigen::VectorXd& gains, const ZoneExclusion& zoneExclusion) const {
    std::vector<bool> excluded = getExcluded(zoneExclusion);
    Eigen::MatrixXd downmix = _zed.downmixForExcluded(excluded);
    Eigen::VectorXd power = gains.array().square();
    return (downmix.transpose() * power).array().sqrt();
  }

}  // namespace ear
