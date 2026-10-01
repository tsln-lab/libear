#include "point_source_panner.hpp"
#include "ear/exceptions.hpp"
#include <boost/math/constants/constants.hpp>

#include <Eigen/Dense>
#include <boost/algorithm/clamp.hpp>
#include <boost/make_unique.hpp>
#include "allocentric.hpp"
#include "convex_hull.hpp"
#include "ear/bs2051.hpp"
#include "ear/helpers/assert.hpp"
#include "facets.hpp"
#include "geom.hpp"
#include "helpers/eigen_helpers.hpp"

namespace ear {

  RegionHandler::RegionHandler(Eigen::VectorXi outputChannels,
                               Eigen::MatrixXd positions)
      : _outputChannels(outputChannels), _positions(positions){};

  boost::optional<Eigen::VectorXd> RegionHandler::handleRemap(
      Eigen::Vector3d position, int numberOfChannels) const {
    boost::optional<Eigen::VectorXd> pv = handle(position);

    if (pv) {
      Eigen::VectorXd _pv = pv.get();
      Eigen::VectorXd out = Eigen::VectorXd::Zero(numberOfChannels);
      for (int i = 0; i < _pv.size(); ++i) {
        out(_outputChannels(i)) = _pv(i);
      }
      return out;
    }
    return boost::none;
  };

  Eigen::VectorXi RegionHandler::outputChannels() { return _outputChannels; }

  Triplet::Triplet(Eigen::Vector3i outputChannels, Eigen::Matrix3d positions)
      : RegionHandler(outputChannels, positions) {
    _basis = _positions.inverse();
  };

  boost::optional<Eigen::VectorXd> Triplet::handle(
      Eigen::Vector3d position) const {
    Eigen::VectorXd pv = position.transpose() * _basis;
    double epsilon = -1e-11;
    if (pv(0) >= epsilon && pv(1) >= epsilon && pv(2) >= epsilon) {
      pv /= pv.norm();
      return pv.cwiseMax(0.0).cwiseMin(1.0).eval();
    }
    return boost::none;
  }

  VirtualNgon::VirtualNgon(Eigen::VectorXi outputChannels,
                           Eigen::MatrixXd positions,
                           Eigen::Vector3d centrePosition,
                           Eigen::VectorXd centreDownmix)
      : RegionHandler(outputChannels, positions),
        _centrePosition(centrePosition),
        _centreDownmix(centreDownmix) {
    int n = static_cast<int>(_outputChannels.size());

    ear_assert(
        n == _positions.rows(),
        "number of downmix coeffs does not match number of output channels");

    ear_assert(
        n == _centreDownmix.size(),
        "number of downmix coeffs does not match number of output channels");

    Eigen::VectorXi order = ngonVertexOrder(positions);

    for (int i = 0; i < n; ++i) {
      int j = (i + 1) % n;
      Eigen::Matrix3d tripletPositions;
      Eigen::RowVector3d position1 = _positions.row(order(i));
      Eigen::RowVector3d position2 = _positions.row(order(j));
      tripletPositions << position1, position2, _centrePosition.transpose();
      Eigen::Vector3i tripletChannels;
      tripletChannels << order(i), order(j), n;
      _regions.push_back(
          boost::make_unique<Triplet>(tripletChannels, tripletPositions));
    }
  }

  boost::optional<Eigen::VectorXd> VirtualNgon::handle(
      Eigen::Vector3d position) const {
    for (const auto& region : _regions) {
      boost::optional<Eigen::VectorXd> pv = region->handleRemap(
          position, static_cast<int>(_centreDownmix.size() + 1));
      if (pv) {
        // downmix the last channel containing the virtual centre
        // speaker into the real speakers, and renormalise
        Eigen::VectorXd _pv = pv.get();
        _pv = _pv.head(_pv.size() - 1) + _centreDownmix * _pv.tail(1);
        _pv /= _pv.norm();
        return _pv;
      }
    }
    return boost::none;
  }

  QuadRegion::QuadRegion(Eigen::VectorXi outputChannels,
                         Eigen::MatrixXd positions)
      : RegionHandler(outputChannels, positions) {
    _order = ngonVertexOrder(positions);
    Eigen::MatrixXd reorderedPositions =
        positions(_order, Eigen::indexing::all);
    Eigen::MatrixXd reorderedAndShiftedPositions =
        reorderedPositions(Eigen::Vector4i{1, 2, 3, 0}, Eigen::indexing::all);
    _polyBasisX = _calcPolyBasis(reorderedPositions);
    _polyBasisY = _calcPolyBasis(reorderedAndShiftedPositions);
  };

  boost::optional<Eigen::VectorXd> QuadRegion::handle(
      Eigen::Vector3d position) const {
    boost::optional<double> x = _pan(position, _polyBasisX);
    boost::optional<double> y = _pan(position, _polyBasisY);

    if (x == boost::none || y == boost::none) {
      return boost::none;
    }

    Eigen::VectorXd pvs = Eigen::Vector4d::Zero();
    pvs(_order) << (1 - x.get()) * (1 - y.get()), x.get() * (1 - y.get()),
        x.get() * y.get(), (1 - x.get()) * y.get();
    if ((pvs.transpose() * _positions) * position <= 0) {
      return boost::none;
    }
    pvs /= pvs.norm();
    return pvs;
  }

  Eigen::Matrix3d QuadRegion::_calcPolyBasis(Eigen::MatrixXd positions) {
    Eigen::Vector3d a = positions.row(0);
    Eigen::Vector3d b = positions.row(1);
    Eigen::Vector3d c = positions.row(2);
    Eigen::Vector3d d = positions.row(3);

    Eigen::Matrix3d polyBasis;
    polyBasis <<  //
        (b - a).cross(c - d),  //
        a.cross(c - d) + (b - a).cross(d),  //
        a.cross(d);

    return polyBasis.transpose();
  }

  /** @brief Calculate the real roots of a quadratic
   *
   * @param a Quadratic term.
   * @param b Linear term.
   * @param c Constant term.
   *
   * @note The equation solved is given by
   *
   *  a x^2 + bx + c = 0
   *
   * @returns The real roots of the quadratic.
   */
  std::vector<double> real_quadratic_roots(double a, double b, double c) {
    double eps = 1e-10;

    if (std::abs(c) < eps) return {0.0};
    if (std::abs(a) < eps) return {-c / b};

    double det = b * b - 4.0 * a * c;
    if (det > eps)
      return {(-b + sqrt(det)) / (2.0 * a), (-b - sqrt(det)) / (2.0 * a)};
    else if (det > -eps)
      return {-b / (2.0 * a)};
    else
      return {};
  }

  boost::optional<double> QuadRegion::_pan(Eigen::Vector3d position,
                                           Eigen::Matrix3d polyBasis) const {
    double epsilon = 1e-10;

    Eigen::Vector3d poly = polyBasis * position;

    for (double root : real_quadratic_roots(poly(0), poly(1), poly(2))) {
      if (-epsilon < root && root < 1.0 + epsilon) {
        return boost::algorithm::clamp(root, 0.0, 1.0);
      }
    }
    return boost::none;
  }

  PolarPointSourcePanner::PolarPointSourcePanner(
      std::vector<std::unique_ptr<RegionHandler>> regions,
      boost::optional<int> numberOfChannels)
      : _regions(std::move(regions)) {
    if (!numberOfChannels) {
      _numberOfOutputChannels = _numberOfRequiredChannels();
    } else {
      _numberOfOutputChannels = numberOfChannels.get();
      ear_assert(_numberOfOutputChannels >= _numberOfRequiredChannels(),
                 "not enough output channels in PolarPointSourcePanner");
    }
  };

  boost::optional<Eigen::VectorXd> PolarPointSourcePanner::handle(
      Eigen::Vector3d position) {
    boost::optional<Eigen::VectorXd> pv;
    for (const auto& region : _regions) {
      pv = region->handleRemap(position, _numberOfOutputChannels);
      if (pv) {
        return pv;
      }
    }
    return boost::none;
  }

  int PolarPointSourcePanner::numberOfOutputChannels() const {
    return _numberOfOutputChannels;
  }

  int PolarPointSourcePanner::_numberOfRequiredChannels() {
    int ret = 0;
    for (const auto& region : _regions) {
      auto maxRegion = region->outputChannels().maxCoeff();
      if (ret < maxRegion) {
        ret = maxRegion;
      }
    }
    return ret + 1;
  }

  PointSourcePanner::~PointSourcePanner() {}

  PointSourcePannerDownmix::PointSourcePannerDownmix(
      std::shared_ptr<PointSourcePanner> psp, Eigen::MatrixXd downmix)
      : _psp(psp), _downmix(downmix){};

  boost::optional<Eigen::VectorXd> PointSourcePannerDownmix::handle(
      Eigen::Vector3d position) {
    boost::optional<Eigen::VectorXd> pv = _psp->handle(position);
    if (pv) {
      Eigen::VectorXd _pv = pv.get();
      _pv = _downmix.transpose() * _pv;
      _pv /= _pv.norm();
      return _pv;
    }
    return boost::none;
  }

  int PointSourcePannerDownmix::numberOfOutputChannels() const {
    return static_cast<int>(_downmix.cols());
  }

  /** @brief Generate extra loudspeaker positions to fill gaps in layers.
   *
   * @param  layout Original layout without the LFE channels
   *
   * @returns
   *   - list of extra channels (layout.Channel).
   *   - downmix matrix to mix the extra channel outputs to the real channels
   */
  std::pair<std::vector<Channel>, Eigen::MatrixXd> extraPosVerticalNominal(
      Layout layout) {
    std::vector<Channel> extraChannels;
    Eigen::MatrixXd downmix = Eigen::MatrixXd::Identity(
        layout.channels().size(), layout.channels().size());

    Layout midLayerLayout;
    std::copy_if(layout.channels().begin(), layout.channels().end(),
                 std::back_inserter(midLayerLayout.channels()), [](Channel c) {
                   return -10 <= c.polarPositionNominal().elevation &&
                          c.polarPositionNominal().elevation <= 10;
                 });

    auto layers = {std::make_tuple(-30.0, -70.0, -10.0),
                   std::make_tuple(30.0, 10.0, 70.0)};

    double layerNominalElevation, layerLowerBound, layerUpperBound;
    for (const auto& layer : layers) {
      std::tie(layerNominalElevation, layerLowerBound, layerUpperBound) = layer;

      Layout currentLayerLayout;
      std::copy_if(
          layout.channels().begin(), layout.channels().end(),
          std::back_inserter(currentLayerLayout.channels()), [&](Channel c) {
            return layerLowerBound <= c.polarPositionNominal().elevation &&
                   c.polarPositionNominal().elevation <= layerUpperBound;
          });

      // for each loudspeaker in the mid layer that has an azimuth greater
      // than az_limit, add a virtual speaker directly above/below it at the
      // elevation of the current layer, which is downmixed directly to the
      // mid layer loudspeaker. az_limit is set to the range of azimuths in
      // the current layer, with some space added to prevent fast vertical
      // source movements when sources move horizontally. If there are no
      // channels on this layer then a copy of all mid layer speakers is
      // made.
      double azimuthLimit = 0.0;
      double layerRealElevation = 0.0;
      if (currentLayerLayout.channels().size() != 0) {
        double azimuthRange = std::numeric_limits<double>::min();
        for (const auto& channel : currentLayerLayout.channels()) {
          if (azimuthRange < std::abs(channel.polarPositionNominal().azimuth)) {
            azimuthRange = std::abs(channel.polarPositionNominal().azimuth);
          }
        }
        azimuthLimit = azimuthRange + 40.0;
        layerRealElevation =
            std::accumulate(currentLayerLayout.channels().begin(),
                            currentLayerLayout.channels().end(), 0.0,
                            [&](double sum, const Channel& c) -> double {
                              return sum + c.polarPosition().elevation;
                            }) /
            static_cast<double>(currentLayerLayout.channels().size());
      } else {
        layerRealElevation = layerNominalElevation;
      }

      double epsilon = 1e-5;
      for (const auto& midChannel : midLayerLayout.channels()) {
        if (std::abs(midChannel.polarPosition().azimuth) >=
            azimuthLimit - epsilon) {
          extraChannels.push_back(Channel(
              "extra",
              PolarPosition(midChannel.polarPosition().azimuth,
                            layerRealElevation, 1.0),
              PolarPosition(midChannel.polarPositionNominal().azimuth,
                            layerNominalElevation, 1.0)

                  ));
          Eigen::VectorXd downmixRow =
              Eigen::VectorXd::Zero(layout.channels().size());
          auto names = layout.channelNames();
          int midChannelIndex = static_cast<int>(std::distance(
              names.begin(),
              std::find(names.begin(), names.end(), midChannel.name())));
          downmixRow(midChannelIndex) = 1.0;
          downmix.conservativeResize(downmix.rows() + 1, Eigen::NoChange);
          downmix.row(downmix.rows() - 1) = downmixRow;
        }
      }
    }
    return std::make_pair(extraChannels, downmix);
  }

  /** @brief Find the adjacent vertices in a hull to the given vertex.
   *
   * @param  facets (list of sets of ints): Convex hull facets, each item
   * represents a facet, with the contents of the set being its vertex indices.
   * @param  vert (int): Vertex index to find vertices adjacent to.
   *
   * @returns Vertices adjacent to `vert`.
   */
  Facet _adjacent_verts(std::vector<Facet> facets, int vert) {
    std::set<int> ret;
    for (const auto& facetVerts : facets) {
      if (std::find(facetVerts.begin(), facetVerts.end(), vert) !=
          facetVerts.end()) {
        ret.insert(facetVerts.begin(), facetVerts.end());
      }
    }
    ret.erase(vert);
    return ret;
  }

  StereoPannerDownmix::StereoPannerDownmix(Eigen::VectorXi outputChannels,
                                           Eigen::MatrixXd positions)
      : RegionHandler(outputChannels, positions) {
    auto layout = getLayout("0+5+0").withoutLfe();
    _psp = configureFullPolarPanner(layout);
  }

  boost::optional<Eigen::VectorXd> StereoPannerDownmix::handle(
      Eigen::Vector3d position) const {
    Eigen::MatrixXd downmix(2, 5);
    downmix << 1.0, 0.0, std::sqrt(3.0) / 3.0, std::sqrt(0.5), 0.0,  //
        0.0, 1.0, std::sqrt(3.0) / 3.0, 0.0, std::sqrt(0.5);

    // pan with 0+5+0, downmix and power normalise
    boost::optional<Eigen::VectorXd> pv = _psp->handle(position);

    if (pv) {
      Eigen::VectorXd _pv = pv.get();
      Eigen::VectorXd pvDownmix = downmix * _pv;
      pvDownmix /= pvDownmix.norm();

      // vary the output level by the balance between the front and rear
      // loudspeakers; 0dB at the front to -3dB at the back
      double front = _pv.head(3).maxCoeff();
      double back = _pv.tail(2).maxCoeff();

      pvDownmix = pvDownmix * std::pow(0.5, (0.5 * back / (front + back)));
      return pvDownmix;
    }
    return boost::none;
  }

  AllocentricPanner::AllocentricPanner(Eigen::MatrixXd positions)
      : _numChannels(static_cast<int>(positions.rows())) {
    // build the speaker tree; see _speaker_tree in the reference
    // implementation. Exact comparisons are used deliberately, as the
    // positions come from a fixed table.
    for (int index = 0; index < _numChannels; index++) {
      double x = positions(index, 0), y = positions(index, 1),
             z = positions(index, 2);

      auto planeIt = _tree.begin();
      for (; planeIt != _tree.end(); ++planeIt) {
        if (planeIt->z == z) break;
        if (planeIt->z > z) {
          planeIt = _tree.end();
          break;
        }
      }
      if (planeIt == _tree.end()) {
        // insert a new plane keeping the planes sorted by z
        auto insertIt = _tree.begin();
        while (insertIt != _tree.end() && insertIt->z < z) ++insertIt;
        _tree.insert(insertIt, Plane{z, {Row{y, {Column{x, index}}}}});
        continue;
      }

      auto& rows = planeIt->rows;
      auto rowIt = rows.begin();
      for (; rowIt != rows.end(); ++rowIt) {
        if (rowIt->y == y) break;
        if (rowIt->y > y) {
          rowIt = rows.end();
          break;
        }
      }
      if (rowIt == rows.end()) {
        auto insertIt = rows.begin();
        while (insertIt != rows.end() && insertIt->y < y) ++insertIt;
        rows.insert(insertIt, Row{y, {Column{x, index}}});
        continue;
      }

      auto& columns = rowIt->columns;
      auto insertIt = columns.begin();
      while (insertIt != columns.end() && insertIt->x < x) ++insertIt;
      if (insertIt != columns.end() && insertIt->x == x)
        throw invalid_argument(
            "two loudspeakers with the same allocentric position");
      columns.insert(insertIt, Column{x, index});
    }
  }

  std::pair<double, double> AllocentricPanner::singleBalancePan(double minimum,
                                                                double maximum,
                                                                double value) {
    if (minimum == maximum) return {1.0, 1.0};
    if (value <= minimum) return {0.0, 1.0};
    if (value >= maximum) return {1.0, 0.0};
    double a = (value - minimum) / (maximum - minimum);
    double aa = a * boost::math::constants::pi<double>() / 2.0;
    return {std::cos(aa), std::sin(aa)};
  }

  std::pair<int, int> AllocentricPanner::findPlanes(double z) const {
    if (z <= _tree[0].z) return {0, 0};
    for (size_t i = 0; i < _tree.size(); i++) {
      if (_tree[i].z == z) return {(int)i, (int)i};
      if (_tree[i].z > z) return {(int)i - 1, (int)i};
    }
    return {(int)_tree.size() - 1, (int)_tree.size() - 1};
  }

  std::pair<int, int> AllocentricPanner::findRows(const Plane& plane,
                                                  double y) {
    if (y <= plane.rows[0].y) return {0, 0};
    for (size_t i = 0; i < plane.rows.size(); i++) {
      if (plane.rows[i].y == y) return {(int)i, (int)i};
      if (plane.rows[i].y > y) return {(int)i - 1, (int)i};
    }
    return {(int)plane.rows.size() - 1, (int)plane.rows.size() - 1};
  }

  std::pair<int, int> AllocentricPanner::findColumns(const Row& row,
                                                     double x) {
    if (x <= row.columns[0].x) return {0, 0};
    for (size_t i = 0; i < row.columns.size(); i++) {
      if (row.columns[i].x == x) return {(int)i, (int)i};
      if (row.columns[i].x > x) return {(int)i - 1, (int)i};
    }
    return {(int)row.columns.size() - 1, (int)row.columns.size() - 1};
  }

  boost::optional<Eigen::VectorXd> AllocentricPanner::handle(
      Eigen::Vector3d position) {
    Eigen::VectorXd ret = Eigen::VectorXd::Zero(_numChannels);

    std::pair<int, int> zPlanes = findPlanes(position(2));
    std::pair<double, double> zGains = singleBalancePan(
        _tree[zPlanes.first].z, _tree[zPlanes.second].z, position(2));

    for (auto zz : {std::make_pair(zGains.first, zPlanes.first),
                    std::make_pair(zGains.second, zPlanes.second)}) {
      const Plane& plane = _tree[zz.second];
      std::pair<int, int> yRows = findRows(plane, position(1));
      std::pair<double, double> yGains = singleBalancePan(
          plane.rows[yRows.first].y, plane.rows[yRows.second].y, position(1));

      for (auto yy : {std::make_pair(yGains.first, yRows.first),
                      std::make_pair(yGains.second, yRows.second)}) {
        const Row& row = plane.rows[yy.second];
        std::pair<int, int> xColumns = findColumns(row, position(0));
        std::pair<double, double> xGains =
            singleBalancePan(row.columns[xColumns.first].x,
                             row.columns[xColumns.second].x, position(0));

        for (auto xx : {std::make_pair(xGains.first, xColumns.first),
                        std::make_pair(xGains.second, xColumns.second)}) {
          ret(row.columns[xx.second].channel) = zz.first * yy.first * xx.first;
        }
      }
    }
    return ret;
  }

  int AllocentricPanner::numberOfOutputChannels() const {
    return _numChannels;
  }

  std::shared_ptr<PointSourcePanner> configureStereoPolarPanner(
      const Layout& layout) {
    auto leftChannel = layout.channelWithName("M+030");
    auto rightChannel = layout.channelWithName("M-030");
    auto channelNames = layout.channelNames();
    auto leftChannelIndex =
        distance(channelNames.begin(),
                 find(channelNames.begin(), channelNames.end(), "M+030"));
    auto rightChannelIndex =
        distance(channelNames.begin(),
                 find(channelNames.begin(), channelNames.end(), "M-030"));

    Eigen::MatrixXd positions(2, 3);
    positions << toCartesianVector3d(leftChannel.polarPosition()).transpose(),
        toCartesianVector3d(rightChannel.polarPosition()).transpose();
    Eigen::Vector2i outputChannels{leftChannelIndex, rightChannelIndex};

    auto panner =
        boost::make_unique<StereoPannerDownmix>(outputChannels, positions);

    std::vector<std::unique_ptr<RegionHandler>> regions;
    regions.push_back(std::move(panner));
    return std::make_shared<PolarPointSourcePanner>(std::move(regions));
  }

  std::tuple<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>,
             std::set<int>, Eigen::MatrixXd>
  getAugmentedLayout(const Layout& layout) {
    // add some extra height speakers that are treated as real speakers until
    // the downmix in PointSourcePannerDownmix
    std::vector<Channel> allChannels(layout.channels());

    std::vector<Channel> extraChannels;
    Eigen::MatrixXd downmix;
    std::tie(extraChannels, downmix) = extraPosVerticalNominal(layout);
    for (const auto& extraChannel : extraChannels) {
      allChannels.push_back(extraChannel);
    }

    // add some virtual speakers above and below that will be used as the
    // centre speaker in a virtual ngon. No upper speaker is added for
    // layouts with UH+180 as this speaker may actually be directly
    // overhead, which may cause a step in the gains wrt the source
    // position.
    auto channelNames = layout.channelNames();

    std::vector<Eigen::Vector3d> virtualPositions;
    virtualPositions.push_back(Eigen::Vector3d{0.0, 0.0, -1.0});
    if (std::find(channelNames.begin(), channelNames.end(), "T+000") ==
            channelNames.end() &&
        std::find(channelNames.begin(), channelNames.end(), "UH+180") ==
            channelNames.end()) {
      virtualPositions.push_back(Eigen::Vector3d{0.0, 0.0, 1.0});
    }

    std::vector<Eigen::Vector3d> positionsReal;
    std::vector<Eigen::Vector3d> positionsNominal;
    std::set<int> virtualVerts;

    for (const Channel& channel : allChannels) {
      positionsReal.push_back(toNormalisedVector3d(channel.polarPosition()));
      positionsNominal.push_back(
          toNormalisedVector3d(channel.polarPositionNominal()));
    }
    for (const Eigen::Vector3d& pos : virtualPositions) {
      virtualVerts.insert(static_cast<int>(positionsReal.size()));
      positionsReal.push_back(pos);
      positionsNominal.push_back(pos);
    }

    return {positionsReal, positionsNominal, virtualVerts, downmix};
  }

  /// Does the layout have the same channels and nominal positions as the
  /// stock BS.2051 layout with the same name?
  bool hasStockNominalPositions(const Layout& layout) {
    const Layout* stock = nullptr;
    for (const Layout& l : loadLayouts())
      if (l.name() == layout.name()) stock = &l;
    if (!stock) return false;

    const auto& channels = layout.channels();
    const auto& stockChannels = stock->channels();
    if (channels.size() != stockChannels.size()) return false;

    for (size_t i = 0; i < channels.size(); i++) {
      if (channels[i].name() != stockChannels[i].name()) return false;
      PolarPosition a = channels[i].polarPositionNominal();
      PolarPosition b = stockChannels[i].polarPositionNominal();
      if (a.azimuth != b.azimuth || a.elevation != b.elevation ||
          a.distance != b.distance)
        return false;
    }
    return true;
  }

  /// Modify layout to set the nominal positions of M+-SC; 45 degrees if they
  /// are wider than 30 degrees, otherwise 15. This ensures that the
  /// triangulation is correct for both the nominal and real loudspeaker
  /// positions when the ordering is changed.
  Layout setScreenSpeakerNominalPositions(const Layout& layout) {
    Layout ret = layout;
    for (Channel& channel : ret.channels()) {
      if (channel.name() == "M+SC" || channel.name() == "M-SC") {
        double oldAz = channel.polarPosition().azimuth;
        double newAz = (oldAz > 0 ? 1.0 : -1.0) * (std::abs(oldAz) > 30.0 ? 45.0 : 15.0);
        channel.polarPositionNominal(PolarPosition(newAz, 0.0, 1.0));
      }
    }
    return ret;
  }

  std::shared_ptr<PointSourcePanner> configureFullPolarPanner(
      const Layout& layoutIn) {
    Layout layout = setScreenSpeakerNominalPositions(layoutIn);

    std::vector<Eigen::Vector3d> positionsReal;
    std::vector<Eigen::Vector3d> positionsNominal;
    std::set<int> virtualVerts;
    Eigen::MatrixXd downmix;
    std::tie(positionsReal, positionsNominal, virtualVerts, downmix) =
        getAugmentedLayout(layout);

    // Facets of the convex hull; each set represents a facet and contains the
    // indices of its corners in positions. Pre-computed facets are only used
    // if the layout has the nominal positions of the stock BS.2051 layout
    // with that name (which is not the case for wide screen loudspeakers).
    auto facets_it = FACETS.find(layout.name());
    std::vector<Facet> facets =
        facets_it != FACETS.end() && hasStockNominalPositions(layout)
            ? facets_it->second
            : convex_hull(positionsNominal);

    // Turn the facets into regions for the point source panner.
    std::vector<std::unique_ptr<RegionHandler>> regions;

    // Facets adjacent to one of the virtual speakers are turned into virtual
    // ngons, with an equal power downmix from the virtual speaker to the real
    // speakers.
    for (int virtualVert : virtualVerts) {
      Facet realVerts = _adjacent_verts(facets, virtualVert);

      ear_assert(!doIntersect(realVerts.begin(), realVerts.end(),
                              virtualVerts.begin(), virtualVerts.end()),
                 "invalid triangulation");

      std::vector<int> realVertsVec(realVerts.begin(), realVerts.end());
      Eigen::Map<Eigen::VectorXi> outputChannels(realVertsVec.data(),
                                                 realVertsVec.size());
      Eigen::MatrixXd positions(outputChannels.size(), 3);
      int rowIndex = 0;
      for (int vert : realVerts) {
        positions.row(rowIndex) = positionsReal[vert];
        ++rowIndex;
      }
      Eigen::Vector3d centrePosition = positionsReal[virtualVert];
      Eigen::VectorXd centreDownmix(outputChannels.size());
      centreDownmix.fill(1.0 /
                         std::sqrt(static_cast<double>(outputChannels.size())));

      regions.push_back(boost::make_unique<VirtualNgon>(
          outputChannels, positions, centrePosition, centreDownmix));
    }
    // Facets not adjacent to virtual speakers are turned into triplets or
    // quads. In the supported layouts there are never facets with more
    // vertices.
    for (const auto& facetVerts : facets) {
      if (doIntersect(facetVerts.begin(), facetVerts.end(),
                      virtualVerts.begin(), virtualVerts.end())) {
        continue;
      }

      if (facetVerts.size() == 3) {
        std::vector<int> facetVertsVec(facetVerts.begin(), facetVerts.end());
        Eigen::Vector3i outputChannels(facetVertsVec.data());
        Eigen::MatrixXd positions(outputChannels.size(), 3);
        int rowIndex = 0;
        for (int vert : facetVerts) {
          positions.row(rowIndex) = positionsReal[vert];
          ++rowIndex;
        }
        regions.push_back(
            boost::make_unique<Triplet>(outputChannels, positions));
      } else if (facetVerts.size() == 4) {
        std::vector<int> facetVertsVec(facetVerts.begin(), facetVerts.end());
        Eigen::Vector4i outputChannels(facetVertsVec.data());
        Eigen::MatrixXd positions(outputChannels.size(), 3);
        int rowIndex = 0;
        for (int vert : facetVerts) {
          positions.row(rowIndex) = positionsReal[vert];
          ++rowIndex;
        }
        regions.push_back(
            boost::make_unique<QuadRegion>(outputChannels, positions));
      } else {
        throw internal_error(
            "facets with more than 4 vertices are not supported");
      }
    }
    return std::make_shared<PointSourcePannerDownmix>(
        std::make_shared<PolarPointSourcePanner>(std::move(regions)), downmix);
  }

  // Check that screen loudspeakers are within allowed ranges.
  void checkScreenSpeakers(const Layout& layout) {
    for (auto& channel : layout.channels()) {
      if (channel.name() == "M+SC" || channel.name() == "M-SC") {
        double abs_az = std::abs(channel.polarPosition().azimuth);
        if (!((5.0 <= abs_az && abs_az <= 25.0) ||
              (35.0 <= abs_az && abs_az <= 60.0))) {
          throw invalid_argument(
              "M+SC or M-SC has azimuth not in the allowed ranges of 5 to 25 "
              "and 35 to 60 degrees");
        }
      }
    }
  }

  std::shared_ptr<PointSourcePanner> configureAllocentricPanner(
      const Layout& layout) {
    checkScreenSpeakers(layout);
    return std::make_shared<AllocentricPanner>(
        allocentric::positionsForLayout(layout));
  }

  std::shared_ptr<PointSourcePanner> configurePolarPanner(
      const Layout& layout) {
    auto isLfe = layout.isLfe();
    if (find(isLfe.begin(), isLfe.end(), true) != isLfe.end()) {
      throw internal_error("lfe channel passed to point source panner");
    }

    checkScreenSpeakers(layout);

    if (layout.name() == std::string("0+2+0")) {
      return configureStereoPolarPanner(layout);
    } else {
      return configureFullPolarPanner(layout);
    }
  }
}  // namespace ear
