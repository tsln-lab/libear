// Compare the Objects gain calculator against gains produced by the reference
// implementation (the EBU ADM Renderer), see
// tools/reference/generate_objects_reference.py.
#include <catch2/catch.hpp>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include "ear/bs2051.hpp"
#include "ear/ear.hpp"
#include "ear/exceptions.hpp"
#include "reference/objects_reference.hpp"

using namespace ear;
using namespace ear::reference;

namespace {
  /// features implemented by libear; cases using other features must throw
  /// not_implemented
  const unsigned implementedFeatures = CARTESIAN | CHANNEL_LOCK | DIVERGENCE |
                                       ZONE_EXCLUSION | SCREEN_REF |
                                       SCREEN_EDGE_LOCK | WIDE_SCREEN_SPEAKERS;

  ObjectsTypeMetadata toMetadata(const ObjectsCase& c) {
    ObjectsTypeMetadata otm;
    otm.cartesian = c.cartesian;
    if (c.cartesian)
      otm.position = CartesianPosition(c.position[0], c.position[1], c.position[2]);
    else
      otm.position = PolarPosition(c.position[0], c.position[1], c.position[2]);
    otm.width = c.width;
    otm.height = c.height;
    otm.depth = c.depth;
    otm.gain = c.gain;
    otm.diffuse = c.diffuse;

    boost::optional<double> maxDistance;
    if (c.channelLockMaxDistance >= 0.0) maxDistance = c.channelLockMaxDistance;
    otm.channelLock = ChannelLock(c.channelLock, maxDistance);

    if (c.divergenceCartesian)
      otm.objectDivergence = CartesianObjectDivergence(
          c.divergence, c.divergenceRange >= 0.0 ? c.divergenceRange : 0.0);
    else
      otm.objectDivergence = PolarObjectDivergence(
          c.divergence, c.divergenceRange >= 0.0 ? c.divergenceRange : 45.0);

    for (const auto& zone : c.zones) {
      if (zone.cartesian)
        otm.zoneExclusion.zones.push_back(CartesianExclusionZone{
            (float)zone.values[0], (float)zone.values[1], (float)zone.values[2],
            (float)zone.values[3], (float)zone.values[4], (float)zone.values[5],
            ""});
      else
        otm.zoneExclusion.zones.push_back(PolarExclusionZone{
            (float)zone.values[0], (float)zone.values[1], (float)zone.values[2],
            (float)zone.values[3], 0.0f, 0.0f, ""});
    }
    otm.screenRef = c.screenRef;
    if (!c.screenEdgeLockHorizontal.empty())
      otm.screenEdgeLock.horizontal = c.screenEdgeLockHorizontal;
    if (!c.screenEdgeLockVertical.empty())
      otm.screenEdgeLock.vertical = c.screenEdgeLockVertical;
    return otm;
  }

  std::string describe(const ObjectsCase& c) {
    return c.name + " " + c.layout + " [" + std::to_string(c.position[0]) + ", " +
           std::to_string(c.position[1]) + ", " + std::to_string(c.position[2]) +
           "] wh=" + std::to_string(c.width) + "/" + std::to_string(c.height) +
           " d=" + std::to_string(c.depth) + " lock=" + std::to_string(c.channelLock) +
           "/" + std::to_string(c.channelLockMaxDistance) +
           " div=" + std::to_string(c.divergence) + "/" + std::to_string(c.divergenceRange) +
           " features=" + std::to_string(c.features);
  }

  void requireVectorApprox(const std::vector<double>& actual,
                           const std::vector<double>& expected, double tol,
                           const std::string& what) {
    REQUIRE(actual.size() == expected.size());
    for (size_t i = 0; i < actual.size(); i++) {
      INFO(what << " channel " << i);
      REQUIRE(actual[i] == Approx(expected[i]).margin(tol));
    }
  }
}  // namespace

TEST_CASE("objects gains match the reference implementation") {
  std::map<std::string, Layout> layouts;
  std::map<std::string, std::unique_ptr<GainCalculatorObjects>> calculators;

  size_t compared = 0, expectedUnsupported = 0;

  for (const auto& c : objectsCases()) {
    INFO(describe(c));

    if (!calculators.count(c.layout)) {
      layouts.emplace(c.layout, layoutForSpec(c.layout));
      calculators.emplace(c.layout, std::make_unique<GainCalculatorObjects>(layouts.at(c.layout)));
    }
    auto& calc = *calculators.at(c.layout);
    const size_t n = layouts.at(c.layout).channels().size();

    std::vector<double> direct(n), diffuse(n);
    ObjectsTypeMetadata otm = toMetadata(c);
    auto warning_cb = [](const Warning&) {};

    if (c.features & ~implementedFeatures) {
      REQUIRE_THROWS_AS(calc.calculate(otm, direct, diffuse, warning_cb), not_implemented);
      expectedUnsupported++;
      continue;
    }

    calc.calculate(otm, direct, diffuse, warning_cb);

    // the extent panner uses single precision internally (observed maximum
    // deviation around 3e-7); everything else matches to double precision
    const double tol = c.usesExtent ? 1e-5 : 1e-9;
    requireVectorApprox(direct, c.direct, tol, "direct");
    requireVectorApprox(diffuse, c.diffuseGains, tol, "diffuse");
    compared++;
  }

  WARN("compared " << compared << " cases with the reference; " << expectedUnsupported
                   << " cases use features that are not implemented yet");
  REQUIRE(compared > 0);
}
