// Compare the DirectSpeakers gain calculator against gains produced by the
// reference implementation (the EBU ADM Renderer), see
// tools/reference/generate_direct_speakers_reference.py.
#include <catch2/catch.hpp>
#include <map>
#include <memory>
#include <string>
#include "ear/bs2051.hpp"
#include "ear/ear.hpp"
#include "reference/direct_speakers_reference.hpp"
#include "reference/objects_reference.hpp"

using namespace ear;
using namespace ear::reference;

namespace {
  DirectSpeakersTypeMetadata toMetadata(const DirectSpeakersCase& c) {
    DirectSpeakersTypeMetadata tm;
    tm.speakerLabels = c.speakerLabels;

    ScreenEdgeLock screenEdgeLock;
    if (!c.screenEdgeLockHorizontal.empty())
      screenEdgeLock.horizontal = c.screenEdgeLockHorizontal;
    if (!c.screenEdgeLockVertical.empty())
      screenEdgeLock.vertical = c.screenEdgeLockVertical;

    if (c.cartesian) {
      CartesianSpeakerPosition pos(c.position[0], c.position[1], c.position[2]);
      if (c.hasMin[0]) pos.XMin = c.minimum[0];
      if (c.hasMin[1]) pos.YMin = c.minimum[1];
      if (c.hasMin[2]) pos.ZMin = c.minimum[2];
      if (c.hasMax[0]) pos.XMax = c.maximum[0];
      if (c.hasMax[1]) pos.YMax = c.maximum[1];
      if (c.hasMax[2]) pos.ZMax = c.maximum[2];
      pos.screenEdgeLock = screenEdgeLock;
      tm.position = pos;
    } else {
      PolarSpeakerPosition pos(c.position[0], c.position[1], c.position[2]);
      if (c.hasMin[0]) pos.azimuthMin = c.minimum[0];
      if (c.hasMin[1]) pos.elevationMin = c.minimum[1];
      if (c.hasMin[2]) pos.distanceMin = c.minimum[2];
      if (c.hasMax[0]) pos.azimuthMax = c.maximum[0];
      if (c.hasMax[1]) pos.elevationMax = c.maximum[1];
      if (c.hasMax[2]) pos.distanceMax = c.maximum[2];
      pos.screenEdgeLock = screenEdgeLock;
      tm.position = pos;
    }

    if (c.lfeFreq) tm.channelFrequency.lowPass = 120.0;
    if (!c.packFormat.empty()) tm.audioPackFormatID = c.packFormat;
    return tm;
  }

  std::string describe(const DirectSpeakersCase& c) {
    std::string labels;
    for (const auto& l : c.speakerLabels) labels += l + " ";
    return c.name + " " + c.layout + " labels=[" + labels + "] cart=" +
           std::to_string(c.cartesian) + " pos=[" + std::to_string(c.position[0]) +
           ", " + std::to_string(c.position[1]) + ", " + std::to_string(c.position[2]) +
           "] min=[" + (c.hasMin[0] ? std::to_string(c.minimum[0]) : "-") + ", " +
           (c.hasMin[1] ? std::to_string(c.minimum[1]) : "-") + ", " +
           (c.hasMin[2] ? std::to_string(c.minimum[2]) : "-") + "] max=[" +
           (c.hasMax[0] ? std::to_string(c.maximum[0]) : "-") + ", " +
           (c.hasMax[1] ? std::to_string(c.maximum[1]) : "-") + ", " +
           (c.hasMax[2] ? std::to_string(c.maximum[2]) : "-") +
           "] lock=" + c.screenEdgeLockHorizontal + "/" + c.screenEdgeLockVertical +
           " lfe=" + std::to_string(c.lfeFreq) + " pack=" + c.packFormat;
  }
}  // namespace

TEST_CASE("DirectSpeakers gains match the reference implementation") {
  std::map<std::string, Layout> layouts;
  std::map<std::string, std::unique_ptr<GainCalculatorDirectSpeakers>> calculators;

  size_t compared = 0;
  for (const auto& c : directSpeakersCases()) {
    INFO(describe(c));

    if (!calculators.count(c.layout)) {
      layouts.emplace(c.layout, layoutForSpec(c.layout));
      calculators.emplace(c.layout, std::make_unique<GainCalculatorDirectSpeakers>(layouts.at(c.layout)));
    }
    auto& calc = *calculators.at(c.layout);
    const size_t n = layouts.at(c.layout).channels().size();

    std::vector<double> gains(n);
    calc.calculate(toMetadata(c), gains, [](const Warning&) {});

    REQUIRE(gains.size() == c.gains.size());
    std::string got, ref;
    for (size_t i = 0; i < n; i++) {
      got += std::to_string(gains[i]) + " ";
      ref += std::to_string(c.gains[i]) + " ";
    }
    INFO("got: " << got);
    INFO("ref: " << ref);
    for (size_t i = 0; i < n; i++) {
      INFO("channel " << i);
      REQUIRE(gains[i] == Approx(c.gains[i]).margin(1e-9));
    }
    compared++;
  }

  WARN("compared " << compared << " DirectSpeakers cases with the reference");
  REQUIRE(compared > 0);
}
