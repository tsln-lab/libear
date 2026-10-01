#include <boost/make_unique.hpp>
#include <cassert>
#include <catch2/catch.hpp>
#include <cmath>
#include <map>
#include <string>
#include <vector>
#include "ear/bs2051.hpp"
#include "ear/ear.hpp"
#include "ear/metadata.hpp"
#include "common/geom.hpp"

using namespace ear;

using Gains = std::vector<float>;

ObjectsTypeMetadata otmWithPos(Position p) {
  ObjectsTypeMetadata otm;
  otm.position = p;
  return otm;
}

using GainsMap = std::map<std::string, double>;
struct GainsMaps {
  GainsMap direct;
  GainsMap diffuse;
};

GainsMap mapGainsToChannelNames(std::vector<std::string> channelNames,
                                Gains gains) {
  assert(channelNames.size() == gains.size());
  GainsMap ret;
  std::transform(channelNames.begin(), channelNames.end(), gains.begin(),
                 std::inserter(ret, ret.end()),
                 std::make_pair<std::string const&, double const&>);
  return ret;
}

class GainCalculatorObjectsTester {
 public:
  GainCalculatorObjectsTester(Layout layout)
      : _layout(layout),
        _gainCalc(boost::make_unique<GainCalculatorObjects>(layout)){};

  GainsMaps run(const ObjectsTypeMetadata& metadata) {
    Gains directGains(_layout.channels().size(), 0.0);
    Gains diffuseGains(_layout.channels().size(), 0.0);
    _gainCalc->calculate(metadata, directGains, diffuseGains);
    auto channelNames = _layout.channelNames();
    auto directGainsMap = mapGainsToChannelNames(channelNames, directGains);
    auto diffuseGainsMap = mapGainsToChannelNames(channelNames, diffuseGains);
    double epsilon = 1e-6;
    for (auto it = directGainsMap.begin(); it != directGainsMap.end();) {
      if (std::abs(it->second) < epsilon) {
        it = directGainsMap.erase(it);
      } else {
        ++it;
      }
    }
    for (auto it = diffuseGainsMap.begin(); it != diffuseGainsMap.end();) {
      if (std::abs(it->second) < epsilon) {
        it = diffuseGainsMap.erase(it);
      } else {
        ++it;
      }
    }
    return GainsMaps{directGainsMap, diffuseGainsMap};
  };

 private:
  Layout _layout;
  std::unique_ptr<GainCalculatorObjects> _gainCalc;
};

TEST_CASE("gain_calculator") {
  auto layout = getLayout("4+7+0").withoutLfe();
  GainCalculatorObjectsTester gainCalc(layout);

  SECTION("basic_centre") {
    auto gainMap = gainCalc.run(otmWithPos(PolarPosition{0.0, 0.0, 1.0}));
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("M+000") == Approx(1.0));
    REQUIRE(gainMap.diffuse.size() == 0);
  }
  SECTION("basic_left") {
    auto gainMap = gainCalc.run(otmWithPos(PolarPosition{30.0, 0.0, 1.0}));
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("M+030") == Approx(1.0));
    REQUIRE(gainMap.diffuse.size() == 0);
  }
  SECTION("basic_left_up") {
    auto gainMap = gainCalc.run(otmWithPos(PolarPosition{45.0, 30.0, 1.0}));
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("U+045") == Approx(1.0));
    REQUIRE(gainMap.diffuse.size() == 0);
  }
}

TEST_CASE("diffuse") {
  auto layout = getLayout("4+7+0").withoutLfe();
  GainCalculatorObjectsTester gainCalc(layout);

  SECTION("half") {
    auto otm = otmWithPos(PolarPosition{0.0, 0.0, 1.0});
    otm.diffuse = 0.5;
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("M+000") == Approx(std::sqrt(0.5)));
    REQUIRE(gainMap.diffuse.size() == 1);
    REQUIRE(gainMap.diffuse.at("M+000") == Approx(std::sqrt(0.5)));
  }

  SECTION("full") {
    auto otm = otmWithPos(PolarPosition{0.0, 0.0, 1.0});
    otm.diffuse = 1.0;
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 0);
    REQUIRE(gainMap.diffuse.size() == 1);
    REQUIRE(gainMap.diffuse.at("M+000") == Approx(1.0));
  }
}

TEST_CASE("gain_value") {
  auto layout = getLayout("4+7+0").withoutLfe();
  GainCalculatorObjectsTester gainCalc(layout);

  auto otm = otmWithPos(PolarPosition{0.0, 0.0, 1.0});
  otm.gain = 0.5;
  auto gainMap = gainCalc.run(otm);
  REQUIRE(gainMap.direct.size() == 1);
  REQUIRE(gainMap.direct.at("M+000") == Approx(0.5));
  REQUIRE(gainMap.diffuse.size() == 0);
}

// channel lock and divergence tests ported from the reference implementation
// (ear/core/objectbased/test/test_gain_calc.py)

TEST_CASE("channel_lock") {
  auto layout = getLayout("4+5+0").withoutLfe();
  GainCalculatorObjectsTester gainCalc(layout);

  SECTION("on_speaker") {
    auto otm = otmWithPos(PolarPosition{0.0, 0.0, 1.0});
    otm.channelLock = ChannelLock(true, 1.0);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("M+000") == Approx(1.0));
  }
  SECTION("close") {
    auto otm = otmWithPos(PolarPosition{14.0, 0.0, 1.0});
    otm.channelLock = ChannelLock(true, 1.0);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("M+000") == Approx(1.0));
  }
  SECTION("not_close_enough") {
    auto otm = otmWithPos(PolarPosition{15.0, 0.0, 1.0});
    double maxDistance = (cart(0.0, 0.0, 1.0) - cart(15.0, 0.0, 1.0)).norm() - 0.01;
    otm.channelLock = ChannelLock(true, maxDistance);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 2);
    REQUIRE(gainMap.direct.at("M+000") == Approx(std::sqrt(0.5)));
    REQUIRE(gainMap.direct.at("M+030") == Approx(std::sqrt(0.5)));
  }
  SECTION("abs_elevation_priority") {
    auto otm = otmWithPos(PolarPosition{30.0, 15.0, 1.0});
    otm.channelLock = ChannelLock(true);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("M+030") == Approx(1.0));
  }
  SECTION("abs_az_priority_left") {
    auto otm = otmWithPos(PolarPosition{15.0, 0.0, 1.0});
    otm.channelLock = ChannelLock(true);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("M+000") == Approx(1.0));
  }
  SECTION("abs_az_priority_right") {
    auto otm = otmWithPos(PolarPosition{-15.0, 0.0, 1.0});
    otm.channelLock = ChannelLock(true);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("M+000") == Approx(1.0));
  }
  SECTION("az_priority_front_top") {
    auto otm = otmWithPos(PolarPosition{0.0, 30.0, 1.0});
    otm.channelLock = ChannelLock(true);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("U-030") == Approx(1.0));
  }
  SECTION("az_priority_rear_top") {
    auto otm = otmWithPos(PolarPosition{180.0, 30.0, 1.0});
    otm.channelLock = ChannelLock(true);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("U-110") == Approx(1.0));
  }
  SECTION("az_priority_rear_mid") {
    auto otm = otmWithPos(PolarPosition{180.0, 0.0, 1.0});
    otm.channelLock = ChannelLock(true);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("M-110") == Approx(1.0));
  }
  SECTION("with LFE in layout") {
    GainCalculatorObjectsTester gainCalcLfe(getLayout("4+5+0"));
    auto otm = otmWithPos(PolarPosition{14.0, 0.0, 1.0});
    otm.channelLock = ChannelLock(true);
    auto gainMap = gainCalcLfe.run(otm);
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("M+000") == Approx(1.0));
  }
}

TEST_CASE("divergence") {
  auto layout = getLayout("4+5+0").withoutLfe();
  GainCalculatorObjectsTester gainCalc(layout);

  SECTION("half") {
    auto otm = otmWithPos(PolarPosition{0.0, 0.0, 1.0});
    otm.objectDivergence = PolarObjectDivergence(0.5, 30.0);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 3);
    REQUIRE(gainMap.direct.at("M+000") == Approx(std::sqrt(1.0 / 3.0)));
    REQUIRE(gainMap.direct.at("M+030") == Approx(std::sqrt(1.0 / 3.0)));
    REQUIRE(gainMap.direct.at("M-030") == Approx(std::sqrt(1.0 / 3.0)));
  }
  SECTION("full") {
    auto otm = otmWithPos(PolarPosition{0.0, 0.0, 1.0});
    otm.objectDivergence = PolarObjectDivergence(1.0, 30.0);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 2);
    REQUIRE(gainMap.direct.at("M+030") == Approx(std::sqrt(0.5)));
    REQUIRE(gainMap.direct.at("M-030") == Approx(std::sqrt(0.5)));
  }
  SECTION("azimuth") {
    auto otm = otmWithPos(PolarPosition{(30.0 + 110.0) / 2.0, 0.0, 1.0});
    otm.objectDivergence = PolarObjectDivergence(1.0, (110.0 - 30.0) / 2.0);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 2);
    REQUIRE(gainMap.direct.at("M+030") == Approx(std::sqrt(0.5)));
    REQUIRE(gainMap.direct.at("M+110") == Approx(std::sqrt(0.5)));
  }
  SECTION("elevation") {
    Eigen::Vector3d p = cart(30.0, 30.0, 1.0);
    p(0) = 0.0;
    auto otm = otmWithPos(PolarPosition{0.0, elevation(p), 1.0});
    otm.objectDivergence = PolarObjectDivergence(
        1.0, degrees(std::asin(cart(-30.0, 30.0, 1.0)(0))));
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 2);
    REQUIRE(gainMap.direct.at("U+030") == Approx(std::sqrt(0.5)));
    REQUIRE(gainMap.direct.at("U-030") == Approx(std::sqrt(0.5)));
  }
  SECTION("azimuth_elevation") {
    Eigen::Vector3d p = cart(40.0, 30.0, 1.0);
    p(0) = 0.0;
    auto otm = otmWithPos(PolarPosition{70.0, elevation(p), 1.0});
    otm.objectDivergence = PolarObjectDivergence(
        1.0, degrees(std::asin(cart(-40.0, 30.0, 1.0)(0))));
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 2);
    REQUIRE(gainMap.direct.at("U+030") == Approx(std::sqrt(0.5)));
    REQUIRE(gainMap.direct.at("U+110") == Approx(std::sqrt(0.5)));
  }
  SECTION("normalised") {
    for (auto params : std::vector<std::pair<double, double>>{
             {1.0, 0.0}, {1.0, 10.0}, {0.5, 10.0}, {1.0, 40.0}, {0.3, 10.0}, {0.7, 10.0}}) {
      auto otm = otmWithPos(PolarPosition{0.0, 0.0, 1.0});
      otm.objectDivergence = PolarObjectDivergence(params.first, params.second);
      auto gainMap = gainCalc.run(otm);
      double sumSquares = 0.0;
      for (auto& g : gainMap.direct) sumSquares += g.second * g.second;
      REQUIRE(std::sqrt(sumSquares) == Approx(1.0));
    }
  }
  SECTION("cartesian divergence in polar mode uses polar divergence with a warning") {
    auto otm = otmWithPos(PolarPosition{0.0, 0.0, 1.0});
    otm.objectDivergence = CartesianObjectDivergence(1.0, 1.0);
    std::vector<Warning> warnings;
    Gains direct(layout.channels().size()), diffuse(layout.channels().size());
    GainCalculatorObjects calc(layout);
    calc.calculate(otm, direct, diffuse,
                   [&](const Warning& w) { warnings.push_back(w); });
    REQUIRE(warnings.size() == 1);
    REQUIRE(warnings[0].code == Warning::Code::DIVERGENCE_POSITIONRANGE_IGNORED);
    // default azimuthRange of 45 degrees: between M+030 and M+110
    auto gainMap = mapGainsToChannelNames(layout.channelNames(), direct);
    REQUIRE(gainMap.at("M+000") == Approx(0.0).margin(1e-6));
    REQUIRE(gainMap.at("M+030") > 0.5);
    REQUIRE(gainMap.at("M+110") > 0.0);
  }
}

// zone exclusion and Cartesian tests ported from the reference implementation
// (ear/core/objectbased/test/test_gain_calc.py)

TEST_CASE("zone_exclusion") {
  auto layout = getLayout("4+5+0").withoutLfe();
  GainCalculatorObjectsTester gainCalc(layout);

  SECTION("front") {
    auto otm = otmWithPos(PolarPosition{0.0, 0.0, 1.0});
    otm.zoneExclusion.zones.push_back(
        PolarExclusionZone{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ""});
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 2);
    REQUIRE(gainMap.direct.at("M+030") == Approx(std::sqrt(0.5)));
    REQUIRE(gainMap.direct.at("M-030") == Approx(std::sqrt(0.5)));
  }
  SECTION("mid_front") {
    auto otm = otmWithPos(PolarPosition{0.0, 0.0, 1.0});
    otm.zoneExclusion.zones.push_back(
        PolarExclusionZone{-180.0, 180.0, 0.0, 0.0, 0.0, 0.0, ""});
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 2);
    REQUIRE(gainMap.direct.at("U+030") == Approx(std::sqrt(0.5)));
    REQUIRE(gainMap.direct.at("U-030") == Approx(std::sqrt(0.5)));
  }
}

TEST_CASE("cartesian") {
  auto layout = getLayout("4+5+0").withoutLfe();
  GainCalculatorObjectsTester gainCalc(layout);

  SECTION("channel_lock_on_speaker_cart") {
    auto otm = otmWithPos(CartesianPosition{1.0, 1.0, 0.0});
    otm.cartesian = true;
    otm.channelLock = ChannelLock(true, 0.01);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("M-030") == Approx(1.0));
  }
  SECTION("channel_lock_not_close_enough_cart") {
    auto otm = otmWithPos(CartesianPosition{0.5, 1.0, 0.0});
    otm.cartesian = true;
    otm.channelLock = ChannelLock(true, 0.01);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 2);
    REQUIRE(gainMap.direct.at("M-030") == Approx(std::sqrt(0.5)));
    REQUIRE(gainMap.direct.at("M+000") == Approx(std::sqrt(0.5)));
  }
  SECTION("channel_lock_no_max_enough_cart") {
    auto otm = otmWithPos(CartesianPosition{0.6, 1.0, 0.0});
    otm.cartesian = true;
    otm.channelLock = ChannelLock(true);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("M-030") == Approx(1.0));
  }
  SECTION("channel_lock_no_max_exclude") {
    auto otm = otmWithPos(CartesianPosition{0.1, 1.0, 0.0});
    otm.cartesian = true;
    otm.channelLock = ChannelLock(true);
    otm.zoneExclusion.zones.push_back(
        PolarExclusionZone{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ""});
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("M-030") == Approx(1.0));
  }
  SECTION("channel_lock_centre_cart") {
    GainCalculatorObjectsTester gainCalc22(getLayout("9+10+3").withoutLfe());
    auto otm = otmWithPos(CartesianPosition{0.0, 0.0, 0.0});
    otm.cartesian = true;
    otm.channelLock = ChannelLock(true);
    auto gainMap = gainCalc22.run(otm);
    REQUIRE(gainMap.direct.size() == 1);
    REQUIRE(gainMap.direct.at("M-090") == Approx(1.0));
  }
  SECTION("diverge_cart") {
    auto otm = otmWithPos(PolarPosition{0.0, 0.0, 1.0});
    otm.cartesian = true;
    otm.objectDivergence = CartesianObjectDivergence(0.5, 1.0);
    auto gainMap = gainCalc.run(otm);
    REQUIRE(gainMap.direct.size() == 3);
    REQUIRE(gainMap.direct.at("M+000") == Approx(std::sqrt(1.0 / 3.0)));
    REQUIRE(gainMap.direct.at("M+030") == Approx(std::sqrt(1.0 / 3.0)));
    REQUIRE(gainMap.direct.at("M-030") == Approx(std::sqrt(1.0 / 3.0)));
  }
  SECTION("custom layout without allocentric positions") {
    Layout custom = layout;
    custom.name("custom");
    GainCalculatorObjectsTester gainCalcCustom(custom);
    auto otm = otmWithPos(CartesianPosition{0.0, 1.0, 0.0});
    otm.cartesian = true;
    REQUIRE_THROWS_AS(gainCalcCustom.run(otm), invalid_argument);
    // polar rendering still works
    REQUIRE(gainCalcCustom.run(otmWithPos(PolarPosition{0.0, 0.0, 1.0})).direct.at("M+000") == Approx(1.0));
  }
}
