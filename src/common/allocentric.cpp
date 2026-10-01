#include "allocentric.hpp"
#include <cmath>
#include <map>
#include <string>
#include "ear/conversion.hpp"
#include "ear/exceptions.hpp"
#include "geom.hpp"

namespace ear {
  namespace allocentric {

    namespace {
      using ChannelPositions = std::map<std::string, Eigen::Vector3d>;

      // from data/allo_positions.yaml in the reference implementation
      const std::map<std::string, ChannelPositions>& alloPositions() {
        static const std::map<std::string, ChannelPositions> positions = {
            {"0+2+0",
             {{"M+030", {-1.0, 1.0, 0.0}}, {"M-030", {1.0, 1.0, 0.0}}}},
            {"0+5+0",
             {{"M+030", {-1.0, 1.0, 0.0}},
              {"M-030", {1.0, 1.0, 0.0}},
              {"M+000", {0.0, 1.0, 0.0}},
              {"M+110", {-1.0, -1.0, 0.0}},
              {"M-110", {1.0, -1.0, 0.0}},
              {"LFE1", {-1.0, 1.0, -1.0}}}},
            {"0+7+0",
             {{"M+030", {-1.0, 1.0, 0.0}},
              {"M-030", {1.0, 1.0, 0.0}},
              {"M+000", {0.0, 1.0, 0.0}},
              {"M+090", {-1.0, 0.0, 0.0}},
              {"M-090", {1.0, 0.0, 0.0}},
              {"M+135", {-1.0, -1.0, 0.0}},
              {"M-135", {1.0, -1.0, 0.0}},
              {"LFE1", {-1.0, 1.0, -1.0}}}},
            {"2+5+0",
             {{"M+030", {-1.0, 1.0, 0.0}},
              {"M-030", {1.0, 1.0, 0.0}},
              {"M+000", {0.0, 1.0, 0.0}},
              {"M+110", {-1.0, -1.0, 0.0}},
              {"M-110", {1.0, -1.0, 0.0}},
              {"U+030", {-1.0, 1.0, 1.0}},
              {"U-030", {1.0, 1.0, 1.0}},
              {"LFE1", {-1.0, 1.0, -1.0}}}},
            {"4+5+0",
             {{"M+030", {-1.0, 1.0, 0.0}},
              {"M-030", {1.0, 1.0, 0.0}},
              {"M+000", {0.0, 1.0, 0.0}},
              {"M+110", {-1.0, -1.0, 0.0}},
              {"M-110", {1.0, -1.0, 0.0}},
              {"U+030", {-1.0, 1.0, 1.0}},
              {"U-030", {1.0, 1.0, 1.0}},
              {"U+110", {-1.0, -1.0, 1.0}},
              {"U-110", {1.0, -1.0, 1.0}},
              {"LFE1", {-1.0, 1.0, -1.0}}}},
            {"4+7+0",
             {{"M+030", {-1.0, 1.0, 0.0}},
              {"M-030", {1.0, 1.0, 0.0}},
              {"M+000", {0.0, 1.0, 0.0}},
              {"M+090", {-1.0, 0.0, 0.0}},
              {"M-090", {1.0, 0.0, 0.0}},
              {"M+135", {-1.0, -1.0, 0.0}},
              {"M-135", {1.0, -1.0, 0.0}},
              {"U+045", {-1.0, 1.0, 1.0}},
              {"U-045", {1.0, 1.0, 1.0}},
              {"U+135", {-1.0, -1.0, 1.0}},
              {"U-135", {1.0, -1.0, 1.0}},
              {"LFE1", {-1.0, 1.0, -1.0}}}},
            {"4+5+1",
             {{"M+030", {-1.0, 1.0, 0.0}},
              {"M-030", {1.0, 1.0, 0.0}},
              {"M+000", {0.0, 1.0, 0.0}},
              {"M+110", {-1.0, -1.0, 0.0}},
              {"M-110", {1.0, -1.0, 0.0}},
              {"U+030", {-1.0, 1.0, 1.0}},
              {"U-030", {1.0, 1.0, 1.0}},
              {"U+110", {-1.0, -1.0, 1.0}},
              {"U-110", {1.0, -1.0, 1.0}},
              {"B+000", {0.0, 1.0, -1.0}},
              {"LFE1", {-1.0, 1.0, -1.0}}}},
            {"3+7+0",
             {{"M+000", {0.0, 1.0, 0.0}},
              {"M+030", {-1.0, 1.0, 0.0}},
              {"M-030", {1.0, 1.0, 0.0}},
              {"U+045", {-1.0, 1.0, 1.0}},
              {"U-045", {1.0, 1.0, 1.0}},
              {"M+090", {-1.0, 0.0, 0.0}},
              {"M-090", {1.0, 0.0, 0.0}},
              {"M+135", {-1.0, -1.0, 0.0}},
              {"M-135", {1.0, -1.0, 0.0}},
              {"UH+180", {0.0, -1.0, 1.0}},
              {"LFE1", {-1.0, 1.0, -1.0}},
              {"LFE2", {1.0, 1.0, -1.0}}}},
            {"4+9+0",
             {{"M+030", {-1.0, 1.0, 0.0}},
              {"M-030", {1.0, 1.0, 0.0}},
              {"M+000", {0.0, 1.0, 0.0}},
              {"M+090", {-1.0, 0.0, 0.0}},
              {"M-090", {1.0, 0.0, 0.0}},
              {"M+135", {-1.0, -1.0, 0.0}},
              {"M-135", {1.0, -1.0, 0.0}},
              {"U+045", {-1.0, 1.0, 1.0}},
              {"U-045", {1.0, 1.0, 1.0}},
              {"U+135", {-1.0, -1.0, 1.0}},
              {"U-135", {1.0, -1.0, 1.0}},
              {"LFE1", {-1.0, 1.0, -1.0}}}},
            {"9+10+3",
             {{"M+060", {-1.0, 0.414214, 0.0}},
              {"M-060", {1.0, 0.414214, 0.0}},
              {"M+000", {0.0, 1.0, 0.0}},
              {"M+135", {-1.0, -1.0, 0.0}},
              {"M-135", {1.0, -1.0, 0.0}},
              {"M+030", {-1.0, 1.0, 0.0}},
              {"M-030", {1.0, 1.0, 0.0}},
              {"M+180", {0.0, -1.0, 0.0}},
              {"M+090", {-1.0, 0.0, 0.0}},
              {"M-090", {1.0, 0.0, 0.0}},
              {"U+045", {-1.0, 1.0, 1.0}},
              {"U-045", {1.0, 1.0, 1.0}},
              {"U+000", {0.0, 1.0, 1.0}},
              {"T+000", {0.0, 0.0, 1.0}},
              {"U+135", {-1.0, -1.0, 1.0}},
              {"U-135", {1.0, -1.0, 1.0}},
              {"U+090", {-1.0, 0.0, 1.0}},
              {"U-090", {1.0, 0.0, 1.0}},
              {"U+180", {0.0, -1.0, 1.0}},
              {"B+000", {0.0, 1.0, -1.0}},
              {"B+045", {-1.0, 1.0, -1.0}},
              {"B-045", {1.0, 1.0, -1.0}},
              {"LFE1", {-1.0, 1.0, -1.0}},
              {"LFE2", {1.0, 1.0, -1.0}}}},
        };
        return positions;
      }

      double sign(double x) { return (x > 0.0) - (x < 0.0); }
    }  // namespace

    Eigen::Vector3d screenSpeakerPositionToCart(const PolarPosition& position) {
      // the y position of these loudspeakers must be identical, and they must
      // be either exactly at the front or the side
      CartesianPosition c = conversion::pointPolarToCart(
          PolarPosition(std::abs(position.azimuth), 0.0, 1.0));
      Eigen::Vector3d posLeft(c.X, c.Y, c.Z);

      bool atFront = std::abs(posLeft(1) - 1.0) < 1e-10;
      bool atSide = std::abs(posLeft(0) + 1.0) < 1e-10;

      if (!(atFront || atSide))
        throw invalid_argument(
            "screen loudspeaker position is neither at the front nor at the "
            "side");

      if (atFront) posLeft(1) = 1.0;
      if (atSide) posLeft(0) = -1.0;

      return posLeft.cwiseProduct(
          Eigen::Vector3d(sign(position.azimuth), 1.0, 1.0));
    }

    boost::optional<Eigen::MatrixXd> positionsForLayoutIfKnown(
        const Layout& layout) {
      auto layoutIt = alloPositions().find(layout.name());
      if (layoutIt == alloPositions().end()) return boost::none;
      const ChannelPositions& layoutPositions = layoutIt->second;

      const auto& channels = layout.channels();
      Eigen::MatrixXd positions(channels.size(), 3);
      for (size_t i = 0; i < channels.size(); i++) {
        const std::string& name = channels[i].name();
        if (name == "M+SC" || name == "M-SC") {
          positions.row(i) =
              screenSpeakerPositionToCart(channels[i].polarPosition())
                  .transpose();
        } else {
          auto channelIt = layoutPositions.find(name);
          if (channelIt == layoutPositions.end()) return boost::none;
          positions.row(i) = channelIt->second.transpose();
        }
      }
      return positions;
    }

    Eigen::MatrixXd positionsForLayout(const Layout& layout) {
      auto positions = positionsForLayoutIfKnown(layout);
      if (!positions)
        throw invalid_argument(
            "no allocentric loudspeaker positions are known for layout '" +
            layout.name() + "' with its channel names");
      return *positions;
    }

    std::vector<bool> getExcluded(const Eigen::MatrixXd& channelPositions,
                                  std::vector<bool> isExcluded) {
      const Eigen::Index n = channelPositions.rows();

      // Remove additional speakers to ensure the layout works well with our
      // panner
      std::vector<bool> original = isExcluded;
      for (Eigen::Index i = 0; i < n; i++) {
        if (original[i] && std::abs(channelPositions(i, 0)) == 1.0 &&
            std::abs(channelPositions(i, 1)) != 1.0) {
          for (Eigen::Index k = 0; k < n; k++) {
            if (channelPositions(k, 1) == channelPositions(i, 1) &&
                channelPositions(k, 2) == channelPositions(i, 2))
              isExcluded[k] = true;
          }
        }
      }

      // Don't do any exclusion if the previous steps result in a layout with
      // no speakers
      bool all = true;
      for (bool ex : isExcluded) all = all && ex;
      if (all) std::fill(isExcluded.begin(), isExcluded.end(), false);

      return isExcluded;
    }

  }  // namespace allocentric
}  // namespace ear
