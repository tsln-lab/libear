#!/usr/bin/env python3
"""Generate reference gains for the DirectSpeakers renderer from the EBU ADM
Renderer; see generate_objects_reference.py for the general idea.

Writes ``tests/reference/direct_speakers_reference_data.cpp``, used by
``direct_speakers_reference_tests``.
"""
import os
import sys
import warnings
from itertools import product

warnings.filterwarnings("ignore")

from ear.core import bs2051  # noqa: E402
from ear.core.direct_speakers.panner import DirectSpeakersPanner  # noqa: E402
from ear.core.metadata_input import DirectSpeakersTypeMetadata, ExtraData  # noqa: E402
from ear.fileio.adm.adm import ADM  # noqa: E402
from ear.fileio.adm.common_definitions import load_common_definitions  # noqa: E402
from ear.fileio.adm.elements import (  # noqa: E402
    AudioBlockFormatDirectSpeakers, BoundCoordinate, DirectSpeakerCartesianPosition,
    DirectSpeakerPolarPosition, Frequency, ScreenEdgeLock)

sys.path.insert(0, os.path.dirname(__file__))
from generate_objects_reference import layout_for_spec, fmt  # noqa: E402

try:
    from importlib.metadata import version
    EAR_VERSION = version("ear")
except Exception:  # pragma: no cover
    EAR_VERSION = "unknown"

CARTESIAN = 1
SCREEN_EDGE_LOCK = 32
WIDE_SCREEN_SPEAKERS = 64

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "reference",
                   "direct_speakers_reference_data.cpp")

adm_common_defs = ADM()
load_common_definitions(adm_common_defs)


class Case(object):
    def __init__(self, name, layout, labels=(), cartesian=False, position=(0.0, 0.0, 1.0),
                 minimum=(None, None, None), maximum=(None, None, None),
                 screenEdgeLockH=None, screenEdgeLockV=None, lfeFreq=False, packFormat=None):
        self.name = name
        self.layout = layout
        self.labels = list(labels)
        self.cartesian = cartesian
        self.position = tuple(float(v) for v in position)
        self.minimum = tuple(minimum)
        self.maximum = tuple(maximum)
        self.screenEdgeLockH = screenEdgeLockH
        self.screenEdgeLockV = screenEdgeLockV
        self.lfeFreq = lfeFreq
        self.packFormat = packFormat

    @property
    def features(self):
        f = 0
        if self.cartesian:
            f |= CARTESIAN
        if self.screenEdgeLockH is not None or self.screenEdgeLockV is not None:
            f |= SCREEN_EDGE_LOCK
        if ":SC=" in self.layout:
            f |= WIDE_SCREEN_SPEAKERS
        return f

    def type_metadata(self):
        screen_edge_lock = ScreenEdgeLock(horizontal=self.screenEdgeLockH, vertical=self.screenEdgeLockV)
        bounds = [BoundCoordinate(v, lo, hi) for v, lo, hi in zip(self.position, self.minimum, self.maximum)]
        if self.cartesian:
            position = DirectSpeakerCartesianPosition(bounded_X=bounds[0], bounded_Y=bounds[1], bounded_Z=bounds[2],
                                                      screenEdgeLock=screen_edge_lock)
        else:
            position = DirectSpeakerPolarPosition(bounded_azimuth=bounds[0], bounded_elevation=bounds[1],
                                                  bounded_distance=bounds[2], screenEdgeLock=screen_edge_lock)
        bf = AudioBlockFormatDirectSpeakers(position=position, speakerLabel=self.labels)
        packs = [adm_common_defs[self.packFormat]] if self.packFormat else None
        return DirectSpeakersTypeMetadata(
            block_format=bf,
            extra_data=ExtraData(channel_frequency=Frequency(lowPass=120.0 if self.lfeFreq else None)),
            audioPackFormats=packs)


def case_from_common_definitions(name, layout, apf_id, acf_id):
    """A case using the block format of a common-definitions channel within
    a common-definitions pack, as the mapping rules require."""
    acf = adm_common_defs[acf_id]
    bf = acf.audioBlockFormats[0]
    pos = bf.position
    assert isinstance(pos, DirectSpeakerPolarPosition)
    bounds = [pos.bounded_azimuth, pos.bounded_elevation, pos.bounded_distance]
    return Case(name, layout, labels=list(bf.speakerLabel),
                position=[b.value for b in bounds],
                minimum=[b.min for b in bounds], maximum=[b.max for b in bounds],
                packFormat=apf_id)


ALL_LABELS = ["M+000", "M+030", "M-030", "M+060", "M-060", "M+090", "M-090", "M+110", "M-110",
              "M+135", "M-135", "M+180", "M+SC", "M-SC", "U+000", "U+030", "U-030", "U+045", "U-045",
              "U+090", "U-090", "U+110", "U-110", "U+135", "U-135", "U+180", "UH+180", "T+000",
              "B+000", "B+045", "B-045", "LFE1", "LFE2", "LFE", "LFEL", "LFER", "X+000"]


def build_cases():
    cases = []
    layouts = ["0+2+0", "0+5+0", "4+5+0", "4+7+0", "4+9+0", "9+10+3", "4+9+0:SC=40"]

    # --- speaker labels, plain and as URNs, with and without LFE frequency
    for layout in layouts:
        for label in ALL_LABELS:
            for lfe in [False, True]:
                cases.append(Case("label", layout, labels=[label], lfeFreq=lfe))
            cases.append(Case("label_urn", layout, labels=["urn:itu:bs:2051:0:speaker:" + label]))
        cases.append(Case("label_multiple", layout, labels=["X+000", "M+030", "M+000"]))
        cases.append(Case("label_multiple", layout, labels=["M+000", "M+030"]))

    # --- polar positions, with and without bounds
    polar_positions = [(0.0, 0.0), (15.0, 0.0), (14.0, 0.0), (30.0, 0.0), (-30.0, 0.0), (0.0, 15.0), (0.0, 14.0),
                       (45.0, 30.0), (110.0, 0.0), (135.0, 0.0), (180.0, 0.0), (0.0, 90.0), (15.0, 90.0),
                       (-100.0, -20.0), (60.0, 45.0)]
    for layout in layouts:
        for az, el in polar_positions:
            cases.append(Case("polar", layout, position=(az, el, 1.0)))
            cases.append(Case("polar_lfe", layout, position=(az, el, 1.0), lfeFreq=True))
        for (az, el), (amin, amax), (emin, emax) in product(
                [(15.0, 0.0), (14.0, 0.0), (0.0, 15.0), (0.0, 14.0), (100.0, 10.0), (15.0, 90.0)],
                [(None, None), (0.0, None), (None, 30.0), (0.0, 30.0), (-45.0, 10.0), (10.0, 20.0)],
                [(None, None), (0.0, 30.0), (-10.0, 10.0)]):
            cases.append(Case("polar_bounds", layout, position=(az, el, 1.0),
                              minimum=(amin, emin, None), maximum=(amax, emax, None)))
        cases.append(Case("polar_distance_bounds", layout, position=(30.0, 0.0, 0.5),
                          minimum=(None, None, 0.0), maximum=(None, None, 1.0)))

    # --- cartesian positions, with and without bounds
    for layout in layouts:
        for x, y, z in product([-1.0, -0.55, -0.5, -0.45, 0.0, 0.5, 1.0], [-1.0, 0.0, 0.1, 1.0], [-1.0, 0.0, 1.0]):
            cases.append(Case("cartesian", layout, cartesian=True, position=(x, y, z)))
        cases.append(Case("cartesian_lfe", layout, cartesian=True, position=(-1.0, 1.0, -1.0), lfeFreq=True))
        for (x, y, z), mn, mx in [
                ((1.0, 0.1, 0.0), (None, 0.0, None), (None, None, None)),
                ((1.0, -0.1, 0.0), (None, None, None), (None, 0.0, None)),
                ((-0.45, 1.0, 0.0), (-1.0, None, None), (1.0, None, None)),
                ((-0.55, 1.0, 0.0), (-1.0, None, None), (1.0, None, None)),
                ((0.0, 1.0, 0.0), (-0.1, 0.5, None), (1.0, None, None)),
                ((0.3, 0.3, 0.3), (-1.0, -1.0, -1.0), (1.0, 1.0, 1.0))]:
            cases.append(Case("cartesian_bounds", layout, cartesian=True, position=(x, y, z), minimum=mn, maximum=mx))

    # --- screen edge lock
    locks = [("left", None), ("right", None), (None, "top"), (None, "bottom"), ("left", "top"), ("right", "bottom")]
    for layout in ["0+5+0", "4+5+0", "9+10+3", "4+9+0"]:
        for h, v in locks:
            cases.append(Case("screenedgelock", layout, position=(-30.0, 0.0, 1.0), screenEdgeLockH=h, screenEdgeLockV=v))
            cases.append(Case("screenedgelock_bounds", layout, position=(0.0, 0.0, 1.0),
                              minimum=(-45.0, -10.0, None), maximum=(10.0, 10.0, None),
                              screenEdgeLockH=h, screenEdgeLockV=v))
            cases.append(Case("screenedgelock_cartesian", layout, cartesian=True, position=(0.0, 1.0, 0.0),
                              screenEdgeLockH=h, screenEdgeLockV=v))
            cases.append(Case("screenedgelock_cartesian_bounds", layout, cartesian=True, position=(0.0, 1.0, 0.0),
                              minimum=(-0.1, 0.5, None), maximum=(1.0, None, None),
                              screenEdgeLockH=h, screenEdgeLockV=v))

    # --- common definitions packs: every channel of a few packs, rendered to several layouts
    packs = {
        "AP_00010002": ["AC_00010001", "AC_00010002"],                                        # stereo
        "AP_00010003": ["AC_00010001", "AC_00010002", "AC_00010003", "AC_00010004", "AC_00010005", "AC_00010006"],  # 5.1
        "AP_00010009": ["AC_00010011", "AC_0001000a"],                                        # 22.2 subset
        "AP_0001000f": ["AC_0001001c"],
        "AP_00010017": ["AC_0001000a"],
    }
    for layout in ["0+2+0", "0+5+0", "4+5+0", "4+7+0", "9+10+3"]:
        for apf, acfs in packs.items():
            for acf in acfs:
                cases.append(case_from_common_definitions("pack_" + apf, layout, apf, acf))

    return cases


def main():
    cases = build_cases()
    panners = {}

    labels = []
    gains_values = []
    case_rows = []
    for case in cases:
        if case.layout not in panners:
            panners[case.layout] = DirectSpeakersPanner(layout_for_spec(case.layout))
        gains = panners[case.layout].handle(case.type_metadata())

        labels_offset = len(labels)
        labels.extend(case.labels)
        gains_offset = len(gains_values)
        gains_values.extend(float(v) for v in gains)

        def opt(values):
            return "{%s}, {%s}" % (", ".join("true" if v is not None else "false" for v in values),
                                   ", ".join(fmt(v if v is not None else 0.0) for v in values))

        case_rows.append("{\"%s\", \"%s\", %du, %d, %d, %s, {%s}, %s, %s, \"%s\", \"%s\", %s, \"%s\", %d, %d}" % (
            case.name, case.layout, case.features,
            labels_offset, len(case.labels),
            "true" if case.cartesian else "false",
            ", ".join(fmt(v) for v in case.position),
            opt(case.minimum), opt(case.maximum),
            case.screenEdgeLockH or "", case.screenEdgeLockV or "",
            "true" if case.lfeFreq else "false",
            case.packFormat or "",
            gains_offset, len(gains)))

    lines = []
    lines.append("// GENERATED by tools/reference/generate_direct_speakers_reference.py -- do not edit")
    lines.append("// reference implementation: ear %s (https://github.com/ebu/ebu_adm_renderer)" % EAR_VERSION)
    lines.append("//")
    lines.append("// The data is stored in plain arrays and converted at runtime, as large")
    lines.append("// braced-init-lists of non-trivial types compile extremely slowly.")
    lines.append('#include "direct_speakers_reference.hpp"')
    lines.append("")
    lines.append("namespace ear {")
    lines.append("  namespace reference {")
    lines.append("    namespace {")
    lines.append("      struct RawCase {")
    lines.append("        const char* name; const char* layout; unsigned features; int labelsOffset; int labelsCount;")
    lines.append("        bool cartesian; double position[3]; bool hasMin[3]; double minimum[3]; bool hasMax[3]; double maximum[3];")
    lines.append("        const char* screenEdgeLockH; const char* screenEdgeLockV; bool lfeFreq; const char* packFormat;")
    lines.append("        int gainsOffset; int numChannels;")
    lines.append("      };")
    lines.append("")
    lines.append("      const char* const k_labels[] = {")
    for i in range(0, len(labels), 8):
        lines.append("        %s," % ", ".join('"%s"' % l for l in labels[i:i + 8]))
    lines.append('        "",  // sentinel')
    lines.append("      };")
    lines.append("")
    lines.append("      const double k_gains[] = {")
    for i in range(0, len(gains_values), 8):
        lines.append("        %s," % ", ".join(fmt(v) for v in gains_values[i:i + 8]))
    lines.append("      };")
    lines.append("")
    lines.append("      const RawCase k_cases[] = {")
    for row in case_rows:
        lines.append("        %s," % row)
    lines.append("      };")
    lines.append("    }  // namespace")
    lines.append("")
    lines.append("    const std::vector<DirectSpeakersCase>& directSpeakersCases() {")
    lines.append("      static const std::vector<DirectSpeakersCase> cases = [] {")
    lines.append("        std::vector<DirectSpeakersCase> ret;")
    lines.append("        for (const RawCase& r : k_cases) {")
    lines.append("          DirectSpeakersCase c;")
    lines.append("          c.name = r.name; c.layout = r.layout; c.features = r.features;")
    lines.append("          for (int i = 0; i < r.labelsCount; i++) c.speakerLabels.push_back(k_labels[r.labelsOffset + i]);")
    lines.append("          c.cartesian = r.cartesian;")
    lines.append("          for (int i = 0; i < 3; i++) {")
    lines.append("            c.position[i] = r.position[i]; c.hasMin[i] = r.hasMin[i]; c.minimum[i] = r.minimum[i];")
    lines.append("            c.hasMax[i] = r.hasMax[i]; c.maximum[i] = r.maximum[i];")
    lines.append("          }")
    lines.append("          c.screenEdgeLockHorizontal = r.screenEdgeLockH; c.screenEdgeLockVertical = r.screenEdgeLockV;")
    lines.append("          c.lfeFreq = r.lfeFreq; c.packFormat = r.packFormat;")
    lines.append("          c.gains.assign(k_gains + r.gainsOffset, k_gains + r.gainsOffset + r.numChannels);")
    lines.append("          ret.push_back(std::move(c));")
    lines.append("        }")
    lines.append("        return ret;")
    lines.append("      }();")
    lines.append("      return cases;")
    lines.append("    }")
    lines.append("")
    lines.append("  }  // namespace reference")
    lines.append("}  // namespace ear")
    lines.append("")

    with open(OUT, "w") as f:
        f.write("\n".join(lines))
    print("wrote %d cases to %s" % (len(cases), os.path.normpath(OUT)))


if __name__ == "__main__":
    sys.exit(main())
