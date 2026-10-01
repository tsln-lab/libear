#!/usr/bin/env python3
"""Generate reference gains for the Objects renderer from the EBU ADM Renderer.

This drives the reference Python implementation (the ``ear`` package,
https://github.com/ebu/ebu_adm_renderer) over a fixed set of layouts and
block formats and writes the resulting direct/diffuse gains as a C++ table,
``tests/reference/objects_reference_data.cpp``. The libear test
``objects_reference_tests`` compares libear against this table.

Usage (from the repository root, with ``ear`` installed, e.g. in a venv)::

    python tools/reference/generate_objects_reference.py

The case list is deterministic, so re-running the script only changes the
output when the reference implementation or the case list changes.
"""
import os
import sys
import warnings
from itertools import product

import numpy as np

warnings.filterwarnings("ignore")

from ear import __name__ as _ear_name  # noqa: E402
from ear.core import bs2051  # noqa: E402
from ear.core.metadata_input import ObjectTypeMetadata  # noqa: E402
from ear.core.objectbased.gain_calc import GainCalc  # noqa: E402
from ear.core.geom import PolarPosition  # noqa: E402
from ear.fileio.adm.elements import (  # noqa: E402
    AudioBlockFormatObjects, CartesianZone, ChannelLock, ObjectCartesianPosition,
    ObjectDivergence, ObjectPolarPosition, PolarZone, ScreenEdgeLock)
from attr import evolve  # noqa: E402

try:
    from importlib.metadata import version
    EAR_VERSION = version("ear")
except Exception:  # pragma: no cover
    EAR_VERSION = "unknown"

# feature flags, must match tests/reference/objects_reference.hpp
CARTESIAN = 1
CHANNEL_LOCK = 2
DIVERGENCE = 4
ZONE_EXCLUSION = 8
SCREEN_REF = 16
SCREEN_EDGE_LOCK = 32
WIDE_SCREEN_SPEAKERS = 64


def layout_for_spec(spec):
    """A BS.2051 layout by name; "<name>:SC=<az>" moves the M+SC/M-SC
    loudspeakers of the layout to +-az degrees."""
    if ":SC=" in spec:
        name, az = spec.split(":SC=")
        az = float(az)
        layout = bs2051.get_layout(name)
        channels = []
        for channel in layout.channels:
            if channel.name == "M+SC":
                channel = evolve(channel, polar_position=PolarPosition(az, 0.0, 1.0))
            elif channel.name == "M-SC":
                channel = evolve(channel, polar_position=PolarPosition(-az, 0.0, 1.0))
            channels.append(channel)
        return evolve(layout, channels=channels)
    return bs2051.get_layout(spec)

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "reference", "objects_reference_data.cpp")


class Case(object):
    def __init__(self, name, layout, cartesian=False, position=(0.0, 0.0, 1.0),
                 width=0.0, height=0.0, depth=0.0, gain=1.0, diffuse=0.0,
                 channelLock=False, channelLockMaxDistance=None,
                 divergence=0.0, divergenceRange=None, divergenceCartesian=False,
                 zones=(), screenRef=False, screenEdgeLockH=None, screenEdgeLockV=None):
        self.name = name
        self.layout = layout
        self.cartesian = cartesian
        self.position = tuple(float(x) for x in position)
        self.width, self.height, self.depth = float(width), float(height), float(depth)
        self.gain, self.diffuse = float(gain), float(diffuse)
        self.channelLock = channelLock
        self.channelLockMaxDistance = channelLockMaxDistance
        self.divergence = float(divergence)
        self.divergenceRange = divergenceRange
        self.divergenceCartesian = divergenceCartesian
        self.zones = list(zones)
        self.screenRef = screenRef
        self.screenEdgeLockH = screenEdgeLockH
        self.screenEdgeLockV = screenEdgeLockV

    @property
    def features(self):
        f = 0
        if self.cartesian:
            f |= CARTESIAN
        if self.channelLock:
            f |= CHANNEL_LOCK
        if self.divergence != 0.0:
            f |= DIVERGENCE
        if self.zones:
            f |= ZONE_EXCLUSION
        if self.screenRef:
            f |= SCREEN_REF
        if self.screenEdgeLockH is not None or self.screenEdgeLockV is not None:
            f |= SCREEN_EDGE_LOCK
        if ":SC=" in self.layout:
            f |= WIDE_SCREEN_SPEAKERS
        return f

    @property
    def uses_extent(self):
        return self.width != 0.0 or self.height != 0.0 or self.depth != 0.0

    def block_format(self):
        screen_edge_lock = ScreenEdgeLock(horizontal=self.screenEdgeLockH, vertical=self.screenEdgeLockV)
        if self.cartesian:
            position = ObjectCartesianPosition(*self.position, screenEdgeLock=screen_edge_lock)
        else:
            position = ObjectPolarPosition(*self.position, screenEdgeLock=screen_edge_lock)

        kwargs = dict(position=position, cartesian=self.cartesian,
                      width=self.width, height=self.height, depth=self.depth,
                      gain=self.gain, diffuse=self.diffuse, screenRef=self.screenRef)

        if self.channelLock:
            kwargs["channelLock"] = ChannelLock(maxDistance=self.channelLockMaxDistance)

        if self.divergence != 0.0:
            if self.divergenceCartesian:
                kwargs["objectDivergence"] = ObjectDivergence(self.divergence, positionRange=self.divergenceRange)
            else:
                kwargs["objectDivergence"] = ObjectDivergence(self.divergence, azimuthRange=self.divergenceRange)

        zones = []
        for zone in self.zones:
            if zone[0] == "polar":
                zones.append(PolarZone(minAzimuth=zone[1], maxAzimuth=zone[2],
                                       minElevation=zone[3], maxElevation=zone[4]))
            else:
                zones.append(CartesianZone(minX=zone[1], maxX=zone[2], minY=zone[3], maxY=zone[4],
                                           minZ=zone[5], maxZ=zone[6]))
        kwargs["zoneExclusion"] = zones

        return AudioBlockFormatObjects(**kwargs)


def build_cases():
    cases = []

    # --- plain point source panning over a coarse grid on every layout
    coarse_az = [-110.0, -30.0, 0.0, 30.0, 110.0]
    coarse_el = [0.0, 30.0]
    for layout in ["0+2+0", "0+5+0", "2+5+0", "4+5+0", "4+7+0", "9+10+3"]:
        for az, el in product(coarse_az, coarse_el):
            cases.append(Case("point", layout, position=(az, el, 1.0)))

    # --- finer grid on two representative layouts
    fine_az = [-180.0, -150.0, -110.0, -90.0, -60.0, -45.0, -30.0, -15.0,
               0.0, 15.0, 30.0, 45.0, 60.0, 90.0, 110.0, 150.0]
    fine_el = [-30.0, 0.0, 30.0, 60.0, 90.0]
    for layout in ["4+5+0", "9+10+3"]:
        for az, el in product(fine_az, fine_el):
            cases.append(Case("point", layout, position=(az, el, 1.0)))

    # --- gain and diffuse
    for layout in ["0+5+0", "4+5+0"]:
        cases.append(Case("gain", layout, position=(20.0, 10.0, 1.0), gain=0.5))
        cases.append(Case("diffuse", layout, position=(20.0, 10.0, 1.0), diffuse=0.5))
        cases.append(Case("diffuse", layout, position=(20.0, 10.0, 1.0), diffuse=1.0, gain=0.25))

    # --- extent
    extents = [(30.0, 0.0, 0.0), (0.0, 30.0, 0.0), (90.0, 45.0, 0.0),
               (0.0, 0.0, 0.5), (360.0, 360.0, 0.0), (60.0, 20.0, 0.3)]
    for layout in ["0+5+0", "4+5+0", "9+10+3"]:
        for (az, el), (w, h, d), dist in product([(0.0, 0.0), (30.0, 0.0), (90.0, 30.0), (180.0, -20.0)],
                                                  extents, [0.5, 1.0, 2.0]):
            cases.append(Case("extent", layout, position=(az, el, dist), width=w, height=h, depth=d))

    # --- channel lock
    for layout in ["4+5+0", "9+10+3"]:
        for az, el in product([-150.0, -100.0, -40.0, -15.0, 0.0, 14.0, 15.0, 30.0, 45.0, 100.0, 180.0],
                              [-20.0, 0.0, 15.0, 30.0, 50.0, 90.0]):
            for max_distance in [None, 0.5, 0.01]:
                cases.append(Case("channellock", layout, position=(az, el, 1.0),
                                  channelLock=True, channelLockMaxDistance=max_distance))
        # channel lock combined with extent and gain
        cases.append(Case("channellock_extent", layout, position=(20.0, 10.0, 1.0),
                          channelLock=True, width=45.0, height=20.0, gain=0.7))

    # --- divergence
    for layout in ["0+5+0", "4+5+0", "9+10+3"]:
        for (az, el), value, rng in product([(0.0, 0.0), (30.0, 0.0), (70.0, 20.0), (180.0, 0.0), (0.0, 45.0)],
                                            [0.3, 0.5, 1.0], [10.0, 30.0, 45.0, 90.0]):
            cases.append(Case("divergence", layout, position=(az, el, 1.0),
                              divergence=value, divergenceRange=rng))
        # default azimuthRange (45)
        cases.append(Case("divergence_default_range", layout, position=(30.0, 0.0, 1.0), divergence=0.5))
        # divergence with extent, distance and diffuse
        cases.append(Case("divergence_extent", layout, position=(0.0, 0.0, 1.0),
                          divergence=0.5, divergenceRange=30.0, width=40.0, height=10.0, diffuse=0.3))
        cases.append(Case("divergence_distance", layout, position=(45.0, 10.0, 0.6),
                          divergence=0.7, divergenceRange=20.0, depth=0.2))
        # cartesian divergence given for a polar block format: reference falls back to polar divergence
        cases.append(Case("divergence_cartesian_in_polar", layout, position=(30.0, 0.0, 1.0),
                          divergence=0.5, divergenceRange=0.3, divergenceCartesian=True))
        # channel lock + divergence
        cases.append(Case("channellock_divergence", layout, position=(10.0, 5.0, 1.0),
                          channelLock=True, divergence=0.5, divergenceRange=30.0))

    # --- cartesian (allocentric) panning
    for layout in ["0+5+0", "4+5+0", "9+10+3"]:
        for x, y, z in product([-1.0, -0.5, 0.0, 0.5, 1.0], [-1.0, -0.3, 0.0, 0.5, 1.0], [-1.0, 0.0, 0.5, 1.0]):
            cases.append(Case("cartesian", layout, cartesian=True, position=(x, y, z)))
        for (x, y, z), (w, h, d) in product([(0.0, 1.0, 0.0), (0.5, 0.5, 0.0), (-0.3, 0.2, 0.6)],
                                            [(0.3, 0.0, 0.0), (0.0, 0.0, 0.5), (0.5, 0.5, 0.5), (1.0, 1.0, 1.0)]):
            cases.append(Case("cartesian_extent", layout, cartesian=True, position=(x, y, z), width=w, height=h, depth=d))
        for x in [-1.0, -0.6, -0.5, 0.0, 0.1, 0.5, 0.6, 1.0]:
            for max_distance in [None, 0.01]:
                cases.append(Case("cartesian_channellock", layout, cartesian=True, position=(x, 1.0, 0.0),
                                  channelLock=True, channelLockMaxDistance=max_distance))
        cases.append(Case("cartesian_divergence", layout, cartesian=True, position=(0.0, 1.0, 0.0),
                          divergence=0.5, divergenceRange=1.0, divergenceCartesian=True))
        cases.append(Case("cartesian_divergence", layout, cartesian=True, position=(0.2, 0.8, 0.0),
                          divergence=1.0, divergenceRange=0.5, divergenceCartesian=True))
        cases.append(Case("cartesian_divergence_polar_in_cartesian", layout, cartesian=True, position=(0.0, 1.0, 0.0),
                          divergence=0.5, divergenceRange=30.0, divergenceCartesian=False))

    # --- zone exclusion
    zones_polar = [
        [("polar", 0.0, 0.0, 0.0, 0.0)],                 # centre speaker only
        [("polar", -180.0, 180.0, 0.0, 0.0)],            # whole middle layer
        [("polar", -180.0, 180.0, 10.0, 90.0)],          # everything above
        [("polar", -30.0, 30.0, -90.0, 90.0)],           # front
        [("polar", 100.0, -100.0, -90.0, 90.0)],         # rear (wrapping range)
        [("polar", 0.0, 0.0, 0.0, 0.0), ("polar", -180.0, 180.0, 10.0, 90.0)],
    ]
    zones_cart = [
        [("cartesian", -1.0, 1.0, -1.0, -0.5, -1.0, 1.0)],   # rear
        [("cartesian", -0.1, 0.1, 0.9, 1.0, -0.1, 0.1)],     # centre
    ]
    for layout in ["0+5+0", "4+5+0", "9+10+3"]:
        for zones in zones_polar + zones_cart:
            for az, el in [(0.0, 0.0), (30.0, 0.0), (-110.0, 0.0), (0.0, 30.0), (180.0, 45.0)]:
                cases.append(Case("zone", layout, position=(az, el, 1.0), zones=zones))
            cases.append(Case("zone_extent", layout, position=(20.0, 10.0, 1.0), zones=zones, width=60.0))
            cases.append(Case("zone_cartesian", layout, cartesian=True, position=(0.3, 0.8, 0.0), zones=zones))

    # --- screen reference (default reference screen, default layout screen)
    for layout in ["0+5+0", "4+5+0", "9+10+3"]:
        for az, el in [(0.0, 0.0), (20.0, 0.0), (-29.0, 10.0), (45.0, 20.0), (120.0, 0.0)]:
            cases.append(Case("screenref", layout, position=(az, el, 1.0), screenRef=True))
        for x, y, z in [(0.0, 1.0, 0.0), (0.5, 1.0, 0.2), (-0.3, 0.9, 0.0)]:
            cases.append(Case("screenref_cartesian", layout, cartesian=True, position=(x, y, z), screenRef=True))

    # --- screen edge lock (default layout screen)
    locks = [("left", None), ("right", None), (None, "top"), (None, "bottom"),
             ("left", "top"), ("right", "bottom")]
    for layout in ["0+5+0", "4+5+0", "9+10+3"]:
        for (h, v), (az, el) in product(locks, [(0.0, 0.0), (45.0, 20.0), (-120.0, 0.0)]):
            cases.append(Case("screenedgelock", layout, position=(az, el, 1.0),
                              screenEdgeLockH=h, screenEdgeLockV=v))
        for (h, v) in locks:
            cases.append(Case("screenedgelock_cartesian", layout, cartesian=True, position=(0.3, 0.8, 0.1),
                              screenEdgeLockH=h, screenEdgeLockV=v))
        cases.append(Case("screenedgelock_extent", layout, position=(10.0, 5.0, 1.0),
                          screenEdgeLockH="left", width=30.0, height=20.0))
        cases.append(Case("screenedgelock_screenref", layout, position=(10.0, 5.0, 1.0),
                          screenEdgeLockV="top", screenRef=True))

    # --- layouts with screen loudspeakers (4+9+0 has M+SC and M-SC at +-15
    # by default); wider than 30 degrees changes their nominal positions
    for layout in ["4+9+0", "4+9+0:SC=20", "4+9+0:SC=40", "4+9+0:SC=55"]:
        for az, el in product([-60.0, -45.0, -40.0, -20.0, -15.0, -10.0, 0.0, 10.0, 15.0, 20.0, 40.0, 45.0, 60.0, 110.0],
                              [0.0, 15.0, 30.0]):
            cases.append(Case("screenspeakers", layout, position=(az, el, 1.0)))
        cases.append(Case("screenspeakers_extent", layout, position=(20.0, 0.0, 1.0), width=45.0))
        cases.append(Case("screenspeakers_channellock", layout, position=(18.0, 0.0, 1.0), channelLock=True))
        for x in [-1.0, -0.5, -0.4, 0.0, 0.4, 0.5, 1.0]:
            cases.append(Case("screenspeakers_cartesian", layout, cartesian=True, position=(x, 1.0, 0.0)))
        cases.append(Case("screenspeakers_zone", layout, position=(0.0, 0.0, 1.0),
                          zones=[("polar", 0.0, 0.0, 0.0, 0.0)]))

    return cases


def fmt(x):
    return repr(float(x))


def main():
    cases = build_cases()
    calcs = {}

    zones_rows = []
    gains_values = []
    case_rows = []
    for case in cases:
        if case.layout not in calcs:
            calcs[case.layout] = GainCalc(layout_for_spec(case.layout))
        gains = calcs[case.layout].render(ObjectTypeMetadata(block_format=case.block_format()))

        zones_offset = len(zones_rows)
        for z in case.zones:
            values = (list(z[1:]) + [0.0, 0.0])[:6]
            zones_rows.append("{%s, {%s}}" % ("true" if z[0] == "cartesian" else "false",
                                              ", ".join(fmt(v) for v in values)))

        gains_offset = len(gains_values)
        gains_values.extend(float(v) for v in gains.direct)
        gains_values.extend(float(v) for v in gains.diffuse)

        case_rows.append("{\"%s\", \"%s\", %du, %s, {%s}, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %d, %d, %s, \"%s\", \"%s\", %s, %d, %d}" % (
            case.name, case.layout, case.features,
            "true" if case.cartesian else "false",
            ", ".join(fmt(v) for v in case.position),
            fmt(case.width), fmt(case.height), fmt(case.depth), fmt(case.gain), fmt(case.diffuse),
            "true" if case.channelLock else "false",
            fmt(case.channelLockMaxDistance if case.channelLockMaxDistance is not None else -1.0),
            fmt(case.divergence),
            fmt(case.divergenceRange if case.divergenceRange is not None else -1.0),
            "true" if case.divergenceCartesian else "false",
            zones_offset, len(case.zones),
            "true" if case.screenRef else "false",
            case.screenEdgeLockH or "", case.screenEdgeLockV or "",
            "true" if case.uses_extent else "false",
            gains_offset, len(gains.direct)))

    lines = []
    lines.append("// GENERATED by tools/reference/generate_objects_reference.py -- do not edit")
    lines.append("// reference implementation: ear %s (https://github.com/ebu/ebu_adm_renderer)" % EAR_VERSION)
    lines.append("//")
    lines.append("// The data is stored in plain arrays and converted at runtime, as large")
    lines.append("// braced-init-lists of non-trivial types compile extremely slowly.")
    lines.append('#include "objects_reference.hpp"')
    lines.append("")
    lines.append("namespace ear {")
    lines.append("  namespace reference {")
    lines.append("    namespace {")
    lines.append("      struct RawCase {")
    lines.append("        const char* name; const char* layout; unsigned features; bool cartesian;")
    lines.append("        double position[3]; double width; double height; double depth; double gain; double diffuse;")
    lines.append("        bool channelLock; double channelLockMaxDistance; double divergence; double divergenceRange;")
    lines.append("        bool divergenceCartesian; int zonesOffset; int zonesCount; bool screenRef;")
    lines.append("        const char* screenEdgeLockH; const char* screenEdgeLockV; bool usesExtent;")
    lines.append("        int gainsOffset; int numChannels;")
    lines.append("      };")
    lines.append("")
    lines.append("      const Zone k_zones[] = {")
    for row in zones_rows:
        lines.append("        %s," % row)
    lines.append("        {false, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},  // sentinel")
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
    lines.append("    const std::vector<ObjectsCase>& objectsCases() {")
    lines.append("      static const std::vector<ObjectsCase> cases = [] {")
    lines.append("        std::vector<ObjectsCase> ret;")
    lines.append("        for (const RawCase& r : k_cases) {")
    lines.append("          ObjectsCase c;")
    lines.append("          c.name = r.name; c.layout = r.layout; c.features = r.features; c.cartesian = r.cartesian;")
    lines.append("          for (int i = 0; i < 3; i++) c.position[i] = r.position[i];")
    lines.append("          c.width = r.width; c.height = r.height; c.depth = r.depth; c.gain = r.gain; c.diffuse = r.diffuse;")
    lines.append("          c.channelLock = r.channelLock; c.channelLockMaxDistance = r.channelLockMaxDistance;")
    lines.append("          c.divergence = r.divergence; c.divergenceRange = r.divergenceRange;")
    lines.append("          c.divergenceCartesian = r.divergenceCartesian;")
    lines.append("          c.zones.assign(k_zones + r.zonesOffset, k_zones + r.zonesOffset + r.zonesCount);")
    lines.append("          c.screenRef = r.screenRef;")
    lines.append("          c.screenEdgeLockHorizontal = r.screenEdgeLockH; c.screenEdgeLockVertical = r.screenEdgeLockV;")
    lines.append("          c.usesExtent = r.usesExtent;")
    lines.append("          c.direct.assign(k_gains + r.gainsOffset, k_gains + r.gainsOffset + r.numChannels);")
    lines.append("          c.diffuseGains.assign(k_gains + r.gainsOffset + r.numChannels, k_gains + r.gainsOffset + 2 * r.numChannels);")
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
