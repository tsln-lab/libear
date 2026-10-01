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
from ear.fileio.adm.elements import (  # noqa: E402
    AudioBlockFormatObjects, CartesianZone, ChannelLock, ObjectCartesianPosition,
    ObjectDivergence, ObjectPolarPosition, PolarZone)

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

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "reference", "objects_reference_data.cpp")


class Case(object):
    def __init__(self, name, layout, cartesian=False, position=(0.0, 0.0, 1.0),
                 width=0.0, height=0.0, depth=0.0, gain=1.0, diffuse=0.0,
                 channelLock=False, channelLockMaxDistance=None,
                 divergence=0.0, divergenceRange=None, divergenceCartesian=False,
                 zones=(), screenRef=False):
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
        return f

    @property
    def uses_extent(self):
        return self.width != 0.0 or self.height != 0.0 or self.depth != 0.0

    def block_format(self):
        if self.cartesian:
            position = ObjectCartesianPosition(*self.position)
        else:
            position = ObjectPolarPosition(*self.position)

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

    return cases


def fmt(x):
    return repr(float(x))


def main():
    cases = build_cases()
    calcs = {}
    lines = []
    lines.append("// GENERATED by tools/reference/generate_objects_reference.py -- do not edit")
    lines.append("// reference implementation: ear %s (https://github.com/ebu/ebu_adm_renderer)" % EAR_VERSION)
    lines.append('#include "objects_reference.hpp"')
    lines.append("")
    lines.append("namespace ear {")
    lines.append("  namespace reference {")
    lines.append("")
    lines.append("    const std::vector<ObjectsCase>& objectsCases() {")
    lines.append("      static const std::vector<ObjectsCase> cases = {")

    for case in cases:
        if case.layout not in calcs:
            calcs[case.layout] = GainCalc(bs2051.get_layout(case.layout))
        gains = calcs[case.layout].render(ObjectTypeMetadata(block_format=case.block_format()))

        zones = ", ".join(
            "{%s, %s}" % ("true" if z[0] == "cartesian" else "false",
                          ", ".join(fmt(v) for v in (list(z[1:]) + [0.0, 0.0])[:6]))
            for z in case.zones)

        lines.append("        {\"%s\", \"%s\", %du, %s, {%s}, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, {%s}, %s, %s," % (
            case.name, case.layout, case.features,
            "true" if case.cartesian else "false",
            ", ".join(fmt(v) for v in case.position),
            fmt(case.width), fmt(case.height), fmt(case.depth), fmt(case.gain), fmt(case.diffuse),
            "true" if case.channelLock else "false",
            fmt(case.channelLockMaxDistance if case.channelLockMaxDistance is not None else -1.0),
            fmt(case.divergence),
            fmt(case.divergenceRange if case.divergenceRange is not None else -1.0),
            "true" if case.divergenceCartesian else "false",
            zones,
            "true" if case.screenRef else "false",
            "true" if case.uses_extent else "false"))
        lines.append("         {%s}," % ", ".join(fmt(v) for v in gains.direct))
        lines.append("         {%s}}," % ", ".join(fmt(v) for v in gains.diffuse))

    lines.append("      };")
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
