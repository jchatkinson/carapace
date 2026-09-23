"""Parses protocol.txt the same way core/examples/cyclic_material.rs does,
so both engines walk through the identical displacement history."""

import pathlib

PROTOCOL_PATH = pathlib.Path(__file__).parent / "protocol.txt"


def load_protocol(path=PROTOCOL_PATH):
    """Returns (tiers, step) where tiers is a list of (amplitude, cycles)."""
    step = None
    tiers = []
    for line in pathlib.Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if parts[0] == "step":
            step = float(parts[1])
        else:
            tiers.append((float(parts[0]), int(parts[1])))
    if step is None:
        raise ValueError("protocol.txt must set `step`")
    return tiers, step


def cyclic_targets(tiers):
    """Expands tiers into the sequence of absolute displacement targets a
    reverse-cyclic protocol visits: for each (amplitude, cycles) tier,
    +amplitude/-amplitude repeated `cycles` times, then back to 0."""
    targets = []
    for amplitude, cycles in tiers:
        for _ in range(cycles):
            targets.append(amplitude)
            targets.append(-amplitude)
        targets.append(0.0)
    return targets
