from __future__ import annotations

import argparse
from pathlib import Path

from load_response_mat import load_response_mat


DEFAULT_RESPONSE_FILES = (
    "Bilayer graphene/simple_strained_BLG_photocurrent_objects_AB_xi1.mat",
    "Twisted bilayer graphene (TBG)/tbg_photocurrent_objects.mat",
    "WTe2/WTe2_ShiSong_6band_deltaZ_response_objects.mat",
    "MBT/MBT_response_objects.mat",
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate response-object .mat files.")
    parser.add_argument("paths", nargs="*", type=Path, help="Response-object .mat files to validate.")
    args = parser.parse_args()

    paths = args.paths or [Path(path) for path in DEFAULT_RESPONSE_FILES]
    ok = True

    for path in paths:
        if not path.exists():
            print(f"SKIP missing: {path}")
            ok = False
            continue

        try:
            data = load_response_mat(path)
        except Exception as exc:
            print(f"FAIL {path}: {exc}")
            ok = False
            continue

        print(f"OK {path}: {data['material_name']} / {data['model_name']}")

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
