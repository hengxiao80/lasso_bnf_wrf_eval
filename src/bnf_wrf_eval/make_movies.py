"""Stitch the batch-rendered PNG frames into per-run / per-case MP4 movies.

Pairs with `bnf_wrf_eval.batch_plot`, which writes one PNG per frame under
`<output-root>/frames/<type>/`. This groups those by everything in the
filename before the trailing `_<YYYYMMDD_HHMM>Z.png` timestamp, sorts each
group by that timestamp, and encodes it to
`<movies-dir>/<that-prefix>.mp4`.

    python -m bnf_wrf_eval.make_movies --which comparison
    python -m bnf_wrf_eval.make_movies --which single --fps 12

Frame rate: 12 fps, except 2 fps for the hourly types (`HOURLY_TYPES`); an
explicit `--fps` applies to every movie.

`--which` takes `comparison` (tb4, tb4_simple, refl3), `single` (the six single-panel
types), `all`, or a comma-separated list of type names. `--task-index` /
`--task-count` (or the SLURM job-array env vars) slice the movie list
`movies[idx::cnt]`. An existing MP4 is skipped unless `--overwrite`.

Encoding uses `imageio` + `imageio-ffmpeg`'s bundled static ffmpeg (H.264,
yuv420p). `bbox_inches="tight"` makes frame dimensions vary by a few
pixels; rather than rescale (which blurs), every frame is padded on white
to the group's maximum width/height before encoding.

STATUS (2026-08-30): new.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

DEFAULT_FRAMES_ROOT = "outputs/frames"
DEFAULT_MOVIES_DIR = "outputs/movies"
DEFAULT_FPS = 12

# Frame types with one frame per hour (HRRR's cadence, or hourly lightning
# windows). DEFAULT_FPS suits the 15-min / 5-min series; at 12 fps an hourly
# movie runs too fast to follow (a 73-frame d1 movie in ~6 s, a d3 one in
# 1-2 s), so these get HOURLY_FPS -- half a second per frame. Includes
# `run_compare`'s hourly `panels_*` types.
HOURLY_FPS = 2
HOURLY_TYPES = {
    "tb4", "tb4_simple", "refl3", "hrrr_tb", "hrrr_refl",
    "d2_tb4", "d2_refl4", "d3_tb4", "d3_refl4", "d2hrrr_tb4", "d2hrrr_refl4",
    "lightning_cg", "lightning_total",
    "panels_tb", "panels_refl",
}

GROUPS = {
    "comparison": ["tb4", "tb4_simple", "refl3"],
    "single": ["wrf_tb", "wrf_refl", "wrf_olr_tb", "hrrr_tb", "hrrr_refl", "goes_tb", "mrms_refl"],
    "d2": ["d2_tb4", "d2_refl4"],
    "d3": ["d3_tb4", "d3_refl4"],
    "d2hrrr": ["d2hrrr_tb4", "d2hrrr_refl4"],
    "lightning": ["lightning_cg", "lightning_total"],
}
GROUPS["all"] = (
    GROUPS["comparison"] + GROUPS["single"] + GROUPS["d2"] + GROUPS["d3"] + GROUPS["d2hrrr"]
    + GROUPS["lightning"]
)

# Frame types rendered on a black background (`plotting.py`'s
# `dark_background=True` default for the lightning comparisons) --
# `encode_movie`'s pad-to-max-size canvas must match, or a size-mismatched
# frame (bbox_inches="tight" makes sizes vary by a few px) gets a visible
# white bar instead of blending in.
DARK_TYPES = {"lightning_cg", "lightning_total"}

_FRAME_RE = re.compile(r"^(?P<base>.+)_(?P<stamp>\d{8}_\d{4})Z\.png$")


def _resolve_types(which: str) -> list[str]:
    names: list[str] = []
    for token in which.split(","):
        token = token.strip()
        if not token:
            continue
        names.extend(GROUPS.get(token, [token]))
    seen: set[str] = set()
    return [n for n in names if not (n in seen or seen.add(n))]


def discover_movies(
    frames_root: str | Path, types: list[str], movies_dir: str | Path
) -> list[tuple[str, list[Path], Path, str]]:
    """`(name, ordered_frames, out_path, ptype)` per movie, across the given
    type subdirs of `frames_root`. Deterministically ordered."""
    frames_root = Path(frames_root)
    movies_dir = Path(movies_dir)
    movies: list[tuple[str, list[Path], Path, str]] = []

    for ptype in types:
        type_dir = frames_root / ptype
        if not type_dir.is_dir():
            continue
        groups: dict[str, list[tuple[str, Path]]] = defaultdict(list)
        for png in type_dir.glob("*.png"):
            m = _FRAME_RE.match(png.name)
            if m is None:
                continue
            groups[m["base"]].append((m["stamp"], png))
        for base in sorted(groups):
            frames = [p for _, p in sorted(groups[base])]
            movies.append((base, frames, movies_dir / f"{base}.mp4", ptype))
    return movies


def encode_movie(
    frames: list[Path], out_path: Path, fps: int, bg_color: str = "white", log=print
) -> None:
    import imageio.v2 as iio
    import numpy as np
    from PIL import Image

    sizes = [Image.open(f).size for f in frames]
    max_w = max(w for w, _ in sizes)
    max_h = max(h for _, h in sizes)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    # Scratch name must still end in a real video extension -- imageio's
    # FFMPEG plugin rejects a URI it can't identify by suffix.
    tmp_path = out_path.with_name(out_path.stem + ".partial.mp4")
    writer = iio.get_writer(
        tmp_path, fps=fps, codec="libx264", quality=8,
        macro_block_size=16, ffmpeg_params=["-pix_fmt", "yuv420p"],
    )
    try:
        for f in frames:
            im = Image.open(f).convert("RGB")
            if im.size != (max_w, max_h):
                canvas = Image.new("RGB", (max_w, max_h), bg_color)
                canvas.paste(im, ((max_w - im.width) // 2, (max_h - im.height) // 2))
                im = canvas
            writer.append_data(np.asarray(im))
    finally:
        writer.close()
    tmp_path.replace(out_path)
    log(f"  {out_path.name}: {len(frames)} frames @ {fps} fps, {max_w}x{max_h}")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m bnf_wrf_eval.make_movies",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--which", default="all",
                   help="comparison | single | all | comma-separated type names")
    p.add_argument("--frames-root", default=DEFAULT_FRAMES_ROOT)
    p.add_argument("--movies-dir", default=DEFAULT_MOVIES_DIR)
    p.add_argument("--fps", type=int, default=None,
                   help=f"frames per second for every movie (default: {DEFAULT_FPS}, "
                        f"or {HOURLY_FPS} for HOURLY_TYPES)")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--task-index", type=int, default=None)
    p.add_argument("--task-count", type=int, default=None)
    p.add_argument("--list", action="store_true",
                   help="print the resolved movie list and exit")
    return p


def _resolve_task_slice(args: argparse.Namespace) -> tuple[int | None, int | None]:
    idx, cnt = args.task_index, args.task_count
    if idx is None and "SLURM_ARRAY_TASK_ID" in os.environ:
        idx = int(os.environ["SLURM_ARRAY_TASK_ID"])
    if cnt is None and "SLURM_ARRAY_TASK_COUNT" in os.environ:
        cnt = int(os.environ["SLURM_ARRAY_TASK_COUNT"])
    if (idx is None) != (cnt is None):
        raise SystemExit("task-index and task-count must be given together")
    return idx, cnt


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    types = _resolve_types(args.which)
    movies = discover_movies(args.frames_root, types, args.movies_dir)

    idx, cnt = _resolve_task_slice(args)
    if idx is not None:
        movies = movies[idx::cnt]

    if args.list:
        for name, frames, out_path, ptype in movies:
            print(f"{name:52s}  {len(frames):5d} frames  -> {out_path}")
        print(f"# {len(movies)} movie(s)  types={','.join(types)}")
        return 0

    if not movies:
        print("no frames found -- run bnf_wrf_eval.batch_plot first", flush=True)
        return 1

    n_ok = n_skip = n_fail = 0
    for name, frames, out_path, ptype in movies:
        if out_path.exists() and not args.overwrite:
            print(f"SKIP  {name} ({out_path.name} exists)", flush=True)
            n_skip += 1
            continue
        print(f"MAKE  {name}  ({len(frames)} frames)", flush=True)
        try:
            bg_color = "black" if ptype in DARK_TYPES else "white"
            fps = args.fps or (HOURLY_FPS if ptype in HOURLY_TYPES else DEFAULT_FPS)
            encode_movie(frames, out_path, fps, bg_color=bg_color)
            n_ok += 1
        except Exception as exc:  # noqa: BLE001
            print(f"FAIL  {name}: {exc!r}", flush=True)
            n_fail += 1

    print(f"make_movies done: ok={n_ok} skipped={n_skip} failed={n_fail}", flush=True)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
