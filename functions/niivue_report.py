#!/usr/bin/env python3
"""Create a local NiiVue volume viewer with bundled, versioned assets."""
import argparse
import html
import json
import shutil
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--volume", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--assets", type=Path,
                        default=Path(__file__).resolve().parents[1] / "assets" / "niivue")
    args = parser.parse_args()
    volume = args.volume.resolve(strict=True)
    if not volume.name.endswith((".nii", ".nii.gz")):
        parser.error("--volume must be a NIfTI file")
    bundle = args.assets / "index.js"
    if not bundle.is_file():
        parser.error(f"NiiVue bundle is missing: {bundle}")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(volume, output / volume.name)
    shutil.copy2(bundle, output / "niivue.js")
    title = html.escape(volume.name)
    config = json.dumps([{"url": volume.name}]).replace("<", "\\u003c")
    page = """<!doctype html>
<html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>NiiVue · TITLE</title>
<style>html,body{margin:0;height:100%;background:#101419;color:#e5e9ef;font:15px system-ui}
header{padding:12px 20px}canvas{display:block;width:100%;height:calc(100% - 76px)}
small{color:#a6b1bf}#error{color:#ffb4a9}</style>
<header>TITLE<br><small>Drag to rotate · Scroll to change slices</small><span id="error"></span></header>
<canvas id="viewer"></canvas>
<script type="module">
import {Niivue} from './niivue.js';
try {
  const viewer = new Niivue({backColor:[0.063,0.078,0.098,1],show3Dcrosshair:true});
  await viewer.attachTo('viewer');
  await viewer.loadVolumes(VOLUMES);
  document.body.dataset.loaded = 'true';
} catch (error) {
  document.getElementById('error').textContent = String(error);
  document.body.dataset.loaded = 'false';
}
</script></html>
"""
    (output / "index.html").write_text(
        page.replace("TITLE", title).replace("VOLUMES", config), encoding="utf-8")
    (output / "README.txt").write_text(
        "Serve this directory with python -m http.server --bind 127.0.0.1 8000\n"
        "then open http://127.0.0.1:8000/ in a WebGL2 browser.\n"
        "The volume and viewer are local; no internet connection is required.\n")
    print(output / "index.html")

if __name__ == "__main__":
    main()
