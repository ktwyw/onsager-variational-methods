"""Gallery of README.md: one figure from each chapter, taken from the executed notebooks.

The figures are not redrawn here: they are copied out of the notebooks, so the gallery
always shows what the chapters show.  Output: docs/images/gallery_chNN.png.

Run from the repository root:  python scripts/make_readme_gallery.py
"""
import base64
import glob
import json
import os

root = os.path.join(os.path.dirname(__file__), "..")
# chapter number -> which figure of the notebook to use (0 = first)
choice = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 1, 7: 1}

for chapter, index in choice.items():
    path = glob.glob(os.path.join(root, "notebooks", "*", f"{chapter:02d}_*.ipynb"))[0]
    notebook = json.load(open(path, encoding="utf-8"))
    images = [output["data"]["image/png"] for cell in notebook["cells"] for output in cell.get("outputs", [])
              if "image/png" in output.get("data", {})]
    out = os.path.join(root, "docs", "images", f"gallery_ch{chapter:02d}.png")
    with open(out, "wb") as handle:
        handle.write(base64.b64decode(images[index]))
    print(f"chapter {chapter}: figure {index + 1} of {len(images)} -> {os.path.normpath(out)}")
