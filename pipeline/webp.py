"""Re-encode the viewer's data PNGs as lossless WebP (about a fifth smaller), checking every file decodes identically.
The PNGs stay the pipeline's output; the viewer loads the .webp copies.

usage: python webp.py <sim data dir>
"""
import glob, io, os, sys
import numpy as np
from PIL import Image

OUT = sys.argv[1]
tp = tw = 0
for f in sorted(glob.glob(os.path.join(OUT, "*.png"))):
    a = np.array(Image.open(f).convert("RGB"))
    buf = io.BytesIO(); Image.fromarray(a).save(buf, "WEBP", lossless=True, quality=100, method=6, exact=True)
    assert np.array_equal(a, np.array(Image.open(io.BytesIO(buf.getvalue())).convert("RGB"))), f
    open(f[:-4] + ".webp", "wb").write(buf.getvalue())
    tp += os.path.getsize(f); tw += len(buf.getvalue())
print(f"png {tp / 1e6:.1f} MB -> webp {tw / 1e6:.1f} MB, all identical")
