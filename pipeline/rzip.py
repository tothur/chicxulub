# fetch single members of a remote zip via HTTP range requests
import sys, struct, urllib.request, zlib
url, pattern, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
def rng(a, b):
    r = urllib.request.Request(url, headers={"Range": f"bytes={a}-{b}"}); return urllib.request.urlopen(r, timeout=120).read()
size = int(urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60).headers["Content-Length"])
tail = rng(size - 65536, size - 1); i = tail.rfind(b"PK\x05\x06")
cd_size, cd_off = struct.unpack("<II", tail[i + 12:i + 20])
cd = rng(cd_off, cd_off + cd_size - 1); p = 0; found = []
while p < len(cd) and cd[p:p + 4] == b"PK\x01\x02":
    meth, = struct.unpack("<H", cd[p + 10:p + 12]); csz, usz = struct.unpack("<II", cd[p + 20:p + 28])
    nl, el, cl = struct.unpack("<HHH", cd[p + 28:p + 34]); lho, = struct.unpack("<I", cd[p + 42:p + 46])
    name = cd[p + 46:p + 46 + nl].decode(); found.append(name)
    if pattern in name:
        h = rng(lho, lho + 29); fnl, fel = struct.unpack("<HH", h[26:30])
        data = rng(lho + 30 + fnl + fel, lho + 30 + fnl + fel + csz - 1)
        out = zlib.decompress(data, -15) if meth == 8 else data
        open(outdir + "/" + name.split("/")[-1], "wb").write(out); print("got", name, len(out))
    p += 46 + nl + el + cl
print([n for n in found if ("_6" in n or "_70Ma" in n) and "MACOSX" not in n])
