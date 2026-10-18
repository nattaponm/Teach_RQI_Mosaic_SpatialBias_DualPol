"""Google Colab helpers: set up the environment, download the teaching data, save results to your own fork.

Typical first cell of every notebook:

    import os, sys, subprocess
    if "google.colab" in sys.modules and not os.path.exists("/content/Teach_RQI_Mosaic_SpatialBias_DualPol"):
        subprocess.run(["git", "clone", "-q", "--depth", "1",
                        "https://github.com/nattaponm/Teach_RQI_Mosaic_SpatialBias_DualPol.git",
                        "/content/Teach_RQI_Mosaic_SpatialBias_DualPol"])
        os.chdir("/content/Teach_RQI_Mosaic_SpatialBias_DualPol")
    sys.path.insert(0, "src")
    from radteach import colab
    P = colab.setup()

Saving results (students):  colab.save_to_fork("6512345678")  — needs two Colab secrets:
    GITHUB_USER   your GitHub user name (the owner of the fork)
    GITHUB_TOKEN  a fine-grained personal access token with access ONLY to your fork,
                  permission "Contents: Read and write"
Results go to the branch results/<student_id> of YOUR fork. The instructor's repository cannot be written:
students are not collaborators, and the function refuses to push to the upstream owner's repository.
"""
from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

REPO = "Teach_RQI_Mosaic_SpatialBias_DualPol"
UPSTREAM_OWNER = "nattaponm"
PIP_PACKAGES = {"pyart": "arm_pyart", "yaml": "pyyaml", "rasterio": "rasterio", "geopandas": "geopandas"}


def in_colab() -> bool:
    return "google.colab" in sys.modules


def _pip(pkgs):
    if pkgs:
        print("installing:", " ".join(pkgs), "(about 1–2 minutes the first time) …")
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", *pkgs], check=False)


def ensure_packages(extra=()):
    import importlib.util
    missing = [pip for mod, pip in PIP_PACKAGES.items() if importlib.util.find_spec(mod) is None]
    _pip(missing + list(extra))


def setup(lang="th", install=True, download=True, packages=("teaching_data",), quiet=False):
    """Prepare paths, packages, fonts and data. Returns a radteach.config.Paths object."""
    from . import config
    root = config.find_root()
    os.chdir(root)
    if str(root / "src") not in sys.path:
        sys.path.insert(0, str(root / "src"))
    if install:
        ensure_packages()
    from . import plotting
    font = plotting.setup(lang)
    P = config.paths(root)
    if download:
        for p in packages:
            download_package(p, P)
    if not quiet:
        print(f"repository : {root}\nColab      : {in_colab()}\nThai font  : {font}\ndata       : {P.raw}")
    return P


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def download_package(name="teaching_data", P=None, force=False):
    """Download and unzip a teaching data package listed in config/data_sources.yaml into data/raw/."""
    from . import config
    P = P or config.paths()
    src = config.load("data_sources")
    pk = src["packages"][name]
    dest_zip = P.raw / pk["file"]
    marker = P.raw / (Path(pk["file"]).stem + ".ok")
    if marker.exists() and not force:
        return P.raw
    url = os.environ.get("RADTEACH_DATA_URL", src["release_base"]) + pk["file"]
    if not dest_zip.exists() or force:
        print(f"downloading {pk['file']} …")
        tmp = str(dest_zip) + ".part"
        urllib.request.urlretrieve(url, tmp)
        os.replace(tmp, dest_zip)
    if pk.get("sha256") and pk["sha256"] != "TO_BE_FILLED":
        h = sha256(dest_zip)
        if h != pk["sha256"]:
            raise RuntimeError(f"checksum mismatch for {pk['file']}: {h}")
    with zipfile.ZipFile(dest_zip) as z:
        z.extractall(P.raw)
    marker.write_text("ok")
    dest_zip.unlink(missing_ok=True)
    print(f"{pk['file']} ready in {P.raw}")
    return P.raw


def _secret(name):
    if name in os.environ:
        return os.environ[name]
    if in_colab():
        try:
            from google.colab import userdata
            return userdata.get(name)
        except Exception:
            return None
    return None


def colab_link(notebook, user=UPSTREAM_OWNER, branch="main"):
    return f"https://colab.research.google.com/github/{user}/{REPO}/blob/{branch}/notebooks/{notebook}"


def save_to_fork(student_id: str, *args, **kwargs):
    """Copy outputs/ into YOUR fork under outputs/<student_id>/ and push to branch results/<student_id>."""
    token = _secret("GITHUB_TOKEN") or "\x00"
    try:
        return _save_to_fork(student_id, *args, **kwargs)
    except subprocess.CalledProcessError as e:   # never show the token in an error message
        msg = ((e.stderr or "") + (e.stdout or "")).replace(token, "***")
        raise RuntimeError("git failed: " + msg.strip()[-500:]) from None


def _save_to_fork(student_id: str, message: str | None = None, folders=("outputs",), max_mb=50.0,
                 branch: str | None = None, workdir="/tmp/radteach_fork"):
    user, token = _secret("GITHUB_USER"), _secret("GITHUB_TOKEN")
    if not user or not token:
        raise RuntimeError("Set the Colab secrets GITHUB_USER and GITHUB_TOKEN first (see docs/student_guide_fork.md)")
    if user.lower() == UPSTREAM_OWNER.lower():
        raise RuntimeError("GITHUB_USER is the instructor account — students push only to their own fork.")
    from . import config
    root = config.find_root()
    branch = branch or f"results/{student_id}"
    url = os.environ.get("RADTEACH_FORK_URL") or f"https://{user}:{token}@github.com/{user}/{REPO}.git"   # override = tests
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    run = lambda *a, **k: subprocess.run(list(a), check=True, cwd=k.get("cwd", workdir), capture_output=True, text=True)
    subprocess.run(["git", "clone", "-q", "--depth", "1", url, workdir], check=True, capture_output=True, text=True)
    exists = subprocess.run(["git", "ls-remote", "--heads", "origin", branch], cwd=workdir,
                            capture_output=True, text=True).stdout.strip()
    if exists:
        run("git", "fetch", "-q", "--depth", "1", "origin", f"{branch}:{branch}")
        run("git", "checkout", "-q", branch)
    else:
        run("git", "checkout", "-q", "-b", branch)
    n, skipped = 0, []
    for fo in folders:
        for f in (root / fo).rglob("*"):
            if not f.is_file() or f.name == ".gitkeep":
                continue
            if f.stat().st_size > max_mb * 1e6:
                skipped.append(f.name); continue
            dst = Path(workdir) / "outputs" / student_id / f.relative_to(root / fo)
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(f, dst); n += 1
    run("git", "config", "user.name", user)
    run("git", "config", "user.email", f"{user}@users.noreply.github.com")
    run("git", "add", "-f", f"outputs/{student_id}")
    st = subprocess.run(["git", "status", "--porcelain"], cwd=workdir, capture_output=True, text=True).stdout
    if not st.strip():
        print("nothing new to save"); return
    run("git", "commit", "-q", "-m", message or f"results of {student_id}")
    run("git", "push", "-q", "origin", branch)
    print(f"saved {n} files to https://github.com/{user}/{REPO}/tree/{branch}/outputs/{student_id}")
    if skipped:
        print("skipped (too large):", skipped)
