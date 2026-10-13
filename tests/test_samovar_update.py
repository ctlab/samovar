"""samovar update fetches a branch in place and then runs install.sh full ask."""

import os
import stat
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
UPDATE = REPO / "bin" / "samovar_update"


def _git(repo, *args):
    subprocess.check_call(
        [
            "git",
            "-c",
            "user.email=samovar-test@example.com",
            "-c",
            "user.name=samovar-test",
            "-C",
            str(repo),
            *args,
        ]
    )


def _write(path: Path, text: str, executable=False):
    path.write_text(text)
    if executable:
        path.chmod(path.stat().st_mode | stat.S_IEXEC)


def _origin(tmp: Path) -> Path:
    origin = tmp / "origin"
    origin.mkdir()
    subprocess.check_call(["git", "init", "-b", "main", str(origin)])
    _write(origin / "environment.yml", "name: samovar\ndependencies: []\n")
    _write(origin / "install.sh", "#!/bin/bash\necho old \"$@\"\n", executable=True)
    _git(origin, "add", "environment.yml", "install.sh")
    _git(origin, "commit", "-m", "base")
    return origin


def _conda(tmp: Path, *, has_env: bool) -> Path:
    base = tmp / "conda-base"
    (base / "etc" / "profile.d").mkdir(parents=True)
    log = tmp / "conda.log"
    marker = tmp / "has-env"
    if has_env:
        marker.write_text("1\n")
    (base / "etc" / "profile.d" / "conda.sh").write_text(
        f"""
conda() {{
  printf '%s\\n' "$*" >> {log}
  if [ "$1" = "env" ] && [ "$2" = "list" ]; then
    if [ -f {marker} ]; then echo "samovar"; fi
    return 0
  fi
  if [ "$1" = "env" ] && [ "$2" = "create" ]; then
    touch {marker}
    return 0
  fi
  if [ "$1" = "activate" ]; then
    export CONDA_DEFAULT_ENV="$2"
    return 0
  fi
}}
"""
    )
    bindir = tmp / "bin"
    bindir.mkdir()
    stub = bindir / "conda"
    _write(
        stub,
        f"""#!/bin/bash
printf '%s\\n' "$*" >> {log}
if [ "$1" = "info" ] && [ "$2" = "--base" ]; then
  echo {base}
fi
""",
        executable=True,
    )
    return bindir


def _run(tmp, clone, bindir, *args, env=None):
    merged = os.environ.copy()
    merged["PATH"] = f"{bindir}:{merged['PATH']}"
    merged["SAMOVAR_ROOT"] = str(clone)
    merged.pop("CONDA_DEFAULT_ENV", None)
    if env:
        merged.update(env)
    return subprocess.run(
        ["bash", str(UPDATE), *args],
        cwd=tmp,
        env=merged,
        text=True,
        capture_output=True,
    )


def test_update_fast_forwards_main_and_creates_env(tmp_path):
    origin = _origin(tmp_path)
    clone = tmp_path / "clone"
    subprocess.check_call(["git", "clone", str(origin), str(clone)])
    _write(
        origin / "install.sh",
        "#!/bin/bash\nprintf '%s\\n' \"$*\" > \"$SAMOVAR_ROOT/ran.txt\"\n",
        executable=True,
    )
    _git(origin, "add", "install.sh")
    _git(origin, "commit", "-m", "new install")
    bindir = _conda(tmp_path, has_env=False)

    result = _run(tmp_path, clone, bindir)
    assert result.returncode == 0, result.stderr
    assert (clone / "ran.txt").read_text().strip() == "full ask"
    head = subprocess.check_output(
        ["git", "-C", str(clone), "rev-parse", "HEAD"], text=True
    ).strip()
    origin_head = subprocess.check_output(
        ["git", "-C", str(origin), "rev-parse", "HEAD"], text=True
    ).strip()
    assert head == origin_head
    log = (tmp_path / "conda.log").read_text()
    assert "env create -f" in log
    assert "activate samovar" in log


def test_update_other_branch_skips_existing_env(tmp_path):
    origin = _origin(tmp_path)
    _git(origin, "checkout", "-b", "feature")
    _write(
        origin / "install.sh",
        "#!/bin/bash\nprintf '%s\\n' \"$*\" > \"$SAMOVAR_ROOT/ran.txt\"\n",
        executable=True,
    )
    _git(origin, "add", "install.sh")
    _git(origin, "commit", "-m", "feature install")
    clone = tmp_path / "clone"
    subprocess.check_call(["git", "clone", "--branch", "main", str(origin), str(clone)])
    bindir = _conda(tmp_path, has_env=True)

    result = _run(tmp_path, clone, bindir, "--branch", "feature")
    assert result.returncode == 0, result.stderr
    branch = subprocess.check_output(
        ["git", "-C", str(clone), "rev-parse", "--abbrev-ref", "HEAD"], text=True
    ).strip()
    assert branch == "feature"
    assert (clone / "ran.txt").read_text().strip() == "full ask"
    log = (tmp_path / "conda.log").read_text()
    assert "env create" not in log
    assert "activate samovar" in log


def test_update_skips_conda_when_env_is_active(tmp_path):
    origin = _origin(tmp_path)
    clone = tmp_path / "clone"
    subprocess.check_call(["git", "clone", str(origin), str(clone)])
    bindir = _conda(tmp_path, has_env=True)
    result = _run(
        tmp_path,
        clone,
        bindir,
        env={"CONDA_DEFAULT_ENV": "samovar"},
    )
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "conda.log").exists()


def test_update_uses_samovar_root_not_the_configured_checkout(tmp_path):
    origin = _origin(tmp_path)
    (origin / "bin").mkdir()
    _write(origin / "bin/samovar_kmer", "built-old\n", executable=True)
    _git(origin, "add", "bin/samovar_kmer")
    _git(origin, "commit", "-m", "track kmer binary")
    _write(
        origin / "install.sh",
        "#!/bin/bash\nprintf '%s\\n' \"$*\" > \"$SAMOVAR_ROOT/ran.txt\"\n",
        executable=True,
    )
    _git(origin, "add", "install.sh")
    _git(origin, "commit", "-m", "new install")

    decoy = tmp_path / "decoy"
    subprocess.check_call(["git", "clone", str(origin), str(decoy)])
    _write(decoy / "bin/samovar_kmer", "dirty-decoy\n", executable=True)
    decoy_head = subprocess.check_output(
        ["git", "-C", str(decoy), "rev-parse", "HEAD"], text=True
    ).strip()
    cfg = tmp_path / "config.json"
    cfg.write_text('{"root": "%s"}\n' % decoy)

    clone = tmp_path / "clone"
    subprocess.check_call(["git", "clone", str(origin), str(clone)])
    subprocess.check_call(["git", "-C", str(clone), "reset", "--hard", "HEAD~1"])
    _write(clone / "bin/samovar_kmer", "dirty-clone\n", executable=True)
    bindir = _conda(tmp_path, has_env=False)

    result = _run(
        tmp_path,
        clone,
        bindir,
        env={"SAMOVAR_CONFIG": str(cfg)},
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "decoy" not in result.stdout
    assert (clone / "ran.txt").read_text().strip() == "full ask"
    assert "dirty-decoy" in (decoy / "bin/samovar_kmer").read_text()
    decoy_after = subprocess.check_output(
        ["git", "-C", str(decoy), "rev-parse", "HEAD"], text=True
    ).strip()
    assert decoy_after == decoy_head
    clone_head = subprocess.check_output(
        ["git", "-C", str(clone), "rev-parse", "HEAD"], text=True
    ).strip()
    origin_head = subprocess.check_output(
        ["git", "-C", str(origin), "rev-parse", "HEAD"], text=True
    ).strip()
    assert clone_head == origin_head
    result = subprocess.run(
        ["bash", str(UPDATE), "--help"],
        text=True,
        capture_output=True,
        check=True,
    )
    assert "--branch" in result.stdout
    assert "main" in result.stdout


def test_update_rejects_unknown_option(tmp_path):
    result = _run(tmp_path, tmp_path, tmp_path, "--nope")
    assert result.returncode != 0
    assert "unknown option" in result.stderr
