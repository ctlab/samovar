from pathlib import Path

import pytest

from samovar import sidecar
from samovar.sidecar import (
    _MEGAHIT_FULL,
    _MEGAHIT_PLAIN,
    _MEGAHIT_POPCNT,
    megahit_core_for_cpu,
    megahit_instruction_violations,
)




@pytest.mark.parametrize("name", ["nanosim", "kraken2", "kaiju"])
def test_sidecar_health_checks_binary_modules_and_versions(tmp_path, monkeypatch, name):
    prefix = Path(tmp_path)
    binary = prefix / "bin" / sidecar.SIDECARS[name]["binary"]
    binary.parent.mkdir()
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o755)

    for helper in sidecar.SIDECARS[name].get("companions", ()):
        path = binary.parent / helper
        path.write_text("#!/bin/sh\n")
        path.chmod(0o755)
    monkeypatch.setattr(sidecar, "env_has_runtime", lambda *_: True)
    monkeypatch.setattr(sidecar, "env_has_versions", lambda *_: True)
    assert sidecar.sidecar_is_healthy(name, prefix)

    monkeypatch.setattr(sidecar, "env_has_versions", lambda *_: False)
    assert not sidecar.sidecar_is_healthy(name, prefix)


def test_motus_install_preserves_independently_configured_companions(tmp_path, monkeypatch):
    from samovar.main_config import lookup_tool_record
    existing = {'exec': {'env': '', 'path': str(tmp_path/'custom-samtools')}, 'type': 'aligner'}
    cfg = {'tools': {'samtools': existing.copy()}}
    monkeypatch.setattr(sidecar, 'load_config', lambda: cfg)
    monkeypatch.setattr(sidecar, 'write_config', lambda value: None)
    prefix = tmp_path/'motus-env'
    sidecar.record_sidecar('motus', prefix, str(prefix/'bin/motus'))
    assert cfg['tools']['samtools'] == existing
    assert lookup_tool_record(cfg, 'motus')['type'] == 'profiler'
    assert lookup_tool_record(cfg, 'bwa')['exec']['path'] == str(prefix/'bin/bwa')
    assert lookup_tool_record(cfg, 'bwa')['type'] == 'aligner'


@pytest.mark.parametrize('valid', [True, False])
def test_source_tool_install_verifies_archive_before_registration(tmp_path, monkeypatch, valid):
    import hashlib
    import io
    import tarfile
    import urllib.request
    from samovar.main_config import lookup_tool_record
    monkeypatch.setenv('SAMOVAR_OFFLINE', '0')  # Transport is a local BytesIO below.

    archive = io.BytesIO()
    with tarfile.open(fileobj=archive, mode='w:gz') as handle:
        member = tarfile.TarInfo('native/run_phanta.py')
        data = b'# upstream entry point\n'
        member.size = len(data)
        handle.addfile(member, io.BytesIO(data))
    content = archive.getvalue()
    spec = {**sidecar.SIDECARS['phanta'], 'source': {
        'revision': 'test', 'directory': 'native', 'url': 'https://invalid.local/source',
        'executable': 'run_phanta.py',
        'sha256': hashlib.sha256(content).hexdigest() if valid else 'invalid',
    }}
    monkeypatch.setitem(sidecar.SIDECARS, 'phanta', spec)
    monkeypatch.setattr('samovar.paths._existing_tool', lambda _: '/configured/tool')
    monkeypatch.setattr(sidecar.subprocess, 'run', lambda *_a, **_kw: __import__('subprocess').CompletedProcess(
        [], 0, stdout='--report-minimizer-data', stderr=''))
    monkeypatch.setattr(sidecar, 'sidecar_envs_dir', lambda: tmp_path / 'envs')
    monkeypatch.setattr(urllib.request, 'urlopen', lambda *_args, **_kw: io.BytesIO(content))
    cfg = {'tools': {}}
    monkeypatch.setattr(sidecar, 'load_config', lambda: cfg)
    written = []
    monkeypatch.setattr(sidecar, 'write_config', written.append)
    if valid:
        dest = sidecar.install_named('phanta')
        record = lookup_tool_record(cfg, 'phanta')
        assert record['type'] == 'profiler'
        assert Path(record['exec']['path']) == dest / 'run_phanta.py'
        assert written and (dest / 'run_phanta.py').read_bytes() == b'# upstream entry point\n'
    else:
        with pytest.raises(RuntimeError, match='checksum'):
            sidecar.install_named('phanta')
        assert not written and not (tmp_path / 'sources/phanta/test').exists()
    assert not list((tmp_path / 'sources/phanta').glob('.install-*'))


def test_griffin_selects_the_popcnt_core_without_bmi2_or_avx2():
    griffin = {'popcnt', 'avx', 'bmi1'}
    assert megahit_core_for_cpu(griffin) == _MEGAHIT_POPCNT
    assert megahit_core_for_cpu({'popcnt', 'bmi2', 'avx2'}) == _MEGAHIT_FULL
    assert megahit_core_for_cpu({'sse2'}) == _MEGAHIT_PLAIN
    assert megahit_core_for_cpu({'bmi2', 'popcnt'}) == _MEGAHIT_POPCNT
    assert megahit_core_for_cpu({'bmi2'}) == _MEGAHIT_PLAIN


def test_restricted_cores_reject_bmi2_avx2_and_plain_rejects_popcnt():
    text = 'bzhi %rax, %rbx, %rcx\nvpermq %ymm0, %ymm1, %ymm2\npopcnt %rax, %rbx\n'
    assert megahit_instruction_violations(text, _MEGAHIT_POPCNT) == ['bmi2', 'avx2']
    assert megahit_instruction_violations(text, _MEGAHIT_PLAIN) == ['bmi2', 'avx2', 'popcnt']
    assert megahit_instruction_violations(text, _MEGAHIT_FULL) == []
    assert megahit_instruction_violations('mov %rax, %rbx\n', _MEGAHIT_PLAIN) == []
