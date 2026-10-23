"""Combine a pinned upstream LDPC checkout with project-owned extensions."""
import io
import json
import re
import shutil
import subprocess
import tarfile
import tempfile
import tokenize
from pathlib import Path

from tqdm.auto import tqdm

PYTHON_SUFFIXES = {'.py', '.pyx', '.pxd', '.pyi'}
CPP_SUFFIXES = {'.h', '.hh', '.hpp', '.hxx', '.c', '.cc', '.cpp', '.cxx', '.inl'}


def rename_python(text):
    """Rename Python/Cython identifiers while preserving strings and spacing."""
    offsets = [0]
    for line in text.splitlines(keepends=True):
        offsets.append(offsets[-1] + len(line))
    replacements = []
    for token in tokenize.generate_tokens(io.StringIO(text).readline):
        if token.type == tokenize.NAME and token.string == 'ldpc':
            start = offsets[token.start[0] - 1] + token.start[1]
            end = offsets[token.end[0] - 1] + token.end[1]
            replacements.append((start, end))
    for start, end in reversed(replacements):
        text = text[:start] + 'ldpc_modified' + text[end:]

    text = re.sub(r'(["\'])ldpc::', r'\1ldpc_modified::', text)
    text = re.sub(r'(\bnamespace\s+["\'])ldpc(["\'])', r'\1ldpc_modified\2', text)
    text = re.sub(
        r'(\b(?:(?:importlib\.)?metadata\.)?version\(\s*)(["\'])ldpc\2(\s*\))',
        r'\1\2ldpc_modified\2\3', text,
    )
    return text


def rename_cpp(text):
    """Rename namespaces and qualified names without changing include paths."""
    text = re.sub(r'\bnamespace\s+ldpc\b', lambda match: match[0].replace('ldpc', 'ldpc_modified'), text)
    return re.sub(r'\bldpc::', 'ldpc_modified::', text)


def prepare():
    custom = Path(__file__).resolve().parent
    root = custom.parents[3]
    upstream = custom.parent / 'outer/ldpc'
    destination = root / 'build/ldpc_modified'
    config = json.loads((custom / 'config.json').read_text())
    assert set(config) == {'upstream_commit', 'print_level'}, config.keys()
    commit, level = config['upstream_commit'], config['print_level']
    assert isinstance(commit, str) and re.fullmatch(r'[0-9a-f]{40}', commit), commit
    assert type(level) is int and level in {1, 2, 3}, level

    actual = subprocess.check_output(['git', '-C', str(upstream), 'rev-parse', 'HEAD'], text=True).strip()
    assert actual == commit, ('Upstream commit mismatch', actual, commit)
    dirty = subprocess.check_output(['git', '-C', str(upstream), 'status', '--porcelain'], text=True)
    assert not dirty, ('Upstream submodule must remain clean', dirty)
    marker = destination / '.ldpc_custom_build.json'
    if destination.exists():
        assert marker.is_file(), ('Refusing to replace an unowned directory', destination)
        assert json.loads(marker.read_text())['owner'] == str(custom), destination

    destination.parent.mkdir(parents=True, exist_ok=True)
    archive = subprocess.check_output(['git', '-C', str(upstream), 'archive', commit])
    with tempfile.TemporaryDirectory(prefix='ldpc-prepare-', dir=destination.parent) as temporary:
        staged = Path(temporary)
        with tarfile.open(fileobj=io.BytesIO(archive)) as contents:
            contents.extractall(staged, filter='data')
        package = staged / 'src_python/ldpc_modified'
        (staged / 'src_python/ldpc').rename(package)

        files = sorted(path for base in [package, staged / 'src_cpp'] for path in base.rglob('*') if path.is_file())
        for path in tqdm(files, desc='Prepare LDPC namespaces', disable=level == 1):
            if path.suffix in PYTHON_SUFFIXES:
                path.write_text(rename_python(path.read_text()))
            elif path.suffix in CPP_SUFFIXES:
                path.write_text(rename_cpp(path.read_text()))
            if level == 3:
                print('Prepared', path.relative_to(staged))

        files = sorted(path for base in [custom / 'src_cpp', custom / 'src_python', custom / 'python_test']
                       for path in base.rglob('*') if path.is_file() and '__pycache__' not in path.parts)
        files += [custom / 'setup.py', custom / 'pyproject.toml']
        for path in tqdm(files, desc='Copy custom LDPC sources', disable=level == 1):
            assert path.suffix not in {'.so', '.pyc', '.o', '.cpp', '.pyi'}, path
            target = staged / path.relative_to(custom)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
            if level == 3:
                print('Copied', path.relative_to(custom))

        (staged / '.ldpc_custom_build.json').write_text(json.dumps({
            'owner': str(custom), 'upstream_commit': commit,
        }, indent=2) + '\n')
        if destination.exists():
            shutil.rmtree(destination)
        staged.rename(destination)

    if level >= 2:
        print('Prepared', destination)
    return destination


if __name__ == '__main__':
    prepare()
