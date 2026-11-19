#!/usr/bin/env python3
"""
pathsafe.py - output-directory validation and bounded, targeted cleanup.

Rationale (reviewer comment C-06): the config schema previously accepted any
string for output_dir, and the workflow performs destructive operations
(rm -rf, GenomicsDB workspace recreation). A malformed value such as "", "/",
"~" or the repository root could, in the worst case, delete unintended files.

This module rejects dangerous output_dir values BEFORE any job or cleanup runs,
and restricts deletion to an explicit, workflow-owned list of RELATIVE file
paths inside the validated run directory. It deletes individual files only -
never a wildcard, a user-supplied list, or a recursive tree.

Safety boundary (unchanged from the design brief): these checks assume the run
directory is controlled by the workflow user and is not swapped/altered by
another process between the check and the deletion. They are not a defence
against an active attacker racing directory swaps on a shared writable path.

Self-test:
    python -I pathsafe.py --selftest
"""

from pathlib import Path


def checked_path(path: Path) -> Path:
    """Resolve `path`, refusing any symlink in the path or its parents."""
    if any(p.is_symlink() for p in (path, *path.parents)):
        raise ValueError("Symlinks are not allowed in run/cleanup paths")
    return path.resolve()


def validate_output(raw, allowed_root, project_root) -> Path:
    """
    Validate and resolve an output_dir value.

    raw          - the configured output_dir (str).
    allowed_root - an existing, dedicated parent directory the output must live
                   strictly below (str or Path).
    project_root - the repository/project root; relative outputs resolve against
                   it (NOT the shell working directory).

    Returns the resolved, validated output Path. Raises ValueError otherwise.
    """
    if not isinstance(raw, str) or not raw.strip():
        raise ValueError("output_dir must be a non-empty string")
    path = Path(raw)
    if ".." in path.parts:
        raise ValueError("Parent traversal ('..') is not allowed in output_dir")

    project = Path(project_root).resolve(strict=True)
    root = checked_path(Path(allowed_root).absolute())
    if not root.is_dir() or root in {Path("/"), Path.home().resolve()}:
        raise ValueError("allowed_root must be an existing dedicated directory")

    path = path if path.is_absolute() else project / path
    out = checked_path(path)

    forbidden = {Path("/"), Path.home().resolve(), project}
    if out in forbidden or out in project.parents:
        raise ValueError(f"Unsafe output_dir: {out}")
    if out == root or root not in out.parents:
        raise ValueError(f"output_dir must be a subdirectory of {root}")
    if out.exists() and not out.is_dir():
        raise ValueError("output_dir exists but is not a directory")
    return out


def resolve_output_dir(raw, project_root) -> Path:
    """
    Convenience wrapper used by the Snakefile. Picks a sensible allowed_root:
      - relative output_dir -> must live strictly under project_root;
      - absolute output_dir -> its parent is the allowed_root (supports HPC
        scratch), while the shared forbidden set still blocks '/', the home
        directory, the project root and any parent of it.
    """
    project = Path(project_root).resolve(strict=True)
    if not isinstance(raw, str) or not raw.strip():
        raise ValueError("output_dir must be a non-empty string")
    candidate = Path(raw)
    if candidate.is_absolute():
        allowed_root = candidate.parent
    else:
        allowed_root = project
    return validate_output(raw, str(allowed_root), str(project))


def cleanup_files(raw, allowed_root, project_root, workflow_files, *, dry_run=True):
    """
    Delete an explicit list of workflow-owned RELATIVE files inside the
    validated output_dir. Returns the list of resolved targets that exist.
    With dry_run=True (default) nothing is deleted - the targets are only
    previewed. A separate explicit call with dry_run=False performs deletion.
    """
    out = validate_output(raw, allowed_root, project_root)
    targets = []
    for name in workflow_files:
        rel = Path(name)
        if rel.is_absolute() or ".." in rel.parts or rel == Path("."):
            raise ValueError("Cleanup entries must be relative file paths")
        target = checked_path(out / rel)
        if out not in target.parents:
            raise ValueError("Cleanup target escaped output_dir")
        if target.exists():
            if not target.is_file():
                raise ValueError("Cleanup accepts individual files only")
            targets.append(target)
    if not dry_run:
        for target in targets:
            target.unlink()
    return targets


# ---------------------------------------------------------------------------
# Self-test (temporary dummy dirs only; no workflow files touched).
# ---------------------------------------------------------------------------
def _selftest():
    import tempfile

    passed = 0

    def expect_ok(fn, label):
        nonlocal passed
        fn()
        passed += 1
        print(f"  [OK]   {label}")

    def expect_fail(fn, label):
        nonlocal passed
        try:
            fn()
        except ValueError:
            passed += 1
            print(f"  [OK]   {label} (rejected as expected)")
        else:
            raise AssertionError(f"expected ValueError: {label}")

    # Resolve the temp base to its physical path first: on many clusters
    # $TMPDIR itself sits under a symlink (e.g. /scratch -> /gpfs/...), which
    # the validator correctly rejects. Using the physical path keeps the test's
    # own scaffolding free of symlink components.
    import shutil
    tmp = str(Path(tempfile.mkdtemp()).resolve())
    try:
        proj = Path(tmp) / "project"
        (proj).mkdir()
        # valid relative output under project
        expect_ok(lambda: resolve_output_dir("results", str(proj)),
                  "relative 'results' under project")
        # valid absolute output under a dedicated scratch parent
        scratch = Path(tmp) / "scratch"
        scratch.mkdir()
        expect_ok(lambda: resolve_output_dir(str(scratch / "run1"), str(proj)),
                  "absolute scratch/run1")
        # rejections
        expect_fail(lambda: resolve_output_dir("", str(proj)), "empty string")
        expect_fail(lambda: resolve_output_dir("   ", str(proj)), "whitespace")
        expect_fail(lambda: resolve_output_dir("/", str(proj)), "root '/'")
        expect_fail(lambda: resolve_output_dir(str(Path.home()), str(proj)),
                    "home directory")
        expect_fail(lambda: resolve_output_dir(".", str(proj)), "project root '.'")
        expect_fail(lambda: resolve_output_dir("../escape", str(proj)),
                    "parent traversal")
        # cleanup: only listed relative files inside output, dry-run previews
        out = resolve_output_dir("results", str(proj))
        out.mkdir(parents=True, exist_ok=True)
        (out / "keep.txt").write_text("keep")
        (out / "drop.bam").write_text("x")
        preview = cleanup_files("results", str(proj), str(proj),
                                ["drop.bam"], dry_run=True)
        assert (out / "drop.bam").exists(), "dry-run must not delete"
        assert len(preview) == 1
        passed += 1
        print("  [OK]   cleanup dry-run previews without deleting")
        cleanup_files("results", str(proj), str(proj), ["drop.bam"], dry_run=False)
        assert not (out / "drop.bam").exists(), "real cleanup must delete target"
        assert (out / "keep.txt").exists(), "cleanup must not touch other files"
        passed += 1
        print("  [OK]   real cleanup deletes only the listed file")
        expect_fail(lambda: cleanup_files("results", str(proj), str(proj),
                                          ["/etc/passwd"], dry_run=False),
                    "absolute cleanup entry")
        expect_fail(lambda: cleanup_files("results", str(proj), str(proj),
                                          ["../outside"], dry_run=False),
                    "traversal cleanup entry")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print(f"\nAll {passed} path-safety self-tests passed.")


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        _selftest()
    else:
        print(__doc__)
