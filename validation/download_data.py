"""Download and verify public inputs used by the external validation."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from typing import Any


def sha256_file(path: Path, *, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _matches(path: Path, expected_size: int, expected_sha256: str) -> bool:
    return (
        path.exists()
        and path.stat().st_size == expected_size
        and sha256_file(path) == expected_sha256.lower()
    )


def download_file(
    url: str,
    destination: Path,
    *,
    expected_size: int,
    expected_sha256: str,
    retries: int = 5,
) -> None:
    """Download with byte-range resume and final size/SHA-256 verification."""

    if _matches(destination, expected_size, expected_sha256):
        print(f"verified {destination}")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    part = destination.with_name(destination.name + ".part")
    for attempt in range(1, retries + 1):
        offset = part.stat().st_size if part.exists() else 0
        headers = {"User-Agent": "celltogwas-validation/1.0"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                append = offset > 0 and getattr(response, "status", 200) == 206
                if offset and not append:
                    offset = 0
                mode = "ab" if append else "wb"
                print(
                    f"downloading {url} -> {destination}"
                    + (f" (resume at {offset:,} bytes)" if append else "")
                )
                with part.open(mode) as handle:
                    shutil.copyfileobj(response, handle, length=8 * 1024 * 1024)
            if part.stat().st_size != expected_size:
                raise ValueError(
                    f"size mismatch for {destination}: expected {expected_size:,}, "
                    f"got {part.stat().st_size:,}"
                )
            digest = sha256_file(part)
            if digest != expected_sha256.lower():
                raise ValueError(
                    f"SHA-256 mismatch for {destination}: expected {expected_sha256}, got {digest}"
                )
            os.replace(part, destination)
            print(f"verified {destination}")
            return
        except (OSError, urllib.error.URLError, ValueError) as error:
            if attempt == retries:
                raise RuntimeError(f"download failed after {retries} attempts: {url}") from error
            print(f"attempt {attempt} failed: {error}; retrying")
            time.sleep(min(2**attempt, 20))


def safe_extract_zip(archive: Path, destination: Path) -> None:
    """Extract a ZIP after rejecting entries outside the destination."""

    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with zipfile.ZipFile(archive) as handle:
        for member in handle.infolist():
            target = (destination / member.filename).resolve()
            if root not in target.parents and target != root:
                raise ValueError(f"unsafe ZIP member: {member.filename}")
        handle.extractall(destination)
    print(f"extracted {archive.name} -> {destination}")


def verify_manifest(
    manifest: dict[str, Any],
    root: Path,
    *,
    verify_only: bool,
    extract: bool,
) -> None:
    entries = manifest.get("datasets")
    if not isinstance(entries, list):
        raise ValueError("Manifest must contain a datasets list.")
    for entry in entries:
        relative_path = str(entry["relative_path"])
        destination = root / relative_path
        expected_size = int(entry["size"])
        expected_sha256 = str(entry["sha256"])
        if verify_only:
            if not _matches(destination, expected_size, expected_sha256):
                raise FileNotFoundError(f"Missing or invalid data file: {destination}")
            print(f"verified {destination}")
            continue
        download_file(
            str(entry["url"]),
            destination,
            expected_size=expected_size,
            expected_sha256=expected_sha256,
        )
        extract_to = entry.get("extract_to")
        if extract and extract_to:
            safe_extract_zip(destination, root / str(extract_to))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(__file__).with_name("data_manifest.json"),
    )
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--no-extract", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    verify_manifest(
        manifest,
        args.data_root.resolve(),
        verify_only=args.verify_only,
        extract=not args.no_extract,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
