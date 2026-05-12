#!/usr/bin/env python3
"""Append BibTeX entries resolved from DOIs into a .bib library."""

from __future__ import annotations

import argparse
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


DOI_PATTERN = re.compile(r"10\.\S+", re.IGNORECASE)
ENTRY_KEY_PATTERN = re.compile(r"@\w+\s*\{\s*([^,\s]+)", re.IGNORECASE)
DOI_FIELD_PATTERN = re.compile(r"\bdoi\s*=\s*\{?([^},\n]+)", re.IGNORECASE)
FIELD_PATTERN = re.compile(
    r"^\s*([A-Za-z][A-Za-z0-9_-]*)\s*=\s*(.+?)\s*,?\s*$",
    re.MULTILINE | re.DOTALL,
)
WORD_PATTERN = re.compile(r"[A-Za-z0-9]+")
STOP_WORDS = {
    "a",
    "an",
    "and",
    "at",
    "by",
    "for",
    "from",
    "in",
    "of",
    "on",
    "or",
    "the",
    "to",
    "via",
    "with",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Resolve a list of DOIs into BibTeX and append them to a .bib library."
    )
    parser.add_argument(
        "dois",
        nargs="*",
        help="DOIs passed directly on the command line.",
    )
    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        help="Text file with one DOI per line. Blank lines and lines starting with # are ignored.",
    )
    parser.add_argument(
        "--bib",
        "-b",
        type=Path,
        required=True,
        help="Path to the BibTeX library file to create or update.",
    )
    parser.add_argument(
        "--email",
        help="Optional contact email for the User-Agent header.",
    )
    parser.add_argument(
        "--pause",
        type=float,
        default=0.2,
        help="Pause in seconds between DOI requests. Default: 0.2",
    )
    return parser.parse_args()


def normalize_doi(raw_value: str) -> str | None:
    value = raw_value.strip()
    if not value or value.startswith("#"):
        return None

    value = value.split("#", 1)[0].strip()
    value = value.replace("https://doi.org/", "").replace("http://doi.org/", "")
    value = value.replace("doi:", "").strip()

    match = DOI_PATTERN.search(value)
    if not match:
        return None
    return match.group(0).rstrip(".,;")


def load_requested_dois(args: argparse.Namespace) -> list[str]:
    ordered: list[str] = []
    seen: set[str] = set()

    sources = list(args.dois)
    if args.input:
        sources.extend(args.input.read_text(encoding="utf-8").splitlines())

    for raw_value in sources:
        doi = normalize_doi(raw_value)
        if doi is None:
            continue
        doi_lower = doi.lower()
        if doi_lower in seen:
            continue
        seen.add(doi_lower)
        ordered.append(doi)

    return ordered


def load_library_state(bib_path: Path) -> tuple[set[str], set[str]]:
    if not bib_path.exists():
        return set(), set()

    content = bib_path.read_text(encoding="utf-8")
    existing_dois = {
        match.group(1).strip().lower()
        for match in DOI_FIELD_PATTERN.finditer(content)
    }
    existing_keys = {
        match.group(1).strip()
        for match in ENTRY_KEY_PATTERN.finditer(content)
    }
    return existing_dois, existing_keys


def build_headers(email: str | None) -> dict[str, str]:
    user_agent = "doi-to-bibtex/1.0"
    if email:
        user_agent = f"{user_agent} ({email})"
    return {
        "Accept": "application/x-bibtex; charset=utf-8",
        "User-Agent": user_agent,
    }


def fetch_bibtex(doi: str, headers: dict[str, str]) -> str:
    encoded_doi = urllib.parse.quote(doi, safe="/")
    request = urllib.request.Request(
        f"https://doi.org/{encoded_doi}",
        headers=headers,
        method="GET",
    )

    try:
        with urllib.request.urlopen(request, timeout=20) as response:
            body = response.read().decode("utf-8", errors="replace").strip()
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"HTTP {exc.code}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"network error: {exc.reason}") from exc

    if not body.startswith("@"):
        raise RuntimeError("response is not BibTeX")
    return body


def replace_entry_key(entry: str, new_key: str) -> str:
    return ENTRY_KEY_PATTERN.sub(lambda match: match.group(0).replace(match.group(1), new_key, 1), entry, count=1)


def unique_key(base_key: str, used_keys: set[str]) -> str:
    if base_key not in used_keys:
        return base_key

    counter = 2
    while f"{base_key}_{counter}" in used_keys:
        counter += 1
    return f"{base_key}_{counter}"


def ensure_trailing_newline(text: str) -> str:
    return text if text.endswith("\n") else text + "\n"


def split_bibtex_fields(body: str) -> list[str]:
    fields: list[str] = []
    current: list[str] = []
    brace_depth = 0
    quote_open = False

    for char in body:
        if char == '"' and brace_depth == 0:
            quote_open = not quote_open
        elif char == "{":
            brace_depth += 1
        elif char == "}":
            brace_depth = max(0, brace_depth - 1)
        elif char == "," and brace_depth == 0 and not quote_open:
            candidate = "".join(current).strip()
            if candidate:
                fields.append(candidate)
            current = []
            continue
        current.append(char)

    trailing = "".join(current).strip()
    if trailing:
        fields.append(trailing)
    return fields


def extract_bibtex_body(entry: str) -> str:
    match = ENTRY_KEY_PATTERN.search(entry)
    if not match:
        return ""

    body = entry[match.end() :].strip()
    if body.startswith(","):
        body = body[1:].strip()
    if body.endswith("}"):
        body = body[:-1].strip()
    return body


def strip_outer_delimiters(value: str) -> str:
    cleaned = value.strip().rstrip(",").strip()
    if len(cleaned) >= 2 and cleaned[0] == "{" and cleaned[-1] == "}":
        return cleaned[1:-1].strip()
    if len(cleaned) >= 2 and cleaned[0] == '"' and cleaned[-1] == '"':
        return cleaned[1:-1].strip()
    return cleaned


def parse_bibtex_fields(entry: str) -> dict[str, str]:
    body = extract_bibtex_body(entry)
    if not body:
        return {}

    fields: dict[str, str] = {}
    for chunk in split_bibtex_fields(body):
        match = FIELD_PATTERN.match(chunk)
        if not match:
            continue
        fields[match.group(1).lower()] = strip_outer_delimiters(match.group(2))
    return fields


def collapse_tex_to_words(value: str) -> list[str]:
    text = re.sub(r"\\[A-Za-z]+\*?(?:\s*\{[^{}]*\})?", " ", value)
    text = re.sub(r"[{}$\\]", " ", text)
    return WORD_PATTERN.findall(text)


def title_case_token(token: str) -> str:
    if token.isdigit():
        return token
    return token[0].upper() + token[1:].lower()


def first_author_token(author_field: str) -> str:
    first_author = author_field.split(" and ", 1)[0].strip()
    if "," in first_author:
        family_name = first_author.split(",", 1)[0].strip()
    else:
        parts = first_author.split()
        family_name = parts[-1] if parts else "Unknown"

    words = collapse_tex_to_words(family_name)
    if not words:
        return "Unknown"
    return "".join(title_case_token(word) for word in words)


def year_token(fields: dict[str, str]) -> str:
    for field_name in ("year", "date"):
        value = fields.get(field_name, "")
        match = re.search(r"\b(19|20)\d{2}\b", value)
        if match:
            return match.group(0)
    return "NoYear"


def is_arxiv_preprint(fields: dict[str, str]) -> bool:
    journal = fields.get("journal", "").lower()
    publisher = fields.get("publisher", "").lower()
    archiveprefix = fields.get("archiveprefix", "").lower()
    return (
        archiveprefix == "arxiv"
        or "arxiv" in journal
        or "arxiv" in publisher
        or "eprint" in fields
    )


def journal_token(fields: dict[str, str]) -> str:
    if is_arxiv_preprint(fields):
        return "Arxiv"

    journal = fields.get("journal") or fields.get("journaltitle") or fields.get("publisher") or ""
    words = [word for word in collapse_tex_to_words(journal) if word.lower() not in STOP_WORDS]
    if not words:
        return "UnknownJournal"
    return "".join(title_case_token(word) for word in words[:2])


def short_title_token(title: str, max_words: int = 3) -> str:
    words = [word for word in collapse_tex_to_words(title) if word.lower() not in STOP_WORDS]
    if not words:
        words = collapse_tex_to_words(title)
    if not words:
        return "Untitled"
    return "".join(title_case_token(word) for word in words[:max_words])


def build_entry_key(entry: str) -> str:
    fields = parse_bibtex_fields(entry)
    author = first_author_token(fields.get("author", ""))
    year = year_token(fields)
    journal = journal_token(fields)
    short_title = short_title_token(fields.get("title", ""))
    return f"{author}{year}{journal}{short_title}"


def append_entries(bib_path: Path, entries: list[str]) -> None:
    bib_path.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if bib_path.exists() else "w"
    with bib_path.open(mode, encoding="utf-8", newline="\n") as handle:
        if mode == "a":
            handle.write("\n")
        handle.write("\n\n".join(ensure_trailing_newline(entry).rstrip() for entry in entries))
        handle.write("\n")


def main() -> int:
    args = parse_args()
    requested_dois = load_requested_dois(args)
    if not requested_dois:
        print("No valid DOIs were provided.", file=sys.stderr)
        return 1

    bib_path = args.bib
    existing_dois, existing_keys = load_library_state(bib_path)
    headers = build_headers(args.email)

    new_entries: list[str] = []
    added = 0
    skipped = 0
    failed = 0

    for doi in requested_dois:
        doi_lower = doi.lower()
        if doi_lower in existing_dois:
            print(f"skip  {doi} (DOI already in library)")
            skipped += 1
            continue

        try:
            entry = fetch_bibtex(doi, headers)
            generated_key = build_entry_key(entry)
            final_key = unique_key(generated_key, existing_keys)
            entry = replace_entry_key(entry, final_key)
            new_entries.append(entry)
            existing_keys.add(final_key)
            existing_dois.add(doi_lower)
            added += 1
            print(f"add   {doi} -> {final_key}")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"fail  {doi} ({exc})", file=sys.stderr)

        if args.pause > 0:
            time.sleep(args.pause)

    if new_entries:
        append_entries(bib_path, new_entries)

    print(
        f"Finished: added={added}, skipped={skipped}, failed={failed}, library={bib_path}"
    )
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
