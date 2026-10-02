"""GWAS Catalog REST API v2 client.

Raw pages are cached exactly as returned. Pagination parameter names and the response
envelope for v2 have NOT been verified against a live response yet (see config.toml),
so this module fails loudly instead of guessing.
"""
import hashlib
import json
import time
from datetime import datetime, timezone

import requests

from ..config import data_path


class GwasApiError(RuntimeError):
    pass


def extract_records(payload):
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict):
        embedded = payload.get("_embedded")
        if isinstance(embedded, dict):
            for value in embedded.values():
                if isinstance(value, list):
                    return value
        for key in ("content", "data", "associations", "results", "items"):
            if isinstance(payload.get(key), list):
                return payload[key]
    keys = list(payload) if isinstance(payload, dict) else type(payload).__name__
    raise GwasApiError(f"Could not locate the record list in the response. Top-level keys: {keys}")


def _get(session, url, params, retries=3):
    for attempt in range(retries):
        try:
            resp = session.get(url, params=params, timeout=60, headers={"Accept": "application/json"})
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code in (429, 500, 502, 503, 504) and attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            raise GwasApiError(f"HTTP {resp.status_code} for {resp.url}: {resp.text[:300]}")
        except requests.RequestException as exc:
            if attempt == retries - 1:
                raise GwasApiError(f"Request failed: {exc}") from exc
            time.sleep(2 ** attempt)


def _request_spec(cfg, efo_trait):
    g = cfg["gwas"]
    return {
        "url": f"{g['api_base'].rstrip('/')}/{g['endpoint']}",
        "trait_param": g["trait_param"],
        "efo_trait": efo_trait,
        "page_size": g["page_size"],
    }


def inspect_first_page(cfg, efo_trait=None):
    """Fetch one page and return (request_spec, raw_payload) so the schema can be checked by eye."""
    g = cfg["gwas"]
    efo_trait = efo_trait or cfg["disease"]["efo_trait"]
    spec = _request_spec(cfg, efo_trait)
    params = {g["trait_param"]: efo_trait, g["size_param"]: g["page_size"], g["page_param"]: g["first_page"]}
    payload = _get(requests.Session(), spec["url"], params)
    path = data_path(cfg, "data", "raw", "gwas", "inspect_first_page.json")
    path.write_text(json.dumps({"request": spec, "payload": payload}, indent=2))
    return spec, payload, path


def fetch_associations(cfg, efo_trait=None, refresh=False, max_pages=1000):
    """Return (records, meta). Raw pages are cached under data/raw/gwas/."""
    g = cfg["gwas"]
    efo_trait = efo_trait or cfg["disease"]["efo_trait"]
    spec = _request_spec(cfg, efo_trait)
    key = hashlib.sha1(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:10]
    cache = data_path(cfg, "data", "raw", "gwas", f"associations_{key}.json")

    if cache.exists() and not refresh:
        blob = json.loads(cache.read_text())
    else:
        session = requests.Session()
        pages, previous = [], None
        for page_no in range(g["first_page"], g["first_page"] + max_pages):
            params = {g["trait_param"]: efo_trait, g["size_param"]: g["page_size"], g["page_param"]: page_no}
            payload = _get(session, spec["url"], params)
            records = extract_records(payload)
            if not records:
                break
            if previous is not None and records == previous:
                raise GwasApiError("Two consecutive pages were identical: the page parameter is probably being ignored.")
            pages.append(payload)
            previous = records
            if len(records) < g["page_size"]:
                break
        else:
            raise GwasApiError(f"Hit max_pages={max_pages}; refusing to silently truncate.")
        blob = {
            "request": spec,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "pages": pages,
        }
        cache.write_text(json.dumps(blob))

    records = [rec for page in blob["pages"] for rec in extract_records(page)]
    meta = {"request": blob["request"], "retrieved_utc": blob["retrieved_utc"],
            "n_pages": len(blob["pages"]), "n_records": len(records), "cache_file": str(cache)}
    return records, meta
