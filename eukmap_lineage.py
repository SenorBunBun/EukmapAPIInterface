#!/usr/bin/env python3
"""
eukmap_lineage.py
=================

Resolve a set of scientific names against the EukMap (UniEuk) taxonomy and, for
every one that matches, return its full taxonomic lineage all the way up to the
root ("Life" / "Eukaryota"), with each ancestor's name and its EukMap rank
("what it is evaluated to be", e.g. species / genus / order / undefined-clade).

Pipeline (the "step 2 -> 3" the request asked for)
--------------------------------------------------
  1. Download the whole UniEuk tree once (the only public way to reconstruct
     ancestry -- EukMap taxa carry no parent pointer), cache it on disk, and
     build id -> (name, parent) and name -> id indexes.
  2. MATCH each input name to an EukMap taxon:
        a. exact, case-insensitive match against a tree node name (offline, fast)
        b. otherwise EukMap's own /search endpoint (covers synonyms / alt names),
           preferring an exact name hit among the candidates.
  3. LINEAGE: walk the parent index from the matched node up to the root, then
     attach each node's rank (fetched from the taxon endpoint, cached & deduped).

Output: JSON (default) or TSV, plus a merged tree of all lineages for a viewer.

Why the EukMap tree is unusual
------------------------------
EukMap is an 18S/phylogeny-based framework focused on protists, so many familiar
species are absent and many deep clades legitimately have rank "undefined"
(they are clades, not Linnaean ranks). Unmatched names are reported, not hidden.

The public API
--------------
  Base            : https://eukmap.unieuk.net/api
  Search  (public): GET /search/taxonomies/{tax}/?query=NAME   Accept: */*
  Taxon   (public): GET /public/taxonomies/{tax}/taxa/{id}     Accept: */*;version=1
  Tree    (public): GET /public/taxonomies/{tax}/taxa/{id}/depth/{n}
                                                                Accept: */*;version=1
No authentication is needed for reads. The versioned Accept header is required
on the /public/* endpoints (a bare Accept yields HTTP 422 "does not specify
version"; a wrong media type yields HTTP 406).

Usage
-----
  python3 eukmap_lineage.py "Acrochaetium homorhizum" "Ciliophora" "Fungi"
  python3 eukmap_lineage.py --names-file species.txt --format tsv -o out.tsv
  echo -e "Fungi\\nMetazoa" | python3 eukmap_lineage.py --stdin
  python3 eukmap_lineage.py "Fungi" --refresh          # force re-download tree
  python3 eukmap_lineage.py "Fungi" --from-domain      # start lineage at Eukaryota

Only depends on the Python 3 standard library.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
DEFAULT_BASE = "https://eukmap.unieuk.net/api"
DEFAULT_TAXONOMY = "unieuk"
ROOT_ID = "1"                    # "Life"
DEFAULT_DOMAIN_NAME = "Eukaryota"  # where --from-domain starts the lineage

CACHE_DIR = os.path.join(os.path.expanduser("~"), ".cache", "eukmap")
TREE_TTL_SECONDS = 24 * 3600     # re-download the full tree at most once a day

# A tree snapshot shipped with the repo, used offline when no fresh cache exists
# (built with: curl .../taxa/1/depth/50 | gzip -9 > data/unieuk_tree.json.gz).
BUNDLED_TREE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "data", "unieuk_tree.json.gz")

# Accept headers the API demands (discovered empirically against the live API).
ACCEPT_SEARCH = "*/*"
ACCEPT_VERSIONED = "*/*;version=1"

HTTP_RETRIES = 4
HTTP_BACKOFF = 2.0               # seconds: 2, 4, 8, 16
RANK_WORKERS = 8                 # parallel rank lookups


# --------------------------------------------------------------------------- #
# HTTP helper
# --------------------------------------------------------------------------- #
def _http_get(url: str, accept: str, timeout: int = 60) -> bytes:
    """GET with retries + exponential backoff. Raises on final failure."""
    last = None
    for attempt in range(HTTP_RETRIES):
        try:
            req = urllib.request.Request(url, headers={"Accept": accept})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            last = exc
            if attempt < HTTP_RETRIES - 1:
                time.sleep(HTTP_BACKOFF * (2 ** attempt))
    raise RuntimeError(f"GET failed after {HTTP_RETRIES} tries: {url}\n  {last}")


# --------------------------------------------------------------------------- #
# NCBI resolver (via the `ncbi-taxonomist` CLI)
# --------------------------------------------------------------------------- #
# EukMap has NO NCBI-taxid crosswalk in its public data, so an NCBI taxid can only
# reach EukMap through a scientific NAME. `ncbi-taxonomist resolve` turns a taxid
# into its authoritative NCBI name AND full NCBI lineage; we use the name to match
# EukMap and, if the exact organism is absent, walk the NCBI lineage upward to the
# nearest ancestor that EukMap does know.
def resolve_ncbi_taxids(taxids: list[str], email: str | None = None,
                        quiet: bool = False) -> dict[str, dict]:
    """taxid -> {name, ncbi_lineage:[{taxid,name,rank}] most-specific-first}.

    Requires the `ncbi-taxonomist` package (pip install ncbi-taxonomist) and
    network access to NCBI Entrez.
    """
    taxids = [str(t) for t in taxids]
    if not taxids:
        return {}
    cmd = ["ncbi-taxonomist", "resolve", "-r", "-t", *taxids]
    if email:
        cmd += ["-e", email]
    if not quiet:
        print(f"[ncbi] resolving {len(taxids)} taxid(s) via ncbi-taxonomist ...",
              file=sys.stderr)
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    except FileNotFoundError:
        raise RuntimeError("ncbi-taxonomist not found; run: pip install ncbi-taxonomist")
    if proc.returncode != 0:
        raise RuntimeError(f"ncbi-taxonomist failed: {proc.stderr.strip()[:400]}")

    out: dict[str, dict] = {}
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        rec = json.loads(line)
        q = str(rec.get("query"))
        taxon = rec.get("taxon") or {}
        lineage = [{"taxid": str(n.get("taxid")), "name": n.get("name"),
                    "rank": n.get("rank")}
                   for n in (rec.get("lineage") or [])]
        out[q] = {"name": taxon.get("name"), "rank": taxon.get("rank"),
                  "ncbi_lineage": lineage}
    return out


# --------------------------------------------------------------------------- #
# EukMap client
# --------------------------------------------------------------------------- #
class EukMap:
    def __init__(self, base: str = DEFAULT_BASE, taxonomy: str = DEFAULT_TAXONOMY,
                 cache_dir: str = CACHE_DIR):
        self.base = base.rstrip("/")
        self.taxonomy = taxonomy
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)

        # Indexes built from the full tree.
        self.name_of: dict[str, str] = {}          # id   -> node name
        self.parent_of: dict[str, str | None] = {}  # id   -> parent id (None at root)
        self.id_by_name: dict[str, list[str]] = {}  # lower(name) -> [ids]

        # Rank cache (persisted): id -> rank string.
        self._rank_cache_path = os.path.join(cache_dir, f"ranks_{taxonomy}.json")
        self._rank_cache: dict[str, str] = self._load_json(self._rank_cache_path, {})

    # ---- tiny json-on-disk helpers ---------------------------------------- #
    @staticmethod
    def _load_json(path: str, default):
        try:
            with open(path, "r", encoding="utf-8") as fh:
                return json.load(fh)
        except (OSError, json.JSONDecodeError):
            return default

    def _save_rank_cache(self) -> None:
        tmp = self._rank_cache_path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(self._rank_cache, fh)
        os.replace(tmp, self._rank_cache_path)

    # ---- full tree + indexes ---------------------------------------------- #
    def _tree_cache_path(self) -> str:
        return os.path.join(self.cache_dir, f"tree_{self.taxonomy}.json")

    @staticmethod
    def _read_maybe_gzip(path: str) -> bytes:
        with open(path, "rb") as fh:
            data = fh.read()
        return gzip.decompress(data) if path.endswith(".gz") else data

    def load_tree(self, refresh: bool = False, quiet: bool = False) -> None:
        """Load the full tree and build id/name/parent indexes.

        Source preference: fresh on-disk cache -> bundled repo snapshot -> live API.
        """
        path = self._tree_cache_path()
        fresh = (not refresh and os.path.exists(path)
                 and (time.time() - os.path.getmtime(path)) < TREE_TTL_SECONDS)
        if fresh:
            if not quiet:
                print(f"[eukmap] using cached tree {path}", file=sys.stderr)
            raw = self._read_maybe_gzip(path)
        elif not refresh and os.path.exists(BUNDLED_TREE):
            if not quiet:
                print(f"[eukmap] using bundled tree snapshot {BUNDLED_TREE}", file=sys.stderr)
            raw = self._read_maybe_gzip(BUNDLED_TREE)
        else:
            url = (f"{self.base}/public/taxonomies/{self.taxonomy}"
                   f"/taxa/{ROOT_ID}/depth/50")
            if not quiet:
                print(f"[eukmap] downloading full tree (~24 MB, once) ...", file=sys.stderr)
            raw = _http_get(url, ACCEPT_VERSIONED, timeout=120)
            tmp = path + ".tmp"
            with open(tmp, "wb") as fh:
                fh.write(raw)
            os.replace(tmp, path)

        tree = json.loads(raw)
        self._index(tree, parent=None)
        if not quiet:
            print(f"[eukmap] indexed {len(self.name_of)} taxa", file=sys.stderr)

    def _index(self, node: dict, parent: str | None) -> None:
        """Iteratively populate name_of / parent_of / id_by_name from the nested tree."""
        stack = [(node, parent)]
        while stack:
            n, par = stack.pop()
            nid = n["id"]["id"]
            name = n.get("name")
            self.name_of[nid] = name
            self.parent_of[nid] = par
            if name:
                self.id_by_name.setdefault(name.strip().lower(), []).append(nid)
            for child in (n.get("children") or []):
                stack.append((child, nid))

    # ---- matching (step 2) ------------------------------------------------ #
    def match(self, query: str) -> dict:
        """
        Match one scientific name to an EukMap taxon.
        Returns a dict: {id, name, rank, matchType, candidates?} or {error}.
        matchType in: exact-tree | exact-search | fuzzy-search
        """
        q = query.strip()
        key = q.lower()

        # (a) exact, offline match against a tree node name
        ids = self.id_by_name.get(key)
        if ids:
            nid = ids[0]
            res = {"query": query, "id": nid, "name": self.name_of[nid],
                   "rank": self.rank_of(nid), "matchType": "exact-tree"}
            if len(ids) > 1:
                res["ambiguous"] = [{"id": i, "name": self.name_of[i]} for i in ids]
            return res

        # (b) EukMap's own search (handles synonyms / alternative names)
        try:
            candidates = self.search(q)
        except RuntimeError as exc:
            return {"query": query, "error": f"search failed: {exc}"}
        if not candidates:
            return {"query": query, "error": "no match in EukMap"}

        exact = [c for c in candidates if (c["name"] or "").strip().lower() == key]
        chosen = exact[0] if exact else candidates[0]
        res = {
            "query": query,
            "id": chosen["id"],
            "name": chosen["name"],
            "rank": chosen.get("rank") or self.rank_of(chosen["id"]),
            "matchType": "exact-search" if exact else "fuzzy-search",
        }
        if not exact or len(candidates) > 1:
            res["candidates"] = [
                {"id": c["id"], "name": c["name"], "rank": c.get("rank")}
                for c in candidates[:10]
            ]
        return res

    def search(self, name: str) -> list[dict]:
        """Call the public /search endpoint; return simplified candidate dicts."""
        url = (f"{self.base}/search/taxonomies/{self.taxonomy}/"
               f"?query={urllib.parse.quote(name)}")
        data = json.loads(_http_get(url, ACCEPT_SEARCH, timeout=45) or b"[]")
        out = []
        for t in data:
            names = t.get("taxonNames") or []
            nm = names[0].get("name") if names else None
            out.append({"id": t["id"]["id"], "name": nm, "rank": t.get("rank")})
        return out

    # ---- ranks ------------------------------------------------------------ #
    def rank_of(self, taxon_id: str) -> str | None:
        """Rank for a single taxon, cached on disk (the tree omits rank)."""
        if taxon_id in self._rank_cache:
            return self._rank_cache[taxon_id]
        url = f"{self.base}/public/taxonomies/{self.taxonomy}/taxa/{taxon_id}"
        try:
            t = json.loads(_http_get(url, ACCEPT_VERSIONED, timeout=45))
        except RuntimeError:
            return None
        rank = t.get("rank")
        self._rank_cache[taxon_id] = rank
        return rank

    def prefetch_ranks(self, ids: list[str]) -> None:
        """Fetch ranks for many ids in parallel, populating the cache."""
        todo = [i for i in dict.fromkeys(ids) if i not in self._rank_cache]
        if not todo:
            return
        with concurrent.futures.ThreadPoolExecutor(max_workers=RANK_WORKERS) as pool:
            pool.map(self.rank_of, todo)
        self._save_rank_cache()

    # ---- lineage (step 3) ------------------------------------------------- #
    def lineage_ids(self, taxon_id: str) -> list[str]:
        """Ordered ids from root down to (and including) taxon_id."""
        chain: list[str] = []
        seen = set()
        cur: str | None = taxon_id
        while cur is not None and cur not in seen:
            seen.add(cur)
            chain.append(cur)
            cur = self.parent_of.get(cur)
        chain.reverse()
        return chain

    def ancestor_names(self, taxon_id: str) -> set[str]:
        """Lowercased names of a taxon and all its EukMap ancestors."""
        return {(self.name_of.get(i) or "").strip().lower()
                for i in self.lineage_ids(taxon_id)}

    def lineage(self, taxon_id: str, from_domain: bool = False) -> list[dict]:
        """Full lineage as [{id, name, rank}], root-first."""
        ids = self.lineage_ids(taxon_id)
        self.prefetch_ranks(ids)
        rows = [{"id": i, "name": self.name_of.get(i), "rank": self.rank_of(i)}
                for i in ids]
        if from_domain:
            for idx, r in enumerate(rows):
                if r["name"] == DEFAULT_DOMAIN_NAME:
                    return rows[idx:]
        return rows


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def resolve_names(client: EukMap, names: list[str], from_domain: bool) -> dict:
    results, unmatched = [], []
    for name in names:
        m = client.match(name)
        if "error" in m:
            unmatched.append({"query": name, "reason": m["error"]})
            continue
        entry = {
            "query": m["query"],
            "matched": {"id": m["id"], "name": m["name"], "rank": m["rank"],
                        "matchType": m["matchType"]},
            "lineage": client.lineage(m["id"], from_domain=from_domain),
        }
        for opt in ("candidates", "ambiguous"):
            if opt in m:
                entry["matched"][opt] = m[opt]
        results.append(entry)
    return {"results": results, "unmatched": unmatched,
            "merged_tree": merged_tree(results)}


# Names so universal they don't corroborate a region (shared by everything).
_UNIVERSAL_NAMES = {"", "life", "biota", "cellular organisms", "eukaryota",
                    "root", "domain"}


def bridge_ncbi(client: EukMap, rec: dict) -> tuple[dict | None, str | None]:
    """Bridge one resolved NCBI record to a EukMap taxon id, homonym-safe.

    Matches ONLY the organism's own scientific name in EukMap -- there is no
    ancestor fallback, so a taxid whose exact organism is absent from EukMap is
    reported unmatched. When several EukMap taxa share the name (homonyms, e.g.
    NCBI animal 'Vertebrata' vs EukMap's red-algal genus 'Vertebrata'), the
    candidate whose EukMap ancestry shares an informative (non-universal) name
    with the NCBI lineage is chosen. A single unambiguous hit is accepted and
    flagged corroborated=False.
    """
    key = (rec.get("name") or "").strip().lower()
    candidates = client.id_by_name.get(key) or []
    if not candidates:
        return None, None

    ncbi_name_set = {(n["name"] or "").strip().lower()
                     for n in (rec.get("ncbi_lineage") or [])}
    ncbi_name_set.add(key)

    for cand in candidates:
        corrob = (client.ancestor_names(cand) & ncbi_name_set) \
            - _UNIVERSAL_NAMES - {key}
        if corrob:
            return ({"ncbi_name": rec["name"], "ncbi_rank": rec.get("rank"),
                     "via": "organism", "corroborated": True,
                     "corroborated_by": sorted(corrob)[:5]}, cand)

    if len(candidates) == 1:
        return ({"ncbi_name": rec["name"], "ncbi_rank": rec.get("rank"),
                 "via": "organism", "corroborated": False}, candidates[0])
    # ambiguous homonyms, none corroborated -> refuse to guess
    return None, None


def resolve_ncbi(client: EukMap, taxids: list[str], from_domain: bool,
                 email: str | None) -> dict:
    """NCBI taxids -> EukMap lineages, bridging on scientific name.

    For each taxid: resolve its NCBI name + lineage, then
      1. match the organism's own name in EukMap; else
      2. walk the NCBI lineage upward and take the nearest ancestor EukMap knows
         (recording which NCBI rank/name the bridge happened at).
    """
    resolved = resolve_ncbi_taxids(taxids, email=email)
    results, unmatched = [], []
    for taxid in taxids:
        rec = resolved.get(str(taxid))
        if not rec or not rec.get("name"):
            unmatched.append({"query": str(taxid), "reason": "unknown NCBI taxid"})
            continue

        bridge, eukmap_id = bridge_ncbi(client, rec)
        if not eukmap_id:
            unmatched.append({"query": str(taxid), "ncbi_name": rec["name"],
                              "reason": "organism not present in EukMap "
                                        "(no direct species/name match)"})
            continue

        results.append({
            "query": str(taxid),
            "ncbi": {"taxid": str(taxid), "name": rec["name"], "rank": rec.get("rank")},
            "bridge": bridge,
            "matched": {"id": eukmap_id, "name": client.name_of.get(eukmap_id),
                        "rank": client.rank_of(eukmap_id)},
            "lineage": client.lineage(eukmap_id, from_domain=from_domain),
        })
    return {"results": results, "unmatched": unmatched,
            "merged_tree": merged_tree(results)}


def merged_tree(results: list[dict]) -> list[dict]:
    """Union all lineages into one nested tree (handy for a taxonomy viewer)."""
    roots: dict[str, dict] = {}
    index: dict[str, dict] = {}
    for entry in results:
        parent_children = roots
        for node in entry["lineage"]:
            nid = node["id"]
            existing = index.get(nid)
            if existing is None:
                existing = {"id": nid, "name": node["name"],
                            "rank": node["rank"], "children": {}}
                index[nid] = existing
                parent_children[nid] = existing
            parent_children = existing["children"]

    def to_list(d: dict) -> list[dict]:
        out = []
        for n in d.values():
            out.append({"id": n["id"], "name": n["name"], "rank": n["rank"],
                        "children": to_list(n["children"])})
        return out

    return to_list(roots)


def to_tsv(payload: dict) -> str:
    lines = ["query\tmatch_type\tdepth\trank\tname\ttaxon_id"]
    for e in payload["results"]:
        mt = e["matched"].get("matchType") or (e.get("bridge") or {}).get("via", "")
        for depth, node in enumerate(e["lineage"]):
            lines.append("\t".join([
                e["query"], mt, str(depth),
                node["rank"] or "", node["name"] or "", node["id"],
            ]))
    for u in payload["unmatched"]:
        lines.append("\t".join([u["query"], "UNMATCHED", "", "", u["reason"], ""]))
    return "\n".join(lines) + "\n"


def collect_names(args: argparse.Namespace) -> list[str]:
    names: list[str] = list(args.names)
    if args.names_file:
        with open(args.names_file, encoding="utf-8") as fh:
            names += [ln.strip() for ln in fh if ln.strip()]
    if args.stdin:
        names += [ln.strip() for ln in sys.stdin if ln.strip()]
    # de-duplicate, preserve order
    seen, uniq = set(), []
    for n in names:
        if n.lower() not in seen:
            seen.add(n.lower())
            uniq.append(n)
    return uniq


def collect_taxids(args: argparse.Namespace) -> list[str]:
    raw: list[str] = []
    if args.ncbi_ids:
        raw += args.ncbi_ids.replace(",", " ").split()
    if args.ncbi_ids_file:
        with open(args.ncbi_ids_file, encoding="utf-8") as fh:
            for ln in fh:
                raw += ln.replace(",", " ").split()
    seen, uniq = set(), []
    for t in raw:
        if t and t not in seen:
            seen.add(t)
            uniq.append(t)
    return uniq


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Resolve scientific names to full EukMap taxonomic lineages.")
    ap.add_argument("names", nargs="*", help="scientific names (quote multi-word)")
    ap.add_argument("--names-file", help="file with one scientific name per line")
    ap.add_argument("--stdin", action="store_true", help="also read names from stdin")
    ap.add_argument("--ncbi-ids", help="comma/space separated NCBI taxids; resolves "
                    "them to names via ncbi-taxonomist, then matches EukMap")
    ap.add_argument("--ncbi-ids-file", help="file with one NCBI taxid per line")
    ap.add_argument("--email", help="contact email for NCBI Entrez (recommended)")
    ap.add_argument("--format", choices=["json", "tsv"], default="json")
    ap.add_argument("-o", "--output", help="write to file instead of stdout")
    ap.add_argument("--from-domain", action="store_true",
                    help="start each lineage at the domain (Eukaryota) rather than root (Life)")
    ap.add_argument("--refresh", action="store_true", help="force re-download of the full tree")
    ap.add_argument("--base", default=DEFAULT_BASE, help="API base URL")
    ap.add_argument("--taxonomy", default=DEFAULT_TAXONOMY, help="taxonomy id (default: unieuk)")
    args = ap.parse_args()

    taxids = collect_taxids(args)
    names = collect_names(args)
    if not names and not taxids:
        ap.error("no input (pass names, --names-file, --stdin, or --ncbi-ids)")

    client = EukMap(base=args.base, taxonomy=args.taxonomy)
    client.load_tree(refresh=args.refresh)

    if taxids:
        payload = resolve_ncbi(client, taxids, from_domain=args.from_domain,
                               email=args.email)
        n_in = len(taxids)
    else:
        payload = resolve_names(client, names, from_domain=args.from_domain)
        n_in = len(names)

    text = (to_tsv(payload) if args.format == "tsv"
            else json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    if args.output:
        with open(args.output, "w", encoding="utf-8") as fh:
            fh.write(text)
        print(f"[eukmap] wrote {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(text)

    matched = len(payload["results"])
    print(f"[eukmap] matched {matched}/{n_in}; "
          f"unmatched {len(payload['unmatched'])}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
