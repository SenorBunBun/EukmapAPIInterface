# EukMap API Interface

Resolve scientific names against the **EukMap / UniEuk** eukaryotic taxonomy and
return, for each match, the **full taxonomic lineage up to the root**, with every
ancestor's **name** and its **rank** ("what it is evaluated to be"). Built to feed
a taxonomic-tree viewer.

## What it does

Given a set of scientific names, `eukmap_lineage.py`:

1. **Downloads the whole UniEuk tree once** and caches it (`~/.cache/eukmap`).
   This is the only public way to reconstruct ancestry — EukMap taxa carry no
   parent pointer and there is no ancestry endpoint.
2. **Matches** each name to an EukMap taxon:
   - exact, case-insensitive match against a tree node name (offline), else
   - EukMap's own `/search` endpoint (covers synonyms / alternative names),
     preferring an exact hit among the candidates.
3. **Builds the lineage** by walking parents to the root, attaching each node's
   rank (from the taxon endpoint, cached and de-duplicated).

Output is JSON (default) or TSV, plus a `merged_tree` (the union of all lineages,
nested) that a viewer can render directly. Unmatched names are reported, never
silently dropped.

## Usage

```bash
python3 eukmap_lineage.py "Acrochaetium homorhizum" "Ciliophora" "Fungi"
python3 eukmap_lineage.py --names-file species.txt --format tsv -o out.tsv
printf 'Fungi\nMetazoa\n' | python3 eukmap_lineage.py --stdin
python3 eukmap_lineage.py "Fungi" --from-domain   # start lineage at Eukaryota
python3 eukmap_lineage.py "Fungi" --refresh       # force re-download of the tree
```

Standard-library only — no dependencies. Python 3.9+.

### JSON shape

```jsonc
{
  "results": [
    {
      "query": "Acrochaetium homorhizum",
      "matched": { "id": "2108", "name": "Acrochaetium homorhizum",
                   "rank": "species", "matchType": "exact-tree" },
      "lineage": [
        { "id": "1",    "name": "Life",       "rank": "life" },
        { "id": "2",    "name": "Eukaryota",  "rank": "undefined" },
        { "id": "836",  "name": "Diaphoretickes", "rank": "undefined" },
        // ... down to ...
        { "id": "2108", "name": "Acrochaetium homorhizum", "rank": "species" }
      ]
    }
  ],
  "unmatched": [ { "query": "Paramecium aurelia", "reason": "no match in EukMap" } ],
  "merged_tree": [ /* nested union of all lineages, each node has children[] */ ]
}
```

`matchType` is one of `exact-tree`, `exact-search`, `fuzzy-search`. On a fuzzy or
ambiguous match the `matched` object also carries a `candidates`/`ambiguous` list
so a caller can disambiguate.

## The EukMap public API (reference)

Base: `https://eukmap.unieuk.net/api` — no auth for reads.

| Purpose | Request |
|---|---|
| Search by name | `GET /search/taxonomies/unieuk/?query=NAME` — header `Accept: */*` |
| Full taxon (has `rank`) | `GET /public/taxonomies/unieuk/taxa/{id}` — header `Accept: */*;version=1` |
| Sub-tree / whole tree | `GET /public/taxonomies/unieuk/taxa/{id}/depth/{n}` — header `Accept: */*;version=1` |

Header gotchas (verified against the live API):
- The `/public/*` endpoints require the **version parameter** in `Accept`
  (`*/*;version=1`); a plain `Accept` gives HTTP 422, a wrong media type HTTP 406.
- **Tree** nodes contain `id`, `name`, `children` but **no rank** — rank lives on
  the full **Taxon** object, hence the separate per-node rank fetch.

## Notes & caveats

- EukMap is an 18S / phylogeny-based framework centred on **protists**, so many
  familiar organisms are absent (e.g. *Paramecium aurelia* above) and many deep
  clades legitimately have rank **`undefined`** — they are clades, not Linnaean
  ranks. That `undefined` is EukMap's honest "evaluated" ranking, not an error.
- **NCBI tax IDs:** EukMap does not index its tree by NCBI taxid — mapping is
  name-based. If you start from NCBI tax IDs, resolve each to a scientific name
  first (e.g. `https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{taxid}`), then feed
  the names to this script.
