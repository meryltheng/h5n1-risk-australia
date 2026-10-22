"""Download eBird Status & Trends relative abundance for a list of species.

Species are given by common name. Each name is looked up in the eBird
taxonomy to get its 6-letter species code. Species with no matching name, or
with no Status & Trends files, are skipped and recorded in the lookup table.
"""

import os
from pathlib import Path

import pandas as pd
import requests

TAXONOMY_URL = "https://api.ebird.org/v2/ref/taxonomy/ebird"
ST_API = "https://st-download.ebird.org/v1"


def get_species_codes(common_names: list[str]) -> dict[str, str | None]:
    """Map each common name to its eBird species code (None if not found).

    Matching ignores upper/lower case, so "gray teal" matches "Gray Teal".
    """
    response = requests.get(TAXONOMY_URL, params={"fmt": "json"}, timeout=120)
    response.raise_for_status()
    lookup = {taxon["comName"].lower(): taxon["speciesCode"] for taxon in response.json()}
    return {name: lookup.get(name.lower()) for name in common_names}


def list_st_files(species_code: str, version: int, key: str) -> list[str]:
    """List a species' Status & Trends files (empty if it has no products)."""
    response = requests.get(
        f"{ST_API}/list-obj/{version}/{species_code}", params={"key": key}, timeout=60
    )
    if response.status_code == 404:
        return []
    response.raise_for_status()
    return response.json()


def run(
    species: list[str], version: int, pattern: str, out_dir: str, lookup_path: str
) -> None:
    key = os.environ["EBIRD_KEY"]
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    records = []
    for name, code in get_species_codes(species).items():
        if code is None:
            print(f"Skipped {name}: no match in the eBird taxonomy")
            records.append({"common_name": name, "species_code": None, "status": "no name match"})
            continue

        wanted = [obj for obj in list_st_files(code, version, key) if pattern in obj]
        if not wanted:
            print(f"Skipped {name} ({code}): no Status & Trends files matching '{pattern}'")
            records.append({"common_name": name, "species_code": code, "status": "no S&T product"})
            continue

        for obj in wanted:
            response = requests.get(
                f"{ST_API}/fetch", params={"objKey": obj, "key": key}, timeout=300
            )
            response.raise_for_status()
            (out_dir / Path(obj).name).write_bytes(response.content)

        print(f"Downloaded {name} ({code}): {len(wanted)} file(s)")
        records.append({"common_name": name, "species_code": code, "status": "downloaded"})

    lookup = pd.DataFrame(records)
    if not (lookup["status"] == "downloaded").any():
        raise RuntimeError("No species were downloaded; check the names and access key.")
    lookup.to_csv(lookup_path, index=False)


if __name__ == "__main__":
    run(
        species=snakemake.params.species,
        version=snakemake.params.version,
        pattern=snakemake.params.pattern,
        out_dir=snakemake.output.tifs,
        lookup_path=snakemake.output.lookup,
    )