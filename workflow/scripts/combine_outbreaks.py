"""Combine WAHIS outbreak records and geocoded DAFF records into one table.

The two sources overlap, so records with the same locality, date and species
in both are flagged as probable duplicates (`dup_candidate`) rather than
dropped. Downstream steps decide what to do with them.
"""

import pandas as pd


def standardise_daff(daff: pd.DataFrame) -> pd.DataFrame:
    """DAFF records in the common column layout; each record is one case."""
    return pd.DataFrame({
        "source": "DAFF",
        "date": pd.to_datetime(daff["date_sampled"]),  # already "YYYY-MM-DD" text
        "locality": daff["locality"],
        "state": daff["jurisdiction"],
        "lga": daff["lga"],
        "species": daff["scientific_name"],
        "latitude": daff["latitude"],
        "longitude": daff["longitude"],
        "id": daff["record_id"].astype(str),
        "cases": 1,
    })


def standardise_wahis(wahis: pd.DataFrame) -> pd.DataFrame:
    """WAHIS records in the common column layout (no state or LGA in WAHIS)."""
    return pd.DataFrame({
        "source": "WAHIS",
        "date": pd.to_datetime(wahis["startDate"]),
        "locality": wahis["locationName"],
        "state": pd.NA,
        "lga": pd.NA,
        "species": wahis["species"],
        "latitude": wahis["latitude"],
        "longitude": wahis["longitude"],
        "id": wahis["outbreakId"].astype(str),
        "cases": wahis["totalCases"],
    })


def flag_duplicates(combined: pd.DataFrame) -> pd.DataFrame:
    """Count the sources reporting each locality + date + species, and flag those with more than one."""
    # lower case and single spaces, so "Port  Macdonnell" and "port macdonnell" match
    combined["loc"] = combined["locality"].str.lower().str.split().str.join(" ")
    combined["n_sources"] = (
        combined.groupby(["loc", "date", "species"], dropna=False)["source"]
        .transform("nunique")
    )
    combined["dup_candidate"] = combined["n_sources"] > 1
    return combined


def run(daff_path: str, wahis_path: str, out_path: str) -> None:
    daff = standardise_daff(pd.read_csv(daff_path))
    wahis = standardise_wahis(pd.read_csv(wahis_path))

    combined = pd.concat([daff, wahis], ignore_index=True)
    combined = flag_duplicates(combined)

    print(
        f"{len(daff)} DAFF + {len(wahis)} WAHIS records; "
        f"{combined['dup_candidate'].sum()} flagged as probable duplicates"
    )
    combined.to_csv(out_path, index=False, date_format="%Y-%m-%d")


if __name__ == "__main__":
    run(
        daff_path=snakemake.input.daff,
        wahis_path=snakemake.input.wahis,
        out_path=snakemake.output[0],
    )