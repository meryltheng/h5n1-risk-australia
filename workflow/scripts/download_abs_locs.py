"""Download ABS Suburbs and Localities (2021) boundaries."""

import tempfile
import zipfile
from pathlib import Path
from urllib.request import urlretrieve

URL = "https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs/edition-3-july-2021-june-2026/access-and-downloads/digital-boundary-files/SAL_2021_AUST_GDA2020_SHP.zip"

def run(shp_path: str, url: str = URL) -> None:
    out_dir = Path(shp_path).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        zip_path = Path(tmp) / "download.zip"
        urlretrieve(url, zip_path)
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(out_dir)

if __name__ == "__main__":
    run(
        shp_path=snakemake.output.shp
        )