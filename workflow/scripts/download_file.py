"""Download a single file from a URL."""

from urllib.request import urlretrieve

def run(url: str, out_path: str) -> None:
    urlretrieve(url, out_path)

if __name__ == "__main__":
    run(
        out_path=snakemake.output[0],        
        url=snakemake.params.url
        )