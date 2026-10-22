rule download_ebird:
    output:
        tifs=directory("data/raw/eBird/abundance"),
        lookup="data/raw/eBird/species_lookup.csv",
    params:
        species=config["ebird"]["species_list"],
        version=config["ebird"]["version"],
        pattern=config["ebird"]["pattern"],
    script:
        "../scripts/download_ebird.py"