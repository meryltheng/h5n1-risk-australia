rule download_boundary:
    output:
        "data/raw/GADM/gadm41_{iso3}_0.json".format(iso3=config["general"]["country_iso3"]),
    params:
        url="https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_{iso3}_0.json".format(
            iso3=config["general"]["country_iso3"]
        ),
    script:
        "../scripts/download_file.py"