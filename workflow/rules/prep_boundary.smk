rule prep_boundary:
    input:
        "data/raw/GADM/gadm41_{iso3}_0.json".format(iso3=config["general"]["country_iso3"]),
    output:
        "data/processed/boundary.gpkg",
    script:
        "../scripts/prep_boundary.py"