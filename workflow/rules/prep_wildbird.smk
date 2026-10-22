rule prep_wildbird:
    input:
        tifs="data/raw/eBird/abundance",
        template="data/processed/template.tif",
    output:
        "data/processed/waterbird_index.tif",
    params:
        temporal_res=config["ebird"]["temporal_res"],
    script:
        "../scripts/prep_wildbird.py"