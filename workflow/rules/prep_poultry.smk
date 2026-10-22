rule prep_poultry:
    input:
        tif="data/raw/FAO/GLW4-2020.D-DA.CHK.tif",
        template="data/processed/template.tif",
    output:
        "data/processed/poultry.tif",
    script:
        "../scripts/prep_poultry.py"