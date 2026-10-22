rule make_template:
    input:
        "data/processed/boundary.gpkg",
    output:
        "data/processed/template.tif",
    script:
        "../scripts/make_template.py"