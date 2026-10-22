rule download_poultry:
    output:
        "data/raw/FAO/GLW4-2020.D-DA.CHK.tif",
    params:
        url="https://storage.googleapis.com/fao-gismgr-glw4-2020-data/DATA/GLW4-2020/MAPSET/D-DA/GLW4-2020.D-DA.CHK.tif",
    script:
        "../scripts/download_file.py"