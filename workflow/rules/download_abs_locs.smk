rule download_abs_locs:
    output:
        shp="data/raw/ABS/SAL_2021_AUST_GDA2020_SHP/SAL_2021_AUST_GDA2020.shp",
    script:
        "../scripts/download_abs_locs.py"