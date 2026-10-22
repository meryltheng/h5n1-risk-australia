rule geocode_daff:
    input:
        daff=f"data/raw/DAFF/events_by_species_{config['outbreak_data']['daff_snapshot']}.csv",
        sal="data/raw/ABS/SAL_2021_AUST_GDA2020_SHP/SAL_2021_AUST_GDA2020.shp",
    output:
        "data/interim/daff_events_geocoded.csv",
    script:
        "../scripts/geocode_daff.py"