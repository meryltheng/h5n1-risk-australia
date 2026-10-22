rule combine_outbreaks:
    input:
        daff="data/interim/daff_events_geocoded.csv",
        wahis=config["outbreak_data"]["wahis_file"],
    output:
        "data/processed/combined_outbreak_data.csv",
    script:
        "../scripts/combine_outbreaks.py"