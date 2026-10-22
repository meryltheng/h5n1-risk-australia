rule scrape_daff_events:
    output:
        csv="data/raw/DAFF/events_by_species_{date}.csv",
    wildcard_constraints:
        date=r"\d{8}",  # exactly 8 digits, e.g. 20261007
    log:
        "logs/scrape_daff_events_{date}.log",
    script:
        "../scripts/scrape_daff_events.py"