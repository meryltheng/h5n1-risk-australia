rule combine_risk:
    input:
        kernel="results/kernel/kernel_layer.tif",
        waterbird="data/processed/waterbird_index.tif",
        poultry="data/processed/poultry.tif",
    output:
        risk="results/risk/risk_layers.tif",
        map="results/risk/risk_maps.png",
    params:
        risk=config["risk"],
    script:
        "../scripts/combine_risk.py"