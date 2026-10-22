rule predict_kernel_layer:
    input:
        fits="results/kernel/kernel_fits.json",
        table="results/kernel/kernel_comparison.csv",
        risk_set="results/kernel/risk_set.csv",
        template="data/processed/template.tif",
    output:
        layer="results/kernel/kernel_layer.tif",
        map="results/kernel/kernel_layer.png",
    params:
        kernel=config["kernel_layer"]["predict_kernel"],
    script:
        "../scripts/predict_kernel_layer.py"