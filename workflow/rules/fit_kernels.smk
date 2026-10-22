rule fit_kernels:
    input:
        outbreaks="data/processed/combined_outbreak_data.csv",
        template="data/processed/template.tif",
    output:
        risk_set="results/kernel/risk_set.csv",
        fits="results/kernel/kernel_fits.json",
        table="results/kernel/kernel_comparison.csv",
        plot="results/kernel/kernel_best.png",
    script:
        "../scripts/fit_kernels.py"