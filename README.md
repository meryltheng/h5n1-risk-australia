# Introduction

This repository contains a workflow to gather the data and run the HPAI H5N1 risk mapping analysis for Australia. 

# Workflow execution

The project requires that [`uv`](https://docs.astral.sh) is installed.

Executing the workflow also requires some environment variables to be set:

`EBIRD_KEY`
: A personal eBird data access key that can be requested at ebird.org/st/request. You can then add this line to your `~/.zshrc`:

```bash
export EBIRD_KEY="your-key-here"
```

The workflow can then be executed by running:

```bash
uv run snakemake
```


# Personal notes

## Create the project and environment

Create the project and environment. In your new, empty folder:

```bash
uv init
uv add snakemake pandas
```

Set up git and `.gitignore`. uv init usually runs `git init` for you and adds a basic `.gitignore`. Add r`esults/`, .`snakemake/` and your raw data to it before your first commit.

Explore your data in a notebook. Install Jupyter as a development-only tool so it doesn't count as a pipeline dependency:

```bash
uv add --dev jupyter
uv run jupyter lab
```

Test your python script and snakemake rule:

```bash
uv run snakemake -n download_abs_locs           # dry run
uv run snakemake --cores 1 download_abs_locs    # real run
uv run snakemake -n data/raw/DAFF/events_by_species_20261007.csv # `scrape_daff_events` has a wildcard constraint (date), so ask for the file its supposed to make instead
```