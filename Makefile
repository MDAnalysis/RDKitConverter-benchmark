# MDAnalysis fork to use
GITHUB_USER ?= cbouy
# Branch of MDAnalysis to install
BRANCH ?= fix-converter
# Use conda or mamba
CONDA ?= conda
# Number of threads to use in parallel
N_WORKERS ?= -1
# Min number of heavy atoms for a molecule
MIN_ATOMS ?= 2
# Max number of heavy atoms for a molecule	 
MAX_ATOMS ?= 50

.ONESHELL:

SHELL := /bin/bash
SET_CONDA_ENV := source $$(conda info --base)/etc/profile.d/conda.sh && conda activate && conda activate rdkitconverter
CHEMBL_SDF := ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_30/chembl_30.sdf.gz

fetch := data/chembl_30.sdf.gz
process := data/chembl_processed_unique.smi.gz
benchmark := data/chembl_failed.smi
report := results/failed_molecules.html

all: report

help:
	@echo 'usage: make <target>'
	@echo
	@echo 'targets:'
	@echo '  help                      Show this help'
	@echo '  install                   Install dependencies'
	@echo '  fetch                     Fetch ChEMBL 30'
	@echo '  process                   Filter, standardize and remove duplicate molecules'
	@echo '  benchmark                 Run the benchmark'
	@echo '  report                    Generate the report'
	@echo '  clean                     Clean all data'
	@echo '  restart                   Clean the benchmark data and report'
	@echo
	@echo 'default behavior: make fetch process benchmark report'

install:
	$(CONDA) env create -f environment.yaml
	@$(SET_CONDA_ENV)
	@pip install git+https://github.com/$(GITHUB_USER)/mdanalysis.git@$(BRANCH)#subdirectory=package

$(fetch):
	@export CHEMBL_SDF=$(CHEMBL_SDF)
	@bash scripts/fetch_chembl.sh

$(process): $(fetch)
	@$(SET_CONDA_ENV)
	@export N_WORKERS=$(N_WORKERS) MIN_ATOMS=$(MIN_ATOMS) MAX_ATOMS=$(MAX_ATOMS)
	@python scripts/process_molecules.py || exit 1
	@bash scripts/drop_duplicates.sh

$(benchmark): $(process)
	@$(SET_CONDA_ENV)
	@export N_WORKERS=$(N_WORKERS)
	@python scripts/benchmark.py

$(report): $(benchmark)
	@$(SET_CONDA_ENV)
	@export N_WORKERS=$(N_WORKERS)
	@python scripts/report.py

clean:
	@rm -rf data/ results/

restart:
	@rm -f data/chembl_failed.smi \
		   data/.failed_count \
		   data/.timing \
		   data/.timestamp
	@rm -rf results/


.PHONY: fetch
fetch: $(fetch)

.PHONY: process
process: $(process)

.PHONY: benchmark
benchmark: $(benchmark)

.PHONY: report
report: $(report)
