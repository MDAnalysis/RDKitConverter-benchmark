GITHUB_USER ?= cbouy
BRANCH ?= fix-converter
CONDA ?= mamba
N_WORKERS ?= -1
MIN_ATOMS ?= 2
MAX_ATOMS ?= 6

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
	@echo '  help          Show this help'
	@echo '  install       Install dependencies'
	@echo '  fetch         Fetch ChEMBL molecules'
	@echo '  process       Standardize and remove duplicate molecules'
	@echo '  benchmark     Run the benchmark'
	@echo '  report        Generate the report'
	@echo '  clean         Clean all data'
	@echo '  restart       Clean all data except the fetched ChEMBL dataset'

install:
	$(CONDA) env create -f environment.yaml
	@$(SET_CONDA_ENV)
	pip install git+https://github.com/$(GITHUB_USER)/mdanalysis.git@$(BRANCH)#subdirectory=package

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
	python scripts/benchmark.py

$(report): $(benchmark)
	@$(SET_CONDA_ENV)
	@export N_WORKERS=$(N_WORKERS)
	python scripts/report.py

restart:
	@rm -f $$(find data/ -name "*" -type f \
			  | grep -v $(fetch) \
			  | grep -v data/.fetched_count \
			  | grep -v data/.timestamp)
	@rm -rf results/

clean:
	@rm -rf data/ results/

.PHONY: fetch
fetch: $(fetch)

.PHONY: process
process: $(process)

.PHONY: benchmark
benchmark: $(benchmark)

.PHONY: report
report: $(report)
