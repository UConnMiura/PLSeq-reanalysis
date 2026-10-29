PYTHON ?= python3
PLSEQ_TMPDIR ?= $(if $(TMPDIR),$(TMPDIR),/tmp)
MPLCONFIGDIR ?= $(PLSEQ_TMPDIR)/plseq-matplotlib-$(shell id -u)
export MPLCONFIGDIR

.PHONY: check figures-core manuscript manuscript-pdf manuscript-biorxiv clean-generated

check:
	$(PYTHON) tests/check_repository_files.py
	$(PYTHON) tests/check_source_checksums.py
	$(PYTHON) tests/verify_figure_checksums.py
	$(PYTHON) -m py_compile scripts/local/*.py tests/*.py
	@for script in scripts/hpc/*.sbatch scripts/install_git_hooks.sh manuscript/build_reader_pdf.sh manuscript/build_biorxiv_pdf.sh; do bash -n "$$script"; done
	bash tests/check_hpc_portability.sh

figures-core:
	$(PYTHON) scripts/local/build_figure1_reanalysis_workflow.py
	$(PYTHON) scripts/local/build_figure3e_replicate_standalone.py
	$(PYTHON) scripts/local/build_nature_communications_figure3e_pilot.py
	$(PYTHON) scripts/local/render_figure6_read_tracks.py
	$(PYTHON) scripts/local/preflight_nature_figure.py
	$(PYTHON) tests/verify_figure_checksums.py

manuscript: figures-core
	bash manuscript/build_reader_pdf.sh

manuscript-pdf:
	bash manuscript/build_reader_pdf.sh

manuscript-biorxiv:
	bash manuscript/build_biorxiv_pdf.sh

clean-generated:
	@echo "Generated outputs are intentionally retained for publication review."
	@echo "Remove ignored build directories manually only when a clean rebuild is required."
