# run this Makefile on your dataset loader script,
# > make check_file=src/<subdir>/<script_name>.py

.PHONY: quality

datasets_dir := src
examples_dir := examples

# Format source code automatically (one file)

quality:
	black --line-length 119 --target-version py38 $(check_file)
	isort $(check_file)
	flake8 $(check_file) --max-line-length 119 --ignore=E203,W503

# Format source code automatically (all files)
quality_all:
	black --check --line-length 119 --target-version py38 $(datasets_dir)
	isort --check-only $(datasets_dir)
	flake8 $(datasets_dir) --max-line-length 119 --ignore=E203,W503