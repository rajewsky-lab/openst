# run this Makefile on your dataset loader script,
# > make check_file=openst/<subdir>/<script_name>.py

.PHONY: quality

source_dir := openst
examples_dir := examples

# Format source code automatically (one file)

quality:
	black --line-length 119 --target-version py38 $(check_file)
	isort $(check_file)
	flake8 $(check_file) --max-line-length 119 --ignore=E203,W503

# Format source code automatically (all files)
quality_all:
	black --check --line-length 119 --target-version py38 $(source_dir)
	isort --check-only $(source_dir)
	flake8 $(source_dir) --max-line-length 119 --ignore=E203,W503

# Format source code automatically (all files)
fix_all:
	black --line-length 119 --target-version py38 $(source_dir)
	isort $(source_dir)
	flake8 $(source_dir) --max-line-length 119 --ignore=E203,W503