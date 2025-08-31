.PHONY: build poetry-install dev-install run-tests run-linter run-linter-fix run-pylint run-mypy run-deptry clean

build:
	cmake -S . -B build -DBUILD_PYTHON_MODULE=ON
	cmake --build build

poetry-install: dev-install
	@echo "Installation complete. Please activate your poetry shell with 'poetry shell'"

dev-install: build
	poetry install
	chmod +x build/bin/mkdssp

run-tests:
	PYTHONPATH=. poetry run pytest --cov=dsspy --cov-report=term-missing

run-linter:
	poetry run ruff check .

run-linter-fix:
	poetry run ruff check . --fix

run-pylint:
	poetry run pylint dsspy tests

run-mypy:
	poetry run mypy --explicit-package-bases dsspy

run-deptry:
	poetry run deptry .

clean:
	rm -rf build
