.PHONY: test test-integration lint fmt docker-build

test:
	uv run pytest

test-integration:
	TEST_DATA_DIR=test_data uv run pytest -m integration -v

lint:
	uv run ruff check src tests
	uv run ruff format --check src tests

fmt:
	uv run ruff format src tests

docker-build:
	docker build -t ladcp .
