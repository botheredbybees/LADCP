"""CLI entry point."""

import click


@click.group()
@click.version_option()
def app() -> None:
    """LADCP processing toolkit."""


@app.command()
@click.argument("cast_file", type=click.Path(exists=True))
def process(cast_file: str) -> None:
    """Process a single LADCP cast file."""
    raise NotImplementedError("ingestion layer not yet implemented")


@app.command()
@click.argument("cast_file", type=click.Path(exists=True))
def check(cast_file: str) -> None:
    """Integrate vertical velocity to estimate zmax and zend."""
    raise NotImplementedError("ingestion layer not yet implemented")
