"""Command-line entry point for cfh."""

from __future__ import annotations

import click

from cfh.genes.registry import available_genes, load_gene_config


@click.group()
def main() -> None:
    """Cancer Fusion Hotspots command-line tools."""


@main.command("list-genes")
def list_genes() -> None:
    """List gene symbols with a registered config."""
    for symbol in available_genes():
        click.echo(symbol)


@main.command("show-gene")
@click.argument("gene_symbol")
def show_gene(gene_symbol: str) -> None:
    """Print the resolved GeneConfig for a gene symbol."""
    config = load_gene_config(gene_symbol)
    click.echo(config.model_dump_json(indent=2))


if __name__ == "__main__":
    main()
