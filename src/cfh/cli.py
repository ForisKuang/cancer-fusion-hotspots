"""Command-line entry point for cfh."""

from __future__ import annotations

from pathlib import Path

import click

from cfh.genes.registry import available_genes, load_gene_config
from cfh.real_benchmark import RealBenchmarkError, run_analysis, run_real_benchmark, write_outputs


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


@main.command("real-benchmark")
@click.argument("gene_symbol")
@click.argument("study_id")
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("runs"),
    show_default=True,
)
@click.option("--n-permutations", type=click.IntRange(min=1), default=1_000, show_default=True)
@click.option("--output-stem", help="Override the output filename stem.")
@click.option(
    "--pdf/--no-pdf",
    default=True,
    show_default=True,
    help="Also render a self-contained report.pdf for human reviewers.",
)
def real_benchmark(
    gene_symbol: str,
    study_id: str,
    output_dir: Path,
    n_permutations: int,
    output_stem: str | None,
    pdf: bool,
) -> None:
    """Run the prototype live cBioPortal/Genome Nexus benchmark."""
    try:
        run = run_real_benchmark(gene_symbol, study_id, n_permutations=n_permutations)
        paths = write_outputs(
            run,
            output_dir,
            output_stem=output_stem,
            pdf=pdf,
            cli_args=[
                "real-benchmark",
                gene_symbol,
                study_id,
                "--output-dir",
                str(output_dir),
                "--n-permutations",
                str(n_permutations),
            ],
        )
    except RealBenchmarkError as exc:
        raise click.ClickException(str(exc)) from None
    except Exception as exc:
        raise click.ClickException(f"Benchmark failed: {type(exc).__name__}: {exc}") from None
    fisher_p_value = run.summary["fisher_p_value"]
    fisher_display = "unavailable" if fisher_p_value is None else f"{fisher_p_value:.6g}"
    click.echo(
        f"Analyzed {run.summary['total_fusions']} {run.gene_symbol} fusions; "
        f"mapped={run.summary['mapped_fusions']}, "
        f"in-frame={run.summary['in_frame_count']}, "
        f"domain-retained={run.summary['kinase_retained_count']}, "
        f"Fisher p={fisher_display}"
    )
    for warning in run.warnings:
        click.echo(f"Warning: {warning}", err=True)
    for kind, path in paths.items():
        click.echo(f"{kind}: {path}")


@main.command("analyze")
@click.argument("gene_symbol")
@click.argument("study_id")
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("runs"),
    show_default=True,
)
@click.option("--n-permutations", type=click.IntRange(min=1), default=1_000, show_default=True)
@click.option(
    "--pdf/--no-pdf",
    default=True,
    show_default=True,
    help="Also render a self-contained report.pdf for human reviewers.",
)
def analyze(
    gene_symbol: str, study_id: str, output_dir: Path, n_permutations: int, pdf: bool
) -> None:
    """Run all registered algorithms for a configured gene and live study."""
    try:
        run = run_analysis(gene_symbol, study_id, n_permutations=n_permutations)
        paths = write_outputs(
            run,
            output_dir,
            pdf=pdf,
            cli_args=[
                "analyze",
                gene_symbol,
                study_id,
                "--output-dir",
                str(output_dir),
                "--n-permutations",
                str(n_permutations),
            ],
        )
    except RealBenchmarkError as exc:
        raise click.ClickException(str(exc)) from None
    except Exception as exc:
        raise click.ClickException(f"Analysis failed: {type(exc).__name__}: {exc}") from None
    click.echo(
        f"Analyzed {run.summary['total_fusions']} {run.gene_symbol} fusions with "
        f"{len(run.results)} registered algorithms"
    )
    for warning in run.warnings:
        click.echo(f"Warning: {warning}", err=True)
    for kind, path in paths.items():
        click.echo(f"{kind}: {path}")


if __name__ == "__main__":
    main()
