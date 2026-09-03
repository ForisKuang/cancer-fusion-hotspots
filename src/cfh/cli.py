"""Command-line entry point for cfh."""

from __future__ import annotations

from pathlib import Path

import click

from cfh.gene_comparison import compare_gene_runs, write_comparison_tsv
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


@main.command("compare-genes")
@click.argument(
    "run_artifacts",
    nargs=-1,
    required=True,
    type=click.Path(path_type=Path, exists=True),
)
@click.option(
    "--output",
    "output_path",
    required=True,
    type=click.Path(path_type=Path, dir_okay=False),
    help="Path for the adjusted-p-value TSV report.",
)
def compare_genes(run_artifacts: tuple[Path, ...], output_path: Path) -> None:
    """BH-adjust p-values from existing run directories or results.json files."""
    try:
        rows = compare_gene_runs(list(run_artifacts))
        if not rows:
            raise click.ClickException("No applicable p-values were found in the run artifacts.")
        write_comparison_tsv(rows, output_path)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from None
    click.echo(f"Adjusted {len(rows)} p-values across {len(run_artifacts)} run artifacts")
    click.echo(f"report: {output_path}")


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
def real_benchmark(
    gene_symbol: str,
    study_id: str,
    output_dir: Path,
    n_permutations: int,
    output_stem: str | None,
) -> None:
    """Run the prototype live cBioPortal/Genome Nexus benchmark."""
    try:
        run = run_real_benchmark(gene_symbol, study_id, n_permutations=n_permutations)
        paths = write_outputs(
            run,
            output_dir,
            output_stem=output_stem,
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
def analyze(gene_symbol: str, study_id: str, output_dir: Path, n_permutations: int) -> None:
    """Run all registered algorithms for a configured gene and live study."""
    try:
        run = run_analysis(gene_symbol, study_id, n_permutations=n_permutations)
        paths = write_outputs(
            run,
            output_dir,
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
