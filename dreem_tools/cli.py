from dreem_tools import run
from dreem_tools import logger
import logging
import click


@click.group()
def cli():
    pass


@cli.command(help="download a run using bs commandline tool")
@click.argument("run_name")
@click.option(
    "-d",
    "--dir",
    default=None,
    help="root directory of where data should be downloaded will default to "
    "$BASESPACE",
)
def download(run_name, dir):
    logger.setup_applevel_logger()
    return run.download(run_name, dir)


@cli.command()
@click.argument("csv")
@click.option("--debug", is_flag=True)
def demultiplex(csv, debug):
    """
    demultiplexes paired fastq files given 3' end barcodes
    """
    logger.setup_applevel_logger(file_name="demultiplex.log")
    return run.demultiplex(csv, debug)


@cli.command()
@click.argument("csv")
@click.argument("data_dir")
@click.argument("seq_path")
@click.option("--hide-dreem-output", is_flag=True)
def runmulti(csv, data_dir, hide_dreem_output, seq_path):
    log = logger.setup_applevel_logger(file_name="run_multi.log")
    log.info("creating processed/ all dreem runs will go here")
    log.info("creating analysis/ all finalized analysis will go here")
    return run.runmulti(csv, data_dir, seq_path, hide_dreem_output)


@cli.command()
@click.argument("json_file")
@click.argument("yml_file")
def parsedata(json_file, yml_file):
    logger.setup_applevel_logger()
    return run.parsedata(json_file, yml_file)


@cli.command()
@click.argument("pickle_file")
def replot(pickle_file):
    return run.replot(pickle_file)


if __name__ == "__main__":
    cli()
