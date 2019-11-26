import json
import logging
import os
import signal

import click
from generate_symlinks import  GenerateSymlink
__author__ = 'vinay kaikala'

def setup_logging(loglevel="INFO"):
    """
    Set up logging
    :param loglevel: loglevel to use, one of ERROR, WARNING, DEBUG, INFO (default INFO)
    :return:
    """
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level,
                        format='%(levelname)s %(asctime)s %(funcName)s - %(message)s')
    logging.info("Started log with loglevel %(loglevel)s" % {"loglevel": loglevel})

@click.command()
@click.option('--outdir',  default='/tmp/autoseq-test', help='output directory', type=click.Path())
def cli(outdir):
    setup_logging()
    create_symlinks = GenerateSymlink(outdir)
    create_symlinks.generateIGVsymlink()
    create_symlinks.create_igv_session_file()
