"""Read Teledyne RDI PD0 binary files. Reference: docs/legacy/loadrdi.m."""

from pathlib import Path


def load_rdi(path: Path) -> dict:
    """Load one RDI PD0 binary (.000) file.

    Returns a dict matching the MATLAB ``d`` struct from loadrdi.m:
    keys will include ``vel`` (velocity, ensembles × bins × beams),
    ``heading``, ``pitch``, ``roll``, ``time``, ``pressure``.
    """
    raise NotImplementedError
