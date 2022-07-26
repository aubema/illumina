__version__ = "2.2.2.20220726.19125213"

try:
    import os as _os

    from . import pytools
    from .batches import batches
    from .init import init

    # from .alternate import alternate
    # from .domain import domain
    # from .extract import extract
    # from .failed import failed
    # from .inputs import inputs
    # from .warp import warp

    path = _os.path.dirname(
        [
            path
            for path in _os.environ["PATH"].split(":")
            if path.endswith("illumina/bin")
        ][0]
    )
except ModuleNotFoundError:
    pass
