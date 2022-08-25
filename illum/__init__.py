__version__ = "2.2.2.20220825.13465213"

try:
    # from .alternate import alternate
    # from .domain import domain
    # from .extract import extract
    # from .failed import failed
    # from .inputs import inputs
    # from .warp import warp
    from . import utils
    from .batches import batches
    from .init import init
except ModuleNotFoundError:
    pass


try:
    import os as _os

    path = _os.path.dirname(
        [
            path
            for path in _os.environ["PATH"].split(":")
            if path.endswith("illumina/bin")
        ][0]
    )
except IndexError:
    raise ValueError(
        "The 'illumina/bin' folder is not in the PATH environment variable."
    )
