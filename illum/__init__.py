__version__ = "2.2.4.20220919.20291447"

try:
    # from .alternate import alternate
    # from .domain import domain
    # from .extract import extract
    # from .failed import failed
    # from .inputs import inputs
    # from .warp import warp
    from . import AngularPowerDistribution as APD
    from . import SpectralPowerDistribution as SPD
    from . import utils
    from .batches import batches
    from .init import init
except ModuleNotFoundError:
    pass


from importlib.resources import files

path = files("illum").as_posix()

del files
