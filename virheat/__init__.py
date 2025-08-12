from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("virheat")
except PackageNotFoundError:
    __version__ = "unknown"
