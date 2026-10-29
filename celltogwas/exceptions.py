"""Package-specific exceptions."""


class CellToGWASError(Exception):
    """Base exception for celltogwas."""


class AnnotationError(CellToGWASError, ValueError):
    """Raised when a gene or regulatory annotation is malformed."""


class InputFormatError(CellToGWASError, ValueError):
    """Raised when an input table does not contain required columns."""


class LinkageError(CellToGWASError, ValueError):
    """Raised when cell-type domains cannot be constructed."""


class EnrichmentError(CellToGWASError, ValueError):
    """Raised when enrichment input or configuration is invalid."""
