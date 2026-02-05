"""Custom exceptions for pyjaspar."""


class PyJasparError(Exception):
    """Base exception for all pyjaspar errors."""


class DatabaseError(PyJasparError):
    """Raised when a database connection or query fails."""


class ReleaseNotFoundError(PyJasparError):
    """Raised when the requested JASPAR release is not available."""


class MotifNotFoundError(PyJasparError):
    """Raised when no motif matches the query criteria."""


class DLNotAvailableError(PyJasparError):
    """Raised when DL collection is accessed on a release that does not support it."""
