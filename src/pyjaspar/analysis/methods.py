"""Shared motif comparison results, validation, and method-specific scoring.

Two methods are available, selected by name:

``"matrix_align"``
    The scoring of JASPAR's "Matrix Align" web tool: semi-global alignment of
    the frequency columns with at most one internal gap. Each aligned column
    contributes between 0 and 2; a gap costs ``open_penalty`` for its first
    column and ``ext_penalty`` for each further column.
``"pearson"``
    The mean column-wise Pearson correlation at the best ungapped offset, with
    at least ``min_overlap`` overlapping columns.

Both try the second motif as given and reverse-complemented. Scores of the two
methods are on different scales and are not comparable with each other.
"""

from __future__ import annotations

import math
import numbers
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

try:
    import numpy as np
except ImportError:
    raise ImportError(
        "The pyjaspar.analysis module requires numpy. Install with: pip install pyjaspar[analysis]"
    ) from None

from ._alignment_search import find_best_offset
from ._matrix_alignment import Coordinates, _frequencies, align_frequencies
from .similarity import pearson_correlation

if TYPE_CHECKING:
    from Bio.motifs.jaspar import Motif

METHODS = ("matrix_align", "pearson")


@dataclass(frozen=True)
class AlignmentPath:
    """The aligned columns of two motifs.

    Coordinates are zero-based column boundaries, as in
    ``Bio.Align.Alignment.coordinates``: ``coordinates[0]`` for the first motif
    (the query) and ``coordinates[1]`` for the second (the target) in the
    orientation that was aligned. Only the aligned core is stored; columns that
    hang off either end are free and not part of the path. A segment that
    advances both rows pairs columns; a segment that advances one row only is an
    internal gap.

    Attributes:
        query_length: Number of columns of the first motif.
        target_length: Number of columns of the second motif.
        coordinates: The path vertices, (query boundaries, target boundaries).
        is_reverse_complement: True if the reverse complement of the second motif
            was aligned; target coordinates then refer to that reverse complement.
            An interval ``[s, e)`` of it is ``[target_length - e, target_length - s)``
            of the motif as given.
    """

    query_length: int
    target_length: int
    coordinates: Coordinates
    is_reverse_complement: bool

    @property
    def offset(self) -> int:
        """Position in the query minus position in the target at the first aligned pair."""
        return self.coordinates[0][0] - self.coordinates[1][0]

    @property
    def overlap(self) -> int:
        """Number of paired columns."""
        query, target = self.coordinates
        return sum(
            q1 - q0
            for q0, q1, t0, t1 in zip(query, query[1:], target, target[1:], strict=False)
            if q1 > q0 and t1 > t0
        )

    @property
    def gaps(self) -> int:
        """Number of internal gap columns (0 if the alignment has no gap)."""
        query, target = self.coordinates
        return sum(
            (q1 - q0) + (t1 - t0)
            for q0, q1, t0, t1 in zip(query, query[1:], target, target[1:], strict=False)
            if (q1 > q0) != (t1 > t0)
        )

    @property
    def alignment_length(self) -> int:
        """Columns of the alignment: paired columns plus gap columns."""
        return self.overlap + self.gaps


@dataclass(frozen=True)
class AlignmentResult:
    """Two motifs compared by ``align_motifs``.

    Attributes:
        method: ``"matrix_align"`` or ``"pearson"``.
        motif1_id: Matrix ID of the first motif.
        motif2_id: Matrix ID of the second motif.
        score: For ``"matrix_align"``, the alignment score: each aligned column
            contributes between 0 and 2, and the gap, if any, subtracts its cost
            (the "Score" column of the JASPAR web tool). For ``"pearson"``, the mean
            column-wise Pearson correlation of the aligned columns, in [-1, 1].
        alignment: The aligned columns (at most one internal gap for
            ``"matrix_align"``, none for ``"pearson"``).
        parameters: The options used (``open_penalty``, ``ext_penalty`` and
            ``both_strands``, or ``min_overlap`` and ``both_strands``).
    """

    method: str
    motif1_id: str | None
    motif2_id: str | None
    score: float
    alignment: AlignmentPath
    parameters: dict[str, Any]

    @property
    def offset(self) -> int:
        return self.alignment.offset

    @property
    def is_reverse_complement(self) -> bool:
        return self.alignment.is_reverse_complement

    @property
    def overlap(self) -> int:
        return self.alignment.overlap

    @property
    def gaps(self) -> int:
        return self.alignment.gaps

    @property
    def alignment_length(self) -> int:
        return self.alignment.alignment_length


@dataclass(frozen=True)
class ProfileHit:
    """A candidate profile found by ``search_profiles``.

    Attributes:
        matrix_id: Matrix ID of the candidate.
        name: Name of the TF the candidate belongs to.
        result: The comparison of the query with this candidate.
        percent_score: For ``"matrix_align"``, the "Percent Score" column of the
            JASPAR web tool, ``100 * score / (2 * m)``, where ``m`` is the narrowest
            profile seen so far. ``m`` starts at the query's width and is lowered by
            each candidate in matrix-ID order (the order of the web tool's table),
            and never raised again. It therefore depends on the set of candidates,
            and it can exceed 100. None for ``"pearson"``.
    """

    matrix_id: str | None
    name: str | None
    result: AlignmentResult
    percent_score: float | None = None

    @property
    def method(self) -> str:
        return self.result.method

    @property
    def score(self) -> float:
        return self.result.score

    @property
    def alignment(self) -> AlignmentPath:
        return self.result.alignment

    @property
    def offset(self) -> int:
        return self.result.offset

    @property
    def is_reverse_complement(self) -> bool:
        return self.result.is_reverse_complement

    @property
    def overlap(self) -> int:
        return self.result.overlap

    @property
    def gaps(self) -> int:
        return self.result.gaps

    @property
    def alignment_length(self) -> int:
        return self.result.alignment_length

    @property
    def width(self) -> int:
        """Number of columns of the candidate."""
        return self.result.alignment.target_length


def _matrix_id(motif: Motif) -> str | None:
    """Matrix ID of a motif; motifs that are not JASPAR records have none."""
    return getattr(motif, "matrix_id", None)


def _label(motif: Motif) -> str:
    """Name of a motif for error messages."""
    return _matrix_id(motif) or getattr(motif, "name", None) or "motif"


def _check_method(method: object) -> str:
    if not isinstance(method, str) or method not in METHODS:
        raise ValueError(f"method must be 'matrix_align' or 'pearson', not {method!r}")
    return method


def _resolve_options(
    method: object,
    *,
    both_strands: bool = True,
    min_overlap: int | None = None,
    open_penalty: float | None = None,
    ext_penalty: float | None = None,
) -> dict[str, Any]:
    """Validate the method name and its options; return the effective options."""
    method = _check_method(method)
    if not isinstance(both_strands, bool):
        raise ValueError(f"both_strands must be True or False, not {both_strands!r}")
    if method == "pearson":
        if open_penalty is not None or ext_penalty is not None:
            raise ValueError("open_penalty and ext_penalty apply only to method='matrix_align'")
        overlap = 4 if min_overlap is None else min_overlap
        if isinstance(overlap, bool) or not isinstance(overlap, numbers.Integral) or overlap < 1:
            raise ValueError(f"min_overlap must be a positive integer, not {overlap!r}")
        return {"min_overlap": int(overlap), "both_strands": both_strands}
    if min_overlap is not None:
        raise ValueError("min_overlap applies only to method='pearson'")
    penalties = {
        "open_penalty": 3.0 if open_penalty is None else open_penalty,
        "ext_penalty": 0.01 if ext_penalty is None else ext_penalty,
    }
    for name, value in penalties.items():
        if (
            isinstance(value, bool)
            or not isinstance(value, numbers.Real)
            or not math.isfinite(value)
            or value < 0
        ):
            raise ValueError(f"{name} must be a finite number >= 0, not {value!r}")
    return {
        "open_penalty": float(penalties["open_penalty"]),
        "ext_penalty": float(penalties["ext_penalty"]),
        "both_strands": both_strands,
    }


def _validate_profile(motif: Motif) -> None:
    """Raise ValueError unless the motif has finite, non-negative A/C/G/T counts."""
    label = _label(motif)
    try:
        rows = [list(motif.counts[base]) for base in "ACGT"]
        counts = np.array(rows, dtype=float)
    except (AttributeError, KeyError, TypeError, ValueError):
        raise ValueError(f"Motif {label} must have A, C, G and T counts") from None
    if counts.ndim != 2 or counts.shape[1] == 0:
        raise ValueError(f"Motif {label} has no columns")
    if not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError(f"Motif {label} has counts that are negative or not finite")
    if (counts.sum(axis=0) <= 0).any():
        raise ValueError(f"Motif {label} has a column with no counts")


def _pearson_path(len1: int, len2: int, offset: int, is_rc: bool) -> AlignmentPath:
    """Path of an ungapped alignment at ``offset``, as ``pearson_correlation`` aligns it."""
    start1 = max(0, offset)
    start2 = max(0, -offset)
    overlap = min(len1 - start1, len2 - start2)
    return AlignmentPath(
        query_length=len1,
        target_length=len2,
        coordinates=((start1, start1 + overlap), (start2, start2 + overlap)),
        is_reverse_complement=is_rc,
    )


def _compare_pearson(a: Motif, b: Motif, options: dict[str, Any]) -> AlignmentResult:
    min_overlap = options["min_overlap"]
    for motif in (a, b):
        if motif.length < min_overlap:
            raise ValueError(
                f"Motif {_label(motif)} has {motif.length} columns, "
                f"fewer than min_overlap={min_overlap}"
            )
    score, offset, is_rc = find_best_offset(
        a, b, pearson_correlation, min_overlap, options["both_strands"]
    )
    return AlignmentResult(
        method="pearson",
        motif1_id=_matrix_id(a),
        motif2_id=_matrix_id(b),
        score=float(score),
        alignment=_pearson_path(a.length, b.length, offset, is_rc),
        parameters=dict(options),
    )


def _compare_matrix_align(a: Motif, b: Motif, options: dict[str, Any]) -> AlignmentResult:
    fa = _frequencies(a)
    fb = _frequencies(b)
    score, is_rc, coordinates = align_frequencies(
        fa, fb, options["open_penalty"], options["ext_penalty"], options["both_strands"]
    )
    return AlignmentResult(
        method="matrix_align",
        motif1_id=_matrix_id(a),
        motif2_id=_matrix_id(b),
        score=score,
        alignment=AlignmentPath(
            query_length=len(fa),
            target_length=len(fb),
            coordinates=coordinates,
            is_reverse_complement=is_rc,
        ),
        parameters=dict(options),
    )


def _compare(a: Motif, b: Motif, method: str, options: dict[str, Any]) -> AlignmentResult:
    """Compare two validated motifs with an already resolved method and options."""
    if method == "pearson":
        return _compare_pearson(a, b, options)
    return _compare_matrix_align(a, b, options)
