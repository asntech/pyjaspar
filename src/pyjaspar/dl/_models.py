"""Data models for the JASPAR DL collection."""

from __future__ import annotations

from dataclasses import dataclass, field

try:
    import numpy as np
    from numpy.typing import NDArray
except ImportError:
    raise ImportError(
        "The pyjaspar.dl module requires numpy. Install with: pip install pyjaspar[dl]"
    ) from None


@dataclass(frozen=True)
class Matrix:
    """A position matrix (PFM or CWM) from the DL collection.

    Attributes:
        kind: Matrix type, either ``"PFM"`` or ``"CWM"``.
        values: 2D array of shape ``(4, L)`` where L is the motif length.
            Rows are ordered A, C, G, T. PFM values are non-negative counts;
            CWM values can be negative floats (contribution weights).
    """

    kind: str
    values: NDArray[np.float64] = field(default_factory=lambda: np.empty((4, 0)))

    @property
    def length(self) -> int:
        """Return the motif length (number of positions)."""
        return int(self.values.shape[1])

    def to_dict(self) -> dict[str, list[float]]:
        """Convert to a dict mapping nucleotide -> list of values."""
        alphabet = ("A", "C", "G", "T")
        return {nt: self.values[i].tolist() for i, nt in enumerate(alphabet)}


@dataclass(frozen=True)
class DLMotifPattern:
    """A motif pattern (cluster) from the DL collection.

    Corresponds to a cluster (MO*) in the JASPAR DL tables. Each cluster
    references a representative motif in DL_MOTIF with associated matrices.

    Attributes:
        motif_id: Composite ID ``"MOXXXXXX.Y"`` (e.g. ``"MO000001.1"``).
        tf_name: Transcription factor name.
        motif_type: One of ``"pmp"`` (primary), ``"amp"`` (alternative), or ``"other"``.
        name: Pattern name from the underlying DL_MOTIF record.
        num_seqlets: Number of seqlets contributing to this pattern.
        num_hits: Number of hits for this pattern.
        matrices: Dict mapping matrix type (``"PFM"``, ``"CWM"``) to Matrix objects.
        jaspar_match: Best matching classic JASPAR motif ID, if available.
    """

    motif_id: str
    tf_name: str
    motif_type: str
    name: str
    num_seqlets: int
    num_hits: int
    matrices: dict[str, Matrix] = field(default_factory=dict)
    jaspar_match: str | None = None

    @property
    def pfm(self) -> Matrix | None:
        """Return the PFM matrix if available."""
        return self.matrices.get("PFM")

    @property
    def cwm(self) -> Matrix | None:
        """Return the CWM matrix if available."""
        return self.matrices.get("CWM")

    def to_biopython(self) -> object:
        """Convert the PFM to a ``Bio.motifs.jaspar.Motif`` object.

        Only works for motif patterns that have a PFM matrix.

        Returns:
            A Bio.motifs.jaspar.Motif with the PFM counts.

        Raises:
            ValueError: If no PFM matrix is available.
        """
        from Bio.motifs import jaspar
        from Bio.motifs.matrix import GenericPositionMatrix

        pfm = self.pfm
        if pfm is None:
            raise ValueError(
                f"Motif pattern {self.motif_id} has no PFM matrix; "
                "cannot convert to BioPython Motif."
            )

        counts = pfm.to_dict()
        return jaspar.Motif(
            matrix_id=self.motif_id,
            name=self.tf_name,
            counts=GenericPositionMatrix("ACGT", counts),
        )


@dataclass(frozen=True)
class DLProfileSummary:
    """Lightweight summary of a DL profile (no matrices loaded).

    Attributes:
        profile_id: Composite ID ``"DLxxxx.Y"`` (e.g. ``"DL0001.1"``).
        base_id: Base ID (e.g. ``"DL0001"``).
        version: Version number.
        tf_name: Transcription factor name.
        tax_group: Taxonomic group (e.g. ``"vertebrates"``).
    """

    profile_id: str
    base_id: str
    version: int
    tf_name: str
    tax_group: str


@dataclass(frozen=True)
class DLProfile:
    """Full DL profile with motif patterns and linked data.

    Attributes:
        profile_id: Composite ID ``"DLxxxx.Y"``.
        base_id: Base ID (e.g. ``"DL0001"``).
        version: Version number.
        tf_name: Transcription factor name.
        tax_group: Taxonomic group.
        primary_motif: The primary motif pattern (PMP) for this profile.
        alt_motifs: Alternative motif patterns (AMPs).
        jaspar_matches: Dict mapping cluster MO* ID to classic JASPAR match ID.
        linked_model_ids: List of model composite IDs (BP*) linked via clusters.
    """

    profile_id: str
    base_id: str
    version: int
    tf_name: str
    tax_group: str
    primary_motif: DLMotifPattern | None
    alt_motifs: list[DLMotifPattern] = field(default_factory=list)
    jaspar_matches: dict[str, str] = field(default_factory=dict)
    linked_model_ids: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class DLModelSummary:
    """Lightweight summary of a DL model (no matrices loaded).

    Attributes:
        model_id: Composite ID ``"BPXXXXXX.Y"`` (e.g. ``"BP000001.1"``).
        base_id: Base ID (e.g. ``"BP000001"``).
        version: Version number.
        tf_name: Transcription factor name.
        cell_line: Cell line used.
        tax_id: Taxonomy ID of the species.
        data_type: Type of data (e.g. ``"ChIP-seq"``).
        model_name: Name of the DL model (e.g. ``"BPNet"``).
        source: Data source (e.g. ``"ENCODE"``).
    """

    model_id: str
    base_id: str
    version: int
    tf_name: str
    cell_line: str
    tax_id: int
    data_type: str
    model_name: str
    source: str


@dataclass(frozen=True)
class DLModel:
    """Full DL model with its motif patterns.

    Attributes:
        model_id: Composite ID ``"BPXXXXXX.Y"``.
        base_id: Base ID (e.g. ``"BP000001"``).
        version: Version number.
        tf_name: Transcription factor name.
        cell_line: Cell line used.
        tax_id: Taxonomy ID.
        data_type: Data type (e.g. ``"ChIP-seq"``).
        model_name: Model architecture name (e.g. ``"BPNet"``).
        source: Data source (e.g. ``"ENCODE"``).
        source_id: Source experiment ID (e.g. ``"ENCSR000ATM"``).
        source_url: URL to the source experiment.
        primary_motif: The primary motif pattern (PMP) for this model, if any.
        motifs: All motif patterns discovered by this model.
    """

    model_id: str
    base_id: str
    version: int
    tf_name: str
    cell_line: str
    tax_id: int
    data_type: str
    model_name: str
    source: str
    source_id: str
    source_url: str
    primary_motif: DLMotifPattern | None
    motifs: list[DLMotifPattern] = field(default_factory=list)
