"""Database interface for JASPAR motif queries."""

from __future__ import annotations

import sqlite3
import warnings
from pathlib import Path

from Bio import BiopythonWarning
from Bio.motifs import jaspar
from Bio.motifs.jaspar import Motif, Record
from Bio.motifs.matrix import GenericPositionMatrix

from ._constants import (
    JASPAR_DFLT_COLLECTION,
    JASPAR_DL_RELEASES,
    JASPAR_LATEST_RELEASE,
    jaspar_releases,
)
from ._exceptions import DatabaseError, DLNotAvailableError, ReleaseNotFoundError
from ._queries import build_motif_query
from .utils import get_jaspardb_path


class JasparDB:
    """Interface to a JASPAR SQLite database.

    This is adapted from the biopython JASPAR5 MYSQL DB.

    Can be used as a context manager::

        with JasparDB(release='JASPAR2026') as jdb:
            motif = jdb.fetch_motif_by_id('MA0001.1')

    Arguments:
        release: JASPAR release name (e.g. 'JASPAR2026', 'JASPAR2020').
            Defaults to the latest available release.
        sqlite_db_path: Path to a custom JASPAR SQLite file.
            If provided, the ``release`` argument is ignored.
    """

    def __init__(
        self,
        release: str = JASPAR_LATEST_RELEASE,
        sqlite_db_path: str | Path | None = None,
    ) -> None:
        self.release = release
        self._conn: sqlite3.Connection

        if sqlite_db_path is not None:
            path = Path(sqlite_db_path)
            if not path.exists():
                raise FileNotFoundError(f"Database file not found: {path}")
            try:
                self._conn = sqlite3.connect(str(path))
            except sqlite3.Error as e:
                raise DatabaseError(f"Failed to connect to {path}: {e}") from e
            self.release = str(path)
        else:
            if release not in jaspar_releases:
                available = ", ".join(jaspar_releases)
                raise ReleaseNotFoundError(
                    f"{release} is not available. Available releases are: {available}"
                )
            db_filename = jaspar_releases[release]
            db_path = get_jaspardb_path(db_filename)
            try:
                self._conn = sqlite3.connect(str(db_path))
            except sqlite3.Error as e:
                raise DatabaseError(f"Failed to connect to database: {e}") from e

    def __str__(self) -> str:
        """Return a string representation of the JASPAR DB connection."""
        return f"JASPAR release:{self.release}:{self._conn}"

    def close(self) -> None:
        """Close the database connection."""
        self._conn.close()

    def __enter__(self) -> JasparDB:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    @property
    def conn(self) -> sqlite3.Connection:
        """Access the underlying SQLite connection (for backward compatibility)."""
        return self._conn

    @property
    def dl(self) -> object:
        """Access the DL collection sub-client.

        Only available for JASPAR releases that include DL tables
        (currently JASPAR2026+).

        Returns:
            A :class:`pyjaspar.dl.DLClient` instance.

        Raises:
            DLNotAvailableError: If the current release does not have DL tables.
        """
        if not hasattr(self, "_dl_client"):
            if self.release not in JASPAR_DL_RELEASES and not self._has_dl_tables():
                raise DLNotAvailableError(
                    f"The DL collection is not available in release '{self.release}'. "
                    f"DL is available in: {', '.join(sorted(JASPAR_DL_RELEASES))}"
                )
            from .dl import DLClient

            self._dl_client = DLClient(self._conn)
        return self._dl_client

    def _has_dl_tables(self) -> bool:
        """Check if the database has DL tables."""
        cur = self._conn.cursor()
        cur.execute(
            "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='DL_PROFILE_SUMMARY'"
        )
        return cur.fetchone()[0] > 0

    def get_releases(self) -> list[str]:
        """Return available JASPAR releases.

        Returns:
            A list of JASPAR release names.
        """
        return list(jaspar_releases.keys())

    def fetch_motif_by_id(self, id: str) -> Motif | None:
        """Fetch a single JASPAR motif by its matrix ID.

        Arguments:
            id: JASPAR matrix ID. May be a fully specified ID including
                the version number (e.g. 'MA0049.2') or just the base ID
                (e.g. 'MA0049'). If only a base ID is provided, the latest
                version is returned.

        Returns:
            A Bio.motifs.jaspar.Motif object, or None if not found.
        """
        (base_id, version) = jaspar.split_jaspar_id(id)
        if not version:
            version = self._fetch_latest_version(base_id)

        int_id = None
        if version:
            int_id = self._fetch_internal_id(base_id, version)

        motif = None
        if int_id:
            motif = self._fetch_motif_by_internal_id(int_id)

        return motif

    def fetch_motifs_by_name(self, name: str | list[str]) -> Record:
        """Fetch JASPAR motifs by TF name(s).

        Arguments:
            name: A single name or list of names.

        Returns:
            A Bio.motifs.jaspar.Record (list) of motifs.
        """
        return self.fetch_motifs(collection=None, tf_name=name)

    def fetch_motifs(
        self,
        collection: str | list[str] | None = JASPAR_DFLT_COLLECTION,
        tf_name: str | list[str] | None = None,
        tf_class: str | list[str] | None = None,
        tf_family: str | list[str] | None = None,
        matrix_id: str | list[str] | None = None,
        tax_group: str | list[str] | None = None,
        species: str | int | list[str | int] | None = None,
        pazar_id: str | list[str] | None = None,
        data_type: str | list[str] | None = None,
        medline: str | list[str] | None = None,
        min_ic: float = 0,
        min_length: int = 0,
        min_sites: int = 0,
        all: bool = False,
        all_versions: bool = False,
    ) -> Record:
        """Fetch a Record (list) of motifs using selection criteria.

        All selection criteria arguments may be specified as a single value
        or a list of values. Motifs must meet ALL the specified criteria.

        Arguments:
            all: If True, return every motif (ignores other criteria).
            matrix_id: JASPAR matrix ID(s). Takes precedence over other criteria
                except ``all``.
            collection: JASPAR collection(s). Defaults to 'CORE'.
                Set to None to search across all collections.
            tf_name: TF name(s).
            tf_class: TF class(es).
            tf_family: TF family(ies).
            tax_group: Taxonomic supergroup(s).
            species: Species taxonomy ID(s).
            data_type: Data type(s) used to compile the matrix.
            pazar_id: PAZAR TF ID(s).
            medline: PubMed ID(s).
            min_ic: Minimum information content (specificity).
            min_length: Minimum motif length.
            min_sites: Minimum number of binding sites.
            all_versions: If True, return all versions of matching motifs.

        Returns:
            A Bio.motifs.jaspar.Record (list) of motifs.
        """
        int_ids = self._fetch_internal_id_list(
            collection=collection,
            tf_name=tf_name,
            tf_class=tf_class,
            tf_family=tf_family,
            matrix_id=matrix_id,
            tax_group=tax_group,
            species=species,
            pazar_id=pazar_id,
            data_type=data_type,
            medline=medline,
            all=all,
            all_versions=all_versions,
        )

        record = jaspar.Record()

        for int_id in int_ids:
            motif = self._fetch_motif_by_internal_id(int_id)
            if motif is None:
                continue

            if min_ic and motif.pssm.mean() < min_ic:
                continue

            if min_length and motif.length < min_length:
                continue

            if min_sites:
                num_sites = sum(motif.counts[nt][0] for nt in "ACGT")
                if num_sites < min_sites:
                    continue

            record.append(motif)

        return record

    # --- Private methods ---

    def _fetch_latest_version(self, base_id: str) -> int | None:
        """Get the latest version number for the given base_id."""
        cur = self._conn.cursor()
        cur.execute(
            "SELECT VERSION FROM MATRIX WHERE BASE_ID = ? ORDER BY VERSION DESC LIMIT 1",
            (base_id,),
        )

        row = cur.fetchone()
        if row:
            return row[0]

        warnings.warn(
            f"Failed to fetch latest version number for JASPAR motif"
            f" with base ID '{base_id}'. No JASPAR motif with this"
            f" base ID appears to exist in the database.",
            BiopythonWarning,
            stacklevel=2,
        )
        return None

    def _fetch_internal_id(self, base_id: str, version: int | str) -> int | None:
        """Fetch the internal id for a base id + version."""
        cur = self._conn.cursor()
        cur.execute(
            "SELECT ID FROM MATRIX WHERE BASE_ID = ? AND VERSION = ? COLLATE NOCASE",
            (base_id, version),
        )

        row = cur.fetchone()
        if row:
            return row[0]

        warnings.warn(
            f"Failed to fetch internal database ID for JASPAR motif"
            f" with matrix ID '{base_id}.{version}'. No JASPAR motif"
            f" with this matrix ID appears to exist.",
            BiopythonWarning,
            stacklevel=2,
        )
        return None

    def _fetch_motif_by_internal_id(self, int_id: int) -> Motif | None:
        """Fetch basic motif information by internal database ID."""
        cur = self._conn.cursor()
        cur.execute(
            "SELECT BASE_ID, VERSION, COLLECTION, NAME FROM MATRIX WHERE ID = ? COLLATE NOCASE",
            (int_id,),
        )

        row = cur.fetchone()
        if not row:
            warnings.warn(
                f"Could not fetch JASPAR motif with internal ID = {int_id}",
                BiopythonWarning,
                stacklevel=2,
            )
            return None

        base_id = row[0]
        version = row[1]
        collection = row[2]
        name = row[3]

        matrix_id = f"{base_id}.{version}"

        counts = self._fetch_counts_matrix(int_id)

        motif = jaspar.Motif(matrix_id, name, collection=collection, counts=counts)

        # Fetch species
        cur.execute("SELECT TAX_ID FROM MATRIX_SPECIES WHERE ID = ?", (int_id,))
        motif.species = [row[0] for row in cur.fetchall()]

        # Fetch protein accession numbers
        cur.execute(
            "SELECT ACC FROM MATRIX_PROTEIN WHERE ID = ? COLLATE NOCASE",
            (int_id,),
        )
        motif.acc = [row[0] for row in cur.fetchall()]

        # Fetch annotations
        cur.execute("SELECT TAG, VAL FROM MATRIX_ANNOTATION WHERE ID = ?", (int_id,))

        tf_family: list[str] = []
        tf_class: list[str] = []
        for row in cur.fetchall():
            attr = row[0]
            val = row[1]
            if attr == "class":
                tf_class.append(val)
            elif attr == "family":
                tf_family.append(val)
            elif attr == "tax_group":
                motif.tax_group = val
            elif attr == "type":
                motif.data_type = val
            elif attr == "pazar_tf_id":
                motif.pazar_id = val
            elif attr == "medline":
                motif.medline = val
            elif attr == "comment":
                motif.comment = val

        motif.tf_family = tf_family
        motif.tf_class = tf_class

        return motif

    def _fetch_counts_matrix(self, int_id: int) -> GenericPositionMatrix:
        """Fetch the counts matrix from the JASPAR DB by internal ID."""
        counts: dict[str, list[float]] = {}
        cur = self._conn.cursor()

        for base in "ACGT":
            cur.execute(
                "SELECT val FROM MATRIX_DATA WHERE ID = ? AND row = ? ORDER BY col",
                (int_id, base),
            )
            counts[base] = [float(row[0]) for row in cur.fetchall()]

        return GenericPositionMatrix("ACGT", counts)

    def _fetch_internal_id_list(
        self,
        collection: str | list[str] | None = JASPAR_DFLT_COLLECTION,
        tf_name: str | list[str] | None = None,
        tf_class: str | list[str] | None = None,
        tf_family: str | list[str] | None = None,
        matrix_id: str | list[str] | None = None,
        tax_group: str | list[str] | None = None,
        species: str | int | list[str | int] | None = None,
        pazar_id: str | list[str] | None = None,
        data_type: str | list[str] | None = None,
        medline: str | list[str] | None = None,
        all: bool = False,
        all_versions: bool = False,
    ) -> list[int]:
        """Fetch list of internal JASPAR motif IDs based on selection criteria."""
        int_ids: list[int] = []
        cur = self._conn.cursor()

        # Special case 1: fetch ALL motifs
        if all:
            cur.execute("SELECT ID FROM MATRIX")
            return [row[0] for row in cur.fetchall()]

        # Special case 2: fetch by specific matrix IDs
        if matrix_id:
            if not isinstance(matrix_id, list):
                matrix_id = [matrix_id]

            if all_versions:
                for mid in matrix_id:
                    (base_id, _version) = jaspar.split_jaspar_id(mid)
                    cur.execute(
                        "SELECT ID FROM MATRIX WHERE BASE_ID = ? COLLATE NOCASE",
                        (base_id,),
                    )
                    int_ids.extend(row[0] for row in cur.fetchall())
            else:
                for mid in matrix_id:
                    (base_id, version) = jaspar.split_jaspar_id(mid)
                    if not version:
                        version = self._fetch_latest_version(base_id)
                    if version:
                        int_id = self._fetch_internal_id(base_id, version)
                        if int_id:
                            int_ids.append(int_id)

            return int_ids

        # General case: build query from criteria
        sql, params = build_motif_query(
            collection=collection,
            tf_name=tf_name,
            tf_class=tf_class,
            tf_family=tf_family,
            species=species,
            tax_group=tax_group,
            pazar_id=pazar_id,
            medline=medline,
            data_type=data_type,
        )

        cur.execute(sql, params)
        rows = cur.fetchall()

        for row in rows:
            id = row[0]
            if all_versions:
                int_ids.append(id)
            else:
                if self._is_latest_version(id):
                    int_ids.append(id)

        if len(int_ids) < 1:
            warnings.warn(
                "Zero motifs returned with current select criteria",
                BiopythonWarning,
                stacklevel=2,
            )

        return int_ids

    def _is_latest_version(self, int_id: int) -> bool:
        """Check if the internal ID represents the latest matrix version."""
        cur = self._conn.cursor()
        cur.execute(
            "SELECT COUNT(*) FROM MATRIX WHERE "
            "BASE_ID = (SELECT BASE_ID FROM MATRIX WHERE ID = ?) "
            "AND VERSION > (SELECT VERSION FROM MATRIX WHERE ID = ?) "
            "COLLATE NOCASE",
            (int_id, int_id),
        )
        row = cur.fetchone()
        return row[0] == 0
