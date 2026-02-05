"""DL sub-client for accessing the JASPAR DL collection."""

from __future__ import annotations

import sqlite3

from pyjaspar._exceptions import MotifNotFoundError

from ._models import (
    DLModel,
    DLModelSummary,
    DLMotifPattern,
    DLProfile,
    DLProfileSummary,
    Matrix,
)
from ._parser import parse_matrix_json
from ._queries import (
    FETCH_CLUSTER,
    FETCH_CLUSTER_JASPAR_MATCH,
    FETCH_CLUSTER_LATEST_VERSION,
    FETCH_CLUSTER_MODELS,
    FETCH_MODEL_BY_ID,
    FETCH_MODEL_LATEST_VERSION,
    FETCH_MODEL_MOTIFS,
    FETCH_MOTIF_BY_ID,
    FETCH_MOTIF_MATRICES,
    FETCH_PROFILE_BY_ID,
    FETCH_PROFILE_CLUSTERS,
    FETCH_PROFILE_LATEST_VERSION,
    build_model_search_query,
    build_profile_search_query,
)


def _split_dl_id(composite_id: str) -> tuple[str, int | None]:
    """Split a composite DL ID like ``"DL0001.1"`` into ``(base_id, version)``.

    Returns ``(base_id, None)`` if no version is present.
    """
    if "." in composite_id:
        parts = composite_id.rsplit(".", 1)
        return parts[0], int(parts[1])
    return composite_id, None


class DLClient:
    """Sub-client for DL collection queries.

    Not instantiated directly; accessed via :attr:`JasparDB.dl`.

    Args:
        conn: An open sqlite3.Connection to a JASPAR database with DL tables.
    """

    def __init__(self, conn: sqlite3.Connection) -> None:
        self._conn = conn

    # --- Public API ---

    def fetch_profile(self, profile_id: str) -> DLProfile:
        """Fetch a full DL profile by its ID.

        Args:
            profile_id: Profile ID, e.g. ``"DL0001.1"`` or ``"DL0001"``
                (latest version).

        Returns:
            A :class:`DLProfile` with primary motif, alt motifs, cluster info,
            and JASPAR matches.

        Raises:
            MotifNotFoundError: If the profile is not found.
        """
        base_id, version = _split_dl_id(profile_id)
        version = self._resolve_version(base_id, version, FETCH_PROFILE_LATEST_VERSION)

        cur = self._conn.cursor()
        cur.execute(FETCH_PROFILE_BY_ID, (base_id, version))
        row = cur.fetchone()
        if row is None:
            raise MotifNotFoundError(f"DL profile not found: {base_id}.{version}")

        _db_id, _base, _ver, primary_motif_id, tf_name, tax_group = row
        composite_id = f"{base_id}.{version}"

        # Fetch primary motif
        primary_motif = self._fetch_motif_with_matrices(primary_motif_id)

        # Fetch linked clusters
        cur.execute(FETCH_PROFILE_CLUSTERS, (base_id, version))
        cluster_rows = cur.fetchall()

        alt_motifs: list[DLMotifPattern] = []
        jaspar_matches: dict[str, str] = {}
        linked_model_ids: list[str] = []

        for cluster_base_id, cluster_version in cluster_rows:
            cluster_composite = f"{cluster_base_id}.{cluster_version}"

            # Fetch cluster details
            cur.execute(FETCH_CLUSTER, (cluster_base_id, cluster_version))
            cluster_row = cur.fetchone()
            if cluster_row is None:
                continue

            _cid, _cbase, _cver, cluster_motif_id, cluster_tf_name = cluster_row

            # Fetch JASPAR match for this cluster
            cur.execute(FETCH_CLUSTER_JASPAR_MATCH, (cluster_base_id, cluster_version))
            match_row = cur.fetchone()
            if match_row:
                jaspar_matches[cluster_composite] = match_row[0]

            # Fetch linked models for this cluster
            cur.execute(FETCH_CLUSTER_MODELS, (cluster_base_id, cluster_version))
            for model_base, model_ver in cur.fetchall():
                linked_model_ids.append(f"{model_base}.{model_ver}")

            # Build motif pattern for this cluster if it's not the primary
            if cluster_motif_id != primary_motif_id:
                motif_data = self._fetch_motif_raw(cluster_motif_id)
                if motif_data is not None:
                    matrices = self._fetch_matrices(cluster_motif_id)
                    alt_motifs.append(
                        DLMotifPattern(
                            motif_id=cluster_composite,
                            tf_name=cluster_tf_name,
                            motif_type=motif_data["motif_type"],
                            name=motif_data["name"],
                            num_seqlets=motif_data["num_seqlets"],
                            num_hits=motif_data["num_hits"],
                            matrices=matrices,
                            jaspar_match=jaspar_matches.get(cluster_composite),
                        )
                    )

        # Build the primary motif pattern with its cluster ID
        primary_motif_pattern: DLMotifPattern | None = None
        if primary_motif is not None:
            # Find which cluster corresponds to the primary motif
            primary_cluster_id: str | None = None
            for cluster_base_id, cluster_version in cluster_rows:
                cur.execute(FETCH_CLUSTER, (cluster_base_id, cluster_version))
                cr = cur.fetchone()
                if cr and cr[3] == primary_motif_id:
                    primary_cluster_id = f"{cluster_base_id}.{cluster_version}"
                    break

            primary_motif_pattern = DLMotifPattern(
                motif_id=primary_cluster_id or composite_id,
                tf_name=tf_name,
                motif_type=primary_motif["motif_type"],
                name=primary_motif["name"],
                num_seqlets=primary_motif["num_seqlets"],
                num_hits=primary_motif["num_hits"],
                matrices=primary_motif["matrices"],
                jaspar_match=jaspar_matches.get(primary_cluster_id or ""),
            )

        # Deduplicate linked model IDs
        seen: set[str] = set()
        unique_model_ids: list[str] = []
        for mid in linked_model_ids:
            if mid not in seen:
                seen.add(mid)
                unique_model_ids.append(mid)

        return DLProfile(
            profile_id=composite_id,
            base_id=base_id,
            version=version,
            tf_name=tf_name,
            tax_group=tax_group,
            primary_motif=primary_motif_pattern,
            alt_motifs=alt_motifs,
            jaspar_matches=jaspar_matches,
            linked_model_ids=unique_model_ids,
        )

    def fetch_model(self, model_id: str) -> DLModel:
        """Fetch a full DL model by its ID.

        Args:
            model_id: Model ID, e.g. ``"BP000001.1"`` or ``"BP000001"``
                (latest version).

        Returns:
            A :class:`DLModel` with its motif patterns.

        Raises:
            MotifNotFoundError: If the model is not found.
        """
        base_id, version = _split_dl_id(model_id)
        version = self._resolve_version(base_id, version, FETCH_MODEL_LATEST_VERSION)

        cur = self._conn.cursor()
        cur.execute(FETCH_MODEL_BY_ID, (base_id, version))
        row = cur.fetchone()
        if row is None:
            raise MotifNotFoundError(f"DL model not found: {base_id}.{version}")

        (
            _db_id,
            _base,
            _ver,
            primary_motif_id,
            tf_name,
            cell_line,
            tax_id,
            data_type,
            model_name,
            source,
            source_id,
            source_url,
        ) = row
        composite_id = f"{base_id}.{version}"

        # Fetch all motifs for this model
        cur.execute(FETCH_MODEL_MOTIFS, (base_id, version))
        motif_rows = cur.fetchall()

        motifs: list[DLMotifPattern] = []
        primary_motif: DLMotifPattern | None = None

        for mrow in motif_rows:
            (
                motif_db_id,
                mname,
                msource,
                msource_id,
                msource_version,
                num_seqlets,
                num_hits,
                mtype,
            ) = mrow
            matrices = self._fetch_matrices(motif_db_id)
            pattern = DLMotifPattern(
                motif_id=f"{msource_id}.{msource_version}:{motif_db_id}",
                tf_name=tf_name,
                motif_type=mtype,
                name=mname,
                num_seqlets=num_seqlets,
                num_hits=num_hits,
                matrices=matrices,
            )
            motifs.append(pattern)
            if motif_db_id == primary_motif_id:
                primary_motif = pattern

        return DLModel(
            model_id=composite_id,
            base_id=base_id,
            version=version,
            tf_name=tf_name,
            cell_line=cell_line,
            tax_id=tax_id,
            data_type=data_type,
            model_name=model_name,
            source=source,
            source_id=source_id,
            source_url=source_url,
            primary_motif=primary_motif,
            motifs=motifs,
        )

    def fetch_motif_pattern(self, motif_id: str) -> DLMotifPattern:
        """Fetch a single DL motif pattern (cluster) by its MO* ID.

        Args:
            motif_id: Motif pattern ID, e.g. ``"MO000001.1"`` or ``"MO000001"``
                (latest version).

        Returns:
            A :class:`DLMotifPattern` with its matrices.

        Raises:
            MotifNotFoundError: If the motif pattern is not found.
        """
        base_id, version = _split_dl_id(motif_id)
        version = self._resolve_version(base_id, version, FETCH_CLUSTER_LATEST_VERSION)

        cur = self._conn.cursor()
        cur.execute(FETCH_CLUSTER, (base_id, version))
        row = cur.fetchone()
        if row is None:
            raise MotifNotFoundError(f"DL motif pattern not found: {base_id}.{version}")

        _db_id, _cbase, _cver, cluster_motif_id, tf_name = row
        composite_id = f"{base_id}.{version}"

        # Fetch the underlying DL_MOTIF
        motif_data = self._fetch_motif_raw(cluster_motif_id)
        if motif_data is None:
            raise MotifNotFoundError(
                f"DL_MOTIF record not found for cluster {composite_id} "
                f"(internal motif ID={cluster_motif_id})"
            )

        matrices = self._fetch_matrices(cluster_motif_id)

        # Fetch JASPAR match
        cur.execute(FETCH_CLUSTER_JASPAR_MATCH, (base_id, version))
        match_row = cur.fetchone()
        jaspar_match = match_row[0] if match_row else None

        return DLMotifPattern(
            motif_id=composite_id,
            tf_name=tf_name,
            motif_type=motif_data["motif_type"],
            name=motif_data["name"],
            num_seqlets=motif_data["num_seqlets"],
            num_hits=motif_data["num_hits"],
            matrices=matrices,
            jaspar_match=jaspar_match,
        )

    def search_profiles(
        self,
        *,
        tf_name: str | list[str] | None = None,
        tax_group: str | list[str] | None = None,
    ) -> list[DLProfileSummary]:
        """Search DL profiles by filter criteria.

        Args:
            tf_name: TF name(s) to filter by.
            tax_group: Taxonomic group(s) to filter by.

        Returns:
            List of :class:`DLProfileSummary` objects.
        """
        sql, params = build_profile_search_query(tf_name=tf_name, tax_group=tax_group)
        cur = self._conn.cursor()
        cur.execute(sql, params)

        results: list[DLProfileSummary] = []
        for row in cur.fetchall():
            _db_id, base_id, version, _pmi, name, tg = row
            results.append(
                DLProfileSummary(
                    profile_id=f"{base_id}.{version}",
                    base_id=base_id,
                    version=version,
                    tf_name=name,
                    tax_group=tg,
                )
            )
        return results

    def search_models(
        self,
        *,
        tf_name: str | list[str] | None = None,
        cell_line: str | list[str] | None = None,
        tax_id: int | list[int] | None = None,
        data_type: str | list[str] | None = None,
        model_name: str | list[str] | None = None,
        source: str | list[str] | None = None,
    ) -> list[DLModelSummary]:
        """Search DL models by filter criteria.

        Args:
            tf_name: TF name(s).
            cell_line: Cell line(s).
            tax_id: Taxonomy ID(s).
            data_type: Data type(s).
            model_name: Model architecture name(s).
            source: Data source(s).

        Returns:
            List of :class:`DLModelSummary` objects.
        """
        sql, params = build_model_search_query(
            tf_name=tf_name,
            cell_line=cell_line,
            tax_id=tax_id,
            data_type=data_type,
            model_name=model_name,
            source=source,
        )
        cur = self._conn.cursor()
        cur.execute(sql, params)

        results: list[DLModelSummary] = []
        for row in cur.fetchall():
            (
                _db_id,
                base_id,
                version,
                _pmi,
                name,
                cl,
                tid,
                dt,
                mn,
                src,
                _sid,
                _surl,
            ) = row
            results.append(
                DLModelSummary(
                    model_id=f"{base_id}.{version}",
                    base_id=base_id,
                    version=version,
                    tf_name=name,
                    cell_line=cl,
                    tax_id=tid,
                    data_type=dt,
                    model_name=mn,
                    source=src,
                )
            )
        return results

    # --- Private helpers ---

    def _resolve_version(self, base_id: str, version: int | None, latest_query: str) -> int:
        """Resolve version: return given version or fetch latest."""
        if version is not None:
            return version
        cur = self._conn.cursor()
        cur.execute(latest_query, (base_id,))
        row = cur.fetchone()
        if row is None:
            raise MotifNotFoundError(f"No DL entity found with base ID '{base_id}'")
        return row[0]

    def _fetch_motif_raw(self, motif_db_id: int) -> dict[str, object] | None:
        """Fetch raw DL_MOTIF data by internal DB ID."""
        cur = self._conn.cursor()
        cur.execute(FETCH_MOTIF_BY_ID, (motif_db_id,))
        row = cur.fetchone()
        if row is None:
            return None

        _mid, name, _source, _source_id, _source_ver, num_seqlets, num_hits, motif_type = row
        return {
            "name": name,
            "num_seqlets": num_seqlets,
            "num_hits": num_hits,
            "motif_type": motif_type,
        }

    def _fetch_matrices(self, motif_db_id: int) -> dict[str, Matrix]:
        """Fetch all matrices for a DL_MOTIF by its internal DB ID."""
        cur = self._conn.cursor()
        cur.execute(FETCH_MOTIF_MATRICES, (motif_db_id,))
        matrices: dict[str, Matrix] = {}
        for mtype, mdata in cur.fetchall():
            values = parse_matrix_json(mdata, mtype)
            matrices[mtype] = Matrix(kind=mtype, values=values)
        return matrices

    def _fetch_motif_with_matrices(self, motif_db_id: int) -> dict[str, object] | None:
        """Fetch a DL_MOTIF row and its matrices by internal DB ID.

        Returns a dict with motif data plus a ``"matrices"`` key,
        or None if the motif doesn't exist.
        """
        data = self._fetch_motif_raw(motif_db_id)
        if data is None:
            return None
        data["matrices"] = self._fetch_matrices(motif_db_id)
        return data
