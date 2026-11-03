#!/usr/bin/env python3

import argparse
import contextlib
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

import freesasa
import gemmi


# Constants

RNA_RESIDUES = frozenset({
    "A", "C", "G", "U",
    "DA", "DC", "DG", "DT",
})

DNA_RESIDUES = frozenset({
    "DA", "DC", "DG", "DT",
})

# Maps CLI --polymer value to a residue filter set.
# None means no filtering (ALL).
POLYMER_FILTER = {
    "RNA": RNA_RESIDUES,
    "DNA": DNA_RESIDUES,
    "ALL": None,
}


# Stderr suppression

@contextlib.contextmanager
def suppress_stderr():
    """
    Redirect FD 2 to /dev/null for the duration of the block.

    Used once around build_freesasa_structure() to silence FreeSASA's
    classifier warnings for non-standard atoms, rather than per atom.
    """
    sys.stderr.flush()
    old_stderr = os.dup(2)
    with open(os.devnull, "w") as devnull:
        try:
            os.dup2(devnull.fileno(), 2)
            yield
        finally:
            sys.stderr.flush()
            os.dup2(old_stderr, 2)
            os.close(old_stderr)


# mmCIF / Gemmi helpers

def read_structure(structure_file):
    """Read a PDB/mmCIF structure with Gemmi."""
    return gemmi.read_structure(
        str(structure_file),
        format=gemmi.CoorFormat.Detect,
    )


def read_cif_block(structure_file):
    """
    Parse the mmCIF file once and return the sole data block.

    Pass the returned block to every helper that needs mmCIF data
    to avoid reading the file multiple times.
    """
    return gemmi.cif.read(str(structure_file)).sole_block()


def get_structure_id(cif_block, structure_file):
    """Return the structure identifier from an already-parsed CIF block."""
    try:
        entry_id = cif_block.find_value("_entry.id")
        if entry_id and entry_id not in (".", "?"):
            return str(entry_id)
    except Exception:
        pass

    return Path(structure_file).stem.split(".")[0]


def get_entity_metadata(cif_block):
    """
    Read entity metadata from the mmCIF _entity category.

    Returns
    -------
    dict
        entity_id -> {"entity_type": str, "entity_description": str}
    """
    table = cif_block.find_mmcif_category("_entity.")
    if not table:
        return {}

    metadata     = {}
    entity_ids   = table.find_column("id")
    entity_types = table.find_column("type")
    descriptions = table.find_column("pdbx_description")

    for i, entity_id in enumerate(entity_ids):
        entity_type = str(entity_types[i]) if i < len(entity_types) else ""
        description = str(descriptions[i]) if i < len(descriptions) else ""

        # Treat mmCIF unknown/missing markers as empty.
        if entity_type in (".", "?"):
            entity_type = ""
        if description in (".", "?"):
            description = ""

        metadata[str(entity_id)] = {
            "entity_type":        entity_type,
            "entity_description": description,
        }

    return metadata


def build_auth_to_entity_map(cif_block):
    """
    Build auth_asym_id -> entity_id from already-parsed mmCIF data.

    Strategy, most to least efficient:

    1. _struct_asym: label_asym_id -> entity_id  [O(C chains)]
    2. _pdbx_poly_seq_scheme: auth_asym_id -> label_asym_id  [O(R residues)]
       Covers all polymer chains in standard wwPDB deposits.
    3. _atom_site fallback for any chains still unresolved after step 2,
       with a seen-set for early termination  [O(N atoms) worst case,
       O(first atom per unresolved chain) in practice].
    """
    # Step 1: label_asym_id -> entity_id
    label_to_entity = {}
    struct_asym = cif_block.find_mmcif_category("_struct_asym.")
    if struct_asym:
        for asym_id, entity_id in zip(
            struct_asym.find_column("id"),
            struct_asym.find_column("entity_id"),
        ):
            label_to_entity[str(asym_id)] = str(entity_id)

    auth_to_entity = {}

    # Step 2: auth_asym_id -> label_asym_id via _pdbx_poly_seq_scheme [O(R)]
    scheme = cif_block.find_mmcif_category("_pdbx_poly_seq_scheme.")
    if scheme:
        auth_col  = scheme.find_column("pdb_strand_id")
        label_col = scheme.find_column("asym_id")
        if auth_col is not None and label_col is not None:
            for auth, label in zip(auth_col, label_col):
                auth = str(auth)
                if auth not in auth_to_entity:
                    entity_id = label_to_entity.get(str(label), "")
                    if entity_id:
                        auth_to_entity[auth] = entity_id

    # Step 3: _atom_site fallback for any chains not covered above [O(N)]
    atom_site = cif_block.find_mmcif_category("_atom_site.")
    if atom_site:
        auth_col  = atom_site.find_column("auth_asym_id")
        label_col = atom_site.find_column("label_asym_id")
        if auth_col is not None and label_col is not None:
            seen = set()
            for auth, label in zip(auth_col, label_col):
                auth = str(auth)
                if auth in seen:
                    continue
                seen.add(auth)
                if auth not in auth_to_entity:
                    entity_id = label_to_entity.get(str(label), "")
                    if entity_id:
                        auth_to_entity[auth] = entity_id

    return auth_to_entity


# Residue metadata

def residue_metadata_from_gemmi(structure, entity_metadata, auth_to_entity):
    """
    Construct residue metadata keyed by
    (model_index, auth_chain, auth_seq_id, insertion_code, auth_comp_id).

    Uses Residue.subchain for the correct label_asym_id so that
    entity lookups are accurate even when auth_asym_id != label_asym_id.
    Entity information is resolved here; no post-processing step needed.
    """
    metadata = {}

    for model_index, model in enumerate(structure, start=1):
        for chain in model:
            auth_chain  = str(chain.name)
            entity_id   = auth_to_entity.get(auth_chain, "")
            entity_info = entity_metadata.get(entity_id, {})

            for residue in chain:
                seqid          = residue.seqid
                auth_seq_id    = str(seqid.num)
                insertion_code = (
                    str(seqid.icode)
                    if seqid.icode not in ("\0", " ")
                    else ""
                )
                auth_comp_id = str(residue.name)

                # residue.subchain holds the label_asym_id from mmCIF.
                # Fall back to auth_chain for PDB files where it is absent.
                label_chain = (
                    str(residue.subchain)
                    if residue.subchain
                    else auth_chain
                )

                key = (
                    model_index,
                    auth_chain,
                    auth_seq_id,
                    insertion_code,
                    auth_comp_id,
                )

                metadata[key] = {
                    "model":              model_index,
                    "auth_chain":         auth_chain,
                    "label_chain":        label_chain,
                    "auth_seq_id":        auth_seq_id,
                    "label_seq_id":       (
                        str(residue.label_seq)
                        if residue.label_seq is not None
                        else ""
                    ),
                    "insertion_code":     insertion_code,
                    "auth_comp_id":       auth_comp_id,
                    "label_comp_id":      auth_comp_id,
                    "entity_id":          entity_id,
                    "entity_type":        entity_info.get("entity_type", ""),
                    "entity_description": entity_info.get("entity_description", ""),
                }

    return metadata



# FreeSASA structure construction

def build_freesasa_structure(model, classifier, include_hydrogen=False):
    """
    Build a FreeSASA Structure directly from Gemmi atoms.

    This function is called inside suppress_stderr() at the call site,
    so classifier/addAtom warnings are silenced without per-atom
    file-descriptor operations.

    resname is stored uppercased in atom_metadata to avoid repeated
    .upper() calls in the SASA aggregation loop.

    Returns
    -------
    fs_structure : freesasa.Structure
    atom_metadata : list of dict
    """
    fs_structure = freesasa.Structure(
        None,
        classifier=classifier,
        options={
            "hetatm":          True,
            "hydrogen":        include_hydrogen,
            "skip-unknown":    False,
            "halt-at-unknown": False,
        },
    )

    atom_metadata = []

    for chain in model:
        auth_chain = str(chain.name)

        # FreeSASA's manual-add API requires a single character chain label;
        # this value is never used for metadata matching.
        fs_chain = auth_chain[-1] if auth_chain else "X"

        for residue in chain:
            seqid          = residue.seqid
            auth_seq_id    = str(seqid.num)
            insertion_code = (
                str(seqid.icode)
                if seqid.icode not in ("\0", " ")
                else ""
            )
            resname_raw   = str(residue.name)
            resname_upper = resname_raw.upper()

            residue_key = (
                auth_chain,
                auth_seq_id,
                insertion_code,
                resname_raw,
            )

            for atom in residue:
                if not include_hydrogen:
                    if str(atom.element.name).upper() == "H":
                        continue

                atom_name = str(atom.name).strip()

                # addAtom() resolves the atomic radius internally and
                # falls back to element-based radii for unknown atoms.
                fs_structure.addAtom(
                    atom_name,
                    resname_raw,
                    auth_seq_id,
                    fs_chain,
                    float(atom.pos.x),
                    float(atom.pos.y),
                    float(atom.pos.z),
                )

                atom_metadata.append({
                    "auth_chain":     auth_chain,
                    "auth_seq_id":    auth_seq_id,
                    "insertion_code": insertion_code,
                    "resname":        resname_raw,
                    "resname_upper":  resname_upper,
                    "atom_name":      atom_name,
                    "residue_key":    residue_key,
                })

    return fs_structure, atom_metadata



# SASA calculation

def make_sasa_parameters(
    algorithm="LeeRichards",
    probe_radius=1.4,
    n_slices=20,
    n_points=1000,
    n_threads=1,
):
    """Construct FreeSASA calculation parameters."""
    if algorithm == "LeeRichards":
        algorithm_value = freesasa.LeeRichards
    elif algorithm == "ShrakeRupley":
        algorithm_value = freesasa.ShrakeRupley
    else:
        raise ValueError(f"Unsupported SASA algorithm: {algorithm}")

    return freesasa.Parameters({
        "algorithm":    algorithm_value,
        "probe-radius": float(probe_radius),
        "n-slices":     int(n_slices),
        "n-points":     int(n_points),
        "n-threads":    int(n_threads),
    })


def calculate_model_sasa(
    model,
    model_index,
    classifier,
    parameters,
    polymer_type="RNA",
    include_hydrogen=False,
    zero_tolerance=1e-8,
):
    """
    Calculate SASA for one structural model.

    All atoms contribute to the calculation so that protein, ligand,
    and solvent atoms correctly occlude the nucleic acid surface.
    Only residues matching polymer_type appear in the returned dict.
    """
    # Suppress FreeSASA classifier warnings once for the whole build,
    # not per atom — avoids O(N) file-descriptor operations.
    with suppress_stderr():
        fs_structure, atom_metadata = build_freesasa_structure(
            model=model,
            classifier=classifier,
            include_hydrogen=include_hydrogen,
        )

    result = freesasa.calc(fs_structure, parameters)

    # Resolve the filter set once before the O(N) loop.
    filter_residues = POLYMER_FILTER[polymer_type]

    residue_sasa = defaultdict(float)

    for atom_index, meta in enumerate(atom_metadata):
        # resname_upper was stored at build time; no .upper() here.
        if filter_residues is not None and meta["resname_upper"] not in filter_residues:
            continue

        sasa = float(result.atomArea(atom_index))

        if abs(sasa) < zero_tolerance:
            sasa = 0.0
        elif sasa < 0:
            raise RuntimeError(
                f"Negative SASA for atom {atom_index}: {sasa}"
            )

        residue_sasa[meta["residue_key"]] += sasa

    return residue_sasa



# Main pipeline

def calculate_rna_sasa(
    structure_file,
    output_file,
    algorithm="LeeRichards",
    probe_radius=1.4,
    n_slices=20,
    n_points=1000,
    n_threads=1,
    model=None,
    polymer_type="RNA",
    include_hydrogen=False,
    zero_tolerance=1e-8,
    precision=4,
):
    """
    Residue-level SASA pipeline.

    mmCIF is the authoritative source of structural identity.
    FreeSASA handles only the geometric SASA calculation.
    The mmCIF file is parsed exactly once.
    """
    structure_file = Path(structure_file)
    output_file    = Path(output_file)

    # Parse structure coordinates and mmCIF metadata — one read each.
    structure = read_structure(structure_file)
    cif_block = read_cif_block(structure_file)

    structure_id   = get_structure_id(cif_block, structure_file)
    entity_meta    = get_entity_metadata(cif_block)
    auth_to_entity = build_auth_to_entity_map(cif_block)
    metadata       = residue_metadata_from_gemmi(
        structure, entity_meta, auth_to_entity
    )

    classifier = freesasa.Classifier()
    parameters = make_sasa_parameters(
        algorithm=algorithm,
        probe_radius=probe_radius,
        n_slices=n_slices,
        n_points=n_points,
        n_threads=n_threads,
    )

    output_rows = []

    for model_index, model_obj in enumerate(structure, start=1):
        if model is not None and model_index != model:
            continue

        residue_sasa = calculate_model_sasa(
            model=model_obj,
            model_index=model_index,
            classifier=classifier,
            parameters=parameters,
            polymer_type=polymer_type,
            include_hydrogen=include_hydrogen,
            zero_tolerance=zero_tolerance,
        )

        for residue_key, sasa in residue_sasa.items():
            auth_chain, auth_seq_id, insertion_code, auth_comp_id = residue_key

            meta = metadata.get(
                (model_index, auth_chain, auth_seq_id, insertion_code, auth_comp_id),
                {},
            )

            output_rows.append({
                "structure_id":       structure_id,
                "model":              model_index,
                "entity_id":          meta.get("entity_id", ""),
                "entity_type":        meta.get("entity_type", ""),
                "entity_description": meta.get("entity_description", ""),
                "auth_chain":         auth_chain,
                "label_chain":        meta.get("label_chain", ""),
                "auth_seq_id":        auth_seq_id,
                "label_seq_id":       meta.get("label_seq_id", ""),
                "insertion_code":     insertion_code,
                "auth_comp_id":       auth_comp_id,
                "label_comp_id":      meta.get("label_comp_id", auth_comp_id),
                "probe_radius_A":     probe_radius,
                "algorithm":          algorithm,
                "n_slices":           n_slices if algorithm == "LeeRichards" else "",
                "n_points":           n_points if algorithm == "ShrakeRupley" else "",
                "sasa_A2":            round(max(0.0, sasa), precision),
            })

    # Write output TSV
    output_path = output_file.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "structure_id",
        "model",
        "entity_id",
        "entity_type",
        "entity_description",
        "auth_chain",
        "label_chain",
        "auth_seq_id",
        "label_seq_id",
        "insertion_code",
        "auth_comp_id",
        "label_comp_id",
        "probe_radius_A",
        "algorithm",
        "n_slices",
        "n_points",
        "sasa_A2",
    ]

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Structure:        {structure_id}")
    print(f"Residues written: {len(output_rows)}")
    print(f"Output:           {output_path}")

    return output_rows



# CLI

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate residue-level SASA directly from "
            "PDB/mmCIF structures using Gemmi + FreeSASA."
        )
    )

    parser.add_argument(
        "structure",
        type=Path,
        help="Input PDB/mmCIF structure file.",
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        type=Path,
        help="Output TSV file.",
    )

    parser.add_argument(
        "--algorithm",
        choices=["LeeRichards", "ShrakeRupley"],
        default="LeeRichards",
        help=(
            "SASA algorithm. "
            "LeeRichards (default) is fast; "
            "ShrakeRupley is more accurate at high point counts."
        ),
    )

    parser.add_argument(
        "--probe-radius",
        type=float,
        default=1.4,
        help="Solvent probe radius in Å. Default: 1.4",
    )

    parser.add_argument(
        "--n-slices",
        type=int,
        default=20,
        help="Lee-Richards slices per atom. Default: 20",
    )

    parser.add_argument(
        "--n-points",
        type=int,
        default=1000,
        help="Shrake-Rupley surface points per atom. Default: 1000",
    )

    parser.add_argument(
        "--n-threads",
        type=int,
        default=1,
        help="Number of threads for the FreeSASA calculation. Default: 1",
    )

    parser.add_argument(
        "--model",
        type=int,
        default=None,
        help="Process only this model number (1-based). Default: all models.",
    )

    parser.add_argument(
        "--polymer",
        choices=["RNA", "DNA", "ALL"],
        default="RNA",
        help=(
            "Polymer type to report SASA for. "
            "All atoms still contribute to surface occlusion. "
            "Default: RNA"
        ),
    )

    parser.add_argument(
        "--include-hydrogen",
        action="store_true",
        help="Include hydrogen atoms in the SASA calculation.",
    )

    parser.add_argument(
        "--zero-tolerance",
        type=float,
        default=1e-8,
        help=(
            "SASA values below this threshold are clamped to zero "
            "to remove floating-point noise. Default: 1e-8"
        ),
    )

    parser.add_argument(
        "--precision",
        type=int,
        default=4,
        help="Decimal places in the output sasa_A2 column. Default: 4",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite the output file if it already exists.",
    )

    args = parser.parse_args()

    # Validation
    if not args.structure.exists():
        parser.error(f"Structure file not found: {args.structure}")

    if args.output.exists() and not args.overwrite:
        parser.error(
            f"Output already exists: {args.output}. Use --overwrite to replace it."
        )

    if args.probe_radius <= 0:
        parser.error("--probe-radius must be greater than 0.")

    if args.n_slices <= 0:
        parser.error("--n-slices must be greater than 0.")

    if args.n_points <= 0:
        parser.error("--n-points must be greater than 0.")

    if args.n_threads <= 0:
        parser.error("--n-threads must be greater than 0.")

    calculate_rna_sasa(
        structure_file=args.structure,
        output_file=args.output,
        algorithm=args.algorithm,
        probe_radius=args.probe_radius,
        n_slices=args.n_slices,
        n_points=args.n_points,
        n_threads=args.n_threads,
        model=args.model,
        polymer_type=args.polymer,
        include_hydrogen=args.include_hydrogen,
        zero_tolerance=args.zero_tolerance,
        precision=args.precision,
    )


if __name__ == "__main__":
    main()