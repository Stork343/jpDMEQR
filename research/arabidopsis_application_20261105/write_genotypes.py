import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml
from bed_reader import create_bed, open_bed

root = Path("data/derived/arabidopsis_20261105")
config = yaml.safe_load(Path("config/discovery_arabidopsis_20261105.yml").read_text())
roles = pd.read_csv(root / "roles.csv", dtype={"accession_id": str})
groups = json.loads((root / "training_sets.json").read_text())
destination = root / "full_genotypes.bed"
assert not destination.exists()
with h5py.File("data/raw/arabidopsis_20261105/1001_SNP_MATRIX/imputed_snps_binary.hdf5", "r") as source:
    source_ids = source["accessions"].asstr()[:]
    lookup = {value: i for i, value in enumerate(source_ids)}
    selected_ids = roles.accession_id.tolist()
    sample_indices = np.array([lookup[value] for value in selected_ids])
    training_indices = [np.array([lookup[value] for value in group]) for group in groups.values()]
    snps = source["snps"]
    keep = np.zeros(snps.shape[0], dtype=bool)
    for start in range(0, snps.shape[0], 10000):
        stop = min(start + 10000, snps.shape[0])
        block = snps[start:stop]
        assert np.isin(block, [0, 1]).all()
        current = np.zeros(len(block), dtype=bool)
        for indices in training_indices:
            frequency = block[:, indices].mean(axis=1)
            current |= np.minimum(frequency, 1 - frequency) >= config["marker_filter"]["training_maf_minimum"]
        keep[start:stop] = current
        if start % 1000000 == 0:
            print("frequency", start, int(keep[:stop].sum()), flush=True)
    retained = np.flatnonzero(keep)
    positions = source["positions"][:]
    regions = source["positions"].attrs["chr_regions"]
    chromosome = np.searchsorted(regions[:, 1], retained, side="right") + 1
    marker_ids = np.array([f"chr{c}_{position}" for c, position in zip(chromosome, positions[retained])])
    assert len(set(marker_ids)) == len(marker_ids)
    properties = {"fid": np.repeat("0", len(selected_ids)), "iid": selected_ids,
                  "chromosome": chromosome, "sid": marker_ids, "bp_position": positions[retained],
                  "allele_1": np.repeat("ONE", len(retained)), "allele_2": np.repeat("ZERO", len(retained))}
    written = 0
    with create_bed(destination, iid_count=len(selected_ids), sid_count=len(retained), properties=properties,
                    count_A1=True) as writer:
        for start in range(0, snps.shape[0], 10000):
            stop = min(start + 10000, snps.shape[0])
            selected = keep[start:stop]
            if not selected.any():
                continue
            block = snps[start:stop][selected][:, sample_indices] * 2
            for vector in block:
                writer.write(vector)
            written += len(block)
            if start % 1000000 == 0:
                print("written", written, flush=True)
    assert written == len(retained)
    reference_rows = np.linspace(0, len(selected_ids) - 1, 12, dtype=int)
    reference_columns = np.linspace(0, len(retained) - 1, 64, dtype=int)
    expected = snps[retained[reference_columns]][:, sample_indices[reference_rows]].T * 2
    actual = open_bed(destination, count_A1=True).read(index=(reference_rows, reference_columns))
    assert np.array_equal(actual, expected)
    reference = pd.DataFrame(expected, columns=marker_ids[reference_columns])
    reference.insert(0, "accession_id", np.array(selected_ids)[reference_rows])
    reference.to_csv(root / "genotype_reader_reference.csv", index=False)
    np.save(root / "retained_source_variant_indices.npy", retained)
audit = Path("results/preflight/empirical_selection_20261101/arabidopsis")
(audit / "genotype_preparation.json").write_text(json.dumps({
    "status": "passed", "accessions": len(selected_ids), "source_variants": len(keep),
    "training_union_variants": len(retained), "checked_calls": int(expected.size),
    "encoding": "PLINK_dosage_equals_twice_author_binary_code",
    "allele_labels": "ONE_and_ZERO_identify_source_binary_states", "outcomes_used_for_variant_filter": False,
    "formal_validation_used": False,
}, indent=2))
