#!/usr/bin/env python3
"""Example: parsing AleRax output with reconparser.

This script demonstrates both run-level and family-level access
to AleRax reconciliation results.
"""

from reconparser import AleRaxRun


def main():
    # ── Point to the AleRax output directory ─────────────────────────
    run = AleRaxRun("data/D1_24genes_uniform/output_UNIFORM_GLOBAL")
    print(run)

    # ── Run-level information ────────────────────────────────────────
    info = run.get_run_info()
    print(f"\nAleRax {info['version']}  |  "
          f"{info['num_families']} families  |  "
          f"{info['num_species']} species  |  "
          f"model: {info['rec_model']}")

    # Species tree
    species_tree = run.get_species_tree()
    print(f"Species tree: {len(species_tree)} leaves")

    # Model parameters (DTL rates)
    params = run.get_model_parameters()
    print(f"\nModel parameters (first 3 genes):\n{params.head(3)}")

    # Per-family likelihoods
    lk = run.get_per_family_likelihoods()
    print(f"\nTotal log-likelihood: {run.get_total_log_likelihood():.2f}")
    print(f"Best family: {lk.iloc[-1]['family']} "
          f"({lk.iloc[-1]['log_likelihood']:.2f})")
    print(f"Worst family: {lk.iloc[0]['family']} "
          f"({lk.iloc[0]['log_likelihood']:.2f})")

    # Global transfers (top 5)
    transfers = run.get_transfers()
    print(f"\nTop 5 global transfers ({len(transfers)} total):")
    print(transfers.head(5).to_string(index=False))

    # ── Family-level analysis ────────────────────────────────────────
    family_name = "K00192"
    fam = run.get_family(family_name)
    print(f"\n{'='*50}")
    print(f"Family: {family_name}  ({fam.get_num_samples()} samples)")
    print(f"{'='*50}")

    # Consensus tree
    consensus = fam.get_consensus_tree()
    print(f"Consensus tree: {len(consensus)} leaves")

    # Summary transfers (top 5)
    fam_transfers = fam.get_transfers()
    print(f"\nTop 5 family transfers ({len(fam_transfers)} total):")
    print(fam_transfers.head(5).to_string(index=False))

    # Event counts across all samples (mean ± std)
    all_counts = fam.get_all_event_counts()
    stats = all_counts[["S", "D", "T"]].agg(["mean", "std"])
    print(f"\nEvent counts (mean ± std across {len(all_counts)} samples):")
    for col in ["S", "D", "T"]:
        print(f"  {col}: {stats.loc['mean', col]:.1f} ± {stats.loc['std', col]:.1f}")

    # Reconciled gene tree with annotations
    rec_tree = fam.get_reconciled_gene_tree()
    t_count = sum(1 for n in rec_tree.traverse() if n.name and ".T@" in n.name)
    d_count = sum(1 for n in rec_tree.traverse() if n.name and ".D@" in n.name)
    print(f"\nReconciled tree annotations: {t_count} transfers, {d_count} duplications")


if __name__ == "__main__":
    main()
