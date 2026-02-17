"""Example: Parsing ALE reconciliation outputs."""

import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from reconparser import ALEParser


def main():
    """Demonstrate ALE parser usage."""

    # Initialize parser with base path
    # This works with "results.ale" or specific files like "results.ale.ucons_tree"
    parser = ALEParser("../path/to/results.ale")

    # Check which files exist
    print("Available files:")
    for file_type, exists in parser.files_exist().items():
        status = "✓" if exists else "✗"
        print(f"  {status} {file_type}")
    print()

    # Parse consensus gene tree
    try:
        gene_tree = parser.get_consensus_tree()
        print(f"Consensus gene tree loaded:")
        print(f"  - {len(gene_tree.get_leaves())} leaves")
        print(f"  - Tree: {gene_tree.write(format=1)[:100]}...")
        print()
    except FileNotFoundError as e:
        print(f"Could not load gene tree: {e}\n")

    # Parse transfer events
    try:
        transfers = parser.get_transfers()
        print(f"Transfer events loaded:")
        print(f"  - {len(transfers)} transfers")
        print(f"  - Frequency range: {transfers['freq'].min():.3f} - {transfers['freq'].max():.3f}")
        print(f"\nTop 5 transfers by frequency:")
        print(transfers.nlargest(5, 'freq').to_string(index=False))
        print()
    except FileNotFoundError as e:
        print(f"Could not load transfers: {e}\n")

    # Parse reconciled tree
    try:
        species_tree = parser.get_reconciled_tree()
        print(f"Reconciled species tree loaded:")
        print(f"  - {len(species_tree.get_leaves())} leaves")
        print(f"  - Tree: {species_tree.write(format=1)[:100]}...")
        print()
    except FileNotFoundError as e:
        print(f"Could not load reconciled tree: {e}\n")

    # Get transfers in Phylustrator format
    try:
        phylustrator_transfers = parser.get_transfers_for_phylustrator()
        print(f"Transfers for Phylustrator: {len(phylustrator_transfers)} events")
        print(f"Example: {phylustrator_transfers[0]}")
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    main()
