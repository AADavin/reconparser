"""Test the ALE parser with real data from the temp folder."""

import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from reconparser import ALEParser


def main():
    """Test parser with the K00192 dataset."""

    # Path to the test data - adjust this to your file location!
    # This assumes the temp folder is in the parent directory (Phylustrator)
    test_data_path = "../../temp/K00192-seqs.ufboot.ale"

    # Convert to absolute path for better error messages
    test_data_path = Path(test_data_path).resolve()

    print("=" * 60)
    print("Testing ALE Parser with K00192 dataset")
    print("=" * 60)
    print(f"Looking for files at: {test_data_path}")
    print()

    # Initialize parser
    parser = ALEParser(test_data_path)

    # Check which files exist
    print("Checking for ALE output files:")
    files_status = parser.files_exist()
    for file_type, exists in files_status.items():
        status = "✓ Found" if exists else "✗ Missing"
        path = parser.get_file_paths()[file_type]
        print(f"  {status:10} {file_type:20} ({path.name})")
    print()

    # Test 1: Parse consensus gene tree
    print("-" * 60)
    print("Test 1: Parsing consensus gene tree (.ucons_tree)")
    print("-" * 60)
    try:
        gene_tree = parser.get_consensus_tree()
        leaves = gene_tree.get_leaves()
        print(f"✓ Successfully parsed gene tree")
        print(f"  - Number of leaves: {len(leaves)}")
        print(f"  - Example leaf names: {', '.join([l.name for l in leaves[:5]])}")
        print(f"  - Tree string (first 150 chars): {gene_tree.write(format=1)[:150]}...")
        print()
    except Exception as e:
        print(f"✗ Failed to parse gene tree: {e}")
        print()

    # Test 2: Parse transfers
    print("-" * 60)
    print("Test 2: Parsing transfer events (.uTs)")
    print("-" * 60)
    try:
        transfers = parser.get_transfers()
        print(f"✓ Successfully parsed transfers")
        print(f"  - Total transfer events: {len(transfers)}")
        print(f"  - Frequency range: {transfers['freq'].min():.4f} to {transfers['freq'].max():.4f}")
        print(f"  - Columns: {list(transfers.columns)}")
        print(f"\n  Top 10 transfers by frequency:")
        top_transfers = transfers.nlargest(10, 'freq')
        for idx, row in top_transfers.iterrows():
            print(f"    {row['from']:15} → {row['to']:15} (freq: {row['freq']:.4f})")
        print()
    except Exception as e:
        print(f"✗ Failed to parse transfers: {e}")
        print()

    # Test 3: Parse reconciled tree
    print("-" * 60)
    print("Test 3: Parsing reconciled species tree (.uml_rec)")
    print("-" * 60)
    try:
        species_tree = parser.get_reconciled_tree()
        leaves = species_tree.get_leaves()
        print(f"✓ Successfully parsed reconciled tree")
        print(f"  - Number of leaves: {len(leaves)}")
        print(f"  - Example leaf names: {', '.join([l.name for l in leaves[:5]])}")
        print(f"  - Tree string (first 150 chars): {species_tree.write(format=1)[:150]}...")
        print()
    except Exception as e:
        print(f"✗ Failed to parse reconciled tree: {e}")
        print()

    # Test 4: Get transfers for Phylustrator
    print("-" * 60)
    print("Test 4: Converting transfers for Phylustrator")
    print("-" * 60)
    try:
        phylustrator_transfers = parser.get_transfers_for_phylustrator()
        print(f"✓ Successfully converted transfers")
        print(f"  - Number of transfer events: {len(phylustrator_transfers)}")
        print(f"  - Format: List of dictionaries")
        print(f"  - Example transfer: {phylustrator_transfers[0]}")
        print()
    except Exception as e:
        print(f"✗ Failed to convert transfers: {e}")
        print()

    print("=" * 60)
    print("All tests completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
