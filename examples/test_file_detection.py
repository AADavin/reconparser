"""Test file detection without requiring ete3."""

import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Test file path detection without importing the full module
from pathlib import Path


def test_file_detection():
    """Test that we can detect the ALE output files."""

    # Path to the test data (use absolute path)
    test_data_path = "/sessions/elegant-gifted-rubin/mnt/Phylustrator/temp/K00192-seqs.ufboot.ale"
    base_path = Path(test_data_path)

    # Define expected file paths
    ucons_tree_path = Path(str(base_path) + '.ucons_tree')
    uts_path = Path(str(base_path) + '.uTs')
    uml_rec_path = Path(str(base_path) + '.uml_rec')

    print("=" * 60)
    print("Testing File Detection")
    print("=" * 60)
    print(f"\nBase path: {base_path}")
    print()

    files = {
        'Consensus tree (.ucons_tree)': ucons_tree_path,
        'Transfers (.uTs)': uts_path,
        'Reconciliation (.uml_rec)': uml_rec_path
    }

    all_exist = True
    for name, path in files.items():
        exists = path.exists()
        status = "✓ Found" if exists else "✗ Missing"
        print(f"{status:10} {name:30} {path}")
        if not exists:
            all_exist = False

    print()

    if all_exist:
        print("✓ All required ALE output files found!")
        print("\nFile sizes:")
        for name, path in files.items():
            size_kb = path.stat().st_size / 1024
            print(f"  - {name:30} {size_kb:>10.2f} KB")

        # Test reading first few lines of each file
        print("\n" + "-" * 60)
        print("Sample content from each file:")
        print("-" * 60)

        # Consensus tree
        print("\n1. Consensus tree (.ucons_tree):")
        with open(ucons_tree_path, 'r') as f:
            lines = f.readlines()[:3]
            for i, line in enumerate(lines, 1):
                preview = line.strip()[:80]
                print(f"   Line {i}: {preview}{'...' if len(line.strip()) > 80 else ''}")

        # Transfers
        print("\n2. Transfers (.uTs):")
        with open(uts_path, 'r') as f:
            lines = f.readlines()[:6]
            for i, line in enumerate(lines, 1):
                print(f"   Line {i}: {line.strip()}")

        # Reconciliation (first line only, since it's huge)
        print("\n3. Reconciliation (.uml_rec):")
        with open(uml_rec_path, 'r') as f:
            lines = f.readlines()[:3]
            for i, line in enumerate(lines, 1):
                preview = line.strip()[:80]
                print(f"   Line {i}: {preview}{'...' if len(line.strip()) > 80 else ''}")

        return True
    else:
        print("✗ Some files are missing!")
        return False


if __name__ == "__main__":
    success = test_file_detection()
    sys.exit(0 if success else 1)
