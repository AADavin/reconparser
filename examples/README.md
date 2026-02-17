# reconparser Examples

This folder contains example scripts and notebooks for testing and using the reconparser library.

## Files

### Test Scripts

- **`test_file_detection.py`** - Simple script to verify ALE files can be found
  - No dependencies required (just pathlib)
  - Good first test to verify file paths are correct

- **`test_with_real_data.py`** - Complete test of all ALE parser functionality
  - Requires: ete3, pandas
  - Tests all three file types (.ucons_tree, .uTs, .uml_rec)
  - Provides detailed output for each parsing step

- **`parse_ale_example.py`** - Basic usage example
  - Template for your own parsing scripts
  - Shows minimal working example

### Notebooks

- **`../notebooks/01_test_ale_parser.ipynb`** - Interactive Jupyter notebook
  - Step-by-step testing and exploration
  - Includes visualizations (if matplotlib installed)
  - Great for debugging path issues
  - **Recommended starting point!**

## Quick Start

### Option 1: Use the Notebook (Recommended)

```bash
cd notebooks
jupyter notebook 01_test_ale_parser.ipynb
```

Then:
1. Update the `ale_base_path` variable in Step 1 to point to your ALE files
2. Run each cell to test the parser step-by-step
3. Debug any path or import issues interactively

### Option 2: Run Test Scripts

First, make sure you're in the examples directory:

```bash
cd reconparser/examples
```

**Step 1: Test file detection**

```bash
python test_file_detection.py
```

If this fails, edit the file and update the path to your ALE files.

**Step 2: Test full parser**

```bash
python test_with_real_data.py
```

## Common Issues

### "File not found" errors

The test scripts use relative paths that assume a specific directory structure:

```
Phylustrator/
├── temp/
│   └── K00192-seqs.ufboot.ale.*    # Your ALE files
└── reconparser/
    └── examples/                    # You are here
```

**Fix**: Edit the test script and update the path:

```python
# Change this line:
test_data_path = "../../temp/K00192-seqs.ufboot.ale"

# To point to your actual file location:
test_data_path = "/path/to/your/data/K00192-seqs.ufboot.ale"
# or
test_data_path = "../your/relative/path/K00192-seqs.ufboot.ale"
```

### "Module not found" errors

Install the required dependencies:

```bash
pip install ete3 pandas numpy
```

Or use conda:

```bash
conda install -c etetoolkit ete3
conda install pandas numpy
```

### Import errors from reconparser

Make sure you're running from the `examples/` directory, or adjust your PYTHONPATH:

```bash
export PYTHONPATH="/path/to/reconparser/src:$PYTHONPATH"
```

## Expected Output

If everything works correctly, you should see:

```
============================================================
Testing ALE Parser with K00192 dataset
============================================================

Checking for ALE output files:
  ✓ Found    consensus_tree       (K00192-seqs.ufboot.ale.ucons_tree)
  ✓ Found    transfers            (K00192-seqs.ufboot.ale.uTs)
  ✓ Found    reconciliation       (K00192-seqs.ufboot.ale.uml_rec)

------------------------------------------------------------
Test 1: Parsing consensus gene tree (.ucons_tree)
------------------------------------------------------------
✓ Successfully parsed gene tree
  - Number of leaves: 101
  - Example leaf names: a2233_67, a2930_87, b1706_2, b65881_64, a839_13
  ...
```

## Next Steps

Once the tests pass, check out:
- `../GETTING_STARTED.md` - Full API documentation
- `parse_ale_example.py` - Template for your own scripts
- Integration examples with Phylustrator in the main README
