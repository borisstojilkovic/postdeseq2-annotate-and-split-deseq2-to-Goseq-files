# postdeseq2-annotate-and-split

Post-processing pipeline for **DESeq2 differential expression results**:  
- Annotates genes using user-supplied annotation files  
- Splits results into multiple significance and fold-change subsets  
- Produces GOseq-compatible True/False expression files  

Supports multiple species (configurable), using either direct `GeneID` matching or locus-based matching.

---

##  Repository Structure

```
postdeseq2-annotate-and-split/
├─ postdeseq2-annotate-and-split.py   # Main script
├─ requirements.txt                   # Python dependencies
├─ annotations.xlsx                    # Annotation index (see below)
├─ annotations/                        # Folder with annotation files
├─ input/                              # Place your DESeq2 result files here
└─ output/                             # Created automatically by the script
```

---

##  Requirements

### 1. Python
- Version **3.9 – 3.12** recommended

### 2. Python packages
Step 1: Open a terminal (Command Prompt on Windows, Terminal on macOS/Linux)
Step 2: Run:
```bash
pip install pandas openpyxl

```

These are needed for:
- **pandas** → data processing
- **openpyxl** → reading Excel annotation files

---

##  Input Files

### 1. DESeq2 Results (in `input/`)
- **Format:** tab- or whitespace-delimited, **7 columns in this order**:
  1. `GeneID`
  2. `Base mean`
  3. `log2FC`
  4. `StdErr`
  5. `Wald-Stats`
  6. `P-value`
  7. `P-adj`

- Decimal **dot** required (`-1.23`, not `-1,23`)
- Leading `gene:` prefixes are automatically removed
- If your files have headers, the script drops the first row in one stage

---

### 2. Annotation Index (`annotations.xlsx`)
- Tells the script which annotation workbook to load
- Minimal columns:

| type | name_file                     |
|------|--------------------------------|
| A    | arabidopsis_annotations.xlsx   |
| S    | tomato_annotations.xlsx        |

---

### 3. Annotation Workbooks (in `annotations/`)
- Listed in `annotations.xlsx`
- Must contain the join column:
  - **For species `A`**: `GeneID`
  - **For species `S`**: `locus` (base ID without version suffix)

Example for `S`:
```
DESeq2 GeneID:   Solyc05g012345.2.1
Annotation locus: Solyc05g012345
```

---

## Species Modes

Set in the script:
```python
species = "A"  # "A" or "S"
```

- `"A"` → join on **GeneID**
- `"S"` → derive `locus` from GeneID (split at `.`) and join on `locus`

---

## ▶️ How to Run

1. **Install requirements** (only needed once):
```bash
pip install pandas openpyxl
```

2. **Prepare folders and files**:
   - Place DESeq2 result files in `input/`
   - Add `annotations.xlsx` listing available annotation files
   - Put actual annotation `.xlsx` files in `annotations/`

3. **Set species** in the script (e.g.,`"A"` or `"S"`)

4. **Run**:
```bash
python postdeseq2-annotate-and-split.py
```

5. **Choose annotation file** when prompted (must match a `name_file` in `annotations.xlsx`)

###  Choosing an Annotation File

When you run the script, it will display a list of available annotation types (from `annotations.xlsx`) and ask:

```
Write the choose for the name_file above for file:
```

This is how you select the annotation to use:  

1. **Look at the table below** to see all available annotation codes and their descriptions.
2. **Type the code** from the **first column** exactly as shown (case-sensitive).
3. Press **Enter**.

| Code | Annotation file | Description |
|------|-----------------|-------------|
| AT   | `annotation_arab.xlsx` | Arabidopsis thaliana |
| AT1  | `annotation_arab_1.xlsx` | Arabidopsis thaliana (alt. set 1) |
| ATG  | `annotation_arabG.xlsx` | Arabidopsis thaliana (G-based IDs) |
| BD   | `annotation_bdello.xlsx` | Bdellovibrio bacteriovorus |
| SL   | `annotation_tom.xlsx` | Solanum lycopersicum (tomato) |
| SLI  | `annotation_tom_with_RKN_Mi-Tomato.xlsx` | Tomato with *Meloidogyne incognita* (RKN) genes |
| SLS  | `annotation_tom_withouext.xlsx` | Tomato without external gene sets |
| OS   | `annotationOS.xlsx` | Oryza sativa (rice) |
| E    | `MG1655_proteins_167_161521.xlsx` | *E. coli* MG1655 (protein-based IDs) |
| LE   | `locus_MG1655_proteins_167_161521.xlsx` | *E. coli* MG1655 (locus-based IDs) |
| SA   | `archaea_S.acidocaldarius.xlsx` | *Sulfolobus acidocaldarius* (archaea) |
| CC   | `Caulobacter crescentus.xlsx` | *Caulobacter crescentus* |

**Example**:
  
If you want to run the script with Arabidopsis thaliana annotations:  
```
Write the choose for the name_file above for file: AT
```

If you want to run the script with Caulobacter crescentus annotations:  
```
Write the choose for the name_file above for file: CC

---

## Output Folders

### `output/annotated/`
- Full DESeq2 table + filter columns
- `annotated_<file>` includes merged annotation

### `output/sigOnly/`
- Significant subsets at multiple thresholds:
  - `padj_low_005`
  - `log2FC_high_0_padj_low_005`
  - `log2FC_low_0_and_padj_low005`
  - `padj_low_001`
  - `log2FC_high_1_and_padj_low_001`
  - `log2FC_low_minus_1_and_padj_low_001`
  - `log2FC_higher_1_or_log2FC_low_minus_1_and_padj_001`

### `output/expression/`
- GOseq-compatible: `GeneID` + `Expression` (True/False)

---

##  Examples

### Arabidopsis
```bash
# In script: species = "A"
python postdeseq2-annotate-and-split.py
# Prompt: arabidopsis_annotations.xlsx
```

### Tomato
```bash
# In script: species = "S"
python postdeseq2-annotate-and-split.py
# Prompt: tomato_annotations.xlsx
```

### Caulobacter crescentus
```bash
# In script: species = "CC"
python postdeseq2-annotate-and-split.py
# Prompt: Caulobacter crescentus.xlsx
```

---

##  Troubleshooting

- **NaN annotation columns** → check join key (`GeneID` or `locus`)
- **ValueError casting floats** → ensure decimal dots, no commas
- **No outputs** → check `input/` contains files and script has write permissions

---

