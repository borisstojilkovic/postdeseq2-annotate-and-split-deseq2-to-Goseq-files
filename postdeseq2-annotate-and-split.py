# ------------------------------------------------------------
# Imports
# ------------------------------------------------------------
import os       # For file and folder handling
import re       # For regular expression operations (used in ID cleaning)
import pandas as pd  # For data loading, processing, and manipulation
import csv      # For CSV file operations
# ------------------------------------------------------------
# Welcome message
# ------------------------------------------------------------
print(
    "#######\nWelcome to postDESeq2\n"
    "Place all files in the folder you want to process.\n"
    "Each file should be tab-delimited and contain the following columns:\n"
    "GeneID, Base mean, log2FC, StdErr, Wald-Stats, P-value, P-adj.\n"
    "You will also need to provide an annotation file.\n"
)
# ------------------------------------------------------------
# Species selection mode
# ------------------------------------------------------------
# This variable controls how the script joins with the annotation file:
# - "A" = join directly on GeneID
# - "S" = derive 'locus' from GeneID by splitting on '.' before joining
species = "A"

# ------------------------------------------------------------
# Load annotation index file
# ------------------------------------------------------------
# annotations.xlsx contains two columns:
#   type      - short code to identify the annotation set
#   name_file - Excel file name inside the 'annotations/' folder
annotations = pd.read_excel("annotations.xlsx")

# Display available annotation sets and prompt user to choose one by 'type' code
annotation = input(
    f'{annotations[["type", "name_file"]].to_string(index=False)} \n'
    "Type the code from the 'type' column above to select an annotation file: "
).upper()

# Retrieve the corresponding annotation file name from the table
annotation_f = annotations.loc[annotations['type'] == annotation, 'name_file'].item()

# Load the selected annotation Excel file from the 'annotations' folder
annotation_file = pd.read_excel(f"annotations/{annotation_f}")

# ------------------------------------------------------------
# Column name constants
# ------------------------------------------------------------
Accession = "GeneID"      # Column name used for joining when species = "A"
Expression = "Expression" # Column name for GOseq-compatible True/False expression

# Column headers for the annotated output files
header_name = [
    "GeneID", "Base mean", "log2FC", "StdErr", "Wald-Stats", "P-value", "P-adj", 'loc', 'Last',
    "padj<0.05", "log2FC>0 and padj<0.05", "log2FC<0 and padj<0.05",
    "padj<0.01", "log2FC>1 and padj<0.01",
    "log2FC<-1 and padj<0.01", "(log2FC>1 or log2FC<-1) and padj<0.01"
]

# Short names for each filtering condition — used in output file naming
lista = [
    "padj_low_005",
    "log2FC_high_0_padj_low_005",
    "log2FC_low_0_and_padj_low005",
    "padj_low_001",
    "log2FC_high_1_and_padj_low_001",
    "log2FC_low_minus_1_and_padj_low_001",
    "log2FC_higher_1_or_log2FC_low_minus_1_and_padj_001"
] 

# Ensure the output directory exists; create it if it doesn't
isExist = os.path.exists("output") 
if isExist!=True:
    os.mkdir("output")

#add solycIDs based on different padj and Log2FC functions
def padj_low_005(row):
    """Genes with padj < 0.05."""
    if row[6] < 0.05:
        return row[0] 
    else:
        return ""
def log2FC_high_0_padj_low_005(row):
    """Genes with padj < 0.05 and log2FC > 0 (upregulated)."""
    if (row[6] < 0.05) and (row[2]>0) :
        return row[0] 
    else:
        return ""
def log2FC_low_0_and_padj_low005(row):
    """Genes with padj < 0.05 and log2FC < 0 (downregulated)."""
    if (row[6] < 0.05) and (row[2]<0):
        return row[0] 
    else:
        return ""
def padj_low_001(row):
    """Genes with padj < 0.01."""
    if row[6] < 0.01:
        return row[0] 
    else:
        return ""        
def log2FC_high_1_and_padj_low_001(row):
    """Genes with padj < 0.01 and log2FC > 1 (strong upregulation)."""
    if (row[6] < 0.01) and (row[2]>1):
        return row[0] 
    else:
        return ""
def log2FC_low_minus_1_and_padj_low_001(row):
    """Genes with padj < 0.01 and log2FC < -1 (strong downregulation)."""
    if (row[6] < 0.01) and (row[2]<-1):
        return row[0] 
    else:
        return ""
def log2FC_higher_1_or_log2FC_low_minus_1_and_padj_001(row):
    """Genes with padj < 0.01 and |log2FC| > 1 (strong up/down regulation)."""
    if (row[6] < 0.01) and ((row[2]<-1) or (row[2]>1) ):
        return row[0] 
    else:
        return ""

# ------------------------------------------------------------
# Create annotated files and classify genes into subsets
# based on padj and log2FC thresholds (UP/DOWN regulated)
# ------------------------------------------------------------

for file_name in os.listdir('input'):
    # --------------------------------------------------------
    # Load DESeq2 result file
    # --------------------------------------------------------
    df = pd.read_csv(
        f"input/{file_name}",
        header=None,        # no header expected
        sep=r'\s+',         # split on any whitespace (tabs or spaces)
        decimal='.'         # decimal separator is a dot
    )

    # Optional lines if input uses commas or has headers:
    # df = df.iloc[1:, :]                      # drop header row
    # df = df.rename(columns=df.iloc[0])       # use first row as header

    # Convert relevant numeric columns to float
    # Column index: 1=Base mean, 2=log2FC, 3=StdErr, 4=Wald-Stats, 5=P-value, 6=P-adj
    df = df.astype({2: 'float', 1: 'float', 3: 'float', 4: 'float', 5: 'float', 6: 'float'})

    # --------------------------------------------------------
    # Ensure annotated output directory exists
    # --------------------------------------------------------
    if not os.path.exists("output/annotated"):
        os.mkdir("output/annotated")

    # --------------------------------------------------------
    # Clean up GeneID column and prepare locus columns
    # --------------------------------------------------------
    # Remove 'gene:' prefix from GeneID strings
    df = df.replace({0: 'gene:'}, {0: ''}, regex=True)

    # Copy GeneID to column index 7
    df[7] = df[0]

    # If species = "S", split into 'locus' (base ID) and 'Last' (version info)
    if species == "S":
        df[[7, 8]] = df[7].str.split(".", expand=True)
    else:
        # For species = "A", keep the same value in column 8
        df[8] = df[7]

    # --------------------------------------------------------
    # Apply padj/log2FC filter functions to classify genes
    # Each function returns GeneID if condition met, else ""
    # --------------------------------------------------------
    df[9]  = df.apply(padj_low_005, axis=1)
    df[10] = df.apply(log2FC_high_0_padj_low_005, axis=1)
    df[11] = df.apply(log2FC_low_0_and_padj_low005, axis=1)
    df[12] = df.apply(padj_low_001, axis=1)
    df[13] = df.apply(log2FC_high_1_and_padj_low_001, axis=1)
    df[14] = df.apply(log2FC_low_minus_1_and_padj_low_001, axis=1)
    df[15] = df.apply(log2FC_higher_1_or_log2FC_low_minus_1_and_padj_001, axis=1)

    # --------------------------------------------------------
    # Save intermediate annotated file (before merging annotation data)
    # --------------------------------------------------------
    df.to_csv(
        f"output/annotated/{file_name}",
        header=header_name,
        sep='\t',
        index=False,
        decimal="."
    )

    # Reload to ensure clean header handling for merging
    df = pd.read_csv(
        f"output/annotated/{file_name}",
        sep='\t',
        low_memory=False,
        decimal="."
    )

    # Add 'locus' column from 'loc'
    df['locus'] = df['loc']

    # --------------------------------------------------------
    # Merge with annotation file based on species mode
    # --------------------------------------------------------
    if species == "S":
        # Merge on 'locus' and drop extra columns not needed anymore
        inner_join = pd.merge(df, annotation_file, on='locus', how='left')
        inner_join = inner_join.drop(['Last', 'loc'], axis=1)
    elif species == "A":
        # Merge directly on 'GeneID'
        inner_join = pd.merge(df, annotation_file, on=Accession, how='left')

    # --------------------------------------------------------
    # Save final annotated file for this input dataset
    # --------------------------------------------------------
    inner_join.to_csv(
        f"output/annotated/annotated_{file_name}",
        sep='\t',
        index=False,
        decimal="."
    )

    print(f"{file_name} is done with annotation")

      
# ------------------------------------------------------------
# GOseq-compatible classification functions
# ------------------------------------------------------------
# These functions are similar to the earlier padj/log2FC filters,
# but instead of returning the GeneID (or ""), they return the
# string "True" or "False". This is the required format for GOseq
# expression category input.
#
# row[2] = log2FC
# row[6] = padj
# ------------------------------------------------------------

def b_padj_low_005(row):
    """Return True if padj < 0.05, else False."""
    return "True" if row[6] < 0.05 else "False"

def b_log2FC_high_0_padj_low_005(row):
    """Return True if padj < 0.05 and log2FC > 0 (upregulated), else False."""
    return "True" if (row[6] < 0.05 and row[2] > 0) else "False"

def b_log2FC_low_0_and_padj_low005(row):
    """Return True if padj < 0.05 and log2FC < 0 (downregulated), else False."""
    return "True" if (row[6] < 0.05 and row[2] < 0) else "False"

def b_padj_low_001(row):
    """Return True if padj < 0.01, else False."""
    return "True" if row[6] < 0.01 else "False"

def b_log2FC_high_1_and_padj_low_001(row):
    """Return True if padj < 0.01 and log2FC > 1 (strong upregulation), else False."""
    return "True" if (row[6] < 0.01 and row[2] > 1) else "False"

def b_log2FC_low_minus_1_and_padj_low_001(row):
    """Return True if padj < 0.01 and log2FC < -1 (strong downregulation), else False."""
    return "True" if (row[6] < 0.01 and row[2] < -1) else "False"

def b_log2FC_higher_1_or_log2FC_low_minus_1_and_padj_001(row):
    """Return True if padj < 0.01 and |log2FC| > 1 (strong up/down regulation), else False."""
    return "True" if (row[6] < 0.01 and (row[2] < -1 or row[2] > 1)) else "False"

# ------------------------------------------------------------
# Compute GOseq-compatible expression tables
# (True/False classification for each filtering condition)
# ------------------------------------------------------------
for file_name in os.listdir('input'):
    # --------------------------------------------------------
    # Load DESeq2 results
    # --------------------------------------------------------
    df = pd.read_csv(
        f"input/{file_name}",
        header=None,        # no header expected
        sep=r'\s+',         # split on any whitespace (tab or space)
        decimal='.'         # decimal separator is a dot
    )

    # Optional adjustments if needed for specific files:
    # df = df.iloc[1:, :]                     # drop header row
    # df = df.rename(columns=df.iloc[0])      # set first row as header

    # Convert log2FC (index 2) and padj (index 6) to floats
    df = df.astype({2: 'float', 6: 'float'})

    # Remove 'gene:' prefix from GeneID (column 0)
    df = df.replace({0: 'gene:'}, {0: ''}, regex=True)

    # Extract file name without extension for naming outputs
    file_nam_without_ext = file_name.split(".")

    # --------------------------------------------------------
    # Ensure output/expression directory exists
    # --------------------------------------------------------
    if not os.path.exists("output/expression"):
        os.mkdir("output/expression")

    # --------------------------------------------------------
    # Apply each GOseq classification function and save result
    # For each filter:
    # - Apply function to return "True"/"False"
    # - Drop statistical columns, keep only GeneID + Expression
    # - Save as tab-delimited file for GOseq
    # --------------------------------------------------------

    # Filter 1: padj < 0.05
    df2 = df.copy()
    df2[9] = df.apply(b_padj_low_005, axis=1)
    df2 = df2.drop([1, 2, 3, 4, 5, 6], axis=1)
    df2.to_csv(
        f"output/expression/expression_{file_nam_without_ext[0]}_{lista[0]}.tab",
        header=[header_name[0], Expression],
        sep='\t',
        index=False,
        decimal="."
    )

    # Filter 2: padj < 0.05 & log2FC > 0 (upregulated)
    df3 = df.copy()
    df3[9] = df.apply(b_log2FC_high_0_padj_low_005, axis=1)
    df3 = df3.drop([1, 2, 3, 4, 5, 6], axis=1)
    df3.to_csv(
        f"output/expression/expression_{file_nam_without_ext[0]}_{lista[1]}.tab",
        header=[header_name[0], Expression],
        sep='\t',
        index=False,
        decimal="."
    )

    # Filter 3: padj < 0.05 & log2FC < 0 (downregulated)
    df4 = df.copy()
    df4[9] = df.apply(b_log2FC_low_0_and_padj_low005, axis=1)
    df4 = df4.drop([1, 2, 3, 4, 5, 6], axis=1)
    df4.to_csv(
        f"output/expression/expression_{file_nam_without_ext[0]}_{lista[2]}.tab",
        header=[header_name[0], Expression],
        sep='\t',
        index=False,
        decimal="."
    )

    # Filter 4: padj < 0.01
    df5 = df.copy()
    df5[9] = df.apply(b_padj_low_001, axis=1)
    df5 = df5.drop([1, 2, 3, 4, 5, 6], axis=1)
    df5.to_csv(
        f"output/expression/expression_{file_nam_without_ext[0]}_{lista[3]}.tab",
        header=[header_name[0], Expression],
        sep='\t',
        index=False,
        decimal="."
    )

    # Filter 5: padj < 0.01 & log2FC > 1 (strong upregulation)
    df6 = df.copy()
    df6[9] = df.apply(b_log2FC_high_1_and_padj_low_001, axis=1)
    df6 = df6.drop([1, 2, 3, 4, 5, 6], axis=1)
    df6.to_csv(
        f"output/expression/expression_{file_nam_without_ext[0]}_{lista[4]}.tab",
        header=[header_name[0], Expression],
        sep='\t',
        index=False,
        decimal="."
    )

    # Filter 6: padj < 0.01 & log2FC < -1 (strong downregulation)
    df7 = df.copy()
    df7[9] = df.apply(b_log2FC_low_minus_1_and_padj_low_001, axis=1)
    df7 = df7.drop([1, 2, 3, 4, 5, 6], axis=1)
    df7.to_csv(
        f"output/expression/expression_{file_nam_without_ext[0]}_{lista[5]}.tab",
        header=[header_name[0], Expression],
        sep='\t',
        index=False,
        decimal="."
    )

    # Filter 7: padj < 0.01 & |log2FC| > 1 (strong up/down regulation)
    df8 = df.copy()
    df8[9] = df.apply(b_log2FC_higher_1_or_log2FC_low_minus_1_and_padj_001, axis=1)
    df8 = df8.drop([1, 2, 3, 4, 5, 6], axis=1)
    df8.to_csv(
        f"output/expression/expression_{file_nam_without_ext[0]}_{lista[6]}.tab",
        header=[header_name[0], Expression],
        sep='\t',
        index=False,
        decimal="."
    )

    # Progress message
    print(f"{file_name} is done with the expression")

# ------------------------------------------------------------
# Create per-threshold "significant only" tables (+ annotations)
#   - Writes UP/DOWN subsets into output/sigOnly/
#   - Joins with the selected annotation workbook
# ------------------------------------------------------------

# Column headers for the temporary standardized TSV we write/read back
name_header = ["GeneID", "Base mean", "log2FC", "StdErr", "Wald-Stats", "P-value", "P-adj", "locus"]

for file_name in os.listdir('input'):
    # --------------------------------------------------------
    # Load DESeq2 results
    # --------------------------------------------------------
    dfs = pd.read_csv(
        f"input/{file_name}",
        header=None,        # no header expected
        sep=r'\s+',         # split on any whitespace (tabs or spaces)
        decimal='.'         # decimal separator is a dot
    )

    # Some input exports include a header row; drop if present
    dfs = dfs.iloc[1:, :]

    # Ensure numeric types for filters: log2FC (2), padj (6)
    dfs = dfs.astype({2: 'float', 6: 'float'})

    # --------------------------------------------------------
    # Make sure the output/sigOnly directory exists
    # --------------------------------------------------------
    if not os.path.exists("output/sigOnly"):
        os.mkdir("output/sigOnly")

    # --------------------------------------------------------
    # Normalize GeneID text and derive locus/version columns
    # --------------------------------------------------------
    # Remove 'gene:' prefix from GeneID (column 0)
    dfs = dfs.replace({0: 'gene:'}, {0: ''}, regex=True)

    # Copy GeneID to column 7 for downstream splitting
    dfs[7] = dfs[0]

    if species == "S":
        # For species S: split versioned IDs: <locus>.<version...>
        dfs[[7, 8]] = dfs[7].str.split(".", expand=True)
    else:
        # For species A: keep a mirror of the ID in col 8 (not used later)
        dfs[8] = dfs[0]

    # Drop the extra column 8 (we keep col 7 as 'locus' when S)
    dfs = dfs.drop([8], axis=1)

    # --------------------------------------------------------
    # Write a standardized TSV and reload with headers
    # (This guarantees consistent column names before merging)
    # --------------------------------------------------------
    dfs.to_csv(
        f"output/{file_name}",
        header=name_header[0:8],
        sep='\t',
        index=False
    )
    dfs = pd.read_csv(
        f"output/{file_name}",
        sep='\t',
        low_memory=False
    )

    # --------------------------------------------------------
    # Merge with annotation data depending on species mode
    # --------------------------------------------------------
    if species == "S":
        # Join on 'locus' (base ID without version suffix)
        dfs = pd.merge(dfs, annotation_file, on='locus', how='left')
    elif species == "A":
        # Join directly on 'GeneID'
        dfs = pd.merge(dfs, annotation_file, on=Accession, how='left')

    # --------------------------------------------------------
    # Generate significance-filtered subsets and save
    # (Names correspond to entries in `lista`)
    # --------------------------------------------------------
    file_nam_without_ext = file_name.split(".")

    # 1) padj < 0.05
    df2 = dfs[dfs["P-adj"] < 0.05]
    df2.to_csv(
        f"output/sigOnly/sig_{file_nam_without_ext[0]}_{lista[0]}.tab",
        sep='\t', index=False, decimal="."
    )

    # 2) padj < 0.05 & log2FC > 0 (up)
    df3 = dfs[(dfs["P-adj"] < 0.05) & (dfs["log2FC"] > 0)]
    df3.to_csv(
        f"output/sigOnly/sig_{file_nam_without_ext[0]}_{lista[1]}.tab",
        sep='\t', index=False, decimal="."
    )

    # 3) padj < 0.05 & log2FC < 0 (down)
    df4 = dfs[(dfs["P-adj"] < 0.05) & (dfs["log2FC"] < 0)]
    df4.to_csv(
        f"output/sigOnly/sig_{file_nam_without_ext[0]}_{lista[2]}.tab",
        sep='\t', index=False, decimal="."
    )

    # 4) padj < 0.01
    df5 = dfs[dfs["P-adj"] < 0.01]
    df5.to_csv(
        f"output/sigOnly/sig_{file_nam_without_ext[0]}_{lista[3]}.tab",
        sep='\t', index=False, decimal="."
    )

    # 5) padj < 0.01 & log2FC > 1 (strong up)
    df6 = dfs[(dfs["P-adj"] < 0.01) & (dfs["log2FC"] > 1)]
    df6.to_csv(
        f"output/sigOnly/sig_{file_nam_without_ext[0]}_{lista[4]}.tab",
        sep='\t', index=False, decimal="."
    )

    # 6) padj < 0.01 & log2FC < -1 (strong down)
    df7 = dfs[(dfs["P-adj"] < 0.01) & (dfs["log2FC"] < -1)]
    df7.to_csv(
        f"output/sigOnly/sig_{file_nam_without_ext[0]}_{lista[5]}.tab",
        sep='\t', index=False, decimal="."
    )

    # 7) padj < 0.01 & |log2FC| > 1 (strong up or down)
    df8 = dfs[(dfs["P-adj"] < 0.01) & ((dfs["log2FC"] > 1) | (dfs["log2FC"] < -1))]
    df8.to_csv(
        f"output/sigOnly/sig_{file_nam_without_ext[0]}_{lista[6]}.tab",
        sep='\t', index=False, decimal="."
    )

    # Progress feedback for this input file
    print(f"{file_name} is done with sig selection")

print("\nFINISHED!")
