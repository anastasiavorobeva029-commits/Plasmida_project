# In Silico Cloning Tool

A tool for automation and digital modeling (in silico) of the molecular cloning process, enabling the insertion of a DNA fragment from a donor plasmid into a vector plasmid via two restriction sites. The script automatically checks the compatibility of sticky or blunt ends, controls the uniqueness of restriction sites, and generates a ready-to-use recombinant sequence.

---

## Key Features

- **Automated FASTA Parsing:** Reads plasmid files (.fasta, .fa, .fna) using Biopython and structures data in a pandas.DataFrame.
- **Circular DNA Support:** Accurate restriction site search using the circular sequence doubling method. The script detects sites even if they are split by the artificial start/end boundary of the sequence file.
- **Interactive Enzyme Selection:** Built-in console menu featuring a list of commonly used commercial restriction enzymes. Selection can be made either by list number or by exact enzyme name.
- **Strict Compatibility Control (Overhang Compatibility):** Analysis of end types (blunt/sticky), overhang orientations, and their mutual complementarity (annealing check prior to ligation).
- **Visual Junction Analysis:** Informative output of left (Vector → Insert) and right (Insert → Vector) junctions in the console with color highlighting (ANSI codes) and precise cut site markers.
- **FASTA Export with Case Marking:** Saves the final plasmid to a file with 60-character line breaks. For visual clarity, the vector backbone is written in **UPPERCASE**, while the insert is written in **lowercase**.

---

## Algorithm Workflow Diagram

```text
       [ START: main() execution ]
                   │
                   ▼
       ┌────────────────────────┐
       │     Tkinter Windows:   │
       │  Select Donor and      │
       │    Vector Files        │
       └───────────┬────────────┘
                   │
                   ▼
       ┌────────────────────────┐
       │   read_fasta_file()    │
       │  Parse FASTA into      │
       │      Pandas DF         │
       └───────────┬────────────┘
                   │
                   ▼
       ┌────────────────────────┐
       │ select_enzyme_pair()   │
       │  Interactive enzyme    │
       │    pair selection      │
       └───────────┬────────────┘
                   │
                   ▼
       ┌────────────────────────┐
       │   perform_cloning()    │
       │  (Main assembly loop)  │
       └───────────┬────────────┘
                   │
                   ▼
       ┌────────────────────────┐
       │find_site_in_circular_dna│
       │ Double the strand      │──► [ Site coordinate search ]
       │    (DNA * 2)           │
       └───────────┬────────────┘
                   │
                   ▼
       ┌────────────────────────┐
       │     validate_site()    │
       │  Uniqueness check      │──► [ ValueError if 0 or >1 sites ]
       └───────────┬────────────┘
                   │
                   ▼
       ┌────────────────────────┐
       │  check_compatibility() │
       │  Sticky/blunt end      │──► [ Abort: incompatible ends ]
       │    compatibility       │
       └───────────┬────────────┘
                   │
                   ├─ Ends compatible (True)
                   │
                   ▼
       ┌────────────────────────┐
       │ Calculate indices      │
       │   (0-based)            │
       │ Cut Insert and         │
       │  Vector Backbone       │
       └───────────┬────────────┘
                   │
                   ▼
       ┌────────────────────────┐
       │     Ligation:          │
       │ final_plasmid =        │
       │ backbone + insert      │
       └───────────┬────────────┘
                   │
                   ▼
       ┌────────────────────────┐
       │    print_results()     │
       │  Color-coded junction  │
       │    analysis in         │
       │     console            │
       └───────────┬────────────┘
                   │
                   ▼
       ┌────────────────────────┐
       │  Save results          │
       │ Vector (UPPERCASE)     │
       │  Insert (lowercase)    │
       └───────────┬────────────┘
                   │
                   ▼
         [ END: FASTA export ]
```

## System Requirements and Dependencies

All third-party dependencies are automatically installed via the `requirements.txt` file.

| Library | Purpose |
| :--- | :--- |
| **biopython** | Biological sequence processing (`Seq`, `SeqIO`) and the built-in restriction enzyme database (`Restriction`) |
| **pandas** | Storage, filtering, and structuring of metadata from FASTA files |
| **numpy** | Vectorized operations for filtering site coordinates and calculating overhang signs |
| **tkinter** | Native file explorer windows for convenient file selection (`filedialog`) |

---

## Installation and Project Structure

### Repository Structure

```text
plasmida_project/
├── data/                    # Folder containing source FASTA plasmid files
├── README.md                # This documentation
├── requirements.txt         # List of external dependencies
├── test_gene.py             # Main cloning script
└── test_lib.py              # Isolated testing script for Biopython functions
```

---

## Description of key functions and logic

### Search for sites in a circular molecule

The circular topology problem is solved in the `find_site_in_circular_dna function`:
```python
doubled_seq = Seq(sequence_str * 2)
sites = np.array(enzyme.search(doubled_seq))
valid_sites = sorted(list(set(sites[sites <= original_length])))
```
The sequence string is doubled, which allows detection of a site at the physical junction between the beginning and end of the text file. The array of coordinates is then filtered, discarding duplicates that fall outside the original plasmid length.

### Compatibility check of overhangs

The `check_compatibility` function implements strict biological logic:

If both enzymes produce blunt ends `(is_blunt())`, they are always compatible.

If the end types differ (one blunt, one sticky) or the overhang lengths/directions `(ovhg)` do not match - `False` is returned.

For sticky ends, strict complementarity of the protruding sequences is checked:
`donor_overhang == reverse_complement(vector_overhang)`

### Validation of site uniqueness
The `validate_site` function ensures that the selected enzyme introduces exactly one cut into the plasmid. If there are 0 sites or more than 1, the script raises a `ValueError`, preventing accidental destruction of the construct.

Note on coordinates: The script accepts and outputs site coordinates in 1-based indexing format (Biopython/GenBank/SnapGene standard), but during the string slicing stage `(perform_cloning)` it automatically translates them into 0-based Python slices `(position - 1)`.


---

## Example session

```bash
python test_gene.py
Selecting plasmid files
[File selection dialog opens]

 Donor: pEGFP-N1_sequence.fasta - 4733 bp
Vector: pBR322.fasta - 4361 bp


 Selection of enzymes for cloning:
 --- Enzymes for the donor ---
 Available enzymes:
 1. EcoRI    9. NotI
 2. HindIII 10. SacI
 3. BamHI   11. XbaI
 ...

 Enter numbers or names separated by spaces
Your choice: 1 3

--- Enzymes for the vector ---
 Available enzymes:
 ...
Your choice: EcoRI BamHI


 Selected enzymes:
  Donor:  EcoRI + BamHI
  Vector: EcoRI + BamHI


Searching for restriction sites:
  Donor: EcoRI (GAATTC) → 1 site. Positions: 1234
  Donor: BamHI (GGATCC) → 1 site. Positions: 5678
  Vector: EcoRI (GAATTC) → 1 site. Positions: 234
  Vector: BamHI (GGATCC) → 1 site. Positions: 3456
Ends are compatible

Insert:   456 bp
Backbone: 3905 bp
Result:   4361 bp

Unique recombinant plasmid created

 Junction analysis
Left junction (Vector -> Insert):
TGCTGGCCGTGATACTCAGCACCATCTCA * GAATTCGTGAGCAAGGGCGAGGAGCTGTTCAC...

Right junction (Insert -> Vector):
...GGCGAGGAGCTGTTCACCGGGGTGGTGCCC * GGATCCTCTAGAGTCGACCTGCAGGCATGCA...
                                        ^ BamHI (GGATCC) site

Result successfully saved to: /path/to/recombinant_plasmid.fasta
```

---

### 3. **Troubleshooting Section**

This will save users from common errors:

## Troubleshooting

| Problem                                                        | Possible solution                                                                                                                                      |
|:----------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Error `ValueError: EcoRI site not found in Donor`**           | The selected enzyme does not cut the plasmid. Check the list of available enzymes and choose a different one.                                          |
| **Error `ValueError: BamHI site occurs 3 times in Vector`** | The enzyme has multiple restriction sites in the plasmid. Select an enzyme with a unique site.                                                         |
| **Message "Ends are incompatible!"**                             | The enzyme overhangs are not complementary. Choose a different pair of enzymes so that the sites on the donor and vector form compatible sticky ends.. |
| **The program does not see `.fasta` files**                           | Make sure the files have `.fasta`, `.fa`, or `.fna` extensions. Check file access permissions.                                                         |
| **The file selection dialog does not open**                           | Ensure the `tkinter` library is installed (on Linux you may need `sudo apt-get install python3-tk`).                                                       |

---

## Output

After successful execution of the script, a `recombinant_plasmid.fasta` file is created in the folder containing the original vector, with the following structure:
```fasta
>recombinant_plasmid_size_4361bp
ATGACCATGATTACGCCAAGCTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAG
CTCGAATTCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCA
