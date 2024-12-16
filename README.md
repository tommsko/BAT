
# BAT - Biomolecule Annotation Toolkit

This tool attempts to identify simulated elements (atoms, molecules) in GROMACS 
input files (.tpr) and associated simulation files.

Please read https://is.muni.cz/th/pf5xb for further reference to methodology of this software.


## Installation

Clone the project
```bash
  git clone https://gitlab.fi.muni.cz/xpavlik8/bat.git
  cd bat
```

Create virtual environment and install requirements. This process was tested on linux and although
it should work in windows, it is likely installation of some components (e.g. networkx) will be complicated.
```bash
  virtualenv venv --python=python3.12
  source venv/bin/activate
  pip install -r requirements.txt
```  

---

## Setup

### Create local BLAST databases (optional, but highly recommended) 

Install blast+ (https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) and make sure it is available in PATH

Build BLAST protein database (PDB)
```bash
  cd <blastp_database_location>
  wget https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
  gunzip pdb_seqres.txt.gz
  makeblastdb -in pdb_seqres.txt -dbtype prot -title pdb -out pdb -parse_seqids
  
  == or alternatively ==
  cd <blastp_database_location>
  wget https://ftp.ncbi.nlm.nih.gov/blast/db/pdbaa.tar.gz
  tar -xvzf pdbaa.tar.gz
```

Build BLAST nucleotide database (PDB)
```bash
  cd <blastn_database_location>
  wget https://ftp.ncbi.nlm.nih.gov/blast/db/pdbnt.tar.gz
  tar -xvzf pdbnt.tar.gz
```

Update `config.ini` such that the following values are set
```ini
[BLAST]
use_local_instance: yes
local_instance_blastp_exec: <blastp_executable_path>
local_instance_blastn_exec: <blastn_executable_path>
local_instance_database_prot: <blastp_database_location>/pdb (or pdbaa)
local_instance_database_nt: <blastn_database_location>/pdbnt
```

### Setup commands for uploading of PDB files

Structural search using RCSB PDB requires PDB file to be remotely accessible using public URL.
Update configuration such the workflow can upload, access, and delete temporary files used as query.
Note that workflow requires ssh keys to be set up, i.e. does not support entering of a password.

Example:
```ini
[PDB_UPLOAD_UPLOAD_COMMAND]
0 = scp
1 = {FILEPATH}
2 = my_username@public_server.cz:/home/my_username/public_html/{FILENAME}
# requires {FILEPATH} (path to source local file) and {FILENAME} {name of the file on remote}
```

```ini
[PDB_UPLOAD_URL]
url: https://www.public_server.cz/~my_username/{FILENAME}
# requires {FILENAME} {name of the file on remote}
```

```ini
[PDB_UPLOAD_DELETE_COMMAND]
0 = ssh
1 = my_username@public_server
2 = rm /home/my_username/public_html/{FILENAME}
# requires {FILENAME} {name of the file on remote}
```

---

## Features/Examples

This toolkit almost always works in two (possibly) identical directories: dataset (containing simulations to analyse)
and dataset_results (containing produced analyses).

In dataset directory (for now `example_dataset`), we have GROMACS input files (.tpr) to analyse.
They may be supplemented by additional files used to generate such input file (e.g. .gro, .itp, ...), 
which will be used to gather additional information that may aid in resolution process. To use them, set
```ini
[input_lookup]
perform_lookup: yes
```
in configuration. If `input_lookup_relaxed_name_requirements` is set to no, such additional files
must have the same name (without extension) as GROMACS input file (.tpr). Otherwise, only prefix must match.
Beware that if wrong additional files are looked-up, i.e. lookup precondition is not true for the data, loading the simulation will fail.

Finally, all commands and their arguments are described in the help of BAT:
```bash
python bat.py --help
```

### 1. Resolution mode
In resolution mode, BAT attempts to resolve all molecules in the simulation to standardized identifiers.
It may require a little fiddling with timeouts in the configuration.

For an example, let's try to resolve all molecules in the example dataset, using 3 threads
```bash
python bat.py -r -i example_dataset -o example_dataset_results -t 3
```

After finishing, we can see the results in the directory *example_dataset_results*, each simulation represented by one JSON file.

**Please note that running the resolution mode is mandatory and must be the first thing done by the workflow**

### 2. Statistics mode
Statistics mode allows for easily seeing performance of the resolution process (works during resolution too).

Let's look at statistics of our example dataset:
```bash
python bat.py -s -i example_dataset -i2 example_dataset_results -o resolution_stats.json
```
For accessing underlying data behind the statistics, we can look at `resolution_stats.json`

### 3. Discover mode
Discover mode allows for generating partial results for molecules from available interpretation paradigms.

Continuing our example, in simulation `5025392_OmpF_SUC`, molecule `seg_4_POPE` was not resolved.
```bash
python bat.py -d -i example_dataset -o export -fs 5025392_OmpF_SUC -fm seg_4_POPE
```

Utilizing data from the simulation, resolvers were not able to determine what kind of molecule it is.
However, after some rudimentary checking, we can see it is [1-Palmitoyl-2-oleoylphosphatidylethanolamine](https://www.ebi.ac.uk/chembl/explore/compound/CHEMBL285376).
In order to make BAT remember this molecule, we will utilize fallback mechanism of fingerprints.

### 4. Fingerprint mode
If molecule cannot be identified, very simple but reliable fingerprint is created. For its location,
please see configuration.

```ini
[fingerprint]
enabled: yes
unknown_fingerprints_export_dir: unknown_hashes
```

We can list all generated fingerprints by:
```bash
python bat.py -lf -i unknown_hashes -o fingerprints.json
```
Which also generates `fingerprints.json` with all the information.

Let's register the fingerprint for POPE above:
```bash
python bat.py -rf -i unknown_hashes
# yes -> yes -> InChI -> WTJKGGKOPKCXLL-MRCUWXFGSA-N -> no
```

Repeating the process in **Resolve** or **Discovery mode**, we can see that `seg_4_POPE` is now resolved.

### 5. Metadata mode
In metadata mode, additional information is fetched to the identifiers resolved previously.
Please see corresponding configuration section for modifying behaviour.

Finally, let's see the example on our test dataset (it does take LONG time):
```bash
python bat.py -m -i example_dataset -o example_dataset_results
```

## Authors

- [@tomsko](tomas.pavlik5055@gmail.com)
