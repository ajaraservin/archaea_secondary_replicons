# Downloading and curating the archaeal dataset
# Overview
This workflow was used to download complete archaeal genomes from the
NCBI RefSeq database and identify assemblies containing secondary 
replicons (secondary chromosomes and/or plasmids). This was dond in 
four steps:

1. Download complete archaeal RefSeq assemblies and sequence reports.
2. Identify assemblies containing more than one chromosome and/or plasmid.
3. Remove single-replicon assemblies from the working dataset.
4. Retrieve and verify genomic and protein FASTA files for all retained
   multi-replicon assemblies.

The resulting dataset was used for the subsequent reclassification of
replicons and downstream functional and taxonomic annotation.

# 1. Download complete archaeal assemblies from NCBI RefSeq
The NCBI Datasets command-line tool was used to download complete
archaeal assemblies from RefSeq, including the sequence reports.

```bash
datasets download genome taxon archaea \
    --assembly-source refseq \
    --assembly-level complete \
    --include seq-report \
    --filename archaea_refseq_complete_reports.zip
```

Create a directory for the archaeal dataset and extract the downloaded
archive. The resulting dataset is organized as:

```text
archaea/
└── ncbi_dataset/
    └── data/
        ├── assembly_data_report.jsonl
        ├── GCF_XXXXXXXXX.X/
        │   └── sequence_report.jsonl
        ├── GCF_XXXXXXXXX.X/
        │   └── sequence_report.jsonl
        └── ...
```

# 2. Identify multi-replicon archaeal genomes
The NCBI assembly and sequence reports were parsed to identify
assemblies containing more than one chromosome and/or plasmid.

The Python script used for this step is:
```bash
python3 identify_multireplicon_genomes.py
```

The script extracts the following information for each chromosome or
plasmid: assembly accession, NCBI taxonomic ID, organism name, BioSample 
accession isolation source, replicon role, assignedMoleculeLocationType, 
replicon name, replicon accession, replicon length, replicon GC content,
number of chromosomes, number of plasmids, total number of replicons,
total genome size, total chromosome size, total plasmid size, fraction 
of the genome represented by plasmids.

The script also creates a list of assemblies containing more than one
replicon and a list of sequence reports that could not be parsed.

# Missing replicon length

When a chromosome or other replicon did not have an explicit length in
the NCBI sequence report, its length was inferred from the reported GC
count and GC percentage:
```text
length = GC_count / (GC_percentage / 100)
```
This reconstructs the sequence length from the GC statistics provided by
NCBI.

# 3. Output files
The script produces the following files:

```text
archaea_multireplicon_replicons.tsv
A replicon-level metadata table with one row per chromosome or plasmid 
belonging to a retained multi-replicon assembly.
archaea_keep_assemblies.txt
This file contains one assembly accession per line for assemblies
containing more than one chromosome and/or plasmid. The accessions in 
this file are used to define the curated multi-replicon archaeal dataset.
archaea_seqreport_unparsed.tsv
This file contains assemblies for which no chromosome or plasmid
replicons could be parsed from the NCBI sequence report.
```

# 4. Curate the dataset by retaining only multi-replicon assemblies
The `retain_multireplicon.py` script removes assembly directories that are not present
in the retained-accession list (`archaea_keep_assemblies.txt`).
After running it, the remaining assembly directories correspond to the
curated multi-replicon archaeal dataset:
```text
archaea/ncbi_dataset/data/GCF_*
```

Check the number of remaining assemblies:
```bash
find archaea/ncbi_dataset/data \
    -maxdepth 1 \
    -type d \
    -name "GCF_*" |
    wc -l
```
```text
243
```

# 5. Retrieve genomic and protein sequences
The retained assemblies were used to retrieve their genomic and protein
FASTA files using the dl_driver.sh script. The retrieval was performed 
from the directory containing the retained `GCF_*` assembly directories. 
The accession list was divided into batches and submitted as SLURM jobs. 
For example:
```bash
TOTAL=$(wc -l < archaea_keep_assemblies.txt)
STEP=500
for ((i=1; i<=TOTAL; i+=STEP)); do
    end=$((i+STEP-1))
    (( end > TOTAL )) && end=$TOTAL
    sbatch \
        --export=ALL,START_LINE=$i,END_LINE=$end \
        dl_driver.sh
done
```
Failed accessions were collected separately in failed_accessions.txt 
and downloaded again. The downloaded archive was then extracted and the 
genomic and protein FASTA files were copied into the corresponding existing 
assembly directory. Recovery loops like the following were used:
```bash
while read -r acc; do
    [ -z "$acc" ] && continue
    datasets download genome accession "$acc" \
        --include genome,protein \
        --filename "${acc}.zip"
done < failed_accessions.txt
```

# 6. Check for incomplete assembly directories
Each retained assembly should contain both:
- a genomic FASTA file (`*.fna`)
- a protein FASTA file (`*.faa`)

Check for incomplete directories:
```bash
for d in GCF_*; do
    if ! compgen -G "$d/*genomic.fna" > /dev/null || \
       ! compgen -G "$d/*protein.faa" > /dev/null; then
        echo "$d"
    fi
done > missing_genomes.txt
```

Count missing genomic and protein FASTA files separately:
```bash
missing_fna=0
missing_faa=0
both_missing=0

for d in GCF_*; do
    has_fna=0
    has_faa=0
    compgen -G "$d/*genomic.fna" > /dev/null && has_fna=1
    compgen -G "$d/*protein.faa" > /dev/null && has_faa=1
    if [[ $has_fna -eq 0 ]]; then
        ((missing_fna++))
    fi
    if [[ $has_faa -eq 0 ]]; then
        ((missing_faa++))
    fi
    if [[ $has_fna -eq 0 && $has_faa -eq 0 ]]; then
        ((both_missing++))
    fi
done

echo "Missing genomic.fna: $missing_fna"
echo "Missing protein.faa: $missing_faa"
echo "Missing both:        $both_missing"
```

A complete dataset should return:
```text
Missing genomic.fna: 0
Missing protein.faa: 0
Missing both:        0
```
# 7. Extract original replicon annotations

The genomic FASTA files retrieved for the curated assemblies contain
the original NCBI replicon annotation in their sequence headers.
Chromosome and plasmid sequences are identified directly in the FASTA
description. For example, a genomic FASTA may contain:

    >NZ_CP120468.1 Halobaculum limi strain YSMS11 chromosome, complete genome
    >NZ_CP120469.1 Halobaculum limi strain YSMS11 plasmid unnamed, complete sequence

These annotations were parsed to generate a genome-specific
`*_genomic_replicons.tsv` file containing the contig accession and its
original replicon classification. The resulting files have the structure:
    contig_accession    replicon_type
For example:
    NZ_CP120468.1    chr1
    NZ_CP120469.1    plasmid

These files preserve the original NCBI replicon annotation and were
used as input for the metadata construction and subsequent replicon
analyses.

# 8. Predict genes and proteins with Prokka

The curated genomic FASTA files were annotated using Prokka with the
Archaea-specific annotation mode. Prokka was run independently for
each retained genome.

The annotation workflow was implemented as:
 ```bash
    #!/bin/bash
    module load prokka
    mkdir -p prokka_output
    for genome in arch_genomes/ncbi_dataset/data/GCF_*/\*.fna
    do
        acc=$(basename "$genome" | sed 's/_ASM.*//')
        if [ -f "prokka_output/${acc}/${acc}.faa" ]; then
            echo "Skipping ${acc} (already completed)"
            continue
        fi
        prokka \
            --kingdom Archaea \
            --cpus 4 \
            --force \
            --outdir prokka_output/${acc} \
            --prefix ${acc} \
            "$genome"
    done
```
Prokka generated genome-specific annotation files, including GFF and
protein FASTA files. The GFF files were subsequently used to associate
predicted proteins with their genomic contigs, while the protein FASTA
files were combined for downstream protein clustering.

# Associated scripts andd datasets
identify_multireplicon_genomes.py   # JSONL parsing and filtering
archaea_keep_assemblies.txt         # Accessions selected for downstream analysis
retain_multireplicon.py             # Selecting those genomes with more than one replicon
dl_driver.sh                        # Sequence downloads