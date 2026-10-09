# Replicon reclassification
# Overview
This step generated replicon-level features to support the identification and reclassification of archaeal replicons. The analysis combined:
1. Replicon sequence statistics.
2. Within-genome GC-content differences.
3. Archaeal GTDB marker counts.
4. Tetranucleotide frequency (TNF) similarity.
5. Curated HMM marker sets for chromosomes, partitioning systems, and plasmids.
6. Marker presence and hit counts from HMM scans.

The main feature table generated in this workflow is `replicon_features_markers.tsv`. The final integrated metadata table, including taxonomy and GC-content features, is `replicon_all_metadata_tax.tsv`.

# 1. Replicon sequence statistics
Individual replicon FASTA files were parsed to retain the original NCBI replicon annotation for each contig. For each genome, the resulting `*_genomic_replicons.tsv` files recorded the contig accession and its original replicon type, for example:
    NZ_CP120468.1    chr1
    NZ_CP120469.1    plasmid

These per-genome files were consolidated into `replicon_metadata.tsv`:
```bash
    (echo -e "accession\tcontig\treplicon_type"
      find ./arch_genomes/ncbi_dataset/data/GC*/*_genomic_replicons.tsv | while read f; do
        accession=$(basename "$f" | sed 's/_ASM.*//')
        awk -v acc="$accession" 'BEGIN{OFS="\t"}{print acc,$1,$2}' "$f"
      done
    ) > replicon_metadata.tsv
```
Replicon sequence length and GC content were then retained together with the genome, contig, and original replicon annotation to form the replicon-level metadata used throughout the reclassification workflow. The original `replicon` annotation was retained as the starting annotation and was not treated as the final reclassification.
The resulting metadata included: genome, contig, replicon, length, gc_content

# 2. Within-genome GC-content comparison
GC content was compared among replicons within each genome. For each genome, the largest replicon was identified based on sequence length, and its GC content was used as the reference value (`chr1_gc`). The difference between the GC content of each replicon and this reference was then calculated, together with its absolute value.
The reference is defined by replicon length rather than by the original `chr1` annotation; the column is retained as `chr1_gc` for consistency with the downstream metadata.
```bash
    awk -F'\t' '
    BEGIN{OFS="\t"}
    FNR==NR{
        if(FNR==1) next
        genome=$1
        len=$4
        gc=$5
        if(!(genome in maxlen) || len>maxlen[genome]){
            maxlen[genome]=len
            chr1gc[genome]=gc
        }
        next
    }
    FNR==1{
        print $0,"chr1_gc","delta_gc","abs_delta_gc"
        next
    }
    {
        genome=$1
        delta=$5-chr1gc[genome]
        absdelta=(delta<0)?-delta:delta
        print $0,chr1gc[genome],delta,absdelta
    }' replicon_metadata_tax.tsv replicon_metadata_tax.tsv > replicon_all_metadata_tax.tsv
```
The resulting variables were:
- `chr1_gc`: GC content of the largest replicon in the genome.
- `delta_gc`: difference between the GC content of a replicon and the reference replicon.
- `abs_delta_gc`: absolute value of `delta_gc`.
These features were retained as replicon-level evidence for downstream reclassification.

# 3. GTDB archaeal marker counts
GTDB-Tk release 232 marker HMMs were combined to create an archaeal marker database containing the Pfam and TIGRFAM markers used by the GTDB archaeal marker set.
```bash
    cat gtdbtk_r232_data/markers/pfam/individual_hmms/*.hmm \
        gtdbtk_r232_data/markers/tigrfam/individual_hmms/*.HMM \
        > archaeal_ar122.hmm

    hmmpress archaeal_ar122.hmm
```
The predicted proteins of each replicon were scanned against this HMM database:
```bash
    hmmscan --tblout output.tbl archaeal_ar122.hmm replicon.faa
```
For each replicon, the number of detected GTDB marker families was summarized as `gtdb_markers`.

# 4. Tetranucleotide frequency similarity
Tetranucleotide frequency (TNF) profiles were calculated for the replicons and compared within each genome.
The TNF analysis produced correlations between replicons and the reference replicon (chr1) used for each genome. The resulting correlation value was incorporated into the replicon feature table as `tnf_corr_best_gtdb`.
The calculation was performed using the `get_tnf.py` script.
The resulting TNF correlations were stored in `TNF_correlations.csv` and subsequently joined to the replicon metadata.

# 5. Replicon marker sets
Three sets of protein markers were selected to provide evidence for different replicon functions.

# Chromosome-associated markers
`Cdc6_lid`, `WHD_Cdc6`, `MCM`, `GINS23_N`, `RNA_pol_Rpb2_1`, `RFC1`, `Gins15_N`

# Partition-system markers
`ParA`, `ParB_N`, `ParB`, `HTH_ParB`, `ParB_C`, `CB_ParB_C`

# Plasmid-associated markers
`RepA_N`, `RepA_C`, `RepB`, `RepB_C`, `RepC`, `Replitron_HUH`

The corresponding Pfam HMMs were retrieved from the Pfam-A database:
```bash
    wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm
```
The individual Pfam HMMs were extracted and combined into three marker databases:
- `archaeal_chr1_markers.hmm`
- `archaeal_partition_markers.hmm`
- `archaeal_plasmid_markers.hmm`

Each database was indexed with `hmmpress`.

# 6. HMM scanning of replicon proteins
Predicted proteins from each replicon were scanned independently against the three marker databases.
# Chromosome markers
```bash
    hmmscan --tblout chr1_scan/${base}.chr1.tbl \
        archaeal_chr1_markers.hmm "$faa" > /dev/null
```
# Partition markers
```bash
    hmmscan --tblout partition_scan/${base}.partition.tbl \
        archaeal_partition_markers.hmm "$faa" > /dev/null
```
# Plasmid markers
```bash
    hmmscan --tblout plasmid_scan/${base}.plasmid.tbl \
        archaeal_plasmid_markers.hmm "$faa" > /dev/null
```
For each replicon, the scans were summarized as:
- **Presence:** number of distinct marker families detected.
- **Hits:** number of protein hits assigned to the marker set.

This produced: `chr1_presence`, `chr1_hits`, `partition_presence`, `partition_hits`, `plasmid_presence`, `plasmid_hits`

# 7. Individual marker presence and hit counts
In addition to the summary counts for each marker set, individual marker families were retained.
The three scan results were combined:
```bash
    markers <- bind_rows(
        read_scan("chr1_scan", "CHR_"),
        read_scan("partition_scan", "PART_"),
        read_scan("plasmid_scan", "PLAS_"))
```
Marker presence was calculated from distinct genome-contig-marker combinations:
```bash
    marker_presence <- markers %>%
        distinct(genome, contig, marker) %>%
        mutate(value = 1L) %>%
        pivot_wider(names_from = marker,
            values_from = value,
            values_fill = 0)
```
Marker hit counts were calculated separately:
```bash
    marker_hits <- markers %>%
        count(genome, contig, marker, name = "hits") %>%
        pivot_wider(names_from = marker,
            values_from = hits,
            values_fill = 0)
```
The resulting individual-marker columns included: CHR_Cdc6_lid, CHR_WHD_Cdc6, CHR_MCM, CHR_GINS23_N, CHR_RNA_pol_Rpb2_1, CHR_RFC1, CHR_Gins15_N,
PART_ParA, PART_ParB_N, PART_ParB, PART_HTH_ParB, PART_ParB_C, PART_CB_ParB_C, PLAS_RepA_N, PLAS_RepA_C, PLAS_RepB, PLAS_RepB_C`, PLAS_RepC, PLAS_Replitron_HUH.
Both marker presence and hit counts were retained.

# 8. Integration of replicon features
The sequence statistics, GTDB marker counts, TNF similarity, marker-set summaries, and individual marker features were integrated into a single replicon-level feature table.
Missing marker values were interpreted as zero when no corresponding HMM hit was detected.

The resulting table was:`replicon_features_markers.tsv`
Its main columns include: genome, contig, replicon, length, gc_content, gtdb_markers, tnf_corr_best_gtdb, chr1_presence, chr1_hits, partition_presence, partition_hits, plasmid_presence, and plasmid_hits followed by the individual marker presence and hit-count columns.

This table provides the marker- and sequence-based evidence used for replicon reclassification. The `replicon` column at this stage retains the original replicon annotation.

# 9. Final integrated replicon metadata
The replicon feature information was combined with the genome-level taxonomic metadata and the within-genome GC-content features.
The resulting table was: `replicon_all_metadata_tax.tsv`
It contains the replicon-level sequence and marker features together with: phylum, class, order, family, genus, species, chr1_gc, delta_gc, abs_delta_gc
Thus, `replicon_all_metadata_tax.tsv` provides the combined replicon-level dataset used for downstream analyses, including:

- sequence length and GC content,
- GTDB marker counts,
- TNF similarity,
- chromosome-associated markers,
- partition-system markers,
- plasmid-associated markers,
- taxonomic assignment,
- GC-content divergence from the largest replicon within the genome.

# Associated scripts and datasets
`get_tnf.py`                     # Calculate TNF
`replicon_features_markers.tsv`  # Replicon-level sequence, GTDB, TNF, marker-set, and individual HMM marker features |
`replicon_all_metadata_tax.tsv`  # Integrated replicon metadata including marker features, taxonomy, and within-genome GC-content differences |
`TNF_correlations.csv`           # TNF similarity between replicons and the reference replicon |
`archaeal_ar122.hmm`             # Combined GTDB archaeal marker HMM database |
`archaeal_chr1_markers.hmm`      # Chromosome-associated marker HMMs |
`archaeal_partition_markers.hmm` # Partition-system marker HMMs |
`archaeal_plasmid_markers.hmm`   # Plasmid-associated marker HMMs |