# Taxonomic and functional annotation
# Overview
The taxonomic and functional annotation workflow consisted of:

1. Assigning taxonomy to the archaeal genomes using GTDB-Tk.
2. Simplifying the GTDB-Tk taxonomy output to the main taxonomic ranks.
3. Integrating genome taxonomy with replicon-level metadata.
4. Calculating GC-content differences between replicons using the largest replicon as the reference.
5. Functionally annotating representative proteins using eggNOG-mapper.
6. Linking functional annotations to protein clusters.
7. Generating individual FASTA files for each genome–replicon combination for downstream analyses.

The representative proteins used for annotation were the 142,855 cluster representatives generated during the clustering workflow described in `02_replicon_and_cluster_tables.md`.

# 1. Genome taxonomy with GTDB-Tk
# 1.1. Generate the genome list and split it into batches
The curated archaeal genomic FASTA files were listed and sorted:
```bash    
    find arch_genomes/ncbi_dataset/data/G* -name "*_genomic.fna" | sort > genome_paths.txt
    split -l 100 genome_paths.txt gtdb_chunk_
    ls gtdb_chunk_* | wc -l
```
The resulting list contained 243 archaeal genomes that were divided in three batches.

# 1.2. Submit GTDB-Tk jobs and merge results
Each batch was submitted independently using the GTDB-Tk execution script:
```bash
    for f in gtdb_chunk_*; do
        sbatch run_gtdbtk.sh $f
    done
```
The taxonomy summary files from the three batches were combined into a single table:
```bash   
    {
        head -n1 gtdb_chunk_aa/output/gtdbtk.ar53.summary.tsv
        tail -n+2 gtdb_chunk_aa/output/gtdbtk.ar53.summary.tsv
        tail -n+2 gtdb_chunk_ab/output/gtdbtk.ar53.summary.tsv
        tail -n+2 gtdb_chunk_ac/output/gtdbtk.ar53.summary.tsv
    } > genome_gtdb_taxonomy.tsv
```
The resulting table contains the GTDB-Tk taxonomic assignment for the 243 archaeal genomes.

# 2. Simplify the GTDB-Tk taxonomy
The GTDB-Tk taxonomy strings were parsed to retain the following taxonomic ranks: Phylum,
Class, Order, Family, Genus, Species
The genome accession was also standardized by removing the assembly suffix:
```bash     
    awk 'BEGIN{FS=OFS="\t"}
    NR==1{
        print "genome","phylum","class","order","family","genus","species"
        next
    }
    {
        genome=$1
        sub("_ASM.*","",genome)
        split($2,tax,";")
        for(i=1;i<=7;i++){
            sub("^[a-z]__","",tax[i])
        }
        print genome,tax[2],tax[3],tax[4],tax[5],tax[6],tax[7]
    }' genome_gtdb_taxonomy.tsv > genome_taxonomy.tsv
```
The resulting `genome_taxonomy.tsv` contains:
    genome    phylum    class    order    family    genus    species

# 3. Integrate taxonomy with replicon metadata
The genome-level taxonomy was integrated with the replicon-level metadata generated during the previous workflow.
Missing taxonomic assignments were represented as `NA`:
```bash
    awk 'BEGIN{FS=OFS="\t"}
    FNR==NR{
        if(FNR>1)
            tax[$1]=$2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7
        next
    }
    FNR==1{
        print $0,"phylum","class","order","family","genus","species"
        next
    }
    {
        if($1 in tax)
            print $0,tax[$1]
        else
            print $0,"NA","NA","NA","NA","NA","NA"
    }' genome_taxonomy_fixed.tsv replicon_features_full.tsv > replicon_metadata_tax.tsv
```
The resulting `replicon_metadata_tax.tsv` combines replicon-level genomic information with the GTDB-Tk taxonomy of the corresponding genome.

# 4. Functional annotation with eggNOG-mapper
The representative protein sequences were functionally annotated using eggNOG-mapper. The annotation workflow was submitted through the `arch_annot_eggnog.sh` script:
The resulting files were stored in `emapper_out/`:
    arch_emapper.emapper.annotations
    arch_emapper.emapper.decorated.gff
    arch_emapper.emapper.orthologs
    arch_emapper.emapper.annotations.xlsx
    arch_emapper.emapper.hits
    arch_emapper.emapper.seed_orthologs

# 5. Extract COG categories and functional descriptions
The eggNOG-mapper annotation table was reduced to the protein identifier, COG category, and functional description:
```bash
    (echo -e "protein\tCOG_category\tDescription"
        grep -v '^#' emapper_out/arch_emapper.emapper.annotations | cut -f1,7,8
    ) > arch_protein_annotation.tsv
```
The resulting `arch_protein_annotation.tsv` provides the functional annotation associated with each representative protein.

# 6. Link functional annotations to protein clusters
A table containing one representative protein per cluster was generated from the cluster master table:
```bash
    awk '!seen[$1]++{print $1"\t"$2}' arch_cluster_master.tsv > arch_cluster_ids.tsv
```
The representative protein identifiers were then matched to their eggNOG-mapper annotations:
```bash
    (echo -e "cluster\trepresentative\tCOG_category\tDescription"
        awk 'BEGIN{FS=OFS="\t"}
        NR==FNR{
            if(FNR==1) next
            cog[$1]=$2
            desc[$1]=$3
            next
        }
        {print $1,$2,cog[$2],desc[$2]}' arch_protein_annotation.tsv arch_cluster_ids.tsv
    ) > arch_cluster_annotation.tsv
```
The resulting `arch_cluster_annotation.tsv` links each protein cluster to:
- Cluster identifier
- Representative protein
- COG category
- Functional description
This table provides the functional annotation of the protein clusters generated during the clustering workflow.

# 7. Functional annotation with FANTASIA
FANTASIA was used to assign Gene Ontology (GO) terms to proteins predicted from the curated archaeal genomes. Protein-coding sequences were first predicted with Prodigal and subsequently used as input for FANTASIA.

# 8.1. Protein prediction
The 243 curated archaeal genomes were submitted to protein prediction using the script `run_prodigal.sh`
The resulting protein FASTA files were stored in `prodigal_out/`. In total, 7,381,413 proteins were predicted.

# 8.2. FANTASIA annotation
Prodigal *.faa outputs were processed using the FANTASIA execution and submission scripts: `submit_fantasia.sh`and `run_fantasia.sh`
To ensure that all genomes had completed successfully, the output directories were compared against the available Prodigal protein files. Missing genomes were identified with `find_missing_fantasia.sh` and reprocessed.
The expected output for each genome was the corresponding:
    fantasia_out/<genome>/<genome>_topgo.txt

# 8.3. Processing GO annotations
The GO annotations from all FANTASIA output files were combined in R. Each comma-separated GO annotation was expanded into an individual row, and the genome identifier was extracted from the output directory name.
```R
    library(tidyverse)
    library(GO.db)
    library(AnnotationDbi)

    files <- list.files(
        "fantasia_out",
        pattern = "_topgo\\.txt$",
        recursive = TRUE,
        full.names = TRUE
    )

    go <- purrr::map_dfr(files, function(f){
        genome <- basename(dirname(f))
        genome <- sub("_ASM.*", "", genome)

        read.delim(
            f,
            header = FALSE,
            sep = "\t",
            stringsAsFactors = FALSE,
            col.names = c("protein", "GO")
        ) %>%
            tidyr::separate_rows(GO, sep = ",\\s*") %>%
            dplyr::mutate(genome = genome)
    })
```
GO identifiers were assigned their ontology and term names using `GO.db`:
```R
    go$ontology <- Ontology(go$GO)
    go$name <- Term(go$GO)

    go <- go %>%
        filter(!is.na(ontology))
```
This produced a table containing: protein, GO, genome, ontology, name
The three GO ontologies were Biological Process (`BP`), Cellular Component (`CC`), and Molecular Function (`MF`).

# 8.4. GO summary tables
The processed FANTASIA annotations were summarized at the protein, genome, and ontology levels. The main output tables were:
- **Complete protein-level GO annotations**
```R
      write.table(go, "GO_analysis/fantasia_GO_annotations.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```
- **GO term frequencies by ontology**
```R      
      go_counts <- go %>% count(ontology, GO, name, sort = TRUE)
      write.table(go_counts, "GO_analysis/GO_term_counts.tsv",sep = "\t", quote = FALSE, row.names = FALSE)
```
- **Genome × GO counts**
```R
      genome_go <- go %>% count(genome, GO)
      write.table(genome_go, "GO_analysis/genome_GO_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```
- **Genome × GO summary including ontology and term name**
```R
      go_summary <- go %>% count(genome, ontology, GO, name)
      write.table(go_summary, "GO_analysis/genome_GO_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```
- **GO richness per genome**
```R
      go_richness <- go %>% distinct(genome, GO) %>% count(genome, name = "GO_richness")
      write.table(go_richness, "GO_analysis/go_richness.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```
- **Genome × GO presence/absence**
```R
      genome_GO_presence <- go %>% distinct(genome, GO, name, ontology)
      write.table(genome_GO_presence, "GO_analysis/genome_GO_presence.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```
- **Genome taxonomy**
```R
      metadata <- read.delim("replicon_metadata_tax.tsv") %>% distinct(genome, phylum, class, order, family, genus, species)
      write.table(metadata, "GO_analysis/genome_taxonomy.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```
- **GO ontology summary**
```R
      ontology_summary <- go %>% count(ontology)
      write.table(ontology_summary, "GO_analysis/ontology_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```
- **Top 100 GO terms per ontology**
```R
      top_GO <- go %>% group_by(ontology, GO, name) %>% summarise(n = n(), .groups = "drop") %>% group_by(ontology) %>% slice_max(n, n = 100) %>% ungroup()
      write.table(top_GO, "GO_analysis/top100_GO_per_ontology.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```
A dataset summary was generated containing the number of GO annotations, unique proteins, unique GO terms, annotated genomes, and annotations belonging to each ontology.
The final FANTASIA dataset contained 4,364,182 GO annotations from 825,277 proteins, representing 7,900 unique GO terms across 228 genomes.

# 9. Generate individual replicon FASTA files
The curated genomic FASTA files were split into individual replicon FASTA files, generating one FASTA file for each genome–contig combination,
using the split_replicons.sh script.
The resulting files were stored in the `replicons/` directory using the genome and contig accessions in the filename.
Each file contains the sequence of an individual replicon from a specific archaeal genome.
The resulting collection contained 833 genome–replicon FASTA files.

# Associated scripts and datasets
`run_gtdbtk.sh`                 # Assign taxonomy using GTDB-Tk database
`arch_annot_eggnog.sh`          # Assign functional annotation using eggNOG-mapper
`genome_gtdb_taxonomy.tsv`      # Complete GTDB-Tk taxonomy output for the archaeal genomes.
`genome_taxonomy.tsv`           # Simplified genome-level taxonomy containing phylum through species.
`replicon_metadata_tax.tsv`     # Replicon metadata integrated with genome-level taxonomy.
`replicon_all_metadata_tax.tsv` # Replicon metadata with GC-content reference, GC difference, and absolute GC difference.
`arch_protein_annotation.tsv`   # COG categories and functional descriptions assigned to representative proteins. 
`arch_cluster_ids.tsv`          # Cluster identifiers and their representative proteins.
`arch_cluster_annotation.tsv`   # Protein cluster annotations linked to representative proteins and eggNOG-mapper functional assignments.
`replicons/GCF_*__NC_*.fna`     # Individual FASTA files containing one replicon per genome–contig combination.
`fantasia_GO_annotations.tsv`   # Complete protein-level GO annotations.

### FANTASIA outputs
| `GO_analysis/fantasia_GO_annotations.tsv` | Complete protein-level GO annotations. |
| `GO_analysis/GO_term_counts.tsv` | GO term frequencies by ontology. |
| `GO_analysis/genome_GO_counts.tsv` | GO counts per genome. |
| `GO_analysis/genome_GO_summary.tsv` | Genome-level GO annotations with ontology and term names. |
| `GO_analysis/genome_GO_presence.tsv` | Non-redundant genome × GO presence table. |
| `GO_analysis/go_richness.tsv` | Number of unique GO terms per genome. |
| `GO_analysis/genome_taxonomy.tsv` | Taxonomic information associated with each genome. |
| `GO_analysis/ontology_summary.tsv` | Number of annotations per GO ontology. |
| `GO_analysis/top100_GO_per_ontology.tsv` | Most frequent GO terms for each ontology. |
| `GO_analysis/dataset_summary.tsv` | Summary statistics for the FANTASIA annotation dataset. |