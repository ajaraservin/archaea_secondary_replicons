# Building metadata and cluster tables
# Overview
This workflow was used to integrate replicon, protein, genome, and protein-cluster information into a set of tables for downstream analyses of multipartite archaeal genomes.
The workflow consisted of the following steps:
1. Build replicon metadata from NCBI replicon annotations.
2. Map predicted proteins to their genomic contigs.
3. Assign replicon classifications to proteins.
4. Create a combined protein FASTA file.
5. Cluster proteins using MMseqs2.
6. Build an integrated cluster provenance table.
7. Generate cluster abundance tables by replicon and genome.
8. Generate a cluster-by-replicon matrix.
9. Generate representative protein IDs for downstream annotation.

# 1. Build replicon metadata
A replicon metadata table was generated from the NCBI `*_genomic_replicons.tsv` files created with the NCBI replicon annotation.
The table was generated with:
```bash    
    (echo -e "accession\tcontig\treplicon_type"
        find ./arch_genomes/ncbi_dataset/data/GC*/*_genomic_replicons.tsv | while read f; do
            accession=$(basename "$f" | sed 's/_ASM.*//')
            awk -v acc="$accession" 'BEGIN{OFS="\t"}{
                print acc, $1, $2
            }' "$f"
        done
    ) > replicon_metadata.tsv
```
The resulting table contains one row per replicon indicating:
- assembly accession
- contig accession
- replicon type

# 2. Build the protein-to-contig mapping
Protein-to-contig relationships were obtained from the Prokka GFF files generated for the curated genomes.
For each CDS, the protein identifier and genomic contig were extracted from the GFF attributes.
```bash    
    (echo -e "protein\taccession\tcontig\treplicon_type"
        find ./arch_genomes/ncbi_dataset/data -path "*/prokka_*/*.gff" | while read gff; do
            accession=$(basename "$gff" .gff)
            awk -F'\t' -v acc="$accession" '
            BEGIN{OFS="\t"}
            /^#/ {next}
            $3=="CDS" {
                contig=$1
                if (match($9, /ID=([^;]+)/, a)) {
                    protein=a[1]
                    print protein, acc, contig
                }
            }' "$gff"
        done
    ) > protein_contig.tmp
```
# 3. Add replicon classification to proteins
The protein-to-contig mapping was joined with the replicon metadata using the genome accession and contig accession.
Proteins whose genome–contig combination could not be matched to the replicon metadata were assigned the category `unknown`.
```bash 
    awk '
    BEGIN{
        FS=OFS="\t"
    }
    NR==FNR {
        if (FNR==1) next
        key=$1"\t"$2
        repl[key]=$3
        next
    }
    FNR==1 {next}
    {
        key=$2"\t"$3
        type=(key in repl ? repl[key] : "unknown")
        print $1, $2, $3, type
    }' replicon_metadata.tsv protein_contig.tmp > protein_metadata.tsv
```
The resulting `protein_metadata.tsv` table contains:
    protein
    accession
    contig
    replicon_type

This table provides the link between individual proteins and their genomic replicons.

# 4. Create the combined protein FASTA
Protein FASTA files from all curated archaeal genomes were concatenated into a single FASTA file for clustering.
```bash
    cat prokka_*/*.faa > arch_all_proteins.faa
```
# 5. Cluster proteins using MMseqs2
The combined protein FASTA was clustered using MMseqs2.
```bash
    mmseqs createdb arch_all_proteins.faa db
    mmseqs cluster db clusters tmp \
        --min-seq-id 0.5 \
        -c 0.8 \
        --cov-mode 1

    mmseqs createtsv db db clusters arch_clusters.tsv
```
The clustering parameters were: minimum sequence identity: 50%, minimum coverage: 80%, coverage mode: 1
The resulting `arch_clusters.tsv` links a protein to its representative protein cluster. It has two columns:
representative_protein and member_protein

# 6. Build the master cluster table
The protein-cluster assignments were combined with the protein metadata to generate an integrated provenance table.
```bash 
    awk 'BEGIN{FS=OFS="\t"}
    NR==FNR{
        meta[$1]=$2"\t"$3"\t"$4
        next
    }
    {
        cluster=$1
        protein=$2
        if(protein in meta){
            print cluster, protein, meta[protein]
        }
    }' protein_metadata.tsv arch_clusters.tsv > arch_cluster_master.tsv
```
The resulting table contains: cluster, protein, genome, contig, replicon
`arch_cluster_master.tsv` provides the relationship between protein clusters and their genomic context.

# 7. Generate cluster × replicon counts
The number of protein copies belonging to each cluster was calculated for each replicon category.
```bash
    awk 'BEGIN{FS=OFS="\t"}
    {
        key=$1"\t"$5
        count[key]++
    }
    END{
        print "cluster","replicon","count"
        for(k in count){
            split(k,a,"\t")
            print a[1], a[2], count[k]
        }
    }' arch_cluster_master.tsv > arch_cluster_replicon_counts.tsv
```
The resulting table contains: cluster, replicon, count

# 8. Generate the cluster × replicon matrix
The cluster × replicon count table was converted into a matrix with one row per protein cluster and separate columns for each replicon category.
```bash    
    awk 'BEGIN{
        FS=OFS="\t"
    }
    NR>1{
        cluster=$1
        repl=$2
        count=$3
        data[cluster,repl]=count
        seen[cluster]=1
    }
    END{
        print "cluster","chr1","chr2plus","plasmid","unknown"
        for(c in seen){
            chr1=((c SUBSEP "chr1") in data ? data[c,"chr1"] : 0)
            chr2=((c SUBSEP "chr2plus") in data ? data[c,"chr2plus"] : 0)
            plasmid=((c SUBSEP "plasmid") in data ? data[c,"plasmid"] : 0)
            unknown=((c SUBSEP "unknown") in data ? data[c,"unknown"] : 0)
            print c, chr1, chr2, plasmid, unknown
        }
    }' arch_cluster_replicon_counts.tsv > arch_cluster_replicon_matrix.tsv
```
The resulting matrix contains: cluster, chr1, chr2plus, plasmid

# 9. Generate cluster × genome × contig counts
A table preserving the complete genomic context of each cluster was generated by counting cluster copies within individual genome-contig combinations.
```bash
    awk 'BEGIN{FS=OFS="\t"}
    NR>1{
        key=$1"\t"$3"\t"$4"\t"$5
        count[key]++
    }
    END{
        print "cluster","genome","contig","replicon","copies"
        for(k in count){
            split(k,a,"\t")
            print a[1],a[2],a[3],a[4],count[k]
        }
    }' arch_cluster_master.tsv > arch_cluster_genome_contig_counts.tsv
```
# 10. Generate cluster × genome × replicon counts
A copy-number table was generated by aggregating cluster occurrences at the genome and replicon level.
```bash
    awk 'BEGIN{FS=OFS="\t"}
    NR>1{
        key=$1"\t"$3"\t"$5
        count[key]++
    }
    END{
        print "cluster","genome","replicon","copies"
        for(k in count){
            split(k,a,"\t")
            print a[1],a[2],a[3],count[k]
        }
    }' arch_cluster_master.tsv > arch_cluster_genome_replicon_counts.tsv
```
This table represents the number of copies of each protein cluster within each genome and replicon category.

# 11. Generate representative protein IDs
The representative proteins from the MMseqs2 clustering results were used to generate the list of proteins for subsequent functional and taxonomic annotation.
```bash
    cut -f1 arch_clusters.tsv | sort -u > representative_ids.txt
```
The resulting file contains one representative protein identifier per cluster.

# Associated scripts and datasets
`replicon_metadata.tsv`                     # Replicon-level metadata linking assemblies to contigs and replicon classifications. |
`protein_metadata.tsv`                      # Protein-to-genome, contig, and replicon mapping. |
`arch_all_proteins.faa`                     # Combined protein FASTA from all curated archaeal genomes. |
`arch_clusters.tsv`                         # Protein-to-cluster assignments generated with MMseqs2. |
`arch_cluster_master.tsv`                   # Integrated cluster, protein, genome, contig, and replicon provenance table. |
`arch_cluster_replicon_counts.tsv`          # Cluster abundance across replicon categories. |
`arch_cluster_replicon_matrix.tsv`          # Cluster-by-replicon abundance matrix. |
`arch_cluster_genome_contig_counts.tsv`     # Cluster copy numbers by genome, contig, and replicon. |
`arch_cluster_genome_replicon_counts.tsv`   # Cluster copy numbers by genome and replicon category. |
`representative_ids.txt`                    # Representative protein IDs used for downstream annotation. |