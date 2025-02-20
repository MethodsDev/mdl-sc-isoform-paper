## Short-read processing

__Short-read__ - run on an Illumina NovaSeq X, `230924_SL-EXC_0072_A2275VKLT3` 
   * `fastq.gz` files: [`gs://mdl-sc-isoform-2025-ms/sequencing_data/pbmc_illumina`](https://console.cloud.google.com/storage/browser/mdl-sc-isoform-2025-ms/sequencing_data/pbmc_illumina) and [`gs://mdl-sc-isoform-2025-ms/sequencing_data/pbmc_illumina`](https://console.cloud.google.com/storage/browser/mdl-sc-isoform-2025-ms/sequencing_data/barnyard_illumina)
   * Contains four samples
     * 10x 3'
     * 10x 5'
     * Fluent PBMC (0.8x SPRI)
     * Fluent Barnyard (K562 + 3T3)

## PBMC processing

#### Making the references

We used the `GRCh38_no_alt` genome reference from [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz), and the [Gencode v39 annotation](https://www.gencodegenes.org/human/release_39.html).

PIPseeker uses a normal STAR index. The sparsity parameter makes it a little smaller.

```bash
STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles GRCh38.fasta --genomeSAsparseD 3 --sjdbGTFfile gencode.v39.annotation.gtf
```

CellRanger has its own format for the reference but in its core it's still a STAR index, so we can reuse the same one from before, or make a fresh one:

```bash
cellranger-7.1.0/cellranger mkref --genome=GRCh38_10x --fasta=GRCh38.fasta --genes=gencode.v39.annotation.gtf
```

### CellRanger

Run CellRanger 7.1.0 with the GRCh38 genome (can tweak additional options for performance): 

```bash
cellranger-7.1.0/cellranger count --id MDL_10x_PBMC_3p --fastqs path/to/pbmc_illumina --sample MDL_10x_PBMC_3p --transcriptome path/to/reference/GRCh38_10x --disable-ui
cellranger-7.1.0/cellranger count --id MDL_10x_PBMC_5p --fastqs path/to/pbmc_illumina --sample MDL_10x_PBMC_5p --transcriptome path/to/reference/GRCh38_10x --disable-ui
```

### PIPseeker

Be warned: this can take a long time and needs a lot of space.
```bash
pipseeker-v2.1.4-linux/pipseeker full \
  --output-path path/to/output/MDL_Fluent_PBMC \
  --fastq path/to/pbmc_illumina/MDL_PIPseq_PBMC \
  --star-index-path path/to/reference/GRCh38
```

## Barnyard

### Reference

We can get the relevant fasta file from the [10x reference](https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-and-mm10-2020-A.tar.gz) and build our own STAR index for PIPseeker using the same fasta and GTF files.

```bash
STAR --runMode genomeGenerate \
  --genomeFastaFiles path/to/refdata-gex-GRCh38-and-mm10-2020-A/fasta/genome.fa \
  --genomeSAsparseD 3 \
  --sjdbGTFfile path/to/refdata-gex-GRCh38-and-mm10-2020-A/genes/genes.gtf
```

### PIPseeker

```bash
pipseeker-v2.1.4-linux/pipseeker full \
  --output-path path/to/output/MDL_Fluent_PBMC \
  --fastq /path/to/pbmc_illumina/MDL_PIPseq_Barnyard \
  --star-index-path path/to/reference/GRCh38-and-mm10
```
