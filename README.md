# Scripts for analysis of Fluent PIPseq and 10x 3' and 10x 5' libraries

This repository serves as the hub for the comparative analysis of Fluent PIPseq and 10x 3' and 5' scRNAseq kits.

 * `notebooks/` contains the Jupyter notebooks needed to perform all analysis.
 * `environment.yaml` has the requirements for the notebooks.
 
To reproduce analysis, clone this repo on to a large machine. We typically used 8 or 16 CPUs with a lot of memory (64GB or more), depending on the processing step. For the initial short-read processing a larger machine might be necessary.

To create the environment, run `conda env create -f environment.yaml`.

Once you have the environment you'll want to register the kernel with `jupyterlab`. Command for this is `ipython kernel install --user --name="mdl-sc-isoform-2025-ms"`, from within the activated environment.

The `conda` package for `graphviz` isn't sufficient for graph layout, you will need to install the system package as well:

```bash
sudo apt install graphviz
```

If you want to reproduce figure panels exactly, you should install the `Open Sans` font from [Google Fonts](https://fonts.google.com/specimen/Open+Sans).

## Code

Besides the notebooks, the analysis relies on the following MDL packages:

 * [`callao`](https://github.com/MethodsDev/callao) (and [`mdl-core`](https://github.com/MethodsDev/mdl-core)) for demultiplexing long reads with artifacts
 * [`bouncer`](https://github.com/MethodsDev/bouncer) and [`barcode-symspell`](https://github.com/MethodsDev/barcode-symspell) for barcode and UMI extraction
 * [`isoscelles`](https://github.com/MethodsDev/isoscelles) for single-cell long-read analysis
 * [`bugzapper`](https://github.com/MethodsDev/bugzapper) is a utility for efficiently subsampling a short-read BAM
 * A package in this repo, `mdl.sc_isoform_paper`, which contains code for this specific analysis.

`callao`, `bouncer`, and `bugzapper` all contain extensions written in Rust, so they need a compiler toolchain to install (at least, until we build a wheel):

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh  # copied from rustup.rs
sudo apt install build-essential
```

## Other tools

To avoid dependency conflicts, I have been installing `minimap2` and `isoquant` in their own environments. `conda` can install all of these from the `bioconda` channel.

 * samtools
 * gffutils
 * minimap2
 * isoquant

The analysis also relies on [`marti`](https://github.com/PopicLab/marti), see the README there for installation details.

## Analysis

There are a few separate types of analysis here.

 * Short read data was processed with `CellRanger` (for 10x 3' and 5') and `PIPseeker` (for Fluent PIPseq)
   * Short-read analysis of PBMC data
   * Short-read barnyard comparison, with ambient and doublet rates
 * Long-read data was tagged with `lima` to find adapters, demultiplexed with `callao`, deconcatenated with `skera` (if MAS-seq), and then run through `marti` for artifact classification. Reads were tagged with `bouncer` and the oriented cDNA was extracted and mapped with `minimap2`. Isoforms were quantified with `IsoQuant`
   * Long-read monomer: focused on artifact comparisons as well as polyA and internal priming
   * Long-read MAS-seq: focused on differences in isoform identification, the how mispriming can lead to unreliable novel IDs

## Data

Data is in Google Cloud Storage: [`gs://mdl-sc-isoform-2025-ms/`](https://console.cloud.google.com/storage/browser/mdl-sc-isoform-2025-ms/)

The analysis was performed on the following datasets:

 * __Reference__ - we used the `GRCh38_no_alt` and Gencode `v39` references.
 * __Metadata__ - the relevant metadata is the barcode lists for 10x and PIPseq, MAS adapters, and our demux adapters (including specialized mas16 adapter files). They are uploaded to `gs://mdl-sc-isoform-2025-ms/metadata.tgz`
 * __Short-read__ - this was run on an Illumina NovaSeq X, `230924_SL-EXC_0072_A2275VKLT3` 
   * `fastq.gz` files: [`gs://mdl-sc-isoform-2025-ms/sequencing_data/pbmc_illumina`](https://console.cloud.google.com/storage/browser/mdl-sc-isoform-2025-ms/sequencing_data/pbmc_illumina) and [`gs://mdl-sc-isoform-2025-ms/sequencing_data/barnyard_illumina`](https://console.cloud.google.com/storage/browser/mdl-sc-isoform-2025-ms/sequencing_data/barnyard_illumina)
   * Contains four samples
     * 10x 3'
     * 10x 5'
     * Fluent PBMC (0.8x SPRI)
     * Fluent Barnyard (K562 + 3T3)
 * __Monomer PBMC__ - run on two Sequel IIe SMRTcells, `m64455e_230922_195123` and `m64455e_230924_064453`
   * Raw data: [`gs://mdl-sc-isoform-2025-ms/sequencing_data/pbmc_monomer`](https://console.cloud.google.com/storage/browser/mdl-sc-isoform-2025-ms/sequencing_data/pbmc_monomer)
   * Contains six samples
     * Fluent PIPseq 0.8x SPRI
     * Fluent PIPseq 0.6x SPRI
     * 10x 3' w/ purification
     * 10x 5' w/ purification
     * 10x 3' w/o purification
     * 10x 5' w/o purification
 * __MAS-seq PBMC__ - run across eight Revio SMRTcells
   * Raw data: [`gs://mdl-sc-isoform-2025-ms/sequencing_data/pbmc_masseq`](https://console.cloud.google.com/storage/browser/mdl-sc-isoform-2025-ms/sequencing_data/pbmc_masseq)
   * Contains four samples
     * Fluent PIPseq 0.8x SPRI
     * Fluent PIPseq 0.6x SPRI
     * 10x 3'
     * 10x 5'
  * __Fix/Fresh MAS-seq__ - A separate lot of PBMCs to reproduce the Methanol+DSP fixation results. We performed two experiments. The "fixation panel" experiment compared five different protocols alongside fresh cells. The "fix-fresh" experiment was a replicate of the Methanol+DSP vs fresh comparison, to achieve higher cells numbers for those conditions. In both experiments the different cells were hashtagged and captured together via PIPseq, then sequenced on both short- and long-read sequencers.
    * __Fix/Fresh PIPseq short-read data__ 
      * `fastq.gz` files: [`gs://mdl-sc-isoform-2025-ms/sequencing_data/fix-fresh_illumina`](https://console.cloud.google.com/storage/browser/mdl-sc-isoform-2025-ms/sequencing_data/fix-fresh_illumina)
      * PIPseeker output for fixation panel: [`gs://mdl-sc-isoform-2025-ms/fixation-panel_pipseq`](https://console.cloud.google.com/storage/browser/mdl-sc-isoform-2025-ms/fixation-panel_pipseq)
      * PIPseeker output for fix-fresh: [`gs://mdl-sc-isoform-2025-ms/fix-fresh_pipseq`](https://console.cloud.google.com/storage/browser/mdl-sc-isoform-2025-ms/fix-fresh_pipseq)
    * __Fix/Fresh PIPseq MAS-seq__
      * Raw data: [`gs://mdl-sc-isoform-2025-ms/sequencing_data/fix-fresh_masseq`](https://console.cloud.google.com/storage/browser/mdl-sc-isoform-2025-ms/sequencing_data/fix-fresh_illumina)
        * Fixation panel data: `m84175_240224_082414_s4.hifi_reads.bcM004.bam`
        * Fix-fresh data: `m84250_240530_044215_s1.hifi_reads.bcM0001.bam` and `m84250_240530_064125_s2.hifi_reads.bcM0001.bam`


The second set of SMRTcells for the MAS-seq PBMC were run at different concentrations as an optimization experiment, but there was no particular relationship between concentration and output, and the data were incorporated normally.

| flowcell | concentration | # reads |
|-|-|-|
| `m84143_230929_195310_s1` | 350 pM | 3,122,400 |
| `m84143_230929_202333_s2` | 400 pM | 4,098,164 |
| `m84143_230929_205439_s3` | 450 pM | 3,361,489 |
| `m84143_230929_212545_s4` | 500 pM | 3,036,894 |
