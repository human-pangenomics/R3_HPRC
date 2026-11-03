# HPRC Assembly QC Pipeline

Practical guide for running the QC workflow applied to HPRC R3 hybrid assemblies and panpatch pilot assemblies. All steps are also available as WDL workflows runnable via Cromwell or Toil — this document covers the underlying commands and expected outputs.

**Assembly type:** Diploid, per-haplotype FASTAs (Verkko hi-c, ref-oriented)  
**Read types:** HiFi and ONT R10  
**Reference:** CHM13v2.0

---

## Pipeline Overview

```
Assembly FASTAs (hap1 + hap2)
    ├── 1. CenSat          centromere + satellite annotation
    ├── 2. Compleasm        BUSCO gene completeness
    ├── 3. Quaak            k-mer QV + CNV
    ├── 4. T2T QC           assembly break annotation vs CHM13
    └── 5. Flagger          read-to-assembly alignment + HMM flagging
              └── 6. NucFlag    secondary-base misassembly detection
```

Steps 1–4 are independent. Flagger (step 5) must finish before NucFlag (step 6), which uses the Flagger output BAM.

---

## Reference Annotation Files

All CHM13v2.0 annotation files are hosted on the HPRC S3 bucket under the T2T annotation prefix.

### CHM13v2.0 genome annotations (used by T2T QC)

```bash
# Genome feature annotation (T2T QC: annotationBed)
aws s3 cp s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_GenomeFeature_v1.0.bed .

# CenSat annotation v2.0 (T2T QC: annotationCENSAT)
aws s3 cp s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.0.bed .

# Segmental duplications — flattened (T2T QC: annotationSD)
# chm13v2.0_SD.flattened.bed is not directly on S3; generate it from the source SD BED:
aws s3 cp s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.bed .
bedtools merge -i chm13v2.0_SD.bed > chm13v2.0_SD.flattened.bed

# CHM13v2.0 reference FASTA (T2T QC alignment)
aws s3 cp s3://human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz .
```

### CenSat HMM profiles (used by CenSat / alphaAnnotation)

```bash
git clone https://github.com/kmiga/alphaAnnotation
# Profiles are in alphaAnnotation/utilities/:
#   AS-HORs-hmmer3.4-071024.hmm
#   AS-SFs-hmmer3.0.290621.hmm
#   rDNA1.0.hmm
#   xy_apes_y_human_newmodels.embl
```

### Quaak reference files

```bash
# sel-500.kmer.fa.gz and chm13v2.ref.path.gz are bundled with the Quaak release:
# https://github.com/Unitig/Quaak/releases
```

### Summary table

| File | Used by | S3 / source |
|---|---|---|
| `chm13v2.0.fa.gz` | T2T QC | `s3://human-pangenomics/T2T/CHM13/assemblies/analysis_set/` |
| `chm13v2.0_GenomeFeature_v1.0.bed` | T2T QC | `s3://human-pangenomics/T2T/CHM13/assemblies/annotation/` |
| `chm13v2.0_SD.flattened.bed` | T2T QC | generated from `chm13v2.0_SD.bed` via `bedtools merge` (see above) |
| `chm13v2.0_censat_v2.0.bed` | T2T QC | `s3://human-pangenomics/T2T/CHM13/assemblies/annotation/` |
| `AS-HORs-hmmer3.4-071024.hmm` | CenSat | alphaAnnotation repo `utilities/` |
| `AS-SFs-hmmer3.0.290621.hmm` | CenSat | alphaAnnotation repo `utilities/` |
| `rDNA1.0.hmm` | CenSat | alphaAnnotation repo `utilities/` |
| `xy_apes_y_human_newmodels.embl` | CenSat | alphaAnnotation repo `utilities/` |
| `primatesnoY_odb10_mb_download.tar.gz` | Compleasm | BUSCO lineage downloads |
| `sel-500.kmer.fa.gz` | Quaak | Quaak release bundle |
| `chm13v2.ref.path.gz` | Quaak | Quaak release bundle |

---

## Step 1: CenSat — Centromere & Satellite Annotation

**Repo:** https://github.com/kmiga/alphaAnnotation  
**Commit:** `35e53eb`  
**Run per haplotype.**

The pipeline is orchestrated by `centromereAnnotation.wdl`. Key underlying tools and versions:

| Tool | Version | Docker |
|---|---|---|
| HMMER (`nhmmer`) | 3.4 | `juklucas/alphasat_hmmer@sha256:ccaedbf4a53f9386017019bd54671b42123c9c6157a289b04f325ffc8a0145da` |
| RepeatMasker | — | same image |
| seqtk | 1.4 | `quay.io/biocontainers/seqtk:1.4--he4a0461_2` |
| bedtools | 2.28.0 | `biocontainers/bedtools:v2.28.0_cv2` |

**Key commands run by the pipeline:**

```bash
# Alpha satellite HOR annotation
nhmmer --cpu 8 \
    --tblout ${SAMPLE}.hap1.AS-HORs.tbl \
    AS-HORs-hmmer3.4-071024.hmm \
    ${SAMPLE}.hap1.fa

# Alpha satellite SF annotation
nhmmer --cpu 8 \
    --tblout ${SAMPLE}.hap1.AS-SFs.tbl \
    AS-SFs-hmmer3.0.290621.hmm \
    ${SAMPLE}.hap1.fa

# rDNA annotation
nhmmer --cpu 8 \
    --tblout ${SAMPLE}.hap1.rDNA.tbl \
    rDNA1.0.hmm \
    ${SAMPLE}.hap1.fa

# RepeatMasker (custom Y/ape models)
RepeatMasker -pa 8 \
    -lib xy_apes_y_human_newmodels.embl \
    -xsmall \
    ${SAMPLE}.hap1.fa
```

**Resources:** 16 CPUs, 64 GB RAM, ~1–2 days

**Output files** (prefix = `{SAMPLE}.hap1`):

```
{SAMPLE}.hap1.cenSat.bed                   main cenSat annotation (used as input to Flagger)
{SAMPLE}.hap1.active.centromeres.bed
{SAMPLE}.hap1.as_hor.bed
{SAMPLE}.hap1.as_hor_sf.bed
{SAMPLE}.hap1.as_sf.bed
{SAMPLE}.hap1.as_strand.bed
{SAMPLE}.hap1.SatelliteStrand.bed
{SAMPLE}.hap1.RepeatMasker.bed
{SAMPLE}.hap1.RepeatMasker.rmsk.bed
{SAMPLE}.hap1.RepeatMasker.out
{SAMPLE}.hap1.RepeatMasker.masked.fasta.gz
{SAMPLE}.hap1_repeat_masker.tar.gz
```

**Note:** Concatenate hap1 and hap2 cenSat BEDs to make the diploid BED required by Flagger:
```bash
cat ${SAMPLE}.hap1.cenSat.bed ${SAMPLE}.hap2.cenSat.bed > ${SAMPLE}.cenSat.dip.bed
```

---

## Step 2: Compleasm — BUSCO Gene Completeness

**Tool:** https://github.com/huangnengcsu/compleasm  
**Docker:** `huangnengcsu/compleasm@sha256:34643e3ebbc9a2dde9d648434e1b79cc576b38967cad86341622de5701da8ede`  
**Lineage:** `primatesnoY_odb10`  
**Run per haplotype.**

```bash
# Extract lineage database first
tar -xzf primatesnoY_odb10_mb_download.tar.gz

docker run --rm \
    -v $(pwd):/data \
    huangnengcsu/compleasm@sha256:34643e3ebbc9a2dde9d648434e1b79cc576b38967cad86341622de5701da8ede \
    compleasm run \
    -a /data/${SAMPLE}.hap1.fa.gz \
    -o /data/compleasm_out/ \
    -l primatesnoY_odb10 \
    -t 24
```

**Resources:** 24 CPUs, 48 GB RAM, ~4–8h

**Output files:**

```
compleasm_out/primatesnoY_odb10/full_table.tsv    per-gene completeness table
compleasm_out/summary.txt                          S/D/F/M/C counts
compleasm_out/primatesnoY_odb10/*.faa              gene sequences
```

---

## Step 3: Quaak — K-mer QV and CNV

**Repo:** https://github.com/YingZhou001/Quaak  
**Commit:** `711991a`  
**Run per haplotype. Not containerized — requires Quaak installed locally.**

```bash
bash quaak.sh \
    -k sel-500.kmer.fa.gz \
    -r chm13v2.ref.path.gz \
    -q ${SAMPLE}.hap1.fa \
    -o outdir/${SAMPLE}.hap1
```

**Resources:** 32 GB RAM, 2 CPUs, ~12h

**Output files** (prefix = `{SAMPLE}.hap1`):

```
{SAMPLE}.hap1.annot.tsv.gz     per-contig k-mer annotation
{SAMPLE}.hap1.block.tsv.gz     haplotype block assignments
{SAMPLE}.hap1.cnv.tsv.gz       copy-number variation
{SAMPLE}.hap1.conf-ctgs.txt    confident contigs
{SAMPLE}.hap1.ctg-strand.tsv   contig strand orientation
{SAMPLE}.hap1.kmer-summary.tsv k-mer QV summary
{SAMPLE}.hap1.que.path.gz      query path
{SAMPLE}.hap1.rescue-ctgs.txt  rescued contigs
{SAMPLE}.hap1.sel-ctgs.txt     selected contigs
{SAMPLE}.hap1.summary.txt      overall QV summary
{SAMPLE}.hap1.suspect-ctgs.txt contigs flagged as suspect
```

---

## Step 4: T2T QC — Assembly Break Annotation

**Repo:** https://github.com/human-pangenomics/hpp_production_workflows  
**Commit:** `afaa8fe`  
**WDL:** `QC/wdl/tasks/find_assembly_breakpoints.wdl`  
**Run per haplotype.**

The pipeline uses mashmap to align the assembly to CHM13, then annotates break points against feature BEDs. Key underlying commands:

```bash
# Align assembly to reference
mashmap \
    -r chm13v2.0.fa.gz \
    -q ${SAMPLE}.hap1.fa.gz \
    --pi 95 -s 50000 \
    -t 32 \
    -o ${SAMPLE}.hap1.mashmap.txt

# Intersect breaks with annotation BEDs
bedtools intersect \
    -a ${SAMPLE}.hap1.breaks.bed \
    -b chm13v2.0_GenomeFeature_v1.0.bed \
    -wao > ${SAMPLE}.hap1.region.breakAnnotation.bed
```

**Docker images:**
- `quay.io/biocontainers/mashmap:2.0--h543ed2d_4`
- `quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6`
- `quay.io/biocontainers/seqtk:1.3--hed695b0_2`
- `quay.io/biocontainers/bioawk:1.0--hed695b0_5`
- `quay.io/biocontainers/datamash:1.1.0--0`
- `rocker/tidyverse`

**Resources:** 32 CPUs, 64 GB RAM, ~6–12h

**Output files** (prefix = `{SAMPLE}.hap1`):

```
{SAMPLE}.hap1.CENSAT.breakAnnotation.bed
{SAMPLE}.hap1.region.breakAnnotation.bed
{SAMPLE}.hap1.SD.breakAnnotation.bed
{SAMPLE}.hap1.filteredFlanks.bed
{SAMPLE}.hap1.formattedBreakAnnotation_CENSAT.txt
{SAMPLE}.hap1.formattedBreakAnnotation_region.txt
{SAMPLE}.hap1.formattedBreakAnnotation_SD.txt
{SAMPLE}.hap1.mashmap.txt
{SAMPLE}.hap1.T2T.contigs.txt      contigs with T2T-quality mapping
{SAMPLE}.hap1.T2T.scaffolds.txt
{SAMPLE}.hap1.unifiedAssembly.fa
assembly.statistics.txt
```

---

## Step 5: Flagger — Read Alignment and HMM Flagging

**Repo:** https://github.com/mobinasri/flagger  
**Commit:** `f30ff97`  
**WDL:** `wdls/workflows/hmm_flagger_end_to_end_with_mapping.wdl`

Run **once per sample per read type** (HiFi and ONT R10 separately). Aligns reads to the diploid assembly (both haplotypes concatenated), then runs HMM-based flagging per haplotype. Requires the diploid cenSat BED from step 1.

**Key steps run internally:**

```bash
# 1. Concatenate haplotypes for diploid alignment
cat ${SAMPLE}.hap1.fa.gz ${SAMPLE}.hap2.fa.gz > ${SAMPLE}.dip.fa.gz

# 2. Align HiFi reads
minimap2 \
    -ax lr:hqae \
    -k 25 \
    --secondary=no \
    -t 32 \
    ${SAMPLE}.dip.fa.gz \
    ${SAMPLE}.HiFi.bam \
    | samtools sort -o ${SAMPLE}.hifi_minimap2.bam
samtools index ${SAMPLE}.hifi_minimap2.bam

# For ONT R10: -ax lr:hq -k 15

# 3. HMM flagging
flagger \
    --bam ${SAMPLE}.hifi_minimap2.bam \
    --cntr-bed ${SAMPLE}.cenSat.dip.bed \
    --preset hifi \
    --output-dir flagger_out/
# For ONT R10: --preset ont_r10
```

**Docker images:**
- `mobinasri/long_read_aligner:v1.1.0` (minimap2 alignment)
- `mobinasri/flagger:v1.2.0` (HMM flagging)

**Resources:** 32 CPUs, 128 GB RAM, ~5–7 days

**Output files** (prefix = `{SAMPLE}.hifi_minimap2`; replace `hifi` with `ont_r10` for ONT):

```
{SAMPLE}.hifi_minimap2.bam                                                              read-to-assembly BAM (input to NucFlag)
{SAMPLE}.hifi_minimap2.bam.bai

{SAMPLE}.hifi_minimap2.downsample_1.0.hifi_flagger_prediction.bed
{SAMPLE}.hifi_minimap2.downsample_1.0.hifi_flagger_prediction.canonical.no_Hap.bed
{SAMPLE}.hifi_minimap2.downsample_1.0.hifi_flagger_prediction.canonical.no_Hap.hap1.bed
{SAMPLE}.hifi_minimap2.downsample_1.0.hifi_flagger_prediction.canonical.no_Hap.hap2.bed
{SAMPLE}.hifi_minimap2.downsample_1.0.hifi_flagger_prediction.conservative.bed
{SAMPLE}.hifi_minimap2.downsample_1.0.hifi_flagger_prediction.conservative.canonical.no_Hap.bed
{SAMPLE}.hifi_minimap2.downsample_1.0.hifi_flagger_prediction.conservative.canonical.no_Hap.hap1.bed
{SAMPLE}.hifi_minimap2.downsample_1.0.hifi_flagger_prediction.conservative.canonical.no_Hap.hap2.bed

{SAMPLE}.hifi_minimap2.downsample_1.0.augmented.cov.gz
{SAMPLE}.hifi_minimap2.downsample_1.0.augmented.stats.tsv
{SAMPLE}.hifi_minimap2.downsample_1.0.augmented.conservative.stats.tsv
{SAMPLE}.hifi_minimap2.downsample_1.0.bias_detection_table.tsv
{SAMPLE}.hifi_minimap2.downsample_1.0.bigwig
{SAMPLE}.hifi_minimap2.downsample_1.0.high_mapq.bigwig
{SAMPLE}.hifi_minimap2.downsample_1.0.high_clip.bigwig
{SAMPLE}.hifi_minimap2.downsample_1.0.region_idx.bigwig
{SAMPLE}.hifi_minimap2.downsample_1.0.mapq_ge_10.cov_ge_4.mappable.hap1.bed
{SAMPLE}.hifi_minimap2.downsample_1.0.mapq_ge_10.cov_ge_4.mappable.hap2.bed
{SAMPLE}.hifi_minimap2.loglikelihood.tsv
```

---

## Step 6: NucFlag — Secondary-Base Misassembly Detection

**Tool:** https://github.com/logsdon-lab/NucFlag  
**Version:** 1.0.0  
**Docker:** `iviolich/nucflag:1.0.0`  
**samtools:** 1.19.2 (`quay.io/biocontainers/samtools@sha256:9cd15e719101ae8808e4c3f152cca2bf06f9e1ad8551ed43c1e626cb6afdaa02`)  
**Run per haplotype per read type.** Input BAM is from Flagger (step 5).

```bash
# 1. Build per-haplotype regions BED from assembly FASTA
samtools faidx ${SAMPLE}.hap1.fa
awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' ${SAMPLE}.hap1.fa.fai > ${SAMPLE}.hap1.genome.bed

# 2. Filter diploid BAM to primary alignments on hap1 contigs only
#    -F 2308 drops unmapped (4), secondary (256), supplementary (2048)
samtools view \
    -F 2308 \
    --bam --with-header \
    --regions-file ${SAMPLE}.hap1.genome.bed \
    --threads 8 \
    ${SAMPLE}.hifi_minimap2.bam \
    -o ${SAMPLE}.hifi_minimap2.hap1.primary.bam
samtools index ${SAMPLE}.hifi_minimap2.hap1.primary.bam

# 3. Run NucFlag
mkdir -p plots/
docker run --rm \
    -v $(pwd):/data \
    iviolich/nucflag:1.0.0 \
    nucflag call \
    -i /data/${SAMPLE}.hifi_minimap2.hap1.primary.bam \
    -f /data/${SAMPLE}.hap1.fa \
    -d /data/plots/ \
    -o /data/${SAMPLE}.hifi.hap1.nucflag.bed \
    -s /data/${SAMPLE}.hifi.hap1.nucflag_status.bed \
    -t 4 \
    -p 4

# For ONT R10, add: -x ont_r10
# Note: NucFlag must be reindexed inside the container; the WDL handles this automatically.

tar -czf ${SAMPLE}.hifi.hap1.nucflag_plots.tar.gz plots/
```

**Important:** Do **not** pass `-b` (regions BED) to `nucflag call`. Without it, NucFlag uses default 10 MB non-overlapping windows — the correct mode for genome-wide QC. Passing a whole-contig BED causes O(N²) self-identity computation on chromosome-scale contigs and runtimes of several days.

**Resources:** 4 CPUs, 64 GB RAM, ~12–24h per haplotype

**Output files** (prefix = `{SAMPLE}.hifi.hap1`):

```
{SAMPLE}.hifi.hap1.nucflag.bed          flagged regions: MISJOIN, HET, COLLAPSE_VAR, COLLAPSE_DEL
{SAMPLE}.hifi.hap1.nucflag_status.bed   per-region pass/fail
{SAMPLE}.hifi.hap1.nucflag_plots.tar.gz per-region coverage plots
```

---

## Output Naming Summary

| Level | Pattern |
|---|---|
| Per-haplotype | `{SAMPLE}.hap{1\|2}.*` |
| Per-sample, per-read-type (Flagger BAM) | `{SAMPLE}.{hifi\|ont_r10}_minimap2.*` |
| Per-haplotype, per-read-type (NucFlag) | `{SAMPLE}.{hifi\|ont_r10}.hap{1\|2}.*` |

`{SAMPLE}` = HPRC sample ID (e.g. `HG00097`)
