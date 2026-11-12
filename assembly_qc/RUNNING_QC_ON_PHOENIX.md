# Running CenSat → Flagger → NucFlag on Phoenix

Step-by-step guide for running the centromere annotation, read-to-assembly alignment/flagging, and misassembly detection pipeline on a single diploid assembly on the Phoenix HPC cluster.

The pipeline runs each step as a SLURM job. All tools execute inside Docker containers via Singularity — no manual software installation needed.

---

## Prerequisites

**Toil virtualenv** — activate before submitting any job:
```bash
source ~/toil_env/bin/activate
```

**AWS credentials** — required only if your assembly FASTAs or read files are on S3:
```bash
# Check whether credentials are already configured
aws sts get-caller-identity
# If not, ask iviolich@ucsc.edu to share the credentials file
```

**Reference files** — already installed on Phoenix at:
```
/private/nanopore/tools/alphaAnnotation/utilities/
    AS-HORs-hmmer3.4-071024.hmm
    AS-SFs-hmmer3.0.290621.hmm
    rDNA1.0.hmm
    xy_apes_y_human_newmodels.embl
```

---

## Inputs

| File | Notes |
|------|-------|
| `{SAMPLE}.haplotype1.fasta[.gz]` | Per-haplotype assembly FASTA |
| `{SAMPLE}.haplotype2.fasta[.gz]` | Per-haplotype assembly FASTA |
| HiFi read BAMs | PacBio HiFi uBAMs or unmapped BAMs (list, one per file) |
| ONT R10 read BAMs | Optional; same format |

Paths can be local (`/path/to/file`) or S3 URIs (`s3://bucket/key`).

---

## Setup

Create a working directory and set a variable for your sample:

```bash
SAMPLE=HG01234    # replace with your sample ID
WORKDIR=/private/nanopore2/hprc/my_qc_runs   # replace with your preferred location

mkdir -p ${WORKDIR}/${SAMPLE}
cd ${WORKDIR}/${SAMPLE}
mkdir -p slurm_logs
```

---

## Step 1: CenSat — centromere & satellite annotation

Run once per haplotype. Each job takes ~1–2 days on 16 CPUs.

**Create input JSONs** — repeat for `hap2`, swapping the FASTA path:

```bash
HAP1_FASTA=/path/to/${SAMPLE}.haplotype1.fasta   # or s3://...
HAP2_FASTA=/path/to/${SAMPLE}.haplotype2.fasta

cat > ${SAMPLE}_hap1_centromereAnnotation.json << EOF
{
  "centromereAnnotation.fasta": "${HAP1_FASTA}",
  "centromereAnnotation.rDNAhmm_profile": "/private/nanopore/tools/alphaAnnotation/utilities/rDNA1.0.hmm",
  "centromereAnnotation.AS_hmm_profile": "/private/nanopore/tools/alphaAnnotation/utilities/AS-HORs-hmmer3.4-071024.hmm",
  "centromereAnnotation.AS_hmm_profile_SF": "/private/nanopore/tools/alphaAnnotation/utilities/AS-SFs-hmmer3.0.290621.hmm",
  "centromereAnnotation.additionalRMModels": "/private/nanopore/tools/alphaAnnotation/utilities/xy_apes_y_human_newmodels.embl"
}
EOF

cat > ${SAMPLE}_hap2_centromereAnnotation.json << EOF
{
  "centromereAnnotation.fasta": "${HAP2_FASTA}",
  "centromereAnnotation.rDNAhmm_profile": "/private/nanopore/tools/alphaAnnotation/utilities/rDNA1.0.hmm",
  "centromereAnnotation.AS_hmm_profile": "/private/nanopore/tools/alphaAnnotation/utilities/AS-HORs-hmmer3.4-071024.hmm",
  "centromereAnnotation.AS_hmm_profile_SF": "/private/nanopore/tools/alphaAnnotation/utilities/AS-SFs-hmmer3.0.290621.hmm",
  "centromereAnnotation.additionalRMModels": "/private/nanopore/tools/alphaAnnotation/utilities/xy_apes_y_human_newmodels.embl"
}
EOF
```

**Submit both haplotypes:**

```bash
source ~/toil_env/bin/activate

WDL=/private/nanopore/tools/alphaAnnotation/cenSatAnnotation/centromereAnnotation.wdl

for HAP in hap1 hap2; do
    mkdir -p ${HAP}
    sbatch \
        --job-name=censat_${SAMPLE}_${HAP} \
        --partition=long \
        --time=72:00:00 \
        --cpus-per-task=16 \
        --mem=64G \
        --output=slurm_logs/censat_${SAMPLE}_${HAP}_%j.log \
        --wrap="
            source ~/toil_env/bin/activate
            cd ${WORKDIR}/${SAMPLE}/${HAP}
            mkdir -p analysis/centromereAnnotation_outputs toil_logs

            SCRATCH=/data/scratch/\$(whoami)/${SAMPLE}_censat_${HAP}
            mkdir -p \${SCRATCH}/workdir

            toil-wdl-runner \
                --jobStore file:\${SCRATCH}/jobstore \
                --workDir \${SCRATCH}/workdir \
                --clean=never \
                --batchSystem single_machine \
                --maxCores 16 \
                --retryCount 1 \
                --caching=false \
                --disableProgress \
                ${WDL} \
                ${WORKDIR}/${SAMPLE}/${SAMPLE}_${HAP}_centromereAnnotation.json \
                --outputDirectory analysis/centromereAnnotation_outputs \
                2>&1 | tee censat_${HAP}.log
        "
done
```

**When both haplotypes finish**, combine into the diploid BED required by Flagger:

```bash
HAP1_BED=$(find hap1/analysis/centromereAnnotation_outputs -name "*.cenSat.bed" | head -1)
HAP2_BED=$(find hap2/analysis/centromereAnnotation_outputs -name "*.cenSat.bed" | head -1)

cat ${HAP1_BED} ${HAP2_BED} > ${SAMPLE}.cenSat.dip.bed
echo "Diploid cenSat BED: $(pwd)/${SAMPLE}.cenSat.dip.bed"
```

---

## Step 2: Flagger — read alignment and HMM flagging

Run once per read type (HiFi and/or ONT R10). Takes ~5–7 days per run on 32 CPUs.

**Create input JSONs:**

Replace the `readFiles` arrays with your actual read file paths. Files can be uBAMs, mapped BAMs, or FASTQs.

```bash
# HiFi JSON
cat > ${SAMPLE}_hifi_HMMFlaggerEndToEndWithMapping.json << EOF
{
  "HMMFlaggerEndToEndWithMapping.sampleName": "${SAMPLE}",
  "HMMFlaggerEndToEndWithMapping.suffixForMapping": "hifi_minimap2",
  "HMMFlaggerEndToEndWithMapping.suffixForFlagger": "hifi_flagger",
  "HMMFlaggerEndToEndWithMapping.hap1AssemblyFasta": "${HAP1_FASTA}",
  "HMMFlaggerEndToEndWithMapping.hap2AssemblyFasta": "${HAP2_FASTA}",
  "HMMFlaggerEndToEndWithMapping.readFiles": [
    "/path/to/${SAMPLE}.hifi.reads.part1.bam",
    "/path/to/${SAMPLE}.hifi.reads.part2.bam"
  ],
  "HMMFlaggerEndToEndWithMapping.presetForFlagger": "hifi",
  "HMMFlaggerEndToEndWithMapping.presetForMapping": "lr:hqae",
  "HMMFlaggerEndToEndWithMapping.kmerSize": 25,
  "HMMFlaggerEndToEndWithMapping.cntrBed": "${WORKDIR}/${SAMPLE}/${SAMPLE}.cenSat.dip.bed",
  "HMMFlaggerEndToEndWithMapping.enableDecomposingCntrBed": true
}
EOF

# ONT R10 JSON (optional)
cat > ${SAMPLE}_ont_r10_HMMFlaggerEndToEndWithMapping.json << EOF
{
  "HMMFlaggerEndToEndWithMapping.sampleName": "${SAMPLE}",
  "HMMFlaggerEndToEndWithMapping.suffixForMapping": "ont_r10_minimap2",
  "HMMFlaggerEndToEndWithMapping.suffixForFlagger": "ont_r10_flagger",
  "HMMFlaggerEndToEndWithMapping.hap1AssemblyFasta": "${HAP1_FASTA}",
  "HMMFlaggerEndToEndWithMapping.hap2AssemblyFasta": "${HAP2_FASTA}",
  "HMMFlaggerEndToEndWithMapping.readFiles": [
    "/path/to/${SAMPLE}.ont_r10.reads.part1.bam",
    "/path/to/${SAMPLE}.ont_r10.reads.part2.bam"
  ],
  "HMMFlaggerEndToEndWithMapping.presetForFlagger": "ont-r10",
  "HMMFlaggerEndToEndWithMapping.presetForMapping": "lr:hqae",
  "HMMFlaggerEndToEndWithMapping.kmerSize": 25,
  "HMMFlaggerEndToEndWithMapping.cntrBed": "${WORKDIR}/${SAMPLE}/${SAMPLE}.cenSat.dip.bed",
  "HMMFlaggerEndToEndWithMapping.enableDecomposingCntrBed": true
}
EOF
```

**Submit** (run after cenSat is complete):

```bash
WDL=/private/nanopore/tools/hpp_production_workflows/QC/wdl/workflows/hmm_flagger_end_to_end_with_mapping.wdl

for READTYPE in hifi ont_r10; do
    mkdir -p ${READTYPE}
    sbatch \
        --job-name=flagger_${SAMPLE}_${READTYPE} \
        --partition=long \
        --time=7-00:00:00 \
        --cpus-per-task=32 \
        --mem=200G \
        --output=slurm_logs/flagger_${SAMPLE}_${READTYPE}_%j.log \
        --wrap="
            source ~/toil_env/bin/activate
            cd ${WORKDIR}/${SAMPLE}/${READTYPE}
            mkdir -p analysis/hmm_flagger_end_to_end_with_mapping_outputs toil_logs

            SCRATCH=/data/scratch/\$(whoami)/${SAMPLE}_flagger_${READTYPE}
            mkdir -p \${SCRATCH}/workdir

            toil-wdl-runner \
                --jobStore file:\${SCRATCH}/jobstore \
                --workDir \${SCRATCH}/workdir \
                --clean=never \
                --batchSystem single_machine \
                --maxCores 32 \
                --retryCount 1 \
                --caching=false \
                --disableProgress \
                ${WDL} \
                ${WORKDIR}/${SAMPLE}/${SAMPLE}_${READTYPE}_HMMFlaggerEndToEndWithMapping.json \
                --outputDirectory analysis/hmm_flagger_end_to_end_with_mapping_outputs \
                2>&1 | tee flagger_${READTYPE}.log
        "
done
```

**Key outputs** (inside `{hifi,ont_r10}/analysis/hmm_flagger_end_to_end_with_mapping_outputs/*/`):

```
{SAMPLE}.hifi_minimap2.bam          # diploid read-to-assembly BAM — input to NucFlag
{SAMPLE}.hifi_minimap2.bam.bai
{SAMPLE}.hifi_minimap2.downsample_1.0.hifi_flagger_prediction.bed
{SAMPLE}.hifi_minimap2.downsample_1.0.hifi_flagger_prediction.canonical.no_Hap.hap1.bed
{SAMPLE}.hifi_minimap2.downsample_1.0.hifi_flagger_prediction.canonical.no_Hap.hap2.bed
```

---

## Step 3: NucFlag — misassembly detection

Run once per haplotype per read type (up to 4 jobs: hap1+hap2 × hifi+ont_r10). Each takes ~12–24h on 4 CPUs.

**Locate the Flagger BAM** from step 2:

```bash
# Find the BAM — it's inside a UUID subdirectory
HIFI_BAM=$(find hifi/analysis/hmm_flagger_end_to_end_with_mapping_outputs \
    -name "${SAMPLE}.hifi_minimap2.bam" | head -1)
HIFI_BAI=${HIFI_BAM}.bai

ONT_BAM=$(find ont_r10/analysis/hmm_flagger_end_to_end_with_mapping_outputs \
    -name "${SAMPLE}.ont_r10_minimap2.bam" | head -1)
ONT_BAI=${ONT_BAM}.bai

echo "HiFi BAM: ${HIFI_BAM}"
echo "ONT BAM:  ${ONT_BAM}"
```

**Create input JSONs** — one per haplotype per read type:

```bash
WDL=/private/nanopore/hprc_qc/wdls/nucflag.wdl

for HAP in hap1 hap2; do
    if [[ "${HAP}" == "hap1" ]]; then
        ASM_FASTA=${HAP1_FASTA}
    else
        ASM_FASTA=${HAP2_FASTA}
    fi

    # HiFi
    cat > ${SAMPLE}.${HAP}_hifi_NucFlag.json << EOF
{
  "NucFlag.inputBam": "${HIFI_BAM}",
  "NucFlag.inputBamBai": "${HIFI_BAI}",
  "NucFlag.assemblyFasta": "${ASM_FASTA}",
  "NucFlag.sampleName": "${SAMPLE}.${HAP}"
}
EOF

    # ONT R10 (omit if not running ONT)
    cat > ${SAMPLE}.${HAP}_ont_r10_NucFlag.json << EOF
{
  "NucFlag.inputBam": "${ONT_BAM}",
  "NucFlag.inputBamBai": "${ONT_BAI}",
  "NucFlag.assemblyFasta": "${ASM_FASTA}",
  "NucFlag.sampleName": "${SAMPLE}.${HAP}",
  "NucFlag.preset": "ont_r10"
}
EOF

done
```

**Submit:**

```bash
for HAP in hap1 hap2; do
    for READTYPE in hifi ont_r10; do
        mkdir -p nucflag_${HAP}_${READTYPE}
        sbatch \
            --job-name=nucflag_${SAMPLE}_${HAP}_${READTYPE} \
            --partition=long \
            --time=2-00:00:00 \
            --cpus-per-task=4 \
            --mem=64G \
            --output=slurm_logs/nucflag_${SAMPLE}_${HAP}_${READTYPE}_%j.log \
            --wrap="
                source ~/toil_env/bin/activate
                cd ${WORKDIR}/${SAMPLE}/nucflag_${HAP}_${READTYPE}
                mkdir -p analysis/NucFlag_outputs

                SCRATCH=/data/scratch/\$(whoami)/${SAMPLE}_nucflag_${HAP}_${READTYPE}
                mkdir -p \${SCRATCH}/workdir

                toil-wdl-runner \
                    --jobStore file:\${SCRATCH}/jobstore \
                    --workDir \${SCRATCH}/workdir \
                    --clean=never \
                    --batchSystem single_machine \
                    --maxCores 4 \
                    --retryCount 1 \
                    --caching=false \
                    --disableProgress \
                    ${WDL} \
                    ${WORKDIR}/${SAMPLE}/${SAMPLE}.${HAP}_${READTYPE}_NucFlag.json \
                    --outputDirectory analysis/NucFlag_outputs \
                    2>&1 | tee nucflag_${HAP}_${READTYPE}.log
            "
    done
done
```

**Key outputs** (inside `nucflag_{hap}_{readtype}/analysis/NucFlag_outputs/*/`):

```
{SAMPLE}.{hap}.nucflag.bed          # flagged regions: MISJOIN, HET, COLLAPSE_VAR, COLLAPSE_DEL
{SAMPLE}.{hap}.nucflag_status.bed   # per-region pass/fail status
{SAMPLE}.{hap}.nucflag_plots.tar.gz # per-window coverage plots
```

---

## Troubleshooting

**Job failed — how to retry**

Toil saves its state in the jobstore (`/data/scratch/{user}/{sample}_{step}/jobstore`). To restart from the last checkpoint rather than from scratch, add `--restart` to the `toil-wdl-runner` call in the failed job:

```bash
# Edit the wrap command to add --restart after toil-wdl-runner, then resubmit
toil-wdl-runner --restart \
    --jobStore file:/data/scratch/$(whoami)/${SAMPLE}_flagger_hifi/jobstore \
    ...
```

Do **not** run `--clean` or delete the jobstore directory before restarting — Toil uses it to resume cached task results.

**Check job status**

```bash
squeue -u $(whoami) --format="%.18i %.30j %.8T %.10M %R"
```

**Read the log**

Each step writes a log alongside the output directory:
```bash
tail -f slurm_logs/flagger_${SAMPLE}_hifi_<jobid>.log
```

**Stale jobstore from a previous run**

If a fresh job fails immediately with "job store already exists", the scratch directory has a leftover jobstore from a prior attempt. Remove it:
```bash
rm -rf /data/scratch/$(whoami)/${SAMPLE}_flagger_hifi/
```
Then resubmit without `--restart`.

---

## Expected runtimes

| Step | CPUs | Memory | Wall time |
|------|------|--------|-----------|
| CenSat (per haplotype) | 16 | 64 GB | 1–2 days |
| Flagger HiFi | 32 | 200 GB | 5–7 days |
| Flagger ONT R10 | 32 | 200 GB | 5–7 days |
| NucFlag (per haplotype per read type) | 4 | 64 GB | 12–24h |

Questions: contact iviolich@ucsc.edu.
