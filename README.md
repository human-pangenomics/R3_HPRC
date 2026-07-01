# R3_HPRC

Tracking repository for the HPRC R3 sample cohort. Contains metadata and sequence data manifests for 345 samples assembled as part of the Human Pangenome Reference Consortium.

## Sequence Data

| Path | Description |
|---|---|
| `sample/R3_HPRC_metadata.csv` | Sample metadata (sex, population, data availability) |
| `sample/R3_HPRC_coverage_summary.csv` | Per-sample sequencing coverage summary |
| `sequence_data/R3_HPRC_hifi_no_dc_map.csv` | HiFi BAM paths (S3) |
| `sequence_data/R3_HPRC_ont.csv` | ONT BAM paths (S3) |
| `sequence_data/R3_HPRC_hic.csv` | Hi-C FASTQ paths (S3) |
| `sequence_data/R3_HPRC_ill.csv` | Illumina FASTQ paths (S3) |
| `sequence_data/R3_HPRC_kinnex.csv` | Kinnex read paths (S3) |
| `sequence_data/R3_HPRC_dc.csv` | Duplex/consensus read paths (S3) |

## Assemblies

R3 assemblies use **verkko v2.3.2** in hybrid mode (HiFi + ONT UL + Hi-C phasing). 

Assemblies are deposited to S3 under the submission prefix:[R3_verkko-v2.3.2_hybrid_assembly](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/bed38f3f-cac4-4703-9852-181bc362cea1--R3_verkko-v2.3.2_hybrid_assembly/)

## Assembly QC

Per-sample QC is run for each batch as assemblies complete, covering:

- **Yak QV** — k-mer-based quality value (Hi-C reads vs assembly)
- **Compleasm** — BUSCO-based gene completeness
- **Dipcall / T2T** — telomere-to-telomere contig and scaffold counts
- **Quaak** — CHM13v2 unique k-mer presence assessment
- **NucFlag** — coverage-based misassembly flagging
- **CenSat** — centromere and satellite annotation
- **Flagger** — HiFi and ONT alignment-based error/duplication/collapse flagging

QC pipeline results: [assembly_qc](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/76E3F239-A9A6-4B6B-840F-B49053807394--R3_verkko-v2.3.2_hybrid_assembly_qc/)

QC summary spreadsheet: [R3 Assembly QC Summary](https://docs.google.com/spreadsheets/d/1nUOeVLDGB_MU6VF3FnHT-muadScXVcKUybEX3Q345ic/edit?usp=sharing)

Full Quaak outputs: [R3_verkko-v2.3.2_hybrid_assembly_quaak](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/77F76FE5-3015-4556-8B28-DB5B8C5F04A4--R3_verkko-v2.3.2_hybrid_assembly_quaak/)
