# scLncR

Primary project documentation is maintained in [ReadMe.md](./ReadMe.md).

Current design summary:
- raw FASTQ-first, technology-aware workflow;
- optional raw FASTQ QC with FastQC and MultiQC before prelnc/count;
- prelnc performs candidate lncRNA discovery from raw FASTQ (with internal BAM intermediates);
- count performs lncRNA-aware quantification from raw FASTQ + augmented reference.

## Raw FASTQ quality control

Run optional QC before prelnc/count:

```bash
scLncR qc -c R/confings/config_QC.yaml
```

Equivalent shell wrapper:

```bash
bash scripts/run_fastq_qc.sh -c R/confings/config_QC.yaml
```

This module only runs FastQC and MultiQC. It does not trim, filter, remove adapters, or modify raw FASTQ files. Outputs are written under `output_dir/fastqc`, `output_dir/multiqc`, `output_dir/logs`, and `output_dir/qc_summary.md`.
