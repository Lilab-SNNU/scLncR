# scLncR

Primary project documentation is maintained in [ReadMe.md](./ReadMe.md).

Current design summary:
- raw FASTQ-first, technology-aware workflow;
- prelnc performs candidate lncRNA discovery from raw FASTQ (with internal BAM intermediates);
- count performs lncRNA-aware quantification from raw FASTQ + augmented reference.
