# ANAPIPE

Pipeline for bacterial amplicon TGS sequencing analysis

## Usage

Example of launch

```bash
nextflow run main.nf -params-file ../data_for_anapipe/params.json -profile singularity --input_worklist worklist.csv --input_reads reads_pwd --outdir out_pwd
```

or

```bash
while IFS=$',' read -r wls reads out; do nextflow run main.nf -params-file ../data_for_anapipe/params.json -profile singularity --input_worklist "${wls}" --input_reads "${reads}" --outdir "${out}"; done < runs.csv
```

where the `runs.csv` is something like

```csv
Plate_1.csv,/path/to/temp_storage/flow_cell/fastq_pass/,anapipe_out_Plate_1
Plate_2.csv,/path/to/temp_storage/flow_cell/fastq_pass/,anapipe_out_Plate_2
Plate_3.csv,/path/to/temp_storage/flow_cell/fastq_pass/,anapipe_out_Plate_3
Plate_4.csv,/path/to/temp_storage/flow_cell/fastq_pass/,anapipe_out_Plate_4
Plate_5.csv,/path/to/temp_storage/flow_cell/fastq_pass/,anapipe_out_Plate_5

```

It's important to set up `params.json` before the pipeline launch!

## Important params

TBD...

## Notes

For the `CLAIR3` nf-core module the container was changed from `biocontainers/clair3:1.0.10--py39hd649744_1` to `biocontainers/clair3:1.0.11--py39hd649744_0` 
to enable `--enable_variant_calling_at_sequence_head_and_tail` handling.
medaka:1.4.4--py38h130def0_0 medaka:2.1.0--py38ha0c3a46_0
