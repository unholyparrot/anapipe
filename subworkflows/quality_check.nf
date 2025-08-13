include { NANOPLOT } from '../modules/nf-core/nanoplot/main'
include { FASTQC } from '../modules/nf-core/fastqc/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_RAW} from '../modules/nf-core/minimap2/align/main'

workflow QC_CONTROLS {
    take:
        input_reads // [meta, [short_reads_1, short_reads_2], long_reads]
        suffix     // String to append to output names (e.g., 'pre_qc' or 'post_qc')
    main:
        ch_multiqc_files = Channel.empty()

        // Modify meta to include the suffix
        input_reads_modified = input_reads.map { meta, long_reads ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_${suffix}"
            [new_meta, long_reads]
        }

        ch_long_reads = input_reads_modified.map{meta, long_reads -> [meta, long_reads]}.filter {
            _meta, reads -> reads && reads.size() > 0
        }

        FASTQC(ch_long_reads)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

        NANOPLOT(ch_long_reads)
        ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.collect{it[1]}.ifEmpty([]))

        if (suffix == "pre") {

            ch_reference = Channel.fromPath(params.reference_fasta).map {
                it -> 
                    def ref_meta = [:]
                    ref_meta.id = "reference"
                    return [ref_meta, it]
            }

            MINIMAP2_ALIGN_RAW(input_reads_modified, ch_reference.collect(), true, "bai", false, false)
        }

    emit:
        multiqc_files = ch_multiqc_files
}
