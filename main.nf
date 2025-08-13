include { INPUT } from './subworkflows/input.nf'
include {QC_CONTROLS as QC_CONTROLS_PRE} from './subworkflows/quality_check.nf'
include {CLEAN_READS} from './subworkflows/clean_reads.nf'
include {QC_CONTROLS as QC_CONTROLS_POST} from './subworkflows/quality_check.nf'
include {MAP_AND_MUTATIONS} from './subworkflows/map_and_mutations.nf'
include {MAP_CONSENSUS} from './subworkflows/map_consensus.nf'


include { MULTIQC } from './modules/nf-core/multiqc/main'

workflow {
    ch_multiqc_files = Channel.empty()

    INPUT()

    QC_CONTROLS_PRE(INPUT.out.ch_input_reads, "pre")
    ch_multiqc_files = ch_multiqc_files.mix(QC_CONTROLS_PRE.out.multiqc_files.collect())

    CLEAN_READS(INPUT.out.ch_input_reads)
    ch_multiqc_files = ch_multiqc_files.mix(CLEAN_READS.out.multiqc_files.collect())

    QC_CONTROLS_POST(CLEAN_READS.out.clean_reads, "post")
    ch_multiqc_files = ch_multiqc_files.mix(QC_CONTROLS_POST.out.multiqc_files.collect())

    MAP_AND_MUTATIONS(CLEAN_READS.out.clean_reads)

    MAP_CONSENSUS(
        CLEAN_READS.out.clean_reads.combine(
            MAP_AND_MUTATIONS.out.out_fasta, by: 0
        )
    )

    MULTIQC(ch_multiqc_files.collect(), [], [], [], [], [])
}
