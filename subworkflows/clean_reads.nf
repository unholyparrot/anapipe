include { CHOPPER } from '../modules/nf-core/chopper/main'
include { KRAKEN2_KRAKEN2 } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_EXTRACTKRAKENREADS } from '../modules/nf-core/krakentools/extractkrakenreads/main'

workflow CLEAN_READS {
    take:
        input_reads
    
    main:
        ch_multiqc_files = Channel.empty()

        if (params.chopper_fasta) {
            CHOPPER(input_reads, Channel.fromPath(params.chopper_fasta))
        }
        else {
            CHOPPER(input_reads, [])
        }

        ch_reads_len_qual_filtered = CHOPPER.out.fastq

        KRAKEN2_KRAKEN2(ch_reads_len_qual_filtered, params.db_kraken2, true, true)
        ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))

        ch_merge_kraken = KRAKEN2_KRAKEN2.out.classified_reads_assignment.combine (
            KRAKEN2_KRAKEN2.out.classified_reads_fastq, by: 0
        ).combine (
            KRAKEN2_KRAKEN2.out.report, by: 0
        )

        KRAKENTOOLS_EXTRACTKRAKENREADS(
            params.desired_kraken2_taxid, 
            ch_merge_kraken.map{it -> [it[0], it[1]]}, 
            ch_merge_kraken.map{it -> [it[0], it[2]]}, 
            ch_merge_kraken.map{it -> [it[0], it[3]]}
        )

        nf_clean_reads = KRAKENTOOLS_EXTRACTKRAKENREADS.out.extracted_kraken2_reads
        multiqc_files = ch_multiqc_files

        clean_reads = nf_clean_reads.filter{
            it -> it[1].countFastq() > params.min_reads_per_sample
        }
    
    emit:
        clean_reads
        multiqc_files
}
