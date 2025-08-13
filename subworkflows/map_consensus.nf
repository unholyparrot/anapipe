include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_CONSENSUS} from '../modules/nf-core/minimap2/align/main'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_CONSENSUS } from '../modules/nf-core/bedtools/genomecov/main'

workflow MAP_CONSENSUS {
    take:
        input_data
    
    main:
        MINIMAP2_ALIGN_CONSENSUS(
            input_data.map { it -> [it[0], it[1]]},
            input_data.map { it -> [it[0], it[2]]},
            true, "bai", false, false
        )

        BEDTOOLS_GENOMECOV_CONSENSUS(
            MINIMAP2_ALIGN_CONSENSUS.out.bam.map {
                it -> [it[0], it[1], 1]
            },
            [],
            "bed",
            true
        )
}
