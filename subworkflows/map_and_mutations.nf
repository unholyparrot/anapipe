include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_CLEAN} from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_SPLIT_ALIGNMENTS } from '../modules/local/samtools/view/main'
include { MEDAKA as MEDAKA_1} from '../modules/local/medaka/main'
include { MEDAKA as MEDAKA_2} from '../modules/local/medaka/main'
include { BEDTOOLS_GENOMECOV } from '../modules/nf-core/bedtools/genomecov/main'
include { SUBSET_COVERAGE as SUBSET_COVERAGE_AMAP } from '../modules/local/awk/subset_coverage'
include { CLAIR3 } from '../modules/nf-core/clair3/main'
include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_PRE} from '../modules/nf-core/bcftools/annotate/main'
include { SNPEFF_SNPEFF } from '../modules/nf-core/snpeff/snpeff/main'
include { SNPEFF_DOWNLOAD } from '../modules/nf-core/snpeff/download/main'
include { BCFTOOLS_ANNOTATE as BCFTOOLS_ANNOTATE_POST} from '../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_VIEW } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_COVERED } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_INDEX } from '../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_CONSENSUS } from '../modules/nf-core/bcftools/consensus/main'
include { SEQKIT_SUBSEQ } from '../modules/local/seqkit/subseq'
include { SEQKIT_REPLACE } from '../modules/nf-core/seqkit/replace/main'


workflow MAP_AND_MUTATIONS {
    take:
        input_reads
    main:

        // индексируем референсную последовательность
        ch_reference = Channel.fromPath(params.reference_fasta).map {
            it -> 
                def ref_meta = [:]
                ref_meta.id = "reference"
                return [ref_meta, it]
        }

        SAMTOOLS_FAIDX(ch_reference, [[:], []], false)

        ch_fai = SAMTOOLS_FAIDX.out.fai

        // картирование прочтений на референсную последовательность
        MINIMAP2_ALIGN_CLEAN(input_reads, ch_reference.collect(), true, "bai", false, false)

        ch_bam = MINIMAP2_ALIGN_CLEAN.out.bam
        ch_bai = MINIMAP2_ALIGN_CLEAN.out.index

        // // разделяю картирования на регионы для коллинга и прочие регионы
        // SAMTOOLS_VIEW_SPLIT_ALIGNMENTS (
        //     ch_bam.combine(
        //         ch_bai, by: 0
        //     ).map {it -> [it[0], it[1], it[2]]},
        //     [[:], []],
        //     Channel.fromPath(params.amplicones_bed).collect(),
        //     "bai"
        // )

        // ch_mapped_for_P1 = SAMTOOLS_VIEW_SPLIT_ALIGNMENTS.out.unselected //
        // ch_mapped_for_calling = SAMTOOLS_VIEW_SPLIT_ALIGNMENTS.out.bam //
        // ch_mapped_for_calling_index = SAMTOOLS_VIEW_SPLIT_ALIGNMENTS.out.bai //

        // SEQKIT_SUBSEQ_MGPA (
        //     ch_reference.collect(),
        //     params.mgpA_bed
        // )

        // SEQKIT_REPLACE_MGPA (
        //     ch_mapped_for_P1.combine (
        //         SEQKIT_SUBSEQ_MGPA.out.fastx.map {it -> it[1]}
        //     ).map{
        //         it -> [it[0], it[2]]
        //     }
        // )

        // SAMTOOLS_FASTQ (
        //     ch_mapped_for_P1,
        //     false
        // )

        // ch_for_medaka = SAMTOOLS_FASTQ.out.other.combine(
        //         SEQKIT_REPLACE_MGPA.out.fastx, by: 0
        // )

        // ch_for_medaka.view()

        // MEDAKA_1 (
        //     ch_for_medaka
        // )

        // MEDAKA_2 (
        //     SAMTOOLS_FASTQ.out.other.combine (
        //         MEDAKA_1.out.assembly, by: 0
        //     )
        // )

        // объединяем картирование
        ch_mapping = ch_bam.combine(
            ch_bai, by: 0
        ).map {
            it -> [it[0], it[1], it[2], [], params.user_model, params.platform]
        }

        BEDTOOLS_GENOMECOV(
            ch_bam.map {
                it -> [it[0], it[1], 1]
            },
            [],
            "bed",
            true
        )

        SUBSET_COVERAGE_AMAP (
            BEDTOOLS_GENOMECOV.out.genomecov
        )

        CLAIR3(ch_mapping, ch_reference.collect(), ch_fai.collect())

        ch_clair3_vcf = CLAIR3.out.vcf

        ch_reannotation_pre = ch_clair3_vcf.combine (
            CLAIR3.out.tbi, by: 0
        ).map {
            it -> [it[0], it[1], it[2], [], []]
        }

        ch_renamer_pre = Channel.fromPath(params.path_rename_chr_pre)

        BCFTOOLS_ANNOTATE_PRE(
            ch_reannotation_pre,
            [],
            ch_renamer_pre.collect(),
        )

        if (params.snpeff_cache) {
                ch_snpeff_cache = Channel.fromPath(params.snpeff_cache).map {
                    it -> 
                        def cache_meta = [:]
                        cache_meta.id = "cache"
                        return [cache_meta, it]
        } 
        } else {
            SNPEFF_DOWNLOAD(
                Channel.value(params.snpeff_db_name).map{
                    it -> 
                        def cache_meta = [:]
                        cache_meta.id = "cache"
                        return [cache_meta, it]
                }
            )
            ch_snpeff_cache = SNPEFF_DOWNLOAD.out.cache
        }

        SNPEFF_SNPEFF(
            BCFTOOLS_ANNOTATE_PRE.out.vcf,
            params.snpeff_db_name,
            ch_snpeff_cache.collect()
        )

        BCFTOOLS_VIEW(
            SNPEFF_SNPEFF.out.vcf.map{
                it -> [it[0], it[1], []]
            },
            [], [], []
        )

        ch_renamer_post = Channel.fromPath(params.path_rename_chr_post)

        BCFTOOLS_ANNOTATE_POST (
            BCFTOOLS_VIEW.out.vcf.map{
                it -> [it[0], it[1], [], [], []]
            },
            [],
            ch_renamer_post.collect()
        )

        // ch_mutations_for_coverage = BCFTOOLS_ANNOTATE_POST.out.vcf.combine (
        //     SUBSET_COVERAGE_AMAP.out.highcov, by: 0
        // )


        // BCFTOOLS_VIEW_COVERED(
        //     ch_mutations_for_coverage.map {
        //         it -> [it[0], it[1], []]
        //     },
        //     [],
        //     ch_mutations_for_coverage.map {
        //         it -> [it[2]]
        //     }, 
        //     []
        // )

        BCFTOOLS_INDEX (
            BCFTOOLS_ANNOTATE_POST.out.vcf
        )
        
        ch_pre_consensus = BCFTOOLS_ANNOTATE_POST.out.vcf.combine (
            BCFTOOLS_INDEX.out.csi, by: 0
        ).combine (
            ch_reference.map { it -> it[1] }
        ).combine (
            SUBSET_COVERAGE_AMAP.out.lowcov, by: 0
        )

        BCFTOOLS_CONSENSUS(
            ch_pre_consensus
        )

        // BCFTOOLS_CONSENSUS.out.fasta.view()

        SEQKIT_SUBSEQ (
            BCFTOOLS_CONSENSUS.out.fasta,
            params.regions_of_interest
        )

        SEQKIT_REPLACE (
            SEQKIT_SUBSEQ.out.fastx
        )

        MEDAKA_1 (
            input_reads.combine(
                SEQKIT_REPLACE.out.fastx, by: 0
            )
        )

        MEDAKA_2 (
            input_reads.combine(
                MEDAKA_1.out.assembly, by: 0
            )
        )

        out_fasta = MEDAKA_2.out.assembly

    emit:
        out_fasta
}
