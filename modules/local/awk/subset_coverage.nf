process SUBSET_COVERAGE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.lowcov.bed"), emit: lowcov
    tuple val(meta), path("*.highcov.bed"), emit: highcov

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk -v mincov="${params.region_min_cov}" '(\$4 < mincov)' ${bed} | awk -v OFS="\\t" '\$1=\$1' > ${meta.id}.lowcov.bed
    awk -v mincov="${params.region_min_cov}" '(\$4 >= mincov)' ${bed} | awk -v OFS="\\t" '\$1=\$1' > ${meta.id}.highcov.bed
    """
}
