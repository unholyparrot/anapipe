include { CAT_FASTQ } from '../modules/nf-core/cat/fastq/main'


workflow INPUT {
    main:
        if (params.input_reads) {
            if (params.input_worklist) {
                if (file(params.input_worklist).getExtension() == "csv") {
                    ch_input_rows = Channel.fromPath(params.input_worklist).splitCsv(header: true).map {
                        row -> 
                            def alias = row.alias
                            def barcode = row.barcode
                            return [alias, barcode]
                    }

                    ch_raw_reads = ch_input_rows.map {
                        id, barcode -> {
                            def meta = [:]
                            meta.id = id
                            meta.single_end = true
                            def files = file("${params.input_reads}/${barcode}/*.fastq.gz")
                            return [meta, files]
                        }
                    }
                } else {
                    exit 1, "Unknown input file format: ${file(params.input_worklist).getExtension()}"
                }
            } else {
                // Assuming input is a directory with barcodeNN subfolders
                def readsDir = file(params.input_reads)
                if (!readsDir.exists() || !readsDir.isDirectory()) {
                    exit 1, "Error: Input reads directory does not exist or is not a directory: ${params.input_reads}"
                }

                ch_raw_reads = Channel.fromPath("${params.input_reads}/barcode*", type: 'dir').map {
                    folder ->
                        def meta = [:]
                        meta.id = folder.getName()
                        meta.single_end = true 
                        def reads = folder.listFiles()
                        return [meta, reads]
                }

            }

        CAT_FASTQ(ch_raw_reads)

        ch_input_reads = CAT_FASTQ.out.reads.filter{
            it -> it[1].countFastq() > params.min_reads_per_sample
        }

        } else {
                exit 1, "Directory with input reads is required."
        }


    emit:
        ch_input_reads
}
