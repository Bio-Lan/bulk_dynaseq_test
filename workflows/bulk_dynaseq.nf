/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FILTER_GTF             } from '../modules/local/filter_gtf'
include { STAR_GENOME            } from '../modules/local/star_genome'
include { PROTOCOL_CMD           } from '../modules/local/protocol_cmd'
include { STARSOLO               } from '../modules/local/starsolo'
include { BAM_SPLIT              } from '../modules/local/bam_split'
include { CONVERSION             } from '../modules/local/conversion'
include { CONVERSION_SUMMARY     } from '../modules/local/conversion_summary'
include { SUBSTITUTION           } from '../modules/local/substitution'
include { QUANT                  } from '../modules/local/quant'
include { MULTIQC                } from '../modules/local/multiqc_sgr'
include { REPORT_SUMMARY         } from '../modules/local/report_summary'

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_bulk_dynaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BULK_DYNASEQ {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    
    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // fastqc
    if (params.run_fastqc) {
        FASTQC (
            ch_samplesheet
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    } 

    // STAR genome
    def star_genome = null
    genome_name = null
    genome_dir = "${workDir}/temp_dir/genome"
    if (params.star_genome) {
        star_genome = params.star_genome
    } else {
        FILTER_GTF(
            params.gtf,
            params.keep_attributes
        )
        
        ch_gtf = FILTER_GTF.out.filtered_gtf
        if(params.genome_name.contains('/')){
            genome_name = params.genome_name.split('/').last()
            genome_dir = params.genome_name
            genome_dir = genome_dir.replaceAll(genome_name+'$', "")
        }else{
            genome_name = params.genome_name
        }

        STAR_GENOME(
            params.fasta,
            ch_gtf,
            genome_name,
            genome_dir
        )
        ch_versions = ch_versions.mix(STAR_GENOME.out.versions.first())
        star_genome = STAR_GENOME.out.index
    }
    
    // create protocol star cmd
    PROTOCOL_CMD (
        ch_samplesheet,
        "${projectDir}/assets/",
        params.protocol,
    )
    ch_multiqc_files = ch_multiqc_files.mix(PROTOCOL_CMD.out.json.collect{it[1]})
    
    // starsolo
    ch_merge = ch_samplesheet.join(PROTOCOL_CMD.out.protocol_cmd.map{ [it[0], it[1].text] })
    STARSOLO (
        ch_merge,
        star_genome,
        "${projectDir}/assets/"
    )
    ch_versions = ch_versions.mix(STARSOLO.out.versions.first())

    // bam split
    BAM_SPLIT (
        STARSOLO.out.bam_sorted
    )
    ch_versions = ch_versions.mix(BAM_SPLIT.out.versions.first())

    // conversion
    split_bam = BAM_SPLIT.out.well_bam.transpose()
    CONVERSION(
        split_bam,
        params.gtf
    )

    ch_merge = CONVERSION.out.conv_postag.groupTuple().join(CONVERSION.out.conv_snp.groupTuple())
    // conversion summary
    CONVERSION_SUMMARY(
        ch_merge
    )

    // substitution
    ch_merge = CONVERSION_SUMMARY.out.conv_sample.join(CONVERSION.out.conv_wellbam.groupTuple())
    SUBSTITUTION(
        ch_merge
    )

    // quant
    ch_merge = CONVERSION_SUMMARY.out.conv_sample.join(CONVERSION.out.conv_wellbam.groupTuple()).join(CONVERSION.out.conv_wellbam_bai.groupTuple()).join(CONVERSION.out.conv_snp.groupTuple()).join(STARSOLO.out.raw_matrix)
    QUANT(
        ch_merge
    )

    // report summary
    ch_merge = STARSOLO.out.read_stats.join(STARSOLO.out.summary)
    ch_merge = ch_merge.join(SUBSTITUTION.out.substitution_stat)
    ch_merge = ch_merge.join(QUANT.out.sample_raw).join(QUANT.out.sample_filter)
    REPORT_SUMMARY(
        ch_merge
    )
    ch_multiqc_files = ch_multiqc_files.mix(REPORT_SUMMARY.out.json.collect{it[1]})

    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }
    
    // MODULE: MultiQC
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        "${projectDir}/multiqc_sgr/",
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList()
    versions       = ch_versions
}




