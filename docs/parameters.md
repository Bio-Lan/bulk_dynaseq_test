# singleron-RD/bulk_dynaseq pipeline parameters

Processing AccuraCode dynascope sequencing data

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-------------|------|---------|----------|--------|
| `input` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.</small></details> | `string` |  | True |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | True |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details> | `string` |  |  |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  |  |

## Genome

Genome files and parameters.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-------------|------|---------|----------|--------|
| `fasta` | Path to genome fasta. | `string` |  |  |  |
| `gtf` | Path to genome gtf. | `string` |  | True |  |
| `star_genome` | Path to STAR genome directory. Required if fasta and gtf are not provided. | `string` |  |  |  |
| `genome_name` | The generated STAR genome index will be saved under this folder. It can then be used for future pipeline runs, reducing processing times. | `string` | star_genome |  |  |
| `keep_attributes` | Attributes in gtf to keep. | `string` | gene_biotype=protein_coding,lncRNA,antisense,IG_LV_gene,IG_V_gene,IG_V_pseudogene,IG_D_gene,IG_J_gene,IG_J_pseudogene,IG_C_gene,IG_C_pseudogene,TR_V_gene,TR_V_pseudogene,TR_D_gene,TR_J_gene,TR_J_pseudogene,TR_C_gene; |  |  |
| `star_genome_additional_args` | Additional args to use when generate STAR genome directory. | `string` |  |  |  |

## Protocol options

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-------------|------|---------|----------|--------|
| `protocol` | Predefined pattern and whitelist. <details><summary>Help</summary><small> If set to "customized", --pattern and --whitelist are required. </small></details>| `string` | bulk_dynaseq-V1 |  |  |
| `well` | The AccuraCode dynascope wells used (384 or 96). | `integer` | 384 |  |  |
| `pattern` | A string to locate cell barcode and UMI in R1 read. For example "C9L16C9L16C9L1U12". <details><summary>Help</summary><small>C: cell barcode<br>L: Linker sequence between segments<br>U: UMI<br>T: poly T</small></details> | `string` |  |  |  |
| `whitelist` | Barcode whitelist files. Multiple whitelists are seperated by whitespace. | `string` |  |  |  |

## STARSolo options

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-------------|------|---------|----------|--------|
| `soloFeatures` | Quantification of different transcriptomic features. <details><summary>Help</summary><small>https://github.com/alexdobin/STAR/issues/1460  </small></details> | `string` | GeneFull_Ex50pAS |  |  |
| `soloCellFilter` | Cell-calling method. <details><summary>Help</summary><small>https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#cell-filtering-calling</small></details> | `string` | None |  |  |
| `soloStrand` |  strandedness of the solo libraries | `string` | Forward |  |  |
| `outFilterMatchNmin` | Alignment will be output only if the number of matched bases is higher than or equal to this value. <details><summary>Help</summary><small>Use default 50 to filter potential short prime sequences.</small></details> | `integer` | 50 |  |  |
| `outSAMattributes` | Output tags in SAM/BAM. <details><summary>Help</summary><small>https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#bam-tags</small></details> | `string` | MD NH HI AS nM CR UR CB UB GX GN sF |  |  |
| `outReadsUnmapped` |  output of unmapped and partially mapped (i.e. mapped only one mate of a paired end read) reads in separate file(s).(None or Fastx) | `string` | None |  |  |
| `starsolo_extra_args` | Extra STARSolo arguments to use. | `string` | --clip3pAdapterSeq AAAAAAAAAAAA --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 |  |  |

## Conversion options

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-------------|------|---------|----------|--------|
| `conversion_type` | Conversion type, TC for dynaseq. | `string` | TC |  |  |
| `basequalilty` | Min base quality of the read sequence. | `integer` | 20 |  |  |
| `snp_threshold` | Snp threshold filter, greater than snp_threshold will be recognized as snp. | `number` | 0.5 |  |  |
| `snp_min_depth` | Minimum depth to call a variant. | `integer` | 10 |  |  |

## Quant options

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-------------|------|---------|----------|--------|
| `umi_cutoff` | If the UMI number exceeds the threshold, it is considered a valid well and reported. | `integer` | 500 | | |
| `gene_cutoff` | If the gene number exceeds the threshold, it is considered a valid well and reported. | `integer` | 0 | | |
| `snp_file` | Backgroud snp file.Can be multiple files, separated by commas. If this option is set, it is valid for all wells. | `string` | | | |
| `snp_matchfile` | Backgroud snp file for each well.<details><summary>Help</summary><small> backgroud snp file for each well, one well per line, the format is  \"well,snp_file\". Can be multiple files, separated by commas. This parameter takes precedence over snp_file. </small></details>| `string` | | | |


## Optional modules

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-------------|------|---------|----------|--------|
| `run_fastqc` | FastQC of raw reads. | `boolean` |  |  |  |

## Optional pipeline

Split fastq based on the well provided.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-------------|------|---------|----------|--------|
| `run_splitfastq` | Split fastq based on information provided by the user. | `boolean` | false |  |  |
| `split_inf` | The split information file.It mus be full path of the file.<details><summary>Help</summary><small> The file has to be a comma-separated file with 3 columns, and a header row like:<br>`sample,sub_sample,wellBC`<br>  `sample` is the same as `sample` in the samplesheet.<br>  `sub_sample` is the prefix of fastq. <br>  `wellBC` is the full path of bc file. The format is one wellbc per line. </small></details>| `string` |  |  |  |
| `split_to_well` | Split fastq into well level. | `string` |  |  |  |

> [!NOTE]
> The path of `split_inf` must be full path. Relative path are not allowed.

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>| `integer` | 16 |  | True |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details>| `string` | 128.GB |  | True |
| `max_time` | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>| `string` | 240.h |  | True |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` |  |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `help` | Display help text. | `boolean` |  |  | True |
| `version` | Display version and exit. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>| `string` |  |  | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `hook_url` | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.</small></details>| `string` |  |  | True |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |
| `multiqc_logo` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file | `string` |  |  | True |
| `multiqc_methods_description` | Custom MultiQC yaml file containing HTML including a methods description. | `string` |  |  | True |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `validationShowHiddenParams` | Show all params when using `--help` <details><summary>Help</summary><small>By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters.</small></details>| `boolean` |  |  | True |
| `validationFailUnrecognisedParams` | Validation of parameters fails when an unrecognised parameter is found. <details><summary>Help</summary><small>By default, when an unrecognised parameter is found, it returns a warinig.</small></details>| `boolean` |  |  | True |
| `validationLenientMode` | Validation of parameters in lenient more. <details><summary>Help</summary><small>Allows string values that are parseable as numbers or booleans. For further information see [JSONSchema docs](https://github.com/everit-org/json-schema#lenient-mode).</small></details>| `boolean` |  |  | True |
