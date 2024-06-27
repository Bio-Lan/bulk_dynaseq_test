## Genome

Genome files and parameters.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `fasta` | Path to genome fasta. | `string` |  |  |  |
| `gtf` | Path to genome gtf. | `string` |  |  |  |
| `star_genome` | Path to STAR genome directory. Required if fasta and gtf are not provided. | `string` |  |  |  |
| `genome_name` | The generated STAR genome index will be saved under this folder. It can then be used for future pipeline runs, reducing processing times. | `string` | star_genome |  |  |
| `keep_attributes` | Attributes in gtf to keep. | `string` | gene_biotype=protein_coding,lncRNA,antisense,IG_LV_gene,IG_V_gene,IG_V_pseudogene,IG_D_gene,IG_J_gene,IG_J_pseudogene,IG_C_gene,IG_C_pseudogene,TR_V_gene,TR_V_pseudogene,TR_D_gene,TR_J_gene,TR_J_pseudogene,TR_C_gene; |  |  |
| `star_genome_additional_args` | Additional args to use when generate STAR genome directory. | `string` |  |  |  |

## Protocol options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `protocol` | Predefined pattern and whitelist.If set to "customized", --pattern and --whitelist are required.| `string` | bulk_rna |  |  |
| `well` | The AccuraCode wells used (384 or 96).| `integer` | 384 |  |  |
| `pattern` | A string to locate cell barcode and UMI in R1 read. For example "C9L16U12". <details><summary>Help</summary><small>C: cell barcode<br>L: Linker sequence between segments<br>U: UMI<br>T: poly T</small></details>| `string` |  |  |  |
| `whitelist` | Barcode whitelist files. Multiple whitelists are seperated by whitespace. | `string` |  |  |  |
