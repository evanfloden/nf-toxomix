#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                         NF-toxomix
========================================================================================
 NF-toxomix Analysis Pipeline. Started 2018-02-15.
 #### Homepage / Documentation
 https://github.com/evanfloden/nf-toxomix
 #### Authors
 Evan Floden (evanfloden) <evan.floden@gmail.com> - https://github.com/evanfloden>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     NF-toxomix v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run skptic/NF-toxomix --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Hardware config to use. docker / aws

    Options:
      --singleEnd                   Specifies that the input is single end reads

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}


/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = '0.1.0'

// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.transcriptomics_data = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE28nnn/GSE28878/matrix/GSE28878_series_matrix.txt.gz"
params.compound_info_excel = "$baseDir/data/Supplementary_Data_1.xls"
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.reads = 'data/*_{1,2}.fq'
params.fasta = "data/l1000_transcripts.fa"
params.outdir = './results'
params.email = false
params.plaintext_email = false
multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Validate inputs
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}

Channel.from(params.transcriptomics_data)
  .set {transcriptomics_data_url_ch}

Channel.fromPath(params.compound_info_excel)
  .set { compound_info_excel_ch}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * Create a channel for input read files
 */
Channel
    .fromFilePairs( params.reads, size: -1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_fastqc; read_files_quantify }


// Header log info
log.info "========================================="
log.info " NF-toxomix v${version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container']    = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}


/*
 * Download transcriptomics Data
 */
process get_transcriptomics_data {

    input:
        val(transcriptomics_data_url) from transcriptomics_data_url_ch

    output:
        file('transcriptomics_data.txt') into transcriptomics_data_ch
        file('transcriptomics_data_raw') into transcriptomics_data_raw_ch

    shell:
        """
        curl -X GET "${transcriptomics_data_url}" > transcriptomics_data_raw.gz
        gunzip transcriptomics_data_raw.gz
        awk -f ${baseDir}/bin/parse_transcript_data.awk  transcriptomics_data_raw|tr -d '"' > transcriptomics_data.txt
        """
}

/*
 * Process the comound info execel sheet in R with pandas
 */
process process_compound_info {

    input:
        file(compound_info_file) from compound_info_excel_ch

    output:
        file('compound_info.tsv') into compound_info_ch

    script:
    """
    #!/opt/conda/bin/python

    import pandas as pd
    print("input  is:" + str("${compound_info_file}"))
    print("output is:" + str("compound_info.tsv"))
    excel_file = pd.read_excel(io="${compound_info_file}", encoding='utf-16')
    excel_file.to_csv(path_or_buf="compound_info.tsv",  encoding='utf-16', sep="\t")
    """
  }

compound_info_ch
    .into { compound_info_ch1; compound_info_ch2 }

/*
 * Create the compound training data
 */
process training_compound_info {

      input:
      file(compound_info) from compound_info_ch1

      output:
      file("training_data_compound_info.tsv") into compound_info_training

      shell:
      """
      a=\$(tempfile -d .)
      cut -f1,10- ${compound_info} |awk 'BEGIN{{print("compound\\tgenotoxicity");}};NR>3'|head -35 > \$a;
      sed -re 's/\\+\$/GTX/g; s/\\-\$/NGTX/g; s/-//g; s/\\]//g; s/\\[//g' \$a > training_data_compound_info.tsv
      """
  }

process validation_compound_info {

     input:
     file(compound_info) from compound_info_ch2

     output:
     file("validation_data_compound_info.tsv") into compound_info_validation_ch

     shell:
     """
     a=\$(tempfile -d .)
     cut -f1,10- ${compound_info} |awk 'NR>41'|sed -re 's/\\+\$/GTX/g; s/\\-\$/NGTX/g;  s/[[:punct:]]//g' > \$a;
     awk 'BEGIN{{print("compound\\tgenotoxicity");}}{{print}}' \$a| sed -re 's/ppDDT\\t/DDT\\t/g; s/\\s+/\\t/g'> validation_data_compound_info.tsv
     """
}

/*
 *  Each series has different solvent, match to the correct solvent
 */
process map_sovent_to_exposure {

    input:
    file(transcriptomics_data_raw) from transcriptomics_data_raw_ch

    output:
    file("solvent2exposure.tsv") into solvent_to_exposure_ch

    shell:
    """
    a=\$(tempfile -d .)
    b=\$(tempfile -d .)
    c=\$(tempfile -d .)
    d=\$(tempfile -d .)
    e=\$(tempfile -d .)
    paste <(grep Sample_title ${transcriptomics_data_raw}|cut -f2-|tr -d '"'|tr '\\t' '\\n') <(grep Series_sample_id ${transcriptomics_data_raw} > \$a;
    cut -f2- \$a|tr -d '"'|sed -re 's/\\s*\$//'|tr ' ' '\\n')|grep 24h > \$b;
    sed -re 's/^Serie\\s*//g; s/, HepG2 exposed to\\s*/\\t/g; s/for 24h, biological rep\\s*/\\t/g' \$b > \$c;
    awk 'BEGIN{{print("series_id\\tcompound\\treplicate\\tarray_name");}}{{print}}' \$c|sed -re 's/\\s+/\\t/g' > \$d;
    sed -re 's/DEPH/DEHP/g; s/Ethyl\\t/EtAc\\t/g; s/NPD\\t/NDP\\t/g; s/Paracres\\t/pCres\\t/g; s/Phenol\\t/Ph\\t/g; s/Resor/RR/g' \$d > \$e;
    sed -re 's/2-Cl\\t/2Cl\\t/g' \$e> solvent2exposure.tsv
    """
}


/*
 * Create a file that stores a mapping between genotoxicity, compound and array information for validation set
 */
process map_compound_to_array_validation {

    input:
    file(validation_data_compound_info) from compound_info_validation_ch
    file(solvent2exposure) from solvent_to_exposure_ch

    output:
    file("compound_array_genotoxicity_val.tsv") into compound_array_genotoxicity_ch

    shell:
    """
    a=\$(tempfile -d .)
    echo -e "series_id\\tcompound\\treplicate\\tarray_name\\tgenotoxicity" > compound_array_genotoxicity_val.tsv;
    LANG=en_EN join -i -o 2.1,2.2,2.3,2.4,1.2 -t \$'\\t' -1 1 -2 2 <(cat ${validation_data_compound_info} |sort -k1.1i,1.3i  -t \$'\\t' ) <(LANG=en_EN sort -fbi -t \$'\\t' -k 2 ${solvent2exposure}) > \$a;
    grep -iv genotoxicity \$a >> compound_array_genotoxicity_val.tsv;
    """
}



/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[NF-toxomix] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[NF-toxomix] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['software_versions'] = software_versions
    email_fields['software_versions']['Nextflow Build'] = workflow.nextflow.build
    email_fields['software_versions']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[NF-toxomix] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[NF-toxomix] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[NF-toxomix] Pipeline Complete"

}
