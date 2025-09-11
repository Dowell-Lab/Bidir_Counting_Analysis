#!/usr/bin/env nextflow

/*
 * =======================================================================
 * 
 *		Bidirectional Counting Pipeline
 * 		
 * =======================================================================
 * 
 * This Source Code is distributed under the GPL-3.0 License
 */


/* 
 * 'Bidir-Counting-Flow' - A Nextflow pipeline for counting over regions of transcription
 * from Nascent RNA sequencing  
 *
 * The pipeline addresses overlapping transcription when counting bidirectionals and genes
 * 
 * =============
 * Authors
 * =============
 * Hope A. Townsend : hope.townsend@colorado.edu
 */
// HELP MESSAGE
def helpMessage() {
    log.info"""
    =========================================
     BidirCountingFlow v${params.version}
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --cons_file '/project/regions.bed' --bams '/project/bams/ --workdir '/project/tempfiles' --outdir '/project/'
    Required arguments:
        --cons_file                   bed file with consensus bidirectionals/tREs of interest where the midpoint of the 2nd & 3rd columns is closest to mu (often with muMerge)
        --mmfiltbams                  Directory pattern for multimap filtered bams: /project/*.mmfilt.sorted.bam (Required if --bams or --crams not specified).
        --crams                       Directory pattern for cram files: /project/crams/ (Required if --mmfiltbams or --bams not specified).
        --bams                        Directory pattern for bam files: /project/bams/ (Required if --mmfiltbams or --crams not specified).
        --workdir                     Nextflow working directory where all intermediate files are saved.

    Intermediate Files:
        --geneputcounts                 Bedtools coverage counts of the genes

    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.
        --prefix                       Prefix used in output files (default = "count_project")  
        --savemmfiltbams               Saves multimapped read filtered bamfiles used for tfit (default 'TRUE')
        --date                         added to prefix in output files
        --tfea                         If you want the TSS and NonTSS bidirectionals separated and in a format to be provided to TFEA as a combined_file, this should be "TRUE" (default="TRUE").

    Analysis Options:
        --make_names                   If the cons_file does NOT have a column 4 with the names, a unique name will be added for each feature (default="YES")
        --count_win                    Window added from mu to consider for initial counting of bidirectionals (so if 1000 --> 2000bp region) (default=1000)
        --tss_win                      Window added from mu to overlap bidirectionals with gene TSS windows to identify 5' gene bidirectionals (default=25)
        --count_limit_genes            Minimum coverage (fraction of gene covered by reads) required to consider a gene transcribed and possibly convoluting transcription (default=0.7)
        --count_limit_bids             Minimum total counts summed from all samples required to consider a bidirectional possibly convoluting gene transcription (default=30)
        --get_fixed_genecounts         Whether or not you want to get counts for genes after removing strongly transcribed bidirectionals overlapping genes


    Files to Use (with Defaults):
        --genome                       Only needed if --crams is used (default is /scratch/Shares/dowell/genomes/hg38/hg38.fa)
        --tss_1kb_file                 sorted bed file with 500bp region +/- around TSS of all genes (including putatives) (Default in assets of github repo)
        --gene_put_file                sorted bed file of all gene isoforms (Default in assets of github repo)
        --gene_put_10kbdntm_file       sorted bed file of the termination site of gene isoforms to 10kb downstream (Default in assets of github repo)
        --gene_count_file              sorted bed file of 5prime truncated genes over which you want to count
        --gene_order_file              sizes of chromosomes in the order of the gene_put_file (default is in assets of github repo: choromosome sizes)
    """.stripIndent()
        
}

nextflow.enable.dsl=1

// process SetupEnvironment {
//     executor 'slurm'
//     time '1h'
//     memory '5 GB'
//     cpus 1
    
//     script:
//     """
//     # source activate /Users/hoto7260/miniconda3/envs/Rjupyter
//     module load bedtools
//     module load samtools/1.8
//     module load subread/1.6.2
//     module load R/4.4.0
//     """
// }

// Configure Variables
// Parameters that a user can edit on the command line with --
// Inputs
params.cons_file = ''
params.mmfiltbams = ''
params.crams = ''
params.bams = ''
params.workdir = './tempfiles'

// Saving
params.outdir = './'
params.prefix = 'count_project'
params.savemmfiltbams = 'TRUE'
params.date = ''
params.tfea = 'TRUE'
params.get_fixed_genecounts = 'TRUE'

// Parameters
params.make_names = 'TRUE'
params.count_win = 1000
params.tss_win = 25
params.count_limit_genes = 0.7
params.count_limit_bids = 30

// Files
//Get the folder for the R and python scripts


params.genome = '/scratch/Shares/dowell/genomes/hg38/hg38.fa'
params.tss_1kb_file = "${workflow.projectDir}/assets/hg38_TSS1kb_refseq_diff53prime_with_putatives.bed"
params.gene_put_file = "${workflow.projectDir}/assets/hg38_refseq_diff53prime_with_putatives_fixnames_sort2.sorted.bed"
params.gene_put_10kbdntm_file = "${workflow.projectDir}/assets/hg38_refseq_diff53prime_10kbdwnstm_with_putatives_fixnames.sorted.bed"
params.gene_count_file = "${workflow.projectDir}/assets/hg38_refseq_diff53prime_5ptrunc_counting.sorted.bed"
params.geneputcounts = ''
params.gene_order_file = "${workflow.projectDir}/assets/hg38.chrom.sizes"


params.version = '0.1'

import java.text.SimpleDateFormat
def date = new java.util.Date()
def sdf = new SimpleDateFormat("yyMMdd")
output_date =  sdf.format(date)
String output_date = new java.text.SimpleDateFormat("yyMMdd").format(new Date())


// Header log info
log.info """=======================================================
BidirCountingFlow v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'BidirCounting-Flow'
summary['Help Message']     = params.help
summary['Pipeline Version'] = params.version
summary['Run Name']         = workflow.runName
if(params.crams) summary['Crams']            = params.crams
if(params.bams) summary['Bams']              = params.bams
if(params.mmfiltbams) summary['MMfiltBams']              = params.mmfiltbams
summary['Consensus Region File']       = params.cons_file
summary['Window for Regions in Counting']        = params.count_win 
summary['Window for Regions in TSS Overlaps']        = params.tss_win 
summary['Required coverage of gene for convolution']        = params.count_limit_genes 
summary['Required total counts of Bid for convolution']        = params.count_limit_bids 
summary['Genome used']        = params.genome 
summary['Gene File for TSS overlap']        = params.tss_1kb_file 
summary['Bed file of putative genes']        = params.gene_put_file 
summary['Bed file of putative genes 10kb downstream']        = params.gene_put_10kbdntm_file 
summary['Bed file of 5p truncated Genes for counting']        = params.gene_count_file 
summary['Working Directory']        = workflow.workDir
summary['Output Directory']        = params.outdir 
summary['Prefix for Saving']        = params.prefix 
summary['Date for Saving']        = params.date 
summary['Saving MMFILTBAMS?']        = params.savemmfiltbams 
summary['Saving TFEA input?']        = params.tfea 
summary['Getting fixed gene counts?']        = params.get_fixed_genecounts 
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Script dir']       = workflow.projectDir
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "======================================================="


println "\nSTARTING PIPELINE"

// PART 0: Outputting software versions

println "[Log 0]: Getting software versions"

// process get_software_versions {
//     time '1h'

//     output:
//     stdout into software_versions

//     script:
//     """
//     #!/bin/bash
//     printf "flow_version: %s\n" ${params.version}
//     printf "nextflow_version: %s\n" ${workflow.nextflow.version}
//     printf "samtools_version: %s\n" \$(samtools --version | head -1 | awk '{print \$NF}')
//     printf "bedtools_version: %s\n" \$(bedtools --version | head -1 | awk -F " v" '{print \$2}')
//     printf "openmpi_version: %s\n" \$(ompi_info | head -2 | tail -1 | awk '{print \$NF}')
//     printf "gcc_version: %s\n" \$(gcc --version | head -1 | awk '{print \$NF}')
//     printf "r_version: %s\n" \$(R --version | head -1 | awk '{print \$3}')
//     printf "pipeline_hash: %s\n" ${workflow.scriptId}
//     """
// }

// software_versions.collectFile(name: "software_versions_bidircount_${output_date}_${workflow.runName}.yaml", storeDir: "${params.outdir}/pipeline_info")

println "[Log 0]: Software versions complete"


// 1. Get the MMFILT BAMS either from bams, crams, or when provided
if (params.crams) {
    println "[Log 1]: Converting CRAM files to Multi-mapped filtered BAM FILES"
    println "[Log 1]: Genome file being used ..... $params.genome "
    println "[Log 1]: Cram file directory ........ $params.crams"
    println "[Log 1]: Working directory ... $workflow.workDir"
    println "[Log 1]: Output directory ... $params.outdir"
    

    

    cramfiles = Channel
                  .fromPath("${params.crams}/*.sorted.cram")
                  .map { file -> tuple(file.baseName.replace('.sorted', ''), file) }


    //cramfiles.view()

    process cram_to_mmfiltbam {
      cpus 16
      queue 'short'
      memory '5 GB'
      time '1h30m'
      tag "$prefix"

      publishDir "${params.outdir}" , mode: 'copy',
      saveAs: {filename ->
              if ((filename == "${prefix}.mmfilt.sorted.bam") & (params.savemmfiltbams == 'TRUE'))    "mmfiltbams/${prefix}.mmfilt.sorted.bam"
              else null
             }
      input:
        tuple val(prefix), file(cram) from cramfiles

      output:
        tuple val(prefix), file("${prefix}.mmfilt.sorted.bam"), file("${prefix}.mmfilt.sorted.bam.bai") into sorted_mmfilt_bam_file_genecount
        tuple val(prefix), file("${prefix}.mmfilt.sorted.bam"), file("${prefix}.mmfilt.sorted.bam.bai") into sorted_mmfilt_bam_file_bidcount
        tuple val(prefix), file("${prefix}.mmfilt.sorted.bam"), file("${prefix}.mmfilt.sorted.bam.bai") into sorted_mmfilt_bam_file_mucount
        tuple val(prefix), file("${prefix}.mmfilt.sorted.bam"), file("${prefix}.mmfilt.sorted.bam.bai") into sorted_mmfilt_bam_file_fingenecount

      script:
      """
      samtools view -@ 16 -b -1 -T ${params.genome} ${cram} | \
            samtools view -@ 16 -h -q 1 | \
            grep -P '(NH:i:1|^@)' | \
            samtools view -h -b > ${prefix}.mmfilt.sorted.bam
      samtools index ${prefix}.mmfilt.sorted.bam ${prefix}.mmfilt.sorted.bam.bai
      """
    }


    
} else if (params.bams) {
    println "[Log 1]: Converting sorted BAM files to Multi-mapped filtered BAM FILES"
    println "[Log 1]: Bam file directory ........ $params.bams"
    println "[Log 1]: Working directory ... $workflow.workDir"
    println "[Log 1]: Output directory ... $params.outdir"

    cramfiles = Channel
                  .fromPath("${params.bams}/*.sorted.bam")
                  .map { file -> tuple(file.baseName.replace('.sorted', ''), file) }

    process bam_to_mmfiltbam {
        cpus 16
        queue 'short'
        memory '5 GB'
        time '1h30m'
        tag "$prefix"

        publishDir "${params.outdir}" , mode: 'copy',
        saveAs: {filename ->
              if ((filename == "${prefix}.mmfilt.sorted.bam") & (params.savemmfiltbams == 'TRUE'))    "mmfiltbams/${prefix}.mmfilt.sorted.bam"
              else null
             }
        input:
        tuple val(prefix), file(bam) from bamfiles

        output:
        tuple val(prefix), file("${prefix}.mmfilt.sorted.bam"), file("${prefix}.mmfilt.sorted.bam.bai") into 
        sorted_mmfilt_bam_file_genecount
        tuple val(prefix), file("${prefix}.mmfilt.sorted.bam"), file("${prefix}.mmfilt.sorted.bam.bai") into 
        sorted_mmfilt_bam_file_bidcount
        tuple val(prefix), file("${prefix}.mmfilt.sorted.bam"), file("${prefix}.mmfilt.sorted.bam.bai") into 
        sorted_mmfilt_bam_file_mucount
        tuple val(prefix), file("${prefix}.mmfilt.sorted.bam"), file("${prefix}.mmfilt.sorted.bam.bai") into 
        sorted_mmfilt_bam_file_fingenecount

        script:
        """
        samtools view -@ 16 -h -q 1 ${bam} | \
               grep -P '(NH:i:1|^@)' | \
               samtools view -h -b > ${prefix}.mmfilt.sorted.bam
        samtools index ${prefix}.mmfilt.sorted.bam ${prefix}.mmfilt.sorted.bam.bai
      """
    }

} else {
    println "[Log 1]: Using the multi-mapped filtered bams from .... $params.mmfiltbams"
    // If mmfilt bams provided instead, use those
    sorted_mmfilt_bam_file_genecount = Channel .fromPath(params.mmfiltbams) .map { file -> tuple((file.simpleName + '.mmfilt.sorted'), file, (file + '.bai'))}
    sorted_mmfilt_bam_file_bidcount = Channel .fromPath(params.mmfiltbams) .map { file -> tuple((file.simpleName + '.mmfilt.sorted'), file, (file + '.bai'))}
    sorted_mmfilt_bam_file_mucount = Channel .fromPath(params.mmfiltbams) .map { file -> tuple((file.simpleName + '.mmfilt.sorted'), file, (file + '.bai'))}
    sorted_mmfilt_bam_file_fingenecount = Channel .fromPath(params.mmfiltbams) .map { file -> tuple((file.simpleName + '.mmfilt.sorted'), file, (file + '.bai'))}
}

println "[Log 1]: Multi-mapped filtered bams are ready\n"

// 2. Get the CONSENSUS REGION FILES FOR ANALYSIS


println "[Log 2]: Getting consensus regions for overlapping and downstream analysis"
println "[Log 2]: Window used for counting ........ $params.count_win"
println "[Log 2]: Window used for identifying Gene 5' bidirectionals ... $params.tss_win"
println "[Log 2]: These files will be saved in ... $params.outdir/regions"

if (!params.cons_file) {
    exit 1, "Error: --cons_file parameter is required"
}


process GetConsFiles {
  cpus 1
  queue 'short'
  memory '5 GB'
  time '1h'
  tag "$prefix"
    
  
  output:
  file ("*win_count.sorted.bed") into count_win_file
  file ("*win_count.sorted.bed") into count_win_file_getbids
  file ("*win_TSS.sorted.bed") into tss_win_file

  script:
  """
  #!/usr/bin/env Rscript
  
  library(data.table)
  library(stringr)
  ### Read in the consensus region
  cons <- data.frame(fread("${params.cons_file}"))

  # Get the names for each bidir
  if (${params.make_names} != FALSE) {
      cons\$unique_name = paste0(cons\$V1, ":", cons\$V2, "-", cons\$V3)
      } else {
        if (nrow(cons) > 3) {
          cons\$unique_name = cons\$V4
          } else {
            print("WARNING: you don't currently have names available in the Consensus file so we are making some for you.")
            cons\$unique_name = paste0(cons\$V1, ":", cons\$V2, "-", cons\$V3) 
          }
      }
  
  ###############
  ## FILTERING ##
  ###############
  cat("\nFILTERING Poor Confidence Regions=============================")

  ## remove any tfit calls longer than 3.5kb (indicates very low confidence)
  cons\$length <- cons\$V3-cons\$V2
  cons <- cons[cons\$length<3500,]
  cat("\nNumber regions saved after removing calls >3.5kb -- not trustworthy", length(unique(cons\$unique)), length(unique(cons\$unique_dreg)))

  ###############
  ## Save according to windows ##
  ###############
  # get mu
  cons\$mu <- as.integer((cons\$V2+cons\$V3)/2)
  count_win = as.integer(${params.count_win})
  tss_win = as.integer(${params.tss_win})
  cons\$count_start <- cons\$mu - count_win
  cons\$count_stop <- cons\$mu + count_win
  cons\$tss_start <- cons\$mu - tss_win
  cons\$tss_stop <- cons\$mu + tss_win

  write.table(cons[,c("V1", "count_start", "count_stop", "unique_name")], 
                  paste0("${params.prefix}", "_MUMERGE_", count_win, "win_count.sorted.bed"), 
              quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
  write.table(cons[,c("V1", "tss_start", "tss_stop", "unique_name")], 
                  paste0("${params.prefix}", "_MUMERGE_", tss_win, "win_TSS.sorted.bed"), 
              quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)  
  cat("\nFinished writing files")
  """
}

// 2. Get overlaps that would introduce possible colliding transcription
println "[Log 3]: Getting overlaps of regions hinting at colliding transcription"
process GetOverlaps {

  cpus 1
  queue 'short'
  memory '5 GB'
  time '30m'
  tag "$prefix"

  input:
  file count_win_file
  file tss_win_file
    
  output:
  file ("TSS_OUT.bed") into tss_out_getbids
  file ("COUNT_OUT.bed") into count_out_getbids
  file ("COUNT_OUT.bed") into count_out_fix
  file ("DWN_OUT.bed") into dwn_out_getbids
  file ("DWN_OUT.bed") into dwn_out_fix
  file ("COUNT_TRUNC_OUT.bed") into count_trunc_out_getbids
  file ("COUNT_TRUNC_OUT.bed") into count_trunc_out_getgenes
  file ("CLOSE_OUT.bed") into close_out_getbids
    
  script:
  """
  #!/bin/bash
  bedtools intersect -wo -a ${tss_win_file} -b ${params.tss_1kb_file} > TSS_OUT.bed
  wc -l TSS_OUT.bed
  bedtools intersect -wo -a ${params.gene_put_file} -b ${count_win_file} > COUNT_OUT.bed
  bedtools intersect -wo -a ${params.gene_put_10kbdntm_file} -b ${count_win_file} > DWN_OUT.bed
  bedtools intersect -wo -a ${params.gene_count_file} -b ${count_win_file} > COUNT_TRUNC_OUT.bed
  ## Bids possibly overlapping with one another (-N means can't have same name)
  bedtools closest -D ref -k 3 -N -a ${tss_win_file} -b ${tss_win_file} > CLOSE_OUT.bed

  if [ \$(wc -l < TSS_OUT.bed) -le 1 ]; then
        echo "Error: TSS_OUT.bed} has ≤ 1 line. Stopping pipeline." 
        wc -l ${tss_win_file}
        wc -l ${params.tss_1kb_file}
        exit 1
 fi
  """
}

println "[Log 3]: Getting transcription levels of genes to roughly establish which are transcribed"
// Get the preliminary coverage of putative genes to roughly establish which ones are transcribed
if (params.geneputcounts) {
    Channel.fromPath(params.geneputcounts)
        .into { put_gene_counts_getbids; put_gene_counts_fix; put_gene_counts_getgenes }
} else {
process CountPutGenes {

  cpus 1
  queue 'short'
  memory '20 GB'
  time '5h'
  tag "$prefix"

  input:
  // only runs once all mmfilt bams are ready
  file ('*') from sorted_mmfilt_bam_file_genecount.collect()
    
  output:
  file ("put_gene_counts.txt") into put_gene_counts_getbids
  file ("put_gene_counts.txt") into put_gene_counts_fix
  file ("put_gene_counts.txt") into put_gene_counts_getgenes
    
  script:
  """
  # Get the genome file to use (order of chromosomes with the sizes in tab-delimited format)
  #cut --fields=1 ${params.gene_put_file} | uniq > chrom_order.txt
  #bedtools coverage -s -sorted -g chrom_order.txt -a ${params.gene_put_file} -b ./*mmfilt.sorted.bam > put_gene_counts.txt
    if [ -d "mmfiltbams" ]; then
        bedtools coverage -s -sorted -g ${params.gene_order_file} -a ${params.gene_put_file} -b ./mmfiltbams/*.bam > put_gene_counts.txt
    else
        bedtools coverage -s -sorted -g ${params.gene_order_file} -a ${params.gene_put_file} -b ./*mmfilt.sorted.bam > put_gene_counts.txt
    fi
  
  """
}
}

// PART 4: Filter Bidirectionals, Get Regions (TSS and tREs)
println "[Log 4]: Getting TSS bidirectionals and Filtered bidirectionals based on transcriptional convolution"
println "[Log 4]: Using coverage limit of genes .... ${params.count_limit_genes}"


// PART 5: GET THE FIXED BID COUNTS
process GetTSSBidsAndGTFs {
  
  cpus 1
  queue 'short'
  memory '5 GB'
  time '10h'
  tag "$prefix"

  publishDir "${params.outdir}/regions/", mode: 'copy',
    saveAs: { filename ->
        if (filename == "tss_bid.txt") { 
            return "tss_bid_${params.prefix}_${params.date}.txt" 
        } else if (filename == "tss_bid_forTFEA.bed") { 
            return "tss_bid_${params.prefix}_${params.date}_forTFEA.bed" 
        } else if (filename == "nontss_bid_forTFEA.bed") { 
            return "nontss_bid_${params.prefix}_${params.date}_forTFEA.bed" 
        } else { 
            return null 
        }
    }

  input:
  file tss_out_getbids
  file count_out_getbids
  file dwn_out_getbids
  file count_trunc_out_getbids
  file close_out_getbids
  file count_win_file_getbids
  file put_gene_counts_getbids

  output:
  file ("filt.saf") into bidirs_saf
  file ("pos.gtf") into pos_gtf
  file ("neg.gtf") into neg_gtf
  file ("uns.gtf") into uns_gtf
  file ("tss_bid.txt") into tss_file
  file ("tss_bid_forTFEA.bed") into tss_TFEA
  file ("nontss_bid_forTFEA.bed") into nontss_TFEA


  script:
  """
  #!/usr/bin/env Rscript
  library(data.table)
  library(stringr)

  ###########################
#####    FUNCTIONS    #####s
###########################

get_info_for_calls <- function(close_df, fixed_length=1000) {
    # This function gets the information needed in the overlaps file to get the
    #     appropriate leftmost and rightmost calls
    # PARAMETERS:
    #     close_df: a dataframe with the columns chrom, start, stop, name of region 1, the same for region 2, and finally the distance between the two regions. It is assumed that the middle of the start and stop of each region (mu) is the starting point of each transcript.
    #     fixed_length: Fixed length distance from mu to consider the "original" length of the transcripts. This fixed_length + 100 is also the max distance the mus of two features can be from one another to consider potential overlap (since the count region wouldn't get farther anyways)
    diff_call = as.integer(fixed_length)+100
    close_df\$mu1 <- as.integer((close_df\$V2+close_df\$V3)/2)
    close_df\$mu2 <- as.integer((close_df\$V6+close_df\$V7)/2)
    close_df\$MUDIFF <- close_df\$mu1-close_df\$mu2
    
    # get the appropriate distances
    close_df\$V2 = close_df\$mu1 - fixed_length
    close_df\$V3 = close_df\$mu1 + fixed_length
    close_df\$V6 = close_df\$mu2 - fixed_length
    close_df\$V7 = close_df\$mu2 + fixed_length
    orig_close_df = close_df
    close_df <- close_df[abs(close_df\$MUDIFF) < diff_call,]
    
    # get leftmost and rightmost
    close_df\$leftmost_Bidir <- "NA"; close_df\$leftmost_left <- 0; close_df\$leftmost_mu <- 0; close_df\$leftmost_right <- 0
    close_df\$rightmost_Bidir <- "NA"; close_df\$rightmost_right <- 0; close_df\$rightmost_mu <- 0; close_df\$rightmost_left <- 0
    # where first is leftmost one
    close_df[close_df\$MUDIFF < 0,]\$leftmost_Bidir <- close_df[close_df\$MUDIFF < 0,]\$V4
    close_df[close_df\$MUDIFF < 0,]\$leftmost_left <- close_df[close_df\$MUDIFF < 0,]\$V2
    close_df[close_df\$MUDIFF < 0,]\$leftmost_mu <- close_df[close_df\$MUDIFF < 0,]\$mu1
    close_df[close_df\$MUDIFF < 0,]\$rightmost_Bidir <- close_df[close_df\$MUDIFF < 0,]\$V8
    close_df[close_df\$MUDIFF < 0,]\$rightmost_right <- close_df[close_df\$MUDIFF < 0,]\$V7
    close_df[close_df\$MUDIFF < 0,]\$rightmost_mu <- close_df[close_df\$MUDIFF < 0,]\$mu2
    close_df[close_df\$MUDIFF < 0,]\$leftmost_right <- close_df[close_df\$MUDIFF < 0,]\$V3
    close_df[close_df\$MUDIFF < 0,]\$rightmost_left <- close_df[close_df\$MUDIFF < 0,]\$V6
    # where second is leftmost one
    close_df[close_df\$MUDIFF > 0,]\$rightmost_Bidir <- close_df[close_df\$MUDIFF > 0,]\$V4
    close_df[close_df\$MUDIFF > 0,]\$rightmost_right <- close_df[close_df\$MUDIFF > 0,]\$V3
    close_df[close_df\$MUDIFF > 0,]\$rightmost_mu <- close_df[close_df\$MUDIFF > 0,]\$mu1
    close_df[close_df\$MUDIFF > 0,]\$leftmost_Bidir <- close_df[close_df\$MUDIFF > 0,]\$V8
    close_df[close_df\$MUDIFF > 0,]\$leftmost_left <- close_df[close_df\$MUDIFF > 0,]\$V6
    close_df[close_df\$MUDIFF > 0,]\$leftmost_mu <- close_df[close_df\$MUDIFF > 0,]\$mu2
    close_df[close_df\$MUDIFF > 0,]\$leftmost_right <- close_df[close_df\$MUDIFF > 0,]\$V7
    close_df[close_df\$MUDIFF > 0,]\$rightmost_left <- close_df[close_df\$MUDIFF > 0,]\$V2
    # check by new MUDIFF should all be positive
    close_df\$leftmost_mu <- as.integer(close_df\$leftmost_mu)
    close_df\$rightmost_mu <- as.integer(close_df\$rightmost_mu)
    close_df\$newMUDIFF <- close_df\$rightmost_mu - close_df\$leftmost_mu
    if (nrow(close_df[close_df\$newMUDIFF < 0,]) > 0) {
        stop("There is a problem in the overlap files of the Bids") }
    
    return(list(close_df, orig_close_df))
    }

get_new_pos <- function(close_df) {
    # This function gets the new positions for regions BUT assumes the regions are only overlapping
    #      one other region if any.
    # PARAMETERS
    #     close_df: a dataframe with the columns chrom, start, stop, name of region 1, the same for region 2, and finally the distance between the two regions. It is assumed that the middle of the start and stop of each region (mu) is the starting point of each transcript.
    # make easier naming of info
    close_df\$left_neg_stop = close_df\$leftmost_mu 
    close_df\$left_neg_start = close_df\$leftmost_left
    close_df\$right_neg_start = pmax(close_df\$leftmost_mu, close_df\$rightmost_left) # new right neg 'start' should be rightmost of left mu or the original end
    close_df\$right_neg_stop = close_df\$rightmost_mu
    close_df\$left_pos_start = close_df\$leftmost_mu
    close_df\$left_pos_stop = pmin(close_df\$rightmost_mu, close_df\$leftmost_right) # new right pos should end at min of right mu or original call
    close_df\$right_pos_start = close_df\$rightmost_mu
    close_df\$right_pos_stop = close_df\$rightmost_right
    print(dim(close_df[(close_df\$left_neg_stop - close_df\$left_neg_start <= 0) | 
                       (close_df\$left_pos_stop - close_df\$left_pos_start <= 0) | 
                       (close_df\$right_neg_stop - close_df\$right_neg_start <= 0) | 
                       (close_df\$right_pos_stop - close_df\$right_pos_start <= 0),]))
    return(close_df)
    }


get_inbetween_lines <- function(close_df, in_between) {
    # This function gets the GTF lines for cases where a region has overlapping regions
    #     on both sides of it.
    # PARAMETERS
    #     close_df: a dataframe with the columns chrom, start, stop, name of region 1, the same for region 2, and finally the distance between the two regions. It is assumed that the middle of the start and stop of each region (mu) is the starting point of each transcript.
    #     in_between: the list of regions that have overlapping regions on both sides 
    # for each in between
    # Initialize a list to store results (more efficient than growing a vector)
    between_lines_list = vector("list", length(in_between))
    for (i in seq_along(in_between)) {
        between_bid <- in_between[i]
        # Filter once for each between_bid to avoid redundant filtering
        left = close_df[close_df\$leftmost_Bid == between_bid,]
        right = close_df[close_df\$rightmost_Bid == between_bid,]
        # Check if left and right have rows to avoid subscript errors
        if (nrow(left) > 0 && nrow(right) > 0) {
            # Build the gene, pos coord, and neg coord strings
            gene_line <- paste0(left\$V1[1], "\t.\tgene\t", left\$leftmost_left[1], "\t", left\$leftmost_right[1], 
                                '\t.\t.\t.\tgene_id "', between_bid, '";')
            
            pos_coord_line <- paste0(left\$V1[1], "\t.\texon\t", left\$leftmost_mu[1], "\t", 
                                     min(c(left\$rightmost_mu, left\$leftmost_right[1])), 
                                     '\t.\t+\t.\tgene_id "', between_bid, '"; transcript_id "', between_bid, '";')
            
            neg_coord_line <- paste0(left\$V1[1], "\t.\texon\t", max(c(right\$leftmost_mu[1], left\$leftmost_left[1])), "\t", 
                                     right\$rightmost_mu[1],  
                                     '\t.\t-\t.\tgene_id "', between_bid, '"; transcript_id "', between_bid, '";')
            
            # Store results in the list instead of concatenating each time
            between_lines_list[[i]] <- c(gene_line, pos_coord_line, neg_coord_line)
        }
    }
    
    # Combine the list into a single vector after the loop
    between_lines <- unlist(between_lines_list)
    return(between_lines)
}

get_isolated_lines <- function(orig_close_df, isolated) {
    # get non duplicated
    orig_close_df = orig_close_df[orig_close_df\$V4 %in% isolated,]
    orig_close_df = orig_close_df[,c("V1", "V2", "V3", "V4", "mu1")]
    orig_close_df = orig_close_df[!duplicated(orig_close_df),]
    if (nrow(orig_close_df) != length(isolated)) {
        cat(nrow(orig_close_df), length(isolated))
        stop("Not a match in duplications")}
    # for each isolated bidir
    isolated_lines = c()
    new_close_df <- data.frame(data.table("gene_row"=paste0(orig_close_df\$V1, "\t.\tgene\t", orig_close_df\$V2, 
                                                            "\t", orig_close_df\$V3, '\t.\t.\t.\tgene_id "', orig_close_df\$V4, '";'), 
                               "pos_row"=paste0(orig_close_df\$V1, "\t.\texon\t", orig_close_df\$mu1, "\t",
                                                           orig_close_df\$V3,  '\t.\t+\t.\tgene_id "', orig_close_df\$V4, 
                                                  '"; transcript_id "', orig_close_df\$V4, '";'), 
                               "neg_row"=paste0(orig_close_df\$V1, "\t.\texon\t", orig_close_df\$V2, "\t",
                                                           orig_close_df\$mu1,  '\t.\t-\t.\tgene_id "', orig_close_df\$V4, 
                                                  '"; transcript_id "', orig_close_df\$V4, '";')))
    isolated_lines <- as.vector(t(new_close_df))
    return(isolated_lines)
    }

get_GTF_lines <- function(close_df, orig_close_df, bids_keep) {
    # This function gets the annotation lines to use for a GTF that accounts for overlapping regions, including regions that have overlap with other regions on both strands.
    # PARAMETERS
    #       close_df: close_df that has already gone through the functions get_info_for_calls and get_new_pos AND has been filtered according
    #           to max_distance.
    #       orig_close_df: close_df that has already gone through the functions get_info_for_calls and get_new_pos but has not been filtered.
    #       bids_keep: list of strings (bids that are to be considered for the GTF)
    inbetween <- intersect(close_df\$rightmost_Bidir, close_df\$leftmost_Bidir)
    # only keep those that can be deconvoluted
    inbetween_keep <- intersect(inbetween, bids_keep)
    isolated <- setdiff(orig_close_df\$V4, union(close_df\$rightmost_Bidir, close_df\$leftmost_Bidir))
    # only keep those that can be deconvoluted
    isolated <- intersect(isolated, bids_keep)
    cat("\nIsolated", length(isolated), "\n")
    isolated_lines = get_isolated_lines(orig_close_df, isolated)
    cat("\nIn between", length(inbetween), "\n")
    inbetween_lines = get_inbetween_lines(close_df, inbetween_keep)
    
    close_df\$left_gene_row = paste0(close_df\$V1, "\t.\tgene\t", close_df\$leftmost_left, "\t", close_df\$leftmost_right, 
                                   '\t.\t.\t.\tgene_id "', close_df\$leftmost_Bidir, '";')
    close_df\$left_neg_row = paste0(close_df\$V1, "\t.\texon\t", close_df\$left_neg_start, "\t", close_df\$left_neg_stop, 
                                   '\t.\t-\t.\tgene_id "', close_df\$leftmost_Bidir, 
                                  '"; transcript_id "', close_df\$leftmost_Bidir, '";' )
    close_df\$left_pos_row = paste0(close_df\$V1, "\t.\texon\t", close_df\$left_pos_start, "\t", close_df\$left_pos_stop, 
                                   '\t.\t+\t.\tgene_id "', close_df\$leftmost_Bidir, 
                                  '"; transcript_id "', close_df\$leftmost_Bidir, '";' )
    close_df\$right_gene_row = paste0(close_df\$V1, "\t.\tgene\t", close_df\$rightmost_left, "\t", close_df\$rightmost_right, 
                                   '\t.\t.\t.\tgene_id "', close_df\$rightmost_Bidir, '";')
    close_df\$right_neg_row = paste0(close_df\$V1, "\t.\texon\t", close_df\$right_neg_start, "\t", close_df\$right_neg_stop, 
                                   '\t.\t-\t.\tgene_id "', close_df\$rightmost_Bidir, 
                                  '"; transcript_id "', close_df\$rightmost_Bidir, '";' )
    close_df\$right_pos_row = paste0(close_df\$V1, "\t.\texon\t", close_df\$right_pos_start, "\t", close_df\$right_pos_stop, 
                                   '\t.\t+\t.\tgene_id "', close_df\$rightmost_Bidir, 
                                  '"; transcript_id "', close_df\$rightmost_Bidir, '";' )
    # split into right and left
    right = close_df[,c("right_gene_row", "right_neg_row", "right_pos_row", "rightmost_Bidir")]
    left = close_df[,c("left_gene_row", "left_neg_row", "left_pos_row", "leftmost_Bidir")]
    # comebine
    colnames(right) = c("gene_row", "neg_row", "pos_row", "bidir")
    colnames(left) = c("gene_row", "neg_row", "pos_row", "bidir")
    # remove any in between
    df = rbind(right, left)
    df = df[!df\$bidir %in% inbetween,]
    # remove any not being kept due to convolution of counts
    df = df[df\$bidir %in% bids_keep,]
    df = df[!duplicated(df),]
    cat("\tFull Number Considered", length(inbetween)+length(isolated)+nrow(df), "\t")
    lines_vector <- as.vector(t(df[,c("gene_row", "neg_row", "pos_row")]))
    # remove any duplicates
    lines_vector <- lines_vector[!duplicated(lines_vector)]
    # if some in between
    if (length(inbetween) != 0) {
        # add the correct lines for the in between
        lines_vector <- c(lines_vector, inbetween_lines)
        }
    # add the isolated ones back in
    lines_vector <- c(lines_vector, isolated_lines)
    
    return(lines_vector)
    }

  ###################################s
  ## 1. Get the TSS Bidirectionals ##
  ###################################
  # read in the TSS region overlaps
  over <- fread("${tss_out_getbids}")
  # Get the distance between mu and TSS
  over\$mu <- as.integer((over\$V3+over\$V2)/2)
  over\$tss <- as.integer((over\$V6+over\$V7)/2)
  over\$MUDIFF <- abs(over\$mu-over\$tss)

  ## Get the Gene (not transcript) name
  over\$Gene <- str_split_fixed(over\$V8, ":", 2)[,1]
  trans_with_tss <- unique(over\$V8)
  genes_with_tss <- unique(over\$Gene)

  ## Get the TSS bidirectionals (note if hand annotation might be needed)
  over <- data.frame(over)
  over\$unique <- seq(1, nrow(over))
  hand_annotate = c(); TSS_used = c(); keep <- c(); remove <- c()
  for (isoform in trans_with_tss) {
    # get the overlaps
    filt = over[over\$V8 == isoform,]
    # get the TSS with the minimum MUDIFF
    tss = filt[filt\$MUDIFF == min(filt\$MUDIFF),]\$V4
    #cat("\t", isoform, tss)
    # if more than one bidirectional has the same min distance --> hand annotate (for now just keep first)
    if (length(tss) != 1) {
        hand_annotate = c(hand_annotate, isoform)
        tss <- tss[1]
    } 
    if (tss %in% TSS_used) {
        # if the TSS has already been assigned to a gene
        # check the difference in MUDIFF between the two genes
        filt2 = over[over\$V4 == tss,]
        # if the MUDIFF of the two genes are significantly different, only keep one
        if (max(abs(filt2\$MUDIFF))-min(abs(filt2\$MUDIFF)) > 50) {
            keep <- c(keep, filt2[abs(filt2\$MUDIFF) == min(abs(filt2\$MUDIFF)),]\$unique)
            remove <- c(remove, filt2[abs(filt2\$MUDIFF) == max(abs(filt2\$MUDIFF)),]\$unique)
        } else {
            keep <- c(keep, filt\$unique)
        }}else {
            keep <- c(keep, filt[filt\$V4 == tss,]\$unique)
            }
    TSS_used <- union(TSS_used, tss)
        }

  # only address the TSS bidirectionals clarified to keep and remove
  over_filt <- over[over\$unique %in% keep,]
  over_filt <- over_filt[!over_filt\$unique %in% remove,]  

  # Record the numbers changed
  trans_with_tss_filt <- unique(over_filt\$V8)
  cat("\tOriginal # Transcripts captured:", length(trans_with_tss), "and new #:", length(trans_with_tss_filt))
  genes_with_tss_filt <- unique(over_filt\$Gene)
  cat("\tOriginal # Genes captured:", length(genes_with_tss), "and new #:", length(genes_with_tss_filt))

  # Save the TSS bidirectionals
  colnames(over_filt) <- c("Chr", "Bid_Start", "Bid_Stop", "BidID", "Gene_Chr", "Gene_Start", "Gene_Stop", 
                        "TranscriptID", ".", "Strand", "Length", "mu", "TSS", "MUDIFF", "GeneID", "unique")
  over_filt <- over_filt[,c("Chr", "Bid_Start", "Bid_Stop", "BidID", "Gene_Start", "Gene_Stop",
             "TranscriptID", "Strand", "mu", "TSS", "GeneID")]

  write.table(over_filt, paste0("tss_bid.txt"), 
           row.names=FALSE, sep="\t", quote=FALSE)

  ##############################################
  ## 2. Filter the bids according to overlaps ##
  ##############################################
  # read in the overlaps of counted regions
  count_over <- fread("${count_out_getbids}")
  dwn_over <- fread("${dwn_out_getbids}")
  # read in the preliminary gene counts
  counts <- fread("${put_gene_counts_getbids}")
  colnames(counts) <- c("chrom", "start", "stop", "Gene", "dot", "strand", "count", "cov_bases", "length", "cov_frac")
  # ensure that the Genes and overlaps have the same Geneid format
  counts\$Geneid <- str_split_fixed(counts\$Gene, ",", 2)[,1]
  
  if (length(setdiff(count_over\$V4, counts\$Geneid)) > 1) {
    stop("There is an inconsistency in the naming between overlaps and counts")
  }
  # get the genes with high enough coverage to be considered transcribed
  high_counts = counts[counts\$cov_frac > ${params.count_limit_genes},]
  high_counts = high_counts\$Geneid
  cat("\tNUM gene isoforms with coverage >", as.character(${params.count_limit_genes}), length(high_counts), high_counts[1:4])
  # only consider overlaps with genes with OUT low counts
  count_over <- count_over[count_over\$V4 %in% high_counts]
  dwn_over <- dwn_over[dwn_over\$V4 %in% high_counts]
    cat("\tGOT HERE 3")
  # Split the overlaps into positive and negative strands
  count_over_neg = count_over[count_over\$V6 == "-", ]
  count_over_pos = count_over[count_over\$V6 == "+", ]
  dwn_over_neg = dwn_over[dwn_over\$V6 == "-", ]
  dwn_over_pos = dwn_over[dwn_over\$V6 == "+", ]

  # Get the NONTSS overlapping on one strand strands
  bid_pos_conflict <- setdiff(union(count_over_pos\$V10, dwn_over_pos\$V10), tss)
  bid_neg_conflict <- setdiff(union(count_over_neg\$V10, dwn_over_neg\$V10), tss)
  remove_bids = intersect(bid_pos_conflict, bid_neg_conflict)

  cat("\tTotal Number of NonTSS bids removing due to overlapping transcribed genes on both strands", 
    length(remove_bids))
  
  # read in the original bids
  bids <- fread("${count_win_file_getbids}")
  cat("\tTotal Number of NonTSS bids removing due to convolution", 
      length(remove_bids), "(Fraction = ", length(remove_bids)/length(setdiff(bids\$V4, over_filt\$BidID)), ")")
  bids <- bids[!bids\$V4 %in% remove_bids]
  cat("\tTotal Number of Bids remaining", length(unique(bids\$V4)), "Num NonTSS", length(setdiff(bids\$V4, over_filt\$BidID)), "\t")

  
  ######################################################
  ## 3. Get the counting files for the remaining bids ##
  ######################################################
  ### SIMPLE SAF
  # save the remaining bids as a saf file
  colnames(bids) <- c("Chr", "Start", "End", "GeneID")
  bids\$Strand <- "+"
  bids[1:2,]
  write.table(bids[,c("GeneID", "Chr", "Start", "End", "Strand")], 
                  "filt.saf",
                quote=FALSE, sep="\t", row.names=FALSE)
  ### MU based counts
  # read in the closest bids overlapping file
  closest <- fread("${close_out_getbids}")
  # get the info needed for positions
  cat("\tUsing Count window", "${params.count_win}")
  closest_list = get_info_for_calls(closest, ${params.count_win})
  # get the proper positions
  closest = get_new_pos(closest_list[[1]])
  # get the GTF lines and write them
  gtf_lines = get_GTF_lines(closest, closest_list[[2]], bids\$GeneID)
  gtf_file = "uns.gtf"
  writeLines(gtf_lines, gtf_file)
  # now get the positive and negative forms
  gtf <- fread(gtf_file)
  pos_gtf = gtf[gtf\$V7 %in% c(".", "+"),]
  neg_gtf = gtf[gtf\$V7 %in% c(".", "-"),]
  gene_only = gtf[gtf\$V3 == "gene",]
  cat("Number of Bids in GTF", nrow(gene_only))
  if (length(unique(gene_only\$V9[duplicated(gene_only\$V9)])) > 0) {
          stop("There are duplicate regions in the GTF") }
      write.table(pos_gtf, "pos.gtf", 
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
      write.table(neg_gtf, "neg.gtf", 
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

  

  ##############################################
  ## 4. If TFEA output requested, print out ##
  ##############################################
  if ("${params.tfea}" == "TRUE") {
    prefix = "${params.prefix}"
    date = "${params.date}"
    # get the 1500bp wwindows
    bids\$mu <- as.integer((as.numeric(bids\$Start) + as.numeric(bids\$End))/2)
    bids\$Start <- bids\$mu - 1500
    bids\$End <- bids\$mu + 1500
      tss <- bids[bids\$GeneID %in% over_filt\$BidID,]
      nontss <- bids[!bids\$GeneID %in% over_filt\$BidID,]
      cat("\tNumber total Bids after filtering:", nrow(bids), "\t\tNumber TSS:", nrow(tss), "\tNumber NonTSS:", nrow(nontss))
      # save
      write.table(tss[,c("Chr", "Start", "End", "GeneID")], 
              "tss_bid_forTFEA.bed",
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
      write.table(nontss[,c("Chr", "Start", "End", "GeneID")], 
              "nontss_bid_forTFEA.bed",
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
      print("Got TFEA Bids")
      } else {
        print("NOT GETTING TFEA FILES")
      }
  """
}

  process SimpleCounts {
    println "[Log 5]: Getting fixed bidirectional counts using Simple window ${params.count_win}"

    cpus 16
    queue 'short'
    memory '10 GB'
    time '3h'
    tag "$prefix"

    input:
    file ('*') from sorted_mmfilt_bam_file_bidcount.collect()
    file bidirs_saf

    output:
    file ("uns_bidirs_simple.txt") into uns_bidir_counts_SIMPLE
    file ("pos_bidirs_simple.txt") into pos_bidir_counts_SIMPLE
    file ("neg_bidirs_simple.txt") into neg_bidir_counts_SIMPLE

    script:
    """
    #!/bin/bash
    if [ -d "mmfiltbams" ]; then
        # count in an unstranded manner
        featureCounts -T 16 -O -s 0 -a ${bidirs_saf} -F 'SAF' -o uns_bidirs_simple.txt ./mmfiltbams/*.bam 
        # count on positive strand
        featureCounts -T 16 -O -s 1 -a ${bidirs_saf} -F 'SAF' -o pos_bidirs_simple.txt ./mmfiltbams/*.bam 
        # count on negative strand
        featureCounts -T 16 -O -s 2 -a ${bidirs_saf} -F 'SAF' -o neg_bidirs_simple.txt ./mmfiltbams/*.bam
    else
        # count in an unstranded manner
        featureCounts -T 16 -O -s 0 -a ${bidirs_saf} -F 'SAF' -o uns_bidirs_simple.txt ./*mmfilt.sorted.bam  
        # count on positive strand
        featureCounts -T 16 -O -s 1 -a ${bidirs_saf} -F 'SAF' -o pos_bidirs_simple.txt ./*mmfilt.sorted.bam  
        # count on negative strand
        featureCounts -T 16 -O -s 2 -a ${bidirs_saf} -F 'SAF' -o neg_bidirs_simple.txt ./*mmfilt.sorted.bam 
    fi
     
    """
  }

  process MuCounts {
    println "[Log 5]: Getting fixed bidirectional counts using MU_counts and original window ${params.count_win}"

    cpus 16
    queue 'short'
    memory '10 GB'
    time '5h'
    tag "$prefix"

    input:
    file ('*') from sorted_mmfilt_bam_file_mucount.collect()
    file neg_gtf
    file pos_gtf
    file uns_gtf

    output:
    file ("str_bidirs.txt") into uns_bidir_counts_MU
    file ("pos_bidirs.txt") into pos_bidir_counts_MU
    file ("neg_bidirs.txt") into neg_bidir_counts_MU

    script:
    """
    #!/bin/bash

    if [ -d "mmfiltbams" ]; then
        # count on both strands
        featureCounts -T 16 -O -s 1 -a ${uns_gtf} -t "exon" -F "GTF" -o str_bidirs.txt ./mmfiltbams/*.bam 
        # count only on positive strand
        featureCounts -T 16 -O -s 1 -a ${pos_gtf} -t "exon" -F "GTF" -o pos_bidirs.txt ./mmfiltbams/*.bam 
        # count only on negatve strand
        featureCounts -T 16 -O -s 1 -a ${neg_gtf} -t "exon" -F "GTF" -o neg_bidirs.txt ./mmfiltbams/*.bam 
    else
        # count on both strands
        featureCounts -T 16 -O -s 1 -a ${uns_gtf} -t "exon" -F "GTF" -o str_bidirs.txt ./*mmfilt.sorted.bam  
        # count only on positive strand
        featureCounts -T 16 -O -s 1 -a ${pos_gtf} -t "exon" -F "GTF" -o pos_bidirs.txt ./*mmfilt.sorted.bam  
        # count only on negatve strand
        featureCounts -T 16 -O -s 1 -a ${neg_gtf} -t "exon" -F "GTF" -o neg_bidirs.txt ./*mmfilt.sorted.bam 
    fi
        
    
    """
  }

process FixCounts {
  cpus 1
  queue 'short'
  memory '5 GB'
  time '1h'
  tag "$prefix"


  publishDir "${params.outdir}/counts/", mode: 'copy',
    saveAs: { filename ->
        if (filename == "Fixed_tss_bidir_counts.txt") { 
            return "fixed_MU_genetss_${params.prefix}_${params.date}_counts.txt" 
        } else if (filename == "Fixed_nontss_bidir_counts.txt") { 
            return "fixed_MU_nongenetss_${params.prefix}_${params.date}_counts.txt" 
        } else if (filename == "Fixed_tss_bidir_counts_SIMPLE.txt") { 
            return "fixed_genetss_${params.prefix}_${params.date}_counts.txt" 
        } else if (filename == "Fixed_nontss_bidir_counts_SIMPLE.txt") { 
            return "fixed_nongenetss_${params.prefix}_${params.date}_counts.txt" 
        } else { 
            return null 
        }
    }

  input:
  file uns_bidir_counts_MU
  file pos_bidir_counts_MU
  file neg_bidir_counts_MU
  file uns_bidir_counts_SIMPLE
  file pos_bidir_counts_SIMPLE
  file neg_bidir_counts_SIMPLE
  file count_out_fix
  file dwn_out_fix
  file tss_file
  file put_gene_counts_fix

  output:
  file("Fixed_tss_bidir_counts.txt")
  file("Fixed_nontss_bidir_counts.txt")
  file("Fixed_tss_bidir_counts_SIMPLE.txt")
  file("Fixed_nontss_bidir_counts_SIMPLE.txt")
  file("nontssbid_uniqueid_above_countreq.bed") into use_bid




  script:
  """
  #!/usr/bin/env Rscript
  library(data.table)
  library(stringr)
  ########################
  ## Read in the counts ##
  ########################
  # read in counts
  uns_counts <- fread("${uns_bidir_counts_MU}")
  pos_counts <- fread("${pos_bidir_counts_MU}")
  neg_counts <- fread("${neg_bidir_counts_MU}")
  uns_counts_simple <- fread("${uns_bidir_counts_SIMPLE}")
  pos_counts_simple <- fread("${pos_bidir_counts_SIMPLE}")
  neg_counts_simple <- fread("${neg_bidir_counts_SIMPLE}")
  #######################################################################
  ## Get the Bids that need to be fixed due to close/overlapping genes ##
  #######################################################################
  # read in the overlaps of bids with transcripts
  count_over <- fread("${count_out_fix}")
  dwn_over <- fread("${dwn_out_fix}")

  # Get the TSS bidirectionals
  tss <- fread("${tss_file}")
  tss <- unique(tss\$BidID)

  # Read in the preliminary gene counts
  counts <- fread("${put_gene_counts_fix}")
  colnames(counts) <- c("chrom", "start", "stop", "bidir", "dot", "strand", "count", "cov_bases", "length", "cov_frac")
  # ensure that the Genes and overlaps have the same Geneid format
  counts\$Geneid <- str_split_fixed(counts\$bidir, ",", 2)[,1]
  # get the high counts
  high_counts = counts[counts\$cov_frac > ${params.count_limit_genes},]\$Geneid
  # only consider overlaps with genes with high counts
  count_over <- count_over[count_over\$V4 %in% high_counts]
  dwn_over <- dwn_over[dwn_over\$V4 %in% high_counts]

  # Split the overlaps into positive and negative strands
  count_over_neg = count_over[count_over\$V6 == "-", ]
  count_over_pos = count_over[count_over\$V6 == "+", ]
  dwn_over_neg = dwn_over[dwn_over\$V6 == "-", ]
  dwn_over_pos = dwn_over[dwn_over\$V6 == "+", ]

  # Get the NONTSS overlapping on one strand strands
  bid_pos_conflict <- setdiff(union(count_over_pos\$V10, dwn_over_pos\$V10), tss)
  bid_neg_conflict <- setdiff(union(count_over_neg\$V10, dwn_over_neg\$V10), tss)
  remove_bids = intersect(bid_pos_conflict, bid_neg_conflict)

  pos_counts = pos_counts[!pos_counts\$Geneid %in% remove_bids,]
  uns_counts = uns_counts[!uns_counts\$Geneid %in% remove_bids,]
  neg_counts = neg_counts[!neg_counts\$Geneid %in% remove_bids,]

  pos_counts_simple = pos_counts_simple[!pos_counts_simple\$Geneid %in% remove_bids,]
  uns_counts_simple = uns_counts_simple[!uns_counts_simple\$Geneid %in% remove_bids,]
  neg_counts_simple = neg_counts_simple[!neg_counts_simple\$Geneid %in% remove_bids,]

  ##########################
  ## Get the Fixed Counts ##
  ##########################
  # use a substring only found in the samples (not the other headings)
  look_str <- "bam"

  # get the pos counts of nontss bidirectionals w/ gene on neg strand
  pos_counts <- data.frame(pos_counts[pos_counts\$Geneid %in% bid_neg_conflict,])
  colnames(pos_counts) <- colnames(uns_counts)
  pos_counts_simple <- data.frame(pos_counts_simple[pos_counts_simple\$Geneid %in% bid_neg_conflict,])
  colnames(pos_counts_simple) <- colnames(uns_counts_simple)
  # multiply by 2
  print(colnames(pos_counts)[grep(look_str, colnames(pos_counts))])
  sample = pos_counts[,colnames(pos_counts)[grep(look_str, colnames(pos_counts))]]
  sample = sample*2
  pos_counts <- cbind(pos_counts[,c("Geneid", "Chr", "Start", "End", "Strand", "Length")], sample)
  sample = pos_counts_simple[,colnames(pos_counts_simple)[grep(look_str, colnames(pos_counts_simple))]]
  sample = sample*2
  pos_counts_simple <- cbind(pos_counts_simple[,c("Geneid", "Chr", "Start", "End", "Strand", "Length")], sample)

  # get the neg counts of bidirectionals w/ gene on pos strand
  neg_counts <- data.frame(neg_counts[neg_counts\$Geneid %in% bid_pos_conflict,])
  colnames(neg_counts) <- colnames(uns_counts)
  # multiply by 2
  sample = neg_counts[,colnames(neg_counts)[grep(look_str, colnames(neg_counts))]]
  sample = sample*2
  neg_counts <- cbind(neg_counts[,c("Geneid", "Chr", "Start", "End", "Strand", "Length")], sample)
  neg_counts_simple <- data.frame(neg_counts_simple[neg_counts_simple\$Geneid %in% bid_pos_conflict,])
  colnames(neg_counts_simple) <- colnames(uns_counts_simple)
  # multiply by 2
  sample = neg_counts_simple[,colnames(neg_counts_simple)[grep(look_str, colnames(neg_counts_simple))]]
  sample = sample*2
  neg_counts_simple <- cbind(neg_counts_simple[,c("Geneid", "Chr", "Start", "End", "Strand", "Length")], sample)

  # Concatenate fixed counts to original counts of unoverlapping
  # First have the unstranded counts not include the conflicting ones
  uns_counts2 <- uns_counts[!uns_counts\$Geneid %in% bid_neg_conflict,]
  uns_counts2 <- uns_counts2[!uns_counts2\$Geneid %in% bid_pos_conflict,]
  uns_counts_simple2 <- uns_counts_simple[!uns_counts_simple\$Geneid %in% bid_neg_conflict,]
  uns_counts_simple2 <- uns_counts_simple2[!uns_counts_simple2\$Geneid %in% bid_pos_conflict,]

  # Combine and save
  full_counts <- rbind(uns_counts2, pos_counts, neg_counts)
  full_counts <- data.frame(full_counts)
  colnames(full_counts) <- colnames(uns_counts)

  full_counts_simple <- rbind(uns_counts_simple2, pos_counts_simple, neg_counts_simple)
  full_counts_simple <- data.frame(full_counts_simple)
  colnames(full_counts_simple) <- colnames(uns_counts_simple)

  # split into tss and nontss
  full_nontss_counts <- full_counts[!full_counts\$Geneid %in% tss,]
  full_tss_counts = full_counts[full_counts\$Geneid %in% tss,]

  full_nontss_counts_simple <- full_counts_simple[!full_counts_simple\$Geneid %in% tss,]
  full_tss_counts_simple = full_counts_simple[full_counts_simple\$Geneid %in% tss,]

  # save full counts
  write.table(full_nontss_counts, "Fixed_nontss_bidir_counts.txt", 
              sep="\t", quote=FALSE, row.names=FALSE)

  write.table(full_tss_counts, "Fixed_tss_bidir_counts.txt", 
              sep="\t", quote=FALSE, row.names=FALSE)

  write.table(full_nontss_counts_simple, "Fixed_nontss_bidir_counts_SIMPLE.txt", 
              sep="\t", quote=FALSE, row.names=FALSE)

  write.table(full_tss_counts_simple, "Fixed_tss_bidir_counts_SIMPLE.txt", 
              sep="\t", quote=FALSE, row.names=FALSE)

  ###################################################################
  ## Get the Bids that should be considered for fixed gene regions ##
  ###################################################################
   # I can't use this as a binary cuz stupid nextflow doesn't allow optional outputs atm
  #if ("${params.get_fixed_genecounts}" == "TRUE") {
  ### Only keep bidirectionals if they have more than X counts across all samples
  full_nontss_counts_moreX <- full_nontss_counts[rowSums(full_nontss_counts[,colnames(full_nontss_counts)[grep(look_str,colnames(full_nontss_counts))]]) > as.integer(${params.count_limit_bids}),]

  full_nontss_counts_moreX <- full_nontss_counts_moreX[,c("Chr", "Start", "End", "Geneid")]
  write.table(full_nontss_counts_moreX, "nontssbid_uniqueid_above_countreq.bed", 
              quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

 print("Got Bids for Genes")
  #}

  """
  
}



// PART 6 (OPTIONAL): GET THE FIXED GENE COUNTS

if (params.get_fixed_genecounts == "TRUE") {
  println "[Log 6]: Getting fixed gene regions"
  println "[Log 6]: Using minimum counts for Bid to be considered convolution .... ${params.count_limit_bids}"
  process GetFixedGeneRegions {
    
    input:
    file use_bid
    file count_trunc_out_getgenes

    output:
    file("new_regions.gtf") into count_gene_regions

    script:
    """
    # get new regions and stats
    python ${workflow.projectDir}/bin/Get_fixed_gene_regions.py --overlap_file ${count_trunc_out_getgenes} --new_regions_output_file "./new_regions.gtf" \
    --stats_output_file "./stats.txt" --full_gene_bed ${params.gene_count_file} --use_bids_bedfile ${use_bid}
    """
  }


process GetFixedGeneCounts {
    println "[Log 6]: Getting fixed gene counts"

   publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
              if ((filename == "fixed_str_genes.txt"))    "counts/fixed_genes_${params.prefix}_${params.date}_counts.txt"
              else null
             }
    
    input:
    file count_gene_regions
    file ('*') from sorted_mmfilt_bam_file_fingenecount.collect()

    output:
    file("fixed_str_genes.txt") into fixed_gene_counts

    script:
    """
    if [ -d "mmfiltbams" ]; then
        featureCounts \
       -T 16 -O -s 1 -a ${count_gene_regions} -t "exon" -F "GTF" -o fixed_str_genes.txt ./mmfiltbams/*.bam  
    else
        featureCounts \
       -T 16 -O -s 1 -a ${count_gene_regions} -t "exon" -F "GTF" -o fixed_str_genes.txt ./*mmfilt.sorted.bam          
    fi
    
    """
  }
}

// /*
//  * Completion report
//  */
// workflow.onComplete {

//     def report_fields = [:]
//     report_fields['version'] = params.version
//     report_fields['runName'] = custom_runName ?: workflow.runName
//     report_fields['success'] = workflow.success
//     report_fields['dateComplete'] = workflow.complete
//     report_fields['duration'] = workflow.duration
//     report_fields['exitStatus'] = workflow.exitStatus
//     report_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
//     report_fields['errorReport'] = (workflow.errorReport ?: 'None')
//     report_fields['commandLine'] = workflow.commandLine
//     report_fields['projectDir'] = workflow.projectDir
//     report_fields['summary'] = summary
//     report_fields['summary']['Date Started'] = workflow.start
//     report_fields['summary']['Date Completed'] = workflow.complete
//     report_fields['summary']['Pipeline script file path'] = workflow.scriptFile
//     report_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
//     report_fields['summary']['Pipeline repository Git URL'] = workflow.repository ?: 'Not stored'
//     report_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId ?: 'See hash'
//     report_fields['summary']['Pipeline Git branch/tag'] = workflow.revision ?: 'See hash'
//     report_fields['summary']['Nextflow Version'] = workflow.nextflow.version
//     report_fields['summary']['Nextflow Build'] = workflow.nextflow.build
//     report_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

//     // Render the TXT template
//     def engine = new groovy.text.GStringTemplateEngine()
//     def tf = new File("$baseDir/assets/report_template.txt")
//     def txt_template = engine.createTemplate(tf).make(report_fields)
//     def report_txt = txt_template.toString()

//     // Render the HTML template
//     def hf = new File("$baseDir/assets/report_template.html")
//     def html_template = engine.createTemplate(hf).make(report_fields)
//     def report_html = html_template.toString()

//     // Write summary HTML to a file
//     def output_d = new File( "${params.outdir}/pipeline_info/" )
//     if( !output_d.exists() ) {
//       output_d.mkdirs()
//     }
//     def output_hf = new File( output_d, "pipeline_report_bidir_${workflow.runName}.html" )
//     output_hf.withWriter { w -> w << report_html }
//     def output_tf = new File( output_d, "pipeline_report_bidir_${workflow.runName}.txt" )
//     output_tf.withWriter { w -> w << report_txt }

//     log.info "[Bidirectional-Flow] Pipeline Complete"

// }
