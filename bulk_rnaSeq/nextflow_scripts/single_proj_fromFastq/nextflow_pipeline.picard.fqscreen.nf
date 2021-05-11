#!/usr/bin/env nextFlow


 /*  ====================================
  *  singleproj_bulkRNA_fromFastq_v210107
  *  ====================================
*/


// set variables
projectID = params.projectID // param has to be set manually in nextflow.config file
launchdir = params.workdir

// automatically assigned in config file
OUTDIR = params.outDir // top outdir (output from individual projects will be added in separate folders herein)
FQDIR = params.fqDir
DMXSTATDIR = params.demuxStats
CNTDIR = params.quantDir
QCDIR = params.qcDir
STARDIR = params.bamDir
PICARDDIR = params.picardDir
FQSDIR = params.fastqScreenDir

lanes = params.lanes // should be set to lanes=0 if all lanes to be included. otherwise, set "1,3" etc.

// Read and process sample sheet
sheet = file(params.sheet)

// create new file for reading into channels that provide sample info!
newsheet = file("${launchdir}/sample_sheet.nf.csv")

// Set species for reference - Species defined in config
if ( params.species.contains("Human") ) {
    genome = params.human_genome
    gtf = params.human_gtf
    picard_RefFlat = params.picard_refFlat_hs
    picard_rRNA=params.picard_rRNA_hs
} else if ( params.species.contains("Mouse") ) {
    genome = params.mouse_genome
    gtf = params.mouse_gtf
    picard_RefFlat = params.picard_refFlat_mm
    picard_rRNA=params.picard_rRNA_mm
} else { // If not mouse or human, set custom references in config 
  genome = params.custom_genome
  gtf = params.custom_gtf    
}

allLines = sheet.readLines()
writeB = false // if next lines has sample info
readS = false // if next line contain species info
newsheet.text=""
for ( line in allLines ) {
    if ( writeB ) {
        newsheet.append(line + "\n")
    }
    if (line.contains("[Data]")) {
        writeB = true
    }
}


println "======= Info =========="
println ">>> Bulk RNA pipeline - single projects >>> "
println "> Project ID          : $projectID "
println "> Species    	       : $params.species " 
println "> Ref Genome	       : $genome "
println "> Ref GTF	       : $gtf " 
println "> Sample sheet	       : $sheet "
println "> Output dir 	       : $OUTDIR "
println "> Common Fastq dir    : $FQDIR "
println "> STAR dir	       : $STARDIR "
println "> fCount dir	       : $CNTDIR " 
println "======================= "


// sample info
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Name ) }
    .unique()
    .tap{infoSamples}
    .into{ fastqc_ch; fastqscreen_ch; star_ch}

infoSamples.subscribe{ println "Samples: $it" }

// fastqc 
process fastqc {

	input:
	set sid, sname from fastqc_ch

        output:
        val "x" into multiqc_fastqc
	
	"""

	mkdir -p ${OUTDIR}
        mkdir -p ${QCDIR}
	mkdir -p ${QCDIR}/FastQC

    	cd ${QCDIR}/FastQC
      
	read1=\$(echo ${FQDIR}/${sid}/${sname}_*R1*fastq.gz)
   	read2=\$(echo ${FQDIR}/${sid}/${sname}_*R2*fastq.gz)

    	# Check if fastq is not containing wildcards (due to sample fastq are not put in sample id folder
    	if [[ \${read1} == *"*R1*"* ]]; then
       	   read1=\$(echo ${FQDIR}/${sname}_*R1*fastq.gz)
       	   read2=\$(echo ${FQDIR}/${sname}_*R2*fastq.gz)
    	fi

	fastqc -t ${task.cpus} --outdir ${QCDIR}/FastQC \${read1}
	fastqc -t ${task.cpus} --outdir ${QCDIR}/FastQC \${read2}
 	
	"""
    
}


// fastq_screen
process fastqScreen {

    input:
    set sid, sname from fastqscreen_ch

    output:
    val "x" into multiqc_fastqscreen
    val "x" into multiqc_fastqscreen_qconly
    
    
    """

    mkdir -p ${FQSDIR}
    echo "HELLO COW"
    echo "${FQSDIR}"
    read1=\$(echo ${FQDIR}/$sid/${sname}_*R1*fastq.gz)
    read2=\$(echo ${FQDIR}/$sid/${sname}_*R2*fastq.gz)

    # Check if fastq is not containing wildcards (due to sample fastq are not put in sample id folder
    if [[ \${read1} == *"*R1*"* ]]; then
       read1=\$(echo ${FQDIR}/${sname}_*R1*fastq.gz)
       read2=\$(echo ${FQDIR}/${sname}_*R2*fastq.gz)
    fi

    
    echo "READ 1 : \${read1}"
    echo "READ 2 : \${read2}"


    /usr/local/bin/FastQ-Screen-0.14.1/fastq_screen \\
    --conf ${params.fastqScreen_config} \\
    --subset 500000 \\
    --outdir ${FQSDIR} \\
    \${read1} \${read2}

    """

}



// Run STAR
process STAR  {

    publishDir "${STARDIR}/", mode: 'copy', overwrite: true

    input:
    set sid, sname from star_ch

    output:
    set val(sname), file("${sname}*") into postStar
    file "${sname}_Aligned.sortedByCoord.out.bam" into starFeatureCounts
    file "${sname}_Aligned.sortedByCoord.out.bam" into starPicardRNA
    file "${sname}_Aligned.sortedByCoord.out.bam" into starPicardMarkDups
    

    when:
    params.align
   
    """
    mkdir -p ${OUTDIR}
    mkdir -p ${STARDIR}

    read1=\$(echo ${FQDIR}/$sid/${sname}_*R1*fastq.gz)
    read2=\$(echo ${FQDIR}/$sid/${sname}_*R2*fastq.gz)

    # Check if fastq is not containing wildcards (due to sample fastq are not put in sample id folder
    if [[ \${read1} == *"*R1*"* ]]; then
       read1=\$(echo ${FQDIR}/${sname}_*R1*fastq.gz)
       read2=\$(echo ${FQDIR}/${sname}_*R2*fastq.gz)
    fi

    
    echo "READ 1 : \${read1}"
    echo "READ 2 : \${read2}"

    STAR --genomeDir ${genome} \\
    --readFilesIn \${read1} \${read2} \\
    --runThreadN ${task.cpus}  \\
    --outSAMtype BAM SortedByCoordinate \\
    --readFilesCommand zcat \\
    --limitBAMsortRAM 10000000000 \\
    --outFileNamePrefix ${sname}

    """ 

}

process featureCounts {

	input:
	file bams from starFeatureCounts.collect()

	output:
	val "x" into postCount

	when:
	params.quant
	
	"""
 	mkdir -p ${CNTDIR}

        outfile=${CNTDIR}/${projectID}_geneid.featureCounts.txt

        featureCounts -T ${task.cpus} -t ${params.feature} --extraAttributes gene_name,gene_type -a ${gtf} -g gene_id  -o \${outfile} -p -s ${params.stranded} ${bams}
   
	"""	

}



process picardRNAmetrics {

    input:
    file bam from starPicardRNA

    output:
    val "x" into postPicardRNAmetrics

    when:
    params.quant

    """
    mkdir -p ${PICARDDIR}/collectRNAmetrics
    outfile=${PICARDDIR}/collectRNAmetrics/${bam}_rna.metrics

    java -jar /usr/local/bin/picard.jar CollectRnaSeqMetrics \\
      I=${bam} \\
      O=\${outfile} \\
      REF_FLAT=${picard_RefFlat} \\
      STRAND=SECOND_READ_TRANSCRIPTION_STRAND \\
      RIBOSOMAL_INTERVALS=${picard_rRNA}
    """

}


process picardMarkDups {

    input:
    file bam_pmd from starPicardMarkDups

    output:
    val "x" into postPicardMarkdups

    when:
    params.quant

    """
    mkdir -p ${PICARDDIR}/markdups
    # mkdir -p ${PICARDDIR}/markdups/bam_out
    mkdir -p ${OUTDIR}/STAR_markedDups

    # outfile=${OUTDIR}/markdups/bam_out/${bam_pmd}_marked_dups.bam
    outfile=${OUTDIR}/STAR_markedDups/${bam_pmd}_markedDups.bam
    metricsfile=${PICARDDIR}/markdups/${bam_pmd}_markedDup_metrics.txt

    java -jar /usr/local/bin/picard.jar MarkDuplicates \\
      I=${bam_pmd} \\
      O=\${outfile}
      M=\${metricsfile} \\
      REMOVE_DUPLICATES=false \\
      ASSUME_SORTED=true \\
      MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=2000
 
 #    O=/dev/null \\     
  

     """

}



process multiqc_postCount {

    input:
    val x from postCount
    val x from postPicardRNAmetrics.collect()
    val x from postPicardMarkdups.collect()
    val x from multiqc_fastqscreen.collect()
    
    output:
    val "${projectID}_multiqc_report.html" into multiqc_outPostCount

    when:
    params.align
    
    script:
    """
    mkdir -p ${QCDIR}/DemuxStats/
    cp ${DMXSTATDIR}/* ${QCDIR}/DemuxStats/

    cd ${OUTDIR}
    multiqc -n ${projectID}_multiqc_report --interactive -o ${QCDIR} .
    """
}


process multiqc_preCount {

    input:
    val x from multiqc_fastqc.collect()
    val x from multiqc_fastqscreen_qconly.collect()
    
    output:
    val "${projectID}_multiqc_report_fq.html" into multiqc_outPreCount

    when:
    params.qcOnly
    
    script:
    """
    mkdir -p ${QCDIR}/DemuxStats/
    cp ${DMXSTATDIR}/* ${QCDIR}/DemuxStats/

    cd ${OUTDIR}
    multiqc -o ${QCDIR}/ -n ${projectID}_multiqc_report_fq --interactive .
    """
}






