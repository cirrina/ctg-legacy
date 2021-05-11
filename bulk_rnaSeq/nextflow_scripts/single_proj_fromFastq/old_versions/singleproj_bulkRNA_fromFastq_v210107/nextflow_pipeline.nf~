#!/usr/bin/env nextFlow


 /*  ====================================
  *  singleproj_bulkRNA_fromFastq_v201124
  *  ====================================
*/


// set variables
projectID = params.projectID // param has to be set manually in nextflow.config file
launchdir = params.workdir

// automatically assigned in config file
OUTDIR = params.outDir // top outdir (output from individual projects will be added in separate folders herein)
FQDIR = params.fqDir
CNTDIR = params.quantDir
QCDIR = params.qcDir
STARDIR = params.bamDir 

lanes = params.lanes // should be set to lanes=0 if all lanes to be included. otherwise, set "1,3" etc.

// Read and process sample sheet
sheet = file(params.sheet)

// create new file for reading into channels that provide sample info!
newsheet = file("${launchdir}/sample_sheet.nf.csv")

// Set species for reference - Species defined in config
if ( params.species.contains("Human") ) {
   genome = params.human_genome
   gtf = params.human_gtf
} else if ( params.species.contains("Mouse") ) {
  genome = params.mouse_genome
  gtf = params.mouse_gtf 
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
    .into{ fastqc_ch; star_ch}

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


// Run STAR
process STAR  {

    publishDir "${STARDIR}/", mode: 'copy', overwrite: true

    input:
    set sid, sname from star_ch

    output:
    set val(sname), file("${sname}*") into postStar
    file "${sname}Aligned.sortedByCoord.out.bam" into count

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
	file bams from count.collect()

	output:
	val "x" into postCount

	when:
	params.quant
	
	"""
	 
 	mkdir -p  ${CNTDIR}

        outfile=${CNTDIR}/${projectID}_geneid.featureCounts.txt

        featureCounts -T ${task.cpus} -t ${params.feature} --extraAttributes gene_name,gene_type -a ${gtf} -g gene_id  -o \${outfile} -p -s ${params.stranded} ${bams}

   
	"""	

}

process multiqc_postCount {

    input:
    val x from postCount

    output:
    val "${projectID}_multiqc_report.html" into multiqc_outPostCount

    when:
    params.align
    
    script:
    """
    cd ${OUTDIR}
    multiqc -n ${projectID}_multiqc_report --interactive -o ${QCDIR} .

    
    """
}



process multiqc_preCount {

    input:
    val x from multiqc_fastqc.collect()

    output:
    val "${projectID}_multiqc_report_fq.html" into multiqc_outPreCount

    when:
    params.qcOnly
    
    script:
    """
    cd ${OUTDIR}
    multiqc -o ${QCDIR}/ -n ${projectID}_multiqc_report_fq --interactive .
    """
}



