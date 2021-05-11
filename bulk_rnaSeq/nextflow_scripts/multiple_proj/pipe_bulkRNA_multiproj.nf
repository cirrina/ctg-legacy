#!/usr/bin/env nextFlow

// set variables
// param has to be set manually in nextflow.config file
runDir = params.runfolder // raw sequnenser data
workDir = params.workdir // where scripts and outputs are put
metaID = params.metaid // overall ID of project etc

// automatically assigned in config file
OUTDIR = params.outDir // top outdir (output from individual projects will be added in separate folders herein)
FQDIR = params.fqDir

// should be set to lanes=0 if all lanes to be included. otherwise, set "1,3" etc.
lanes = params.lanes 

// Read and process sample sheet
sheet = file(params.sheet)

// create new file for reading into channels that provide sample info!
newsheet = file("${launchdir}/sample_sheet.nf.csv")

allLines = sheet.readLines()
writeB = false // if next lines has sample info
newsheet.text=""     

for ( line in allLines ) {
    if ( writeB ) {
	newsheet.append(line + "\n")
    }
    if (line.contains("[Data]")) {
	writeB = true
    }
}

// Set reference genomes / gtf
human_genome = params.human_genome
human_gtf = params.human_gtf
mouse_genome = params.mouse_genome
mouse_gtf = params.mouse_gtf
custom_genome = params.custom_genome
custom_gtf = params.custom_gtf


println "======= Info =========="
println ">>> Bulk RNA pipeline - multiple projects >>> "
println "> Experiment 	       : $exp "
println "> Sample sheet	       : $sheet "
println "> Experiment ID       : $metaID "
println "> Output dir 	       : $OUTDIR "
println "======================= "

// Projects
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> row.Sample_Project, row.Sample_Species }
    .unique()
    .tap{infoProject}
    .into{mvFastq_ch; quant_ch}

infoProject.subscribe{ println "Projects: $it" }

// sample info
Channel
    .fromPath(newsheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Name, row.Sample_Project, row.Sample_Species ) }
    .unique()
    .tap{infoSamples}
    .into{fastqc_ch; star_ch; fcount_ch }

infoSamples.subscribe{ println "Samples: $it" }

// Run bcl2fastq
process bcl2fastq {

    input:
    val sheet 

    output:
    val "x" into moveFastq
        
    """
    bcl2fastq -R $exp \\
              --sample-sheet $sheet \\
              --no-lane-splitting  \\
              -r $task.cpus \\
              -p $task.cpus  \\
              -w $task.cpus  \\
              --output-dir $FQDIR
   
     """
}

process moveFastq {

    input:
    val x from moveFastq
    set projid,species from mvFastq_ch

    output:
    val "x" into fqc_ch, star_go

    """
    mkdir -p ${OUTDIR}/${projid}
    mkdir -p ${OUTDIR}/${projid}/Fastq_Raw

    mv ${FQDIR}/${projid}/* ${OUTDIR}/${projid}/Fastq_Raw/

    """

}



// fastqc 
process fastqc {

	input:
	val x from fqc_ch.collect()
	set sid, sname, projid from fastqc_ch

        output:
        val projid into multiqc_fastqc
	
	"""

        mkdir -p ${OUTDIR}/$projid/QC/
        mkdir -p ${OUTDIR}/$projid/QC/FastQC

    	cd $OUTDIR/$projid/QC
      
	read1=\$(echo ${OUTDIR}/$projid/Fastq_Raw/$sid/${sname}*R1*fastq.gz)
   	read2=\$(echo ${OUTDIR}/$projid/Fastq_Raw/$sid/${sname}*R2*fastq.gz)

    	# Check if fastq is not containing wildcards (due to sample fastq are not put in sample id folder
    	if [[ \$read1 == *"*R1*"* ]]; then
       	   read1=\$(echo ${OUTDIR}/${projid}/Fastq_Raw/${sname}*R1*fastq.gz)
       	   read2=\$(echo ${OUTDIR}/${projid}/Fastq_Raw/${sname}*R2*fastq.gz)
    	fi

	fastqc -t $task.cpus --outdir $OUTDIR/$projid/QC/FastQC \$read1
	fastqc -t $task.cpus --outdir $OUTDIR/$projid/QC/FastQC \$read2
 	
	"""
    
}


// Run STAR
process STAR  {

    publishDir "${OUTDIR}/$projid/STAR/", mode: 'copy', overwrite: true

    input:
    val projDone from star_go.collect()
    set sid, sname, projid, species from star_ch

    output:
    set val(sname), file("${sname}*") into postStar
    file "${sname}Aligned.sortedByCoord.out.bam" into count

    when:
    params.align
   
    
    """

    read1=\$(echo ${OUTDIR}/$projid/Fastq_Raw/$sid/${sname}*R1*fastq.gz)
    read2=\$(echo ${OUTDIR}/$projid/Fastq_Raw/$sid/${sname}*R2*fastq.gz)

    # Check if fastq is not containing wildcards (due to sample fastq are not put in sample id folder
    if [[ \$read1 == *"*R1*"* ]]; then
       read1=\$(echo ${OUTDIR}/${projid}/Fastq_Raw/${sname}*R1*fastq.gz)
       read2=\$(echo ${OUTDIR}/${projid}/Fastq_Raw/${sname}*R2*fastq.gz)
    fi

    # Species 
    if [ $species == "Human" ] || [ $species == "human" ]; then
       genome=$human_genome
    elif [ $species == "Murine"  ] || [ $species == "murine" ] || [ $species == "mouse" ] || [ $species == "Mouse" ]; then
       genome=$mouse_genome
    else
       echo ">SPECIES NOT RECOGNIZED - using nextflow.config params!"
       genome=$custom_genome
    fi


    STAR --genomeDir \${genome} \\
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
    	set projid, species from quant_ch
	file bams from count.collect()

	output:
	val projid into postCount

	when:
	params.align
	params.quant
	
	"""

    # Species 
    if [ $species == "Human" ] || [ $species == "human" ]; then
       gtf=$human_gtf
    elif [ $species == "Murine"  ] || [ $species == "murine" ] || [ $species == "mouse" ] || [ $species == "Mouse" ]; then
       gtf=$mouse_gtf
    else
       echo ">SPECIES NOT RECOGNIZED - using nextflow.config params!"
       gtf=$custom_gtf
    fi
	 
 	mkdir -p  ${OUTDIR}/$projid/Quantification/

        # Get bam files for project
        bams=\$(echo ${OUTDIR}/$projid/STAR/*.bam)

        outfile=${OUTDIR}/$projid/Quantification/${projid}_genename.featureCounts.txt


        featureCounts -T ${task.cpus} -t ${params.feature} -a \${gtf} -g gene_name -o \${outfile} -p -s ${params.stranded} \${bams}

   
	"""	

}

process multiqc_subCount {

    input:
    val projid from postCount

    output:
    val "multiqc_report.html" into multiqc_outch

    when:
    params.align
    
    script:
    """
    cd $OUTDIR/$projid
    multiqc -n $projid -o ${OUTDIR}/$projid/QC/ .

    
    """
}



process multiqc_preCount {

    input:
    val projid from multiqc_fastqc

    output:
    val "multiqc_report.html" into multiqc_outPre

    when:
    params.demuxOnly
    
    script:
    """
    cd $OUTDIR/$projid
    multiqc -o ${OUTDIR}/$projid/QC/ -n $projid .
    """
}



