#!/usr/bin/env nextflow

nextflow.enable.dsl=2

IONICE = 'ionice -c2 -n7'

def get_star_index (genome) {
    return(params.star_index[genome])
}

def get_gtf (genome) {
    return(params.gtf[genome])
}

def get_chrom_sizes (genome) {
    return(get_star_index(genome) + '/chrNameLength.txt')
}

def get_genome (library) {
    return(params.libraries[library].genome)
}

def library_to_readgroups (library) {
    return(params.libraries[library].readgroups.keySet())
}

def library_and_readgroup_to_fastqs (library, readgroup) {
    return(params.libraries[library].readgroups[readgroup])
}


process star {

    publishDir "${params.results}/star/${library}", mode: 'rellink', overwrite: true
    container 'library://porchard/default/star:2.7.9a'
    tag "$library"
    memory '75 GB'
    cpus 10

    input:
    tuple val(library), val(genome), path(first_fastq), path(second_fastq)

    output:
    tuple val(library), path("${library}.Aligned.sortedByCoord.out.bam"), path("${library}.Log.final.out"), path("${library}.Log.out"), path("${library}.Log.progress.out")
    tuple val(library), val(genome), path("${library}.Aligned.sortedByCoord.out.bam"), emit: bam
    path("${library}.Log.final.out"), emit: for_multiqc

    """
    ${IONICE} STAR --runThreadN 10 --seedPerWindowNmax 30 --genomeLoad NoSharedMemory --outFileNamePrefix ${library}. --runRNGseed 789727 --readFilesCommand gunzip -c --outSAMattributes NH HI nM AS --genomeDir ${get_star_index(genome)} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs --sjdbGTFfile ${get_gtf(genome)} --readFilesIn ${first_fastq.join(',')} ${second_fastq.join(',')}
    """

}


process star_multiqc {

    publishDir "${params.results}/multiqc/star", mode: 'rellink', overwrite: true
    container 'library://porchard/default/general:20220107'

    input:
    path(x)

    output:
    path('multiqc_data')
    path('multiqc_report.html')

    """
    multiqc .
    """

}


process prune {

    publishDir "${params.results}/prune", mode: 'rellink', overwrite: true
    container 'library://porchard/default/general:20220107'
    maxForks 10
    tag "$library"

    input:
    tuple val(library), val(genome), path(bam)

    output:
    tuple path("${library}.pruned.bam"), path("${library}.pruned.bam.bai")

    """
    ${IONICE} samtools view -h -b -q 255 -F 4 -F 256 -F 2048 $bam | samtools sort -m 1g -O bam -T sort_tmp -o ${library}.pruned.bam && samtools index ${library}.pruned.bam
    """

}


process fastqc {

    publishDir "${params.results}/fastqc", mode: 'rellink', overwrite: true
    container 'library://porchard/default/general:20220107'
    maxForks 6
    tag "${library} ${readgroup}"

    input:
    tuple val(library), val(readgroup), path(fastq)

    output:
    tuple path(outfile_1), path(outfile_2)

    script:
    outfile_1 = fastq.getName().replaceAll('.fastq.gz', '_fastqc.html')
    outfile_2 = fastq.getName().replaceAll('.fastq.gz', '_fastqc.zip')

    """
    fastqc $fastq
    """

}


process fastq_multiqc {

    publishDir "${params.results}/multiqc/fastq", mode: 'rellink', overwrite: true
    container 'library://porchard/default/general:20220107'

    input:
    path(x)

    output:
    path('multiqc_data')
    path('multiqc_report.html')

    """
    multiqc .
    """

}


process qorts {

    publishDir "${params.results}/qorts"
    memory '50 GB'
    container 'library://porchard/default/r-general:20220112'
    tag "$library"

    input:
    tuple val(library), val(genome), path(bam)

    output:
    path("${library}")

    """
    mkdir -p $library
    ionice -c 2 -n 7 java -Xmx25g -jar \$QORTS_JAR QC --stranded --generatePlots --title $library --chromSizes ${get_chrom_sizes(genome)} $bam ${get_gtf(genome)} $library
    """

}


process qorts_multi {

    publishDir "${params.results}/qorts-multi"
    memory '50 GB'
    container 'library://porchard/default/r-general:20220112'

    input:
    path(qorts)

    output:
    path('multiqc')

    """
    echo "sample.ID qc.data.dir" > decoder.txt
    echo ${qorts.join(' ')} | perl -pe 's/ /\\n/g' | awk '{print(\$1, \$1)}' | perl -pe 's/ /\\t/' >> decoder.txt
    mkdir -p multiqc
    Rscript /usr/local/lib/R/site-library/QoRTs/extdata/scripts/qortsGenMultiQC.R ./ decoder.txt multiqc/
    """

}


workflow {

    libraries = params.libraries.keySet()
    fastq_in = []
    fastqc_in = []

    for (library in libraries) {
        for (readgroup in library_to_readgroups(library)) {
            fastqs = library_and_readgroup_to_fastqs(library, readgroup)
            genome = get_genome(library)
            first_read = fastqs['1']
            second_read = fastqs['2']
            fastqc_in << [library, readgroup, file(first_read)]
            fastqc_in << [library, readgroup, file(second_read)]
            fastq_in << [library, genome, file(first_read), file(second_read)]
        }
    }

    fastqc(Channel.from(fastqc_in)).flatten().toSortedList() | fastq_multiqc
    star_in = Channel.from(fastq_in).groupTuple(by: [0, 1])
    star_out = star(star_in)
    star_multiqc(star_out.for_multiqc.toSortedList())
    prune(star_out.bam)
    qorts(star_out.bam).toSortedList() | qorts_multi

}