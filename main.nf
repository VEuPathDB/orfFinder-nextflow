#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------

workflow {
  if(!params.fastaSubsetSize) {
    error "Missing params.fastaSubsetSize"
  }

  if(params.inputFilePath) {
    seqs = channel.fromPath( params.inputFilePath )
             .splitFasta( by:params.fastaSubsetSize, file:true  )
  }
  else {
    error "Missing params.inputFilePath"
  }

  gff = orfFinder(seqs, params.minPepLength)
  index = indexResults(gff.collectFile(),params.outputFileName)
}

process orfFinder {
  container 'bioperl/bioperl:stable'

  input:
    path subsetFasta
    val minPepLength

  output:
    path 'orf_subset.gff'

  script:
  """
  orfFinder --inputFile $subsetFasta \
    --minPepLength $minPepLength \
    --outFile orf_subset.gff
  """
}

process indexResults {
  container 'biocontainers/tabix:v1.9-11-deb_cv1'

  publishDir params.outputDir, mode: 'copy'

  input:
    path gff
    val outputFileName

  output:
    path '*.gz'
    path '*.gz.tbi'

  script:
  """
  sort -k1,1 -k4,4n $gff > ${outputFileName}
  bgzip ${outputFileName}
  tabix -p gff ${outputFileName}.gz
  """
}
