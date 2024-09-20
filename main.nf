#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//--------------------------------------------------------------------------
// Param Checking
//--------------------------------------------------------------------------

if(!params.fastaSubsetSize) {
  throw new Exception("Missing params.fastaSubsetSize")
}

if(params.inputFilePath) {
  seqs = Channel.fromPath( params.inputFilePath )
           .splitFasta( by:params.fastaSubsetSize, file:true  )
}
else {
  throw new Exception("Missing params.inputFilePath")
}

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------

workflow {
  gff = orfFinder(seqs, params.minPepLength)
  index = indexResults(gff.collectFile(),params.outputFileName)
}

process orfFinder {
  container = 'bioperl/bioperl:stable'

  input:
    path subsetFasta
    val minPepLength

  output:
    path 'orf_subset.gff'

  script:
  """
  orfFinder --dataset $subsetFasta \
    --minPepLength $minPepLength \
    --outFile orf_subset.gff
  """
}

process indexResults {
  container = 'biocontainers/tabix:v1.9-11-deb_cv1'

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
