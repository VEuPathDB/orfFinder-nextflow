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
  index = indexResults(gff.collectFile())
  index.gff.collectFile(storeDir: params.outputDir, name: params.outputFileName)  
  index.indexed.collectFile(storeDir: params.outputDir, name: params.outputFileName + ".gz.tbi")
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

  input:
    path gff

  output:
    path gff, emit: gff
    path 'sorted_input.gff.gz.tbi', emit: indexed

  script:
  """
  sort -k1,1 -k4,4n $gff > sorted_input.gff
  bgzip sorted_input.gff
  tabix -p gff sorted_input.gff.gz
  """
}
