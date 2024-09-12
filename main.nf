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
  orfFinder(seqs, params.minPepLength) | collectFile(storeDir: params.outputDir, name: params.outputFileName)
}


process orfFinder {
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
