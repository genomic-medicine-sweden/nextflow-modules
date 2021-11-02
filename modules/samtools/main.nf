// Samtools functions

process view {
  tag "$inputFile.simpleName"
  label "process_medium"

  input:
    path inputFile

  output:
    path outputFile

  script:
    def defaultFileType = input.getExtension()
    """
    samtools view $inputFile > $outputFile
    """
}

process sort {
  tag "$inputFile.simpleName"
  label "process_medium"

}
