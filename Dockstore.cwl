#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "ADTEx"
label: "ADTEx workflow"
cwlVersion: v1.0
description: |
    A Docker container for the ADTEx workflow. See the [github repo](https://github.com/Jeltje/adtex) for more information.
    ```
    Usage:
    # fetch CWL
    $> dockstore cwl --entry quay.io/jeltje/adtex:v1.0.1 > Dockstore.cwl
    # make a runtime JSON template and edit it (or use the content of sample_configs.json in this git repo)
    $> dockstore convert cwl2json --cwl Dockstore.cwl > Dockstore.json
    # run it locally with the Dockstore CLI
    $> dockstore launch --entry quay.io/jeltje/adtex:v1.0.1 \
        --json Dockstore.json
    ```

dct:creator:
  "@id": "jeltje"
  foaf:name: Jeltje van Baren
  foaf:mbox: "mailto:jeltje.van.baren@gmail.com"

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/jeltje/adtex:v1.0.1"

hints:
  - class: ResourceRequirement
    coresMin: 1

inputs:

  - id: "#centromeres"
    type: File
    doc: "Centromere bed file"
    format: "http://edamontology.org/format_3003"
    inputBinding:
      prefix: -c

  - id: "#targets"
    type: File
    doc: "Exome Targets bed file"
    format: "http://edamontology.org/format_3003"
    inputBinding:
      prefix: -b

  - id: "#control_bam_input"
    type: File
    doc: "The control exome BAM file used as input, it must be sorted."
    format: "http://edamontology.org/format_2572"
    inputBinding:
      prefix: -n 

  - id: "#tumor_bam_input"
    type: File
    doc: "The tumor exome BAM file used as input, it must be sorted."
    format: "http://edamontology.org/format_2572"
    inputBinding:
      prefix: -t 

  - id: "#sample_id"
    type: string
    doc: "sample ID to use in output"
    inputBinding:
      prefix: -s 


stdout: output.cnv

outputs: 
  - id: output
    type: stdout

baseCommand: ["-o", "/var/spool/cwl/workdir", "--tostdout"]

