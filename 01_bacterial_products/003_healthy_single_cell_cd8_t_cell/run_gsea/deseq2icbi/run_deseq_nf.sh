#!/bin/bash
nextflow run main.nf -profile icbi -resume "$@"
