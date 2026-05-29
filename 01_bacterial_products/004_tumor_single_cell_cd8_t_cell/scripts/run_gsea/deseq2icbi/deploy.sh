#!/bin/bash
# Specifies bash as the interpreter.

# Function to display help message
display_help() {
    echo "Usage: $0 [options]"
    echo
    echo "   Deploys the DESeq2ICBI Nextflow pipeline by creating necessary"
    echo "   directories and symbolic links, and copying the configuration file."
    echo
    echo "   Options:"
    echo "     -h, --help    Display this help message and exit."
    echo
    echo "   The script will create a 'bin' directory (if it doesn't exist) one level"
    echo "   above its own location and link R scripts into it. It will also link"
    echo "   the main Nextflow script (main.nf) and a run script (run_deseq_nf.sh)"
    echo "   to the directory one level above its own location (project root)."
    echo "   The nextflow.config file will be copied to the project root for editing."
    exit 0
}

# Check for help flag
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    display_help
fi

# Indicates the start of deployment.
echo "Deploying DESeq2ICBI"

# Gets the absolute path of the script's directory.
BASE_DIR=$(realpath $(dirname "$0"))

# Defines a 'bin' directory one level above the script's location.
NF_BIN_DIR=${BASE_DIR}/../bin

# Extracts the name of the script's directory (e.g., "DESeq2ICBI_source").
DESEQ2_DIR=${BASE_DIR##*/}

# Prints the script's base directory path.
echo "Base directory: $BASE_DIR"
echo "Pipeline source directory name: $DESEQ2_DIR"
echo "Nextflow bin directory: $NF_BIN_DIR"


# Creates the NF_BIN_DIR if it doesn't exist.
echo "Creating ${NF_BIN_DIR} if it doesn't exist..."
mkdir -p ${NF_BIN_DIR}

# Changes current directory to NF_BIN_DIR.
echo "Changing directory to ${NF_BIN_DIR}..."
cd ${NF_BIN_DIR}

# Links Report.R from the source DESeq2 directory into NF_BIN_DIR.
# -s: create a symbolic link
# -f: force (remove existing destination files)
echo "Linking Report.R..."
ln -sf ../${DESEQ2_DIR}/Report.R .

# Links runDESeq2_ICBI.R from the source DESeq2 directory into NF_BIN_DIR.
echo "Linking runDESeq2_ICBI.R..."
ln -sf ../${DESEQ2_DIR}/runDESeq2_ICBI.R .

# Moves up one directory level (to the project root).
echo "Changing directory to project root (../)..."
cd ..

# Links main.nf from the source DESeq2 directory to the project root.
echo "Linking main.nf..."
ln -sf ${DESEQ2_DIR}/main.nf .

# Links run_deseq_nf.sh from the source DESeq2 directory to the project root.
echo "Linking run_deseq_nf.sh..."
ln -sf ${DESEQ2_DIR}/run_deseq_nf.sh .

# Copies nextflow.config from the source DESeq2 directory to the project root.
# -f: force overwrite if the destination file already exists.
echo "Copying nextflow.config..."
cp -f ${DESEQ2_DIR}/example_nextflow.config nextflow.config

# Informs the user to customize the copied nextflow.config.
echo
echo
echo "DONE: Please edit nextflow.config in the current directory to fit your needs."
echo
echo

