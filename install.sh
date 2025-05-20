#!/bin/bash

# Exit on error
set -e

echo "Installing SamovaR..."

# Check if R is installed
if ! command -v R &> /dev/null; then
    echo "R is not installed. Please install R first."
    exit 1
fi

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Python 3 is not installed. Please install Python 3 first."
    exit 1
fi

# Get R version
R_VERSION=$(R --version | grep "R version" | awk '{print $3}')
echo "Detected R version: $R_VERSION"

# Check for custom R library path
if [ -f .Rprofile ]; then
    R_LIB_PATH=$(grep "libPaths" .Rprofile | grep -o "'.*'" | tr -d "'")
elif [ -f ~/.Renviron ]; then
    R_LIB_PATH=$(grep "R_LIBS.*=" ~/.Renviron | cut -d'=' -f2 | tr -d '"')
elif [ -f ~/.Rprofile ]; then
    R_LIB_PATH=$(grep "libPaths" ~/.Rprofile | grep -o "'.*'" | tr -d "'")
fi

# If no custom path found, use default
if [ -z "$R_LIB_PATH" ]; then
    R_LIB_PATH=$(Rscript -e "cat(.libPaths()[1])")
fi

echo "Using R library path: $R_LIB_PATH"

# Create config.json with R settings
cat > config.json << EOF
{
    "r_lib_path": "$R_LIB_PATH",
    "r_version": "$R_VERSION"
}
EOF

# Install R dependencies
echo "Installing R dependencies..."
R -e "if (!require('devtools')) install.packages('devtools', repos='https://cloud.r-project.org/', lib='$R_LIB_PATH')"
R -e "devtools::install_deps(dependencies = TRUE, lib='$R_LIB_PATH')"

# Install R package
echo "Installing R package..."
R -e "devtools::install_local('.', lib='$R_LIB_PATH')"

# Install Python dependencies
echo "Installing Python dependencies..."
python3 -m pip install --upgrade pip
python3 -m pip install -e .

# Install Snakemake if not present
if ! command -v snakemake &> /dev/null; then
    echo "Installing Snakemake..."
    python3 -m pip install snakemake
fi

# Install ISS if not present
if ! command -v iss &> /dev/null; then
    echo "Installing InSilicoSeq..."
    python3 -m pip install in_silico_seq
fi

echo "Installation complete!" 