#!/bin/bash

# Exit on error
set -e

echo "Installing SamovaR..."

# Set default values
DEFAULT_R_PATH="R"
DEFAULT_PYTHON_PATH="python"

# Get paths from user if specified in config.json
if [ -f build/config.json ]; then
    PYTHON_PATH=$(grep -o '"python_path": *"[^"]*"' build/config.json | sed 's/"python_path": *"\(.*\)"/\1/')
    R_PATH=$(grep -o '"r_path": *"[^"]*"' build/config.json | sed 's/"r_path": *"\(.*\)"/\1/')
    R_LIB_PATH=$(grep -o '"r_lib_path": *"[^"]*"' build/config.json | sed 's/"r_lib_path": *"\(.*\)"/\1/')
fi

# Check for custom R library path
if [ -z "$R_LIB_PATH" ]; then
    if [ -f .Rprofile ]; then
        R_LIB_PATH=$(grep "libPaths" .Rprofile | grep -o "'.*'" | tr -d "'")
    elif [ -f ~/.Renviron ]; then
        R_LIB_PATH=$(grep "R_LIBS*=" ~/.Renviron | cut -d'=' -f2 | tr -d '"')
    elif [ -f ~/.Rprofile ]; then
        R_LIB_PATH=$(grep "libPaths" ~/.Rprofile | grep -o "'.*'" | tr -d "'")
    fi
fi

# Set default values if not provided
PYTHON_PATH=${PYTHON_PATH:-$DEFAULT_PYTHON_PATH}
R_PATH=${R_PATH:-$DEFAULT_R_PATH}

# If R_LIB_PATH is not provided, get it from R
if [ -z "$R_LIB_PATH" ]; then
    R_LIB_PATH=$($R_PATH --quiet -e "cat(.libPaths()[1])" | head -2 | tail +2 | tr -d ' >')
fi

# Create or update config.json with all settings
cat > build/config.json << EOF
{
    "r_path": "$R_PATH",
    "python_path": "$PYTHON_PATH",
    "r_lib_path": "$R_LIB_PATH"
}
EOF

# Check if R is installed
if ! command -v $R_PATH &> /dev/null; then
    echo "R is not installed. Please install R first."
    exit 1
fi

# Check if Python is installed
if ! command -v $PYTHON_PATH &> /dev/null; then
    echo "Python 3 is not installed. Please install Python 3 first."
    exit 1
fi

# Get python version
PYTHON_VERSION=$($PYTHON_PATH --version | grep "Python" | awk '{print $2}')
echo "Detected Python version: $PYTHON_VERSION"

# Get R version
R_VERSION=$($R_PATH --version | grep "R version" | awk '{print $3}')
echo "Detected R version: $R_VERSION"

# If no custom path found, use default
if [ -z "$R_LIB_PATH" ]; then
    R_LIB_PATH=$($R_PATH --quiet -e "cat(.libPaths()[1])" | head -2 | tail +2)
fi

echo "Detected R library path: $R_LIB_PATH"

# Install R dependencies
echo "Installing R dependencies..."
$R_PATH --quiet -e "if (!require('remotes', lib='$R_LIB_PATH')) install.packages('remotes', repos='https://cloud.r-project.org/', lib='$R_LIB_PATH')"&> /dev/null
$R_PATH --quiet -e "if (!require('samovaR', lib='$R_LIB_PATH')) library(remotes, lib='$R_LIB_PATH'); remotes::install_deps(dependencies = TRUE, lib='$R_LIB_PATH')"&> /dev/null

# Install R package
echo "Installing R package..."
$R_PATH --quiet -e "library(remotes, lib='$R_LIB_PATH'); remotes::install_local('.', lib='$R_LIB_PATH')"&> /dev/null

# Install Python dependencies
echo "Installing Python dependencies..."
$PYTHON_PATH -m pip install --upgrade pip&> /dev/null
$PYTHON_PATH -m pip install -e .

# Create necessary directories
mkdir -p build
mkdir -p bin

# Make scripts executable
chmod +x bin/samovar
chmod +x bin/samovar_generate
chmod +x bin/samovar_preprocess
chmod +x bin/samovar_exec

# Add bin directory to PATH in .bashrc if not already present
if ! grep -q "export PATH=\$PATH:$(pwd)/bin" ~/.bashrc; then
    echo "export PATH=\$PATH:$(pwd)/bin" >> ~/.bashrc
    echo "Added SamovaR bin directory to PATH in ~/.bashrc"
    echo "Please run 'source ~/.bashrc' or restart your terminal to update PATH"
fi

# Create config.json if it doesn't exist
if [ ! -f build/config.json ]; then
    echo "Creating build/config.json..."
    cat > build/config.json << EOF
{
    "python_path": "$(which python3)",
    "r_path": "$(which R)",
    "r_lib_path": "$(R -e 'cat(.libPaths()[1])' 2>/dev/null)"
}
EOF
fi

echo "Installation completed successfully!" 