#!/bin/bash
# This shell script is running the figure script of the article with the right trailling argument

echo "Starting script..."

mkdir Figures
savePath="./Figures/"
pythonPath=$(which python3)

Rscript figureGenerator_article.R "$savePath" "$pythonPath"

echo "Script finished."
