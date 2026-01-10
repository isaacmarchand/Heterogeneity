#!/bin/bash
# This shell script is running the figure script of the article with the right trailling argument

echo "Starting script..."

mkdir Figures
savePath="./Figures/"

Rscript figureGenerator_article.R "$savePath"

echo "Script finished."
