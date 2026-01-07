#!/bin/bash
# This shell script is running the figure script of the article with the right trailling argument

echo "Starting script..."

savePath="/Users/macbook/Library/Mobile Documents/com~apple~CloudDocs/School/SFU/Research/Coding/Plots/testFigScript/"
pythonPath="/opt/anaconda3/bin/python3"

Rscript figureGenerator_article.R "$savePath" "$pythonPath"

echo "Script finished."
