#!/bin/bash

# Loop through all files matching the pattern NAME-draft.pdf
for file in *-draft.pdf; do
  # Extract the base name (everything before -draft.pdf)
  base_name="${file%-draft.pdf}"

  # Use pdfcrop to create a cropped version with the suffix -crop.pdf
  pdfcrop "$file" "${base_name}-draft-crop.pdf"

  # Rename the cropped file to NAME.pdf
  mv "${base_name}-draft-crop.pdf" "${base_name}.pdf"
done
rm *-draft.pdf

echo "Cropping and renaming completed."