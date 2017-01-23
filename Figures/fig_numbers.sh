#!/usr/bin/env bash

declare -a FIGURES=(
  "site_map.pdf"
  "distance_decay.pdf"
  "pam_sil.pdf"
  "pam_in_space.pdf"
  "diversity_distance.pdf"
  "diversity.pdf"
  "otu_in_space_select.pdf"
)

for (( i=1; i<="${#FIGURES[@]}"; ++i)) ; do
    ID=$((i - 1 )) 
    cp "${FIGURES[ID]}"  Fig"${i}".pdf
done



