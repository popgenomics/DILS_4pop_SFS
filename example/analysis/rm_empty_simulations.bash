#!/bin/bash

# Chemin vers le fichier contenant la liste des répertoires
file="unfinished_sim.txt"

# Lire chaque ligne du fichier
while IFS= read -r directory
do
  # Vérifier si le répertoire existe
  if [ -d "$directory" ]; then
    echo "Contenu de $directory:"
    ls -l "$directory" # Liste détaillée des fichiers dans le répertoire
    rm -rf "$directory"
    echo "" # Ajoute une ligne vide pour la lisibilité
  else
    echo "Le répertoire $directory n'existe pas."
  fi
done < "$file"

