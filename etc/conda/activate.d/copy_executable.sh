#!/bin/sh
# etc/conda/activate.d/copy_executable.sh

echo "Inizio esecuzione dello script."

# Naviga alla directory del repository per copiare gli eseguibili
SOURCE_DIR="$(dirname "$0")/../../.."
echo "SOURCE_DIR: $SOURCE_DIR"

# Variabile per verificare se i file esistono
FILES="$SOURCE_DIR/bin/k223d*"

# Controlla se ci sono file corrispondenti
if ls $FILES 1> /dev/null 2>&1; then
   echo "File trovati, procedo con la copia."
   # Copia tutti i file che iniziano con "executable" nella directory bin dell'ambiente
   cp $FILES "$CONDA_PREFIX/bin/"
   echo "Eseguibili copiati in $CONDA_PREFIX/bin/"
else
   echo "Errore: Nessun eseguibile trovato."
fi
#
