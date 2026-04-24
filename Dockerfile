# Base image con Miniconda
FROM continuumio/miniconda3

# Copia dei file di configurazione
COPY ANTIFASc.yml /tmp/ANTIFASc.yml
#COPY requirements.txt /tmp/requirements.txt

# Creazione dell'ambiente Conda
RUN conda env create -f /tmp/ANTIFASc.yml
ENV MPLCONFIGDIR=/tmp/matplotlib
RUN mkdir -p /tmp/matplotlib

# Attivazione dell'ambiente
SHELL ["conda", "run", "-n", "antifasc", "/bin/bash", "-c"]

# Installa Jupyter nel tuo env conda
RUN conda run -n antifasc pip install jupyterlab ipykernel

# Registra il kernel
RUN conda run -n antifasc python -m ipykernel install \
    --name antifasc \
    --display-name "Python (pyANTI-FASc)"

# Installazione delle dipendenze Python
# RUN pip install -r /tmp/requirements.txt

# Copia del codice sorgente
COPY . /app
WORKDIR /app

# Copia dell'eseguibile Fortran
RUN cp bin/k223d.x /opt/conda/envs/antifasc/bin/ && chmod +x /opt/conda/envs/antifasc/bin/k223d.x

# Comando di avvio (da personalizzare)
WORKDIR /app/bin
CMD ["conda", "run", "--no-capture-output", "-n", "antifasc", "python", "antifasc_main.py"]

