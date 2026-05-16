# Base image con Miniconda
FROM anaconda/miniconda

# Copia dei file di configurazione
COPY ANTIFASc.yml /tmp/ANTIFASc.yml
#COPY requirements.txt /tmp/requirements.txt

# Creazione dell'ambiente Conda
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
RUN conda env create -f /tmp/ANTIFASc.yml
RUN mkdir -p /tmp/matplotlib && chmod 777 /tmp/matplotlib
ENV MPLCONFIGDIR=/tmp/matplotlib

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
RUN mkdir -p "$CONDA_PREFIX/bin" && \
    cp bin/k223d.x "$CONDA_PREFIX/bin/k223d.x" && \
    chmod +x "$CONDA_PREFIX/bin/k223d.x"
#RUN cp bin/k223d.x /opt/conda/envs/antifasc/bin/ && chmod +x /opt/conda/envs/antifasc/bin/k223d.x

# Comando di avvio (da personalizzare)
WORKDIR /app/bin
CMD ["conda", "run", "--no-capture-output", "-n", "antifasc", "python", "antifasc_main.py"]

