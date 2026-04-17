# Usa una base image con Miniconda
FROM continuumio/miniconda3:latest

# Setta la directory di lavoro
WORKDIR /app

# Copia il file YAML di Conda nel container
# COPY ANTIFASc.yml /app/ANTIFASc.yml

# Copia tutto il repository dentro il container
COPY . /app/

# Crea l'ambiente Conda
RUN conda env create -f ANTIFASc.yml

# Attiva l'ambiente Conda
RUN echo "conda activate antifasc" >> ~/.bashrc

# Installa Jupyter Notebook
RUN conda install -n antifasc jupyter

# Copia dell'eseguibile Fortran
RUN cp bin/k223d.x /opt/conda/envs/antifasc/bin/ && chmod +x /opt/conda/envs/antifasc/bin/k223d.x

# Esponi la porta per Jupyter
EXPOSE 8888

# Comando per avviare Jupyter Notebook nel container
CMD ["bash", "-c", "source activate antifasc && jupyter notebook --ip=0.0.0.0 --port=8888 --allow-root"]

