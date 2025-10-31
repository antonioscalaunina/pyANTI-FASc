# PyANTI-FASc Docker Guide

This guide provides all instructions to install Docker Desktop (or Docker on Linux), build the Docker image, and run the PyANTI-FASc Jupyter Notebook, with outputs immediately available in your local folder.

---

## 1️⃣ Install Docker

### Windows + WSL
1. Download Docker Desktop: [https://www.docker.com/products/docker-desktop](https://www.docker.com/products/docker-desktop)  
2. During installation, **enable WSL2 integration**  
3. Launch Docker Desktop and ensure it is running

### Mac
1. Download Docker Desktop: [https://www.docker.com/products/docker-desktop](https://www.docker.com/products/docker-desktop)  
2. Install and launch Docker Desktop

### Linux (Ubuntu)
```bash
sudo apt update
sudo apt install -y docker.io
sudo usermod -aG docker $USER
# Either restart your terminal or run:
newgrp docker
This allows you to run Docker commands without sudo.

## 2️⃣ Clone the repository
bash
Copy code
git clone https://github.com/antonioscalaunina/pyANTI-FASc.git
cd pyANTI-FASc
3️⃣ Build the Docker image
bash
Copy code
docker build -t pyantifasc-jupyter .
Includes Python, Conda, dependencies from ANTIFASc.yml, Fortran executable, and Jupyter Notebook.

4️⃣ Run the container with Jupyter Notebook
The -v $(pwd):/app (Linux/WSL/Mac) or ${PWD}:/app (Windows/PowerShell) mounts your local folder inside the container so all outputs are available on your host system.

Linux / WSL
bash
Copy code
docker run -it --rm -p 8888:8888 -v $(pwd):/app pyantifasc-jupyter
Windows PowerShell
powershell
Copy code
docker run -it --rm -p 8888:8888 -v ${PWD}:/app pyantifasc-jupyter
Mac
bash
Copy code
docker run -it --rm -p 8888:8888 -v $(pwd):/app pyantifasc-jupyter
Make sure you are in the folder where you want outputs to be saved.

5️⃣ Access the Jupyter Notebook
The terminal will show a URL with a token, for example:

ruby
Copy code
http://127.0.0.1:8888/?token=abc123...
Copy the URL including the token into your browser to access the notebook.

The token ensures a secure connection and is required unless disabled.

