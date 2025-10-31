# PyANTI-FASc Docker Guide

This guide provides all instructions to install Docker Desktop (or Docker on Linux), build the Docker image, and run the PyANTI-FASc Jupyter Notebook, with outputs immediately available in your local folder.

---

## 1️⃣ Install Docker

### Windows + WSL
1. Download Docker Desktop: [https://www.docker.com/products/docker-desktop](https://www.docker.com/products/docker-desktop)  
2. During installation, **enable WSL2 integration** (if you want to run it on WSL2 distribution on Windows) 
3. Launch Docker Desktop and ensure **it is running**.

### Mac
1. Download Docker Desktop: [https://www.docker.com/products/docker-desktop](https://www.docker.com/products/docker-desktop)  
2. Install and **launch Docker Desktop**

### Linux (Ubuntu)
```bash
sudo apt update
sudo apt install -y docker.io
sudo usermod -aG docker $USER
```
Either restart your terminal or run:
```
newgrp docker
```
This allows you to run Docker commands without sudo.


## 2️⃣ Clone the repository
```
git clone https://github.com/antonioscalaunina/pyANTI-FASc.git
cd pyANTI-FASc
```
If you are using Windows PowerShell you have to be sure that git is installed. Otherwise you can download the zipped code [here](https://github.com/antonioscalaunina/pyANTI-FASc/archive/refs/heads/main.zip)

## 3️⃣ Build the Docker image

From the root of the project (where the Dockerfile is located):
```
docker build -t pyantifasc-jupyter .
```

The image includes Python, Conda, all dependencies from ANTIFASc.yml, the Fortran executable, and Jupyter Notebook.

## 4️⃣ Run the container with Jupyter Notebook
In the following commands, the strings ```-v $(pwd):/app``` (Linux/WSL/Mac) or ```${PWD}:/app``` (Windows/PowerShell) mounts your local folder inside the container so all outputs are available on your host system.

### Linux / WSL
```
docker run -it --rm -p 8888:8888 -v $(pwd):/app pyantifasc-jupyter
```

### Windows PowerShell
```
docker run -it --rm -p 8888:8888 -v ${PWD}:/app pyantifasc-jupyter
```

### Mac
```
docker run -it --rm -p 8888:8888 -v $(pwd):/app pyantifasc-jupyter
```

Make sure you are in the main folder of the software. After building the Docker image, you can always run it again without rebuilding, as long as Docker Desktop is running.

In this case, you may verify that the image is available with the command:

```
docker images
```

## 5️⃣ Access the Jupyter Notebook
The terminal will show a URL with a token, for example:

```
http://127.0.0.1:8888/?token=abc123...
```

Copy the URL including the token into your browser to access the notebook.

The token ensures a secure connection and is required unless disabled.

All repository files will be accessible and editable from the Notebook.

Any generated outputs (GeoJSON, PDFs, reports) will be immediately available in your local folder, without manual copying.

This ensures 1:1 synchronization between container and host.

## 6️⃣  Stopping the container

Once your run is finished, Press ```CTRL + C``` in the terminal where Jupyter is running.
The container will be automatically removed thanks to the ```--rm``` flag.

