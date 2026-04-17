# 🐳 PyANTI-FASc Docker Guide (CLI Version)

This guide provides all instructions to install Docker Desktop (or Docker on Linux), build the Docker image, and run PyANTI-FASc directly from the terminal, with outputs immediately available in your local folder.

---

# ✅ Quick Start (if you already have installed and run Docker)

## Linux & Mac

``` bash
git clone https://github.com/antonioscalaunina/pyANTI-FASc.git  
cd pyANTI-FASc  
docker build -t pyantifasc .  
chmod +x antifasc  
./antifasc
```

## Windows Power Shell

``` bash
git clone https://github.com/antonioscalaunina/pyANTI-FASc.git  
cd pyANTI-FASc  
docker build -t pyantifasc .  
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned 
./antifasc.ps1
```
**WARNING**

`Set-ExecutionPolicy` might be not necessary

If `git` is not installed, you can download the ZIP version of the repository.




# 1️⃣ Install Docker

## Windows + WSL
Download Docker Desktop: https://www.docker.com/products/docker-desktop  
Run the installer. During installation, enable **WSL2 integration** (if you want to use WSL).  
Launch Docker Desktop and ensure it is running.

---

## Mac
Download Docker Desktop: https://www.docker.com/products/docker-desktop  
Install and launch Docker Desktop.

---

## Linux (Ubuntu)

``` bash
sudo apt update  
sudo apt install -y docker.io  
sudo usermod -aG docker $USER  
```

Then either restart your terminal or run:

``` bash
newgrp docker
```  

This allows you to run Docker commands without `sudo`.

---

# 2️⃣ Clone the repository

``` bash
git clone https://github.com/antonioscalaunina/pyANTI-FASc.git  
cd pyANTI-FASc
```  

If you are using Windows PowerShell, make sure that `git` is installed. Otherwise, you can download the ZIP version of the repository.

---

# 3️⃣ Build the Docker image

From the root of the project (where the Dockerfile is located):

``` bash
docker build -t pyantifasc .
```  

The image includes:
- Python and Conda
- All dependencies from `ANTIFASc.yml`
- The Fortran executable
- The runtime environment to execute the software

---

# 4️⃣ Create the `antifasc` command

To simplify usage, a launcher script is provided.

## Linux / Mac

Make the script executable:

``` bash
chmod +x antifasc
```  

Run the software:

``` bash
./antifasc
```  

Optional (to use it as a global command):

``` bash
sudo cp antifasc /usr/local/bin/
```  

Then you can simply run:

``` bash
antifasc
```  

---

## Windows PowerShell

Run:

``` bash
.\antifasc.ps1
```

If you get an error concerning the abilitation to script execution run this command

``` bash
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```
It is enough to run this command just once. From that moment on you can run

``` bash
.\antifasc.ps1
```

without any issue

---

# 5️⃣ How it works

The `antifasc` command runs:

``` bash
docker run --rm -it -v $(pwd):/app pyantifasc
```  

This means:

- Your local project folder is mounted inside the container at `/app`
- All generated outputs are written directly to your local filesystem
- No manual file copying is required

This ensures **1:1 synchronization between container and host**.

---

# 6️⃣ Customizing Inputs and Configuration

Before running the software, you may want to adapt the inputs or configuration parameters.

You can do this in two ways:

### 📂 Modify used input files
Edit the following lines where the input files are defined in the main script (`bin/antifasc_main.py`).  

``` bash
#Specify input json file and scaling_relationship file paths
input_file='../config_files/Parameters/input.json'
scaling_file='../config_files/Parameters/scaling_relationship.json'
``` 

### ⚙️ Modify configuration directly in the input files
Edit the parameter files located in:

``` bash
config_files/Parameters/
``` 

These control:
- Model parameters  
- Thresholds and settings  
- Processing options  

👉 This is the recommended approach to tune the model or run different scenarios.
refer to the Example.md files for more details
---

### 🔄 Important note

Because the project folder is mounted into the container:

``` bash
-v $(pwd):/app
```  

any changes you make locally are **immediately available inside Docker**.

👉 No rebuild is required.

---

### 📘 Examples

Refer to the example notebooks and sample configurations included in the repository to understand:
- Input data structure  
- Parameter usage  
- Typical workflows  

---

# 7️⃣ Run the software

Make sure you are in the main project directory:

``` bash
cd pyANTI-FASc
```   

Then execute:

``` bash
./antifasc
```   

The program will run directly inside the container.

---

# 8️⃣ Outputs

All generated outputs (GeoJSON, ASCII, etc.) will be immediately available in your local project folder:

``` bash
pyANTI-FASc/output/
```   

---

# 9️⃣ Re-running the container

After building the Docker image, you can run the software again at any time without rebuilding:

``` bash
./antifasc
```    

As long as Docker is running, the image will be reused.

---

# 🔟 Stopping the container

- The container will be automatically removed thanks to the `--rm` flag  

---
