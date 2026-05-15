# pyANTI-FASc 
## Recommended version - Docker 🐳
Python implementation of **ANTI-FASc**  
(*Automatic Numerical Tsunami Initial conditions: on-the-Fly rupture Areas and earthquake Scenarios*)

---

## 🧠 General description

pyANTI-FASc is a software for the fast generation of large ensembles of stochastic slip distributions on complex non-planar fault geometries (Herrero & Murphy 2018; Maesano et al. 2017; Tonini et al. 2021).

These slip models can be used as initial conditions for tsunami simulations in:

- Probabilistic Tsunami Hazard Assessment (PTHA, see Basili et al. 2021)
- Real-time Probabilistic Tsunami Forecasting (PTF, see Selva et al. 2021)

👉 Please cite the software as indicated on Zenodo:  
https://doi.org/10.5281/zenodo.13614657

---
## 🐳 Install Docker

Before running pyANTI-FASc, make sure Docker is installed and running.

### Windows (also for WSL) / Mac

Download and install Docker Desktop:  
https://www.docker.com/products/docker-desktop  

- Launch Docker Desktop after installation  
- Make sure it is running before using the software

If you are using pyANTI-FASc on WSL on Windows make sure to *enable wsl integration* on Docker Desktop. Note that if you are on Windows 11 Home, WSL is the only mode to install Docker pyANTI-FASc image

---

### Linux (Ubuntu)

    sudo apt update
    sudo apt install -y docker.io
    sudo usermod -aG docker $USER


Then restart your terminal or run:

    newgrp docker


This allows you to run Docker without `sudo`.

## 🚀 Quick Start 🐳

### Linux / Mac

    git clone https://github.com/antonioscalaunina/pyANTI-FASc.git

    cd pyANTI-FASc
    docker build -t pyantifasc .
    chmod +x antifasc
    ./antifasc


### Windows PowerShell


    git clone https://github.com/antonioscalaunina/pyANTI-FASc.git

    cd pyANTI-FASc
    docker build -t pyantifasc .
    
    Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
    .\antifasc.ps1


---

## ⚙️ How to run

pyANTI-FASc can be launched either in CLI mode or through the interactive Jupyter notebook interface.

### <img width="20" height="20" alt="image" src="https://github.com/user-attachments/assets/75177071-5d70-4223-9026-541e322d5ddb" >  CLI execution

The command

    ./antifasc

or, on Windows PowerShell,

    .ntifasc.ps1

runs the full pipeline using the default configuration file:

    config_files/Parameters/input.json

A different input file can be selected with:

    ./antifasc --input input_test.json
    ./antifasc --input input_test.yaml

or, on Windows PowerShell,

    .ntifasc.ps1 --input input_test.json
    .ntifasc.ps1 --input input_test.yaml

If only a filename is provided, the file is searched in:

    config_files/Parameters/

Both JSON and YAML input files are supported. YAML files are automatically converted to the corresponding JSON file before execution, so the internal ANTIFASc workflow still reads a JSON configuration.

---

### <img width="30" height="30" alt="image" src="https://github.com/user-attachments/assets/ed54452d-a63f-48ec-8db2-080ed147353d" /> Jupyter Notebook (interactive)

The command

    ./antifasc notebook

or, on Windows PowerShell,

    .ntifasc.ps1 notebook

starts the JupyterLab environment.

The main notebook interface is available at:

    bin/antifasc_main.ipynb

The notebook provides interactive widgets to create or update `input.json`, including mesh selection, scaling laws, magnitude bins, application mode, event parameters, and ensemble configuration.

### 🔑 Accessing JupyterLab

After running the notebook command, a URL will appear in the terminal, similar to:


http://127.0.0.1:8888/lab?token=xxxxxxxxxxxxxxxx


👉 To open JupyterLab:

- **Linux (also WSL) / Mac**: press `Ctrl + Click` on the link (or copy and paste the URL into your browser)
- **Windows PowerShell**: copy and paste the URL into your browser  

⚠️ This token is required for authentication and changes every time you start the container.

You can then navigate within the file system of the repository and run the different Jupyter Notebook through the Jupyter interface

the main pipeline is available in the notebook:

    bin/antifasc_main.ipynb

## 📊 Interactive visualization

Interactive notebook available at:

     bin/interactive_slip_maps.ipynb

it can be run after the running of either the CLI pipeline or the notebook.

Running it the user can:

- select the folder via interactive widgets 👆
- visualize slip and rake maps through the GeoJSON output files
- optionally generate, visualize and export HTML interactive maps

*Note that to access to the notebooks you should enter in the folder `bin`*

---

## 🧩 Workflow structure

The software is composed of three main modules:

### 1. Preprocess
- Mesh reading from slab geometries (soon available interaction with [EFSM20 services](https://seismofaults.eu/services/efsm20-services))
- Area computation
- Connectivity and distances

### 2. Rupture areas
- Hazard mode → full coverage of mesh
- PTF mode → event-based scenarios

### 3. k223d module
- Refined distance computation
- Stochastic slip generation
- Variable rigidity support

---

## 📂 Inputs and configuration

The main configuration file is:

     config_files/Parameters/input.json

YAML input files are also supported for CLI runs, for example:

     config_files/Parameters/input.yaml

When a YAML file is used with `--input`, it is automatically converted to the corresponding JSON file before the run. This allows users to keep a commented, human-readable YAML file while preserving the internal JSON-based workflow.

Example:

    ./antifasc --input input.yaml

will read:

    config_files/Parameters/input.yaml

and create/use:

    config_files/Parameters/input.json

or, for a custom name:

    ./antifasc --input input_test.yaml

will create/use:

    config_files/Parameters/input_test.json

The input configuration controls:

- mesh input mode: existing mesh, GeoJSON mesh, or rectangular fault generation
- scaling laws and magnitude bins
- application mode: Hazard or Event-based PTF
- event information for PTF runs
- ensemble settings, including stochastic realizations, coupling limits, rupture-area filters, optional sub-boundaries, and variable-rigidity settings

In the notebook version, the input file can be created or updated through interactive widgets. If `input.json` is not present, the notebook starts from default values and creates the file when the user saves the input.

External optional CSV files are expected in:

     config_files/Rigidity/     # rigidity profiles: Depth,Rigidity
     config_files/Mesh/         # sub-boundaries: Longitude,Latitude

👉 Changes to input files are immediately available inside Docker because the repository is mounted into the container. No rebuild is required after editing input files.

**IMPORTANT: please have a look at the [Example1](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/Example1_Tohoku.md) and [Example2](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/Example2_Mediterranean.md) documentation for more details about the input files.**

---

## 📁 Outputs

Generated in:

     output/

Includes:
  - ASCII slip files
  - GeoJSON
  - HTML maps (optionally, look at `bin/plot_interactive_maps.ipynb`)

**IMPORTANT: please have a look at the [Example1](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/Example1_Tohoku.md) and [Example2](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/Example2_Mediterranean.md) documentation for
more details about the output folder and its folder tree**

---

## 🐳 How Docker works

The command:

    ./antifasc

runs pyANTI-FASc inside the Docker container in CLI mode.

The command:

    ./antifasc notebook

starts JupyterLab inside the same Docker image.

In both cases, the repository is mounted inside the container at `/app`. This ensures:

- full synchronization between container and host
- input changes immediately visible inside Docker
- outputs directly available locally
- reproducibility

---

## ⚠️ Notes

- On Windows, containers run as root → this does not affect file usability
- Port `8888` is used for Jupyter → change if already in use
- Docker must be running

---

## 🧪 Examples

See:

- [bin/antifasc_main.ipynb](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/bin/antifasc_main.ipynb)
- [bin/antifasc_main_Ex2.ipynb](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/bin/antifasc_main_Ex2.ipynb)

And the corresponding example descriptions:

- [Example1_Tohoku.md](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/Example1_Tohoku.md)
- [Example2_Mediterranean.md](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/Example2_Mediterranean.md)
   

---

## 🧱 Manual installation (Conda)

If you prefer a manual setup without Docker, see:

👉 [MANUAL_INSTALL.md](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/MANUAL_INSTALL.md)

⚠️ This approach requires manual dependency management and is not recommended for new users.

---

## 📚 Complementary Documentation

- Manual installation → `MANUAL_INSTALL.md`
- Interactive notebooks → included in `bin/`
- Wiki → under construction

---

## 🙏 Acknowledgements

A special thanks to Stefano Lorito, Fabrizio Romano, Manuela Volpe, Hafize Basak Bayraktar, Jacopo Selva, Gaetano Festa and Antonio Giovanni Iaccarino for participating at the different phases of conceiving, revising, developing and testing of the current version of the platform.

Thanks to Roberto Basili, Francesco Emanuele Maesano, Mara Monica Tiberti and Gareth Davies for their contribution in providing slab geometries.

The k223d module used in pyANTI-FASc was originally developed by Shane Murphy and Andre Herrero; See the original version [here](https://github.com/s-murfy/k223d). We gratefully acknowledge their fundamental contribution to this component of the workflow.


---



## 4 BIBLIOGRAPHY

Basili R. et al. (2021), The Making of the NEAM Tsunami Hazard Model 2018 (NEAMTHM18), Frontiers in Earth Science, DOI: 10.3389/feart.2020.616594 

Herrero and Murphy (2018), 	Self-similar slip distributions on irregular shaped faults, Geophysical Journal International, DOI: 10.1093/gji/ggy104

Maesano, F.E., Tiberti, M.M. and Basili, R. (2017), The Calabrian Arc: Three-dimensional modelling of the subduction interface, Scientific Reports DOI: 10.1038/s41598-017-09074-8

Murotani, S., Satake, K., and Fujii, Y. (2013), Scaling relations of seismic moment, rupture area, average slip, and asperity size for M~9 subduction-zone earthquakes, Geophysical Research Letters, 40, DOI: 10.1002/grl.50976.

Scala A. et al. (2020), Effect of Shallow Slip Amplification Uncertainty on Probabilistic Tsunami Hazard Analysis in Subduction Zones: Use of Long-Term Balanced Stochastic Slip Models, Pure and Applied Geophysics, DOI: 10.1007/s00024-019-02260-x

Selva, J., Lorito, S., Volpe, M. et al. (2021). Probabilistic tsunami forecasting for early warning. Nat Commun 12, 5677 (2021). DOI: 10.1038/s41467-021-25815-w

Strasser, F. O., Arango, M. C., & Bommer, J. J. (2010), Scaling of the source dimensions of interface and intraslab subduction-zone earthquakes with moment magnitude. Seismological Research Letters, DOI: 10.1785/gssrl.81.6.941.

Tonini et al. (2020), Importance of earthquake rupture geometry on tsunami modelling: The Calabrian Arc subduction interface (Italy) case study, Geophysical Journal International, DOI: 10.1093/gji/ggaa409
