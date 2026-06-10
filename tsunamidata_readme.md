# Running pyANTI-FASc on the TsunamiData Jupyter platform



## 1. Clone the TsunamiData branch

Open a JupyterLab terminal and run:

```bash
git clone --branch tsunamidata-platform https://github.com/antonioscalaunina/pyANTI-FASc.git
cd pyANTI-FASc
```


## 2. Create the micromamba environment

Create the `antifasc` environment from the provided environment file:

```bash
micromamba env create -f ANTIFASc.yml
```

Copy the `k223d.x` executable into the `antifasc` environment:

```bash
ENVPATH=$(micromamba run -n antifasc python -c "import sys; print(sys.prefix)")
cp bin/k223d.x "$ENVPATH/bin/k223d.x"
chmod +x "$ENVPATH/bin/k223d.x"
```


## 3. Register the Jupyter kernel

Register the environment as a Jupyter kernel:

```bash
micromamba run -n antifasc python -m ipykernel install --user --name antifasc --display-name "Python (pyANTI-FASc)"
```

Then reload JupyterLab.

## 4. Open the main notebook

Open the following notebook in JupyterLab:

```text
pyANTI-FASc/bin/antifasc_main.ipynb
```

Then select the correct kernel:

```text
Kernel → Change Kernel → Python (pyANTI-FASc)
```

If `Python (pyANTI-FASc)` is available, you can run the notebook.

## 5. Run pyANTI-FASc from the terminal

The code can also be launched directly from the terminal.

From the repository `bin` folder:

```bash
cd pyANTI-FASc/bin
```

run:

```bash
micromamba run -n antifasc python antifasc_main.py --input input.json
```

The input file must be available in the configuration folder and can be either in JSON or YAML format.
Please refer to the main documentation and examples for more details:

```text
config_files/Parameters/input.json
config_files/Parameters/input.yaml
```

## 6. Visualize outputs

Outputs are written in:

```text
output/
```

To inspect them, open:

```text
bin/interactive_slip_maps.ipynb
```

Then select:

```text
Kernel → Change Kernel → Python (pyANTI-FASc)
```

Use the notebook widgets to choose the output folder, magnitude, scaling law, background map and slip realization.

## 7. Download outputs

In the final section of the interactive notebook `interactive_slip_maps.ipynb`, you can select an output folder and create a `.tar.gz` archive.

The notebook prints the archive path, for example:

```text
/home/jovyan/work/pyANTI-FASc/bin/download_archives/<output_name>_YYYYMMDD_HHMMSS.tar.gz
```

To download it:

1. Open the JupyterLab file browser.
2. Go to the printed path.
3. Right-click the `.tar.gz` file.
4. Select **Download**.

## 8. Clean outputs from the platform

After checking that the archive has been downloaded, the user can utilize the cleanup button in the notebook.

Cleanup requires a second confirmation and removes both:

```text
output/<selected_output_folder>/
```

and the generated `.tar.gz` archive.
