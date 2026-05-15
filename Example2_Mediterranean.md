# EXAMPLE 2 — Mediterranean Hazard simulation (EFSM20)

This guide shows how to run a practical **Hazard** pyANTI-FASc application for an offshore Sicily fault system in the central Mediterranean.

The goal is to generate a large ensemble of stochastic slip distributions on a fault mesh derived from the EFSM20 mesh service, using all selected magnitude bins and a normal-faulting scaling law.

👉 The final output consists of stochastic slip distributions ready to be used as initial conditions for tsunami simulations.

The same test case can also be configured and launched through the Jupyter notebook:

```text
bin/antifasc_main_Ex2.ipynb
```

The notebook version includes additional widgets that allow the user to visualize and verify intermediate steps, such as barycenter selection and rupture-area definition. For large Hazard applications, the CLI run is usually faster.

---

# 1 — Mesh

pyANTI-FASc can use meshes made available through the European Database of Seismogenic Faults services ([EFSM20](https://seismofaults.eu/component/tags/tag/efsm20)).

In particular, this example uses a mesh provided by the *EFSM20 Meshes* service. Details are available [here](https://seismofaults.eu/services/efsm20-services). The WFS service can also be accessed through the GetCapabilities XML file:

```text
https://services.seismofaults.eu/EFSM20_Meshes/ows?service=WFS&request=GetCapabilities
```

If you use this service, please cite:

```text
https://doi.org/10.13127/efsm20/meshes
```

Figure 1 shows a screenshot of some of the meshed faults available within the EFSM20 database.

![Screenshot of EFSM mesh database imported to QGis](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/Image_qgis.jpg)

*Figure 1 — Screenshot of EFSM20 mesh database imported into QGIS.*

For this example, the mesh is:

```text
ITCF00G
```

and the starting GeoJSON file is expected at:

```text
utils/sz_slabs/ITCF00G_mesh.geojson
```

Since this is a GeoJSON mesh, the input uses:

```text
mesh_gen: 1
```

During the first run, pyANTI-FASc reads the GeoJSON file and creates the corresponding nodes and faces files in:

```text
utils/sz_slabs/ITCF00G/subfaults/
```

After that, the same mesh can also be reused with:

```text
mesh_gen: 0
```

---

# 2 — Input configuration

pyANTI-FASc can be configured using either:

```text
config_files/Parameters/input.json
```

or, for CLI runs, a YAML file such as:

```text
config_files/Parameters/input_Sicily.yaml
```

YAML files are useful because they support comments. When a YAML input is provided to the CLI, pyANTI-FASc automatically converts it to the corresponding JSON file and then runs the code using the JSON input. The internal workflow therefore remains based on `input.json`.

The notebook interface can also create or update `input.json` interactively through widgets.

---

## 2.1 Recommended YAML input for this example

A compact and documented YAML configuration for the Mediterranean Hazard example is shown below.

```yaml
# ============================================================
# Mesh input
# ============================================================

# 0 = use an existing mesh from utils/sz_slabs/<zone_name>/subfaults
# 1 = build mesh from utils/sz_slabs/<zone_name>_mesh.geojson
# 2 = build a rectangular fault from user-defined parameters
mesh_gen: 1

# EFSM20 mesh name.
# The file utils/sz_slabs/ITCF00G_mesh.geojson must exist.
zone_name: ITCF00G

# Short acronym used internally for mesh and output filenames
acronym: ITC

# Mercator / UTM projection zone for the selected mesh
Merc_zone: 33

# Constant rake angle in degrees.
# For this normal-faulting example, rake is set to -90.
rake: -90.0


# ============================================================
# Scaling laws and magnitude bins
# ============================================================

Scaling:

  # Scaling laws must match the names defined in bin/utils/scaling_laws.py.
  # This example uses a Wells & Coppersmith normal-faulting relation.
  laws:
    - WC1994_normal

  # Magnitude bins used to generate rupture dimensions
  magnitude_bins:
    mode: range
    min: 5.5
    max: 7.5
    step: 0.1


# ============================================================
# Application
# ============================================================

# Hazard mode uses all selected magnitude bins and all eligible barycenters.
application: Hazard


# ============================================================
# Event information
# Not used when application: Hazard
# ============================================================

Event:

  # This block is ignored in Hazard mode, but can be left in the file.
  Name: Sicily_test
  Hypo_LonLat:
    - 12.8
    - 38.5
  Magnitude: 7.0


# ============================================================
# Ensemble configuration
# ============================================================

Configure:

  # Rupture shape:
  # Rectangle = preserve the expected aspect ratio from the scaling law
  # Circle = more isotropic rupture area along strike and dip directions
  shape: Rectangle

  # Number of stochastic slip distributions per rupture area
  numb_stoch: 5

  # Coupling limits in km
  coupling_shallow_limit: 1.0
  coupling_deep_limit: 55.0

  # Rupture-area selection parameters.
  # These values reduce the number of very similar rupture areas,
  # especially at the largest magnitude bins.
  minimum_bnd_distance: 0.1
  minimum_interdistance: 0.2
  Fact_area_scaling: 1.0

  # Variable rigidity workflow
  # 0 = homogeneous rigidity only
  # 1 = also compute variable-rigidity slip distributions
  variable_mu: 1

  # 0 = use default rigidity model
  # 1 = use a CSV file from config_files/Rigidity
  Rigidity_file_logic: 0

  # Stress-drop variation flag
  Stress_drop_var: 0

  # Factor used by the default rigidity model
  Fact_rigidity: 0.5

  # PTF-only parameters.
  # They are ignored in Hazard mode but can be left in the file.
  Magnitude_lb: 0.3
  Magnitude_ub: 0.3
  hypo_baryc_distance: 1.0

  # Optional sub-boundary
  # 0 = use the full inferred mesh boundary
  # 1 = use a CSV boundary file from config_files/Mesh
  mesh_sub_boundary: 0
```

---

## 2.2 Equivalent JSON input

The same configuration can also be provided directly as JSON in:

```text
config_files/Parameters/input_Sicily.json
```

In the current workflow, the scaling laws and magnitude bins are defined directly inside the `Scaling` block of the main input file. A separate `scaling_relationship.json` file is no longer required for the standard workflow.

A minimal JSON version of the same configuration is:

```json
{
  "mesh_gen": 1,
  "zone_name": "ITCF00G",
  "acronym": "ITC",
  "Merc_zone": 33,
  "rake": -90.0,

  "Scaling": {
    "laws": [
      "WC1994_normal"
    ],
    "magnitude_bins": {
      "mode": "range",
      "min": 5.5,
      "max": 7.5,
      "step": 0.1
    }
  },

  "application": "Hazard",

  "Event": {
    "Name": "Sicily_test",
    "Hypo_LonLat": [
      12.8,
      38.5
    ],
    "Magnitude": 7.0
  },

  "Configure": {
    "shape": "Rectangle",
    "numb_stoch": 5,
    "coupling_shallow_limit": 1.0,
    "coupling_deep_limit": 55.0,
    "minimum_bnd_distance": 0.1,
    "minimum_interdistance": 0.2,
    "Fact_area_scaling": 1.0,
    "variable_mu": 1,
    "Rigidity_file_logic": 0,
    "Stress_drop_var": 0,
    "Fact_rigidity": 0.5,
    "Magnitude_lb": 0.3,
    "Magnitude_ub": 0.3,
    "hypo_baryc_distance": 1.0,
    "mesh_sub_boundary": 0
  }
}
```

---

## 2.3 Notebook configuration

The same input can be configured interactively using:

```text
bin/antifasc_main_Ex2.ipynb
```

The notebook widget lets the user select:

- mesh source;
- scaling laws;
- magnitude binning;
- application mode;
- ensemble parameters;
- optional rigidity and sub-boundary files.

For this example, the application is **Hazard**, so the event-related and PTF-only options are not used. The notebook can still show additional widgets to inspect intermediate steps, such as:

- barycenter selection;
- rupture-area definition;
- slip-distribution plotting.

These checks are useful for understanding and validating the workflow.

---

# 3 — Run pyANTI-FASc

Once the input file is ready, the run can be launched from the repository root.

## Docker CLI run

Using the default input:

```bash
./antifasc
```

This uses:

```text
config_files/Parameters/input.json
```

Using a YAML input:

```bash
./antifasc --input input_Sicily.yaml
```

This searches for:

```text
config_files/Parameters/input_Sicily.yaml
```

converts it to:

```text
config_files/Parameters/input_Sicily.json
```

and then runs the code.

Using a custom JSON input:

```bash
./antifasc --input input_Sicily.json
```

This searches for:

```text
config_files/Parameters/input_Sicily.json
```

**Look at [README.md](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/README.md) to see how to run the same pipeline on Windows PowerShell.**

## Docker notebook run

```bash
./antifasc notebook
```

Then open the JupyterLab URL printed in the terminal and run:

```text
bin/antifasc_main_Ex2.ipynb
```
**Look at [README.md](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/README.md) to see how to run the same pipeline on Windows PowerShell.**

## Manual installation

If running without Docker:

```bash
conda activate antifasc
cd bin
python antifasc_main.py
```

or with a custom input:

```bash
python antifasc_main.py --input input_Sicily.yaml
```


---

# 4 — Expected screen output

During the run, pyANTI-FASc prints the main processing steps.

First, the input file and scaling laws are read. Since this example starts from a GeoJSON mesh, the first run also creates the corresponding nodes and faces files:

```text
reading input.json file
reading scaling laws from input.json
../utils/sz_slabs/ITCF00G_mesh.geojson
Great! I really love to create mesh files from GeoJSON format

Files ../utils/sz_slabs/ITCF00G/subfaults/ITCF00G_mesh_nodes.dat and ../utils/sz_slabs/ITCF00G/subfaults/ITCF00G_mesh_faces.dat created successfully inside subfaults!
reading mesh
```

This means that from now on the same mesh can also be reused with:

```text
mesh_gen: 0
```

because the mesh is now available in the `utils/sz_slabs` database.

Then the barycenters are selected. Since this is a Hazard example, all selected magnitude bins are used:

```text
Barycenter selection
Magnitude bin # 0 - Mw=5.5000
Magnitude bin # 1 - Mw=5.6000
Magnitude bin # 2 - Mw=5.7000
...
Magnitude bin # 20 - Mw=7.5000
Mw: [5.5 5.6 5.7 ... 7.5]
Scaling names: ['WC1994_normal']
```

After that, rupture areas are computed for each selected magnitude bin and scaling law:

```text
Computing Rupturing areas: 100%|████████████████████████████████| 21/21
Mw=5.5, Name scaling: WC1994_normal, N=...
Mw=5.6, Name scaling: WC1994_normal, N=...
...
Mw=7.5, Name scaling: WC1994_normal, N=...
Writing Output
```

The code then writes the rupture-area inputs and computes the stochastic slip distributions:

```text
Computing slip distributions for the homogeneous and variable rigidity cases
```

If `variable_mu` is set to `1`, both homogeneous and variable-rigidity slip distributions are generated. If it is set to `0`, only the homogeneous case is computed.

---

# 5 — Output structure

When the process is complete, the output is organized under:

```text
output/
```

For this example, the main output folder will have a name similar to:

```text
output/ITCF00G_Hazard_slip_ITC/
```

with a structure like:

```text
output/
└── ITCF00G_Hazard_slip_ITC
    ├── homogeneous_mu
    │   ├── 5_5000
    │   │   └── WC1994_normal
    │   │       ├── output_file.log
    │   │       ├── Slip4HySea00001_000.dat
    │   │       ├── Slip4HySea00001_000.json
    │   │       └── ...
    │   └── ...
    └── variable_mu
        ├── 5_5000
        │   └── WC1994_normal
        │       ├── Slip4HySea00001_000.dat
        │       ├── Slip4HySea00001_000.json
        │       └── ...
        └── ...
```

The number of generated slip distributions is controlled by:

```text
Configure.numb_stoch
```

For each rupture area, `numb_stoch` stochastic slip distributions are produced.

For smaller magnitudes, when there are not enough cells to build fully stochastic slip distributions, homogeneous or rigidity-modulated slip distributions may be produced with file names such as:

```text
Slip4HySeaXXXXX_000.dat
```

For larger magnitudes, multiple stochastic distributions are produced for each rupture area, and the second index indicates the stochastic realization number:

```text
Slip4HySeaXXXXX_001.dat
Slip4HySeaXXXXX_002.dat
...
```

The run also produces corresponding GeoJSON/JSON files for the generated slip distributions. These files can be used directly for visualization.

---

# 6 — Slip file format

The `Slip4HySea*.dat` files are written in the standard format used as input by [Tsunami-HySea](https://edanya.uma.es/hysea/), one of the most widely used tsunami simulators in the community.

Each row describes one triangular subfault element and contains:

```text
LON1 LAT1 DEPTH1(km) LON2 LAT2 DEPTH2(km) LON3 LAT3 DEPTH3(km) RAKE SLIP(m)
```

Example:

```text
12.869199 38.554981 8.857142 12.812249 38.546745 8.857142 12.850728 38.512363 11.142860 -90.000000 9.104819
12.812249 38.546745 8.857142 12.869199 38.554981 8.857142 12.807246 38.586212 6.571429 -90.000000 8.936686
12.812249 38.546745 8.857142 12.794422 38.503761 11.142860 12.850728 38.512363 11.142860 -90.000000 9.739596
```

For each `Slip4HySea*.dat` file, a corresponding ***GeoJSON*** file containing the same distribution is also produced. It can be promptly used for visualization. See next section.

---

# 7 — Post-process

The slip distributions can be plotted using personal scripts, GIS tools, or web services supporting GeoJSON.

The GeoJSON files can be uploaded to QGIS, as shown in Figure 2, or to any web service using the GeoJSON standard, such as [kepler.gl](https://kepler.gl/).

![Screenshot of a slip distributions imported to QGis](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/Image_qgis2.jpg)

*Figure 2 — Screenshot of a slip distribution imported into QGIS.*

An interactive notebook is available at:

```text
bin/interactive_slip_maps.ipynb
```

It allows the user to:

- select one of the output folders;
- visualize slip and rake maps from the GeoJSON files;
- optionally generate and export interactive HTML maps.

![screenshot](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/Screenshot_interactive_plot_JN_Sicily.png)

*Figure 3 — Screenshot showing the interactive slip plotter Jupyter Notebook.*

Within the selected folder, for each `Slip4HySea*.dat` file, the notebook can also generate an interactive HTML map, such as this [example](https://antonioscalaunina.github.io/pyANTI-FASc/utils/Slip4HySea00136_005.html).
