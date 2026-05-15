# EXAMPLE 1 — Tohoku earthquake

This guide shows how to run a practical **Event-based PTF** pyANTI-FASc application for the 11 March 2011 Mw 9.0 Tohoku earthquake.

The goal is to generate stochastic slip distributions on the Kuril–Japan slab mesh, using magnitude bins around the target event magnitude and empirical rupture scaling laws.

👉 The final output consists of stochastic slip distributions ready to be used as initial conditions for tsunami simulations.

The same test case can also be configured and launched through the Jupyter notebook:

```text
bin/antifasc_main.ipynb
```

---

# 1 — Mesh

pyANTI-FASc provides an ensemble of predefined mesh discretizations representing the seismogenic portions of most subducting plates worldwide. These meshes are derived from geometries defined within the Slab 2.0 project and are available at this [webpage](https://www.sciencebase.gov/catalog/item/5aa1b00ee4b0b1c392e86467).

They are discretized using a relatively uniform node spacing, with a minimum inter-node distance of approximately 10–15 km.

All available meshes can be found [here](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main/utils/sz_slabs), while Figure 1 provides a global overview of their distribution.

![Map of the subducting plate meshes available in the current version of pyANTI-FASc](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/map_of_slabs.png)

*Figure 1 — Map of the subducting plate meshes available in the current version of pyANTI-FASc.*

For this example, the predefined mesh is:

```text
kurilsjapan
```

and is read from:

```text
utils/sz_slabs/kurilsjapan/subfaults/
```

---

# 2 — Input configuration

pyANTI-FASc can be configured using either:

```text
config_files/Parameters/input.json
```

or, for CLI runs, a YAML file such as:

```text
config_files/Parameters/input.yaml
```

YAML files are useful because they support comments. When a YAML input is provided to the CLI, pyANTI-FASc automatically converts it to the corresponding JSON file and then runs the code using the JSON input. The internal workflow therefore remains based on `input.json`.

The notebook interface can also create or update `input.json` interactively through widgets.

---

## 2.1 Recommended YAML input for this example

A compact and documented YAML configuration for the Tohoku example is shown below.

```yaml
# ============================================================
# Mesh input
# ============================================================

# 0 = use an existing mesh from utils/sz_slabs/<zone_name>/subfaults
# 1 = build mesh from utils/sz_slabs/<zone_name>_mesh.geojson
# 2 = build a rectangular fault from user-defined parameters
mesh_gen: 0

# Precomputed Kuril–Japan slab mesh
zone_name: kurilsjapan

# Short acronym used internally for mesh and output filenames
acronym: KuJ

# Mercator / UTM projection zone for the selected mesh
Merc_zone: 54

# Constant rake angle in degrees
rake: 90.0


# ============================================================
# Scaling laws and magnitude bins
# ============================================================

Scaling:

  # Scaling laws must match the names defined in bin/utils/scaling_laws.py
  laws:
    - Strasser2010_interface
    - Murotani2013

  # Magnitude bins used to generate rupture dimensions
  magnitude_bins:
    mode: range
    min: 8.8
    max: 9.2
    step: 0.1


# ============================================================
# Application
# ============================================================

# Event-based Probabilistic Tsunami Forecasting
application: PTF


# ============================================================
# Event information
# Required only for application: PTF
# ============================================================

Event:

  # Name used to build output folder names
  Name: Tohoku_test

  # Hypocenter coordinates [longitude, latitude]
  Hypo_LonLat:
    - 142.369
    - 38.322

  # Event magnitude
  Magnitude: 9.0


# ============================================================
# Ensemble configuration
# ============================================================

Configure:

  # Rupture shape: 
  # Rectangle = preserve the expected aspect ratio from scaling law)
  # Circle = pretty isotrope shae along strike and dip directions)
  shape: Rectangle

  # Number of stochastic slip distributions per rupture area
  numb_stoch: 2

  # Coupling limits in km
  coupling_shallow_limit: 1.0
  coupling_deep_limit: 55.0

  # Rupture-area selection parameters
  minimum_bnd_distance: 0.25
  minimum_interdistance: 0.1
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

  # PTF-only parameters
  # Magnitude bins in [Mw - Magnitude_lb, Mw + Magnitude_ub] are selected
  Magnitude_lb: 0.1
  Magnitude_ub: 0.1

  # Barycenters closer than hypo_baryc_distance * rupture length
  # from the hypocenter are selected
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
config_files/Parameters/input.json
```

In the current workflow, the scaling laws and magnitude bins are defined directly inside the `Scaling` block of the main input file. A separate `scaling_relationship.json` file is no longer required for the standard workflow.

A minimal JSON version of the same configuration is:

```json
{
  "mesh_gen": 0,
  "zone_name": "kurilsjapan",
  "acronym": "KuJ",
  "Merc_zone": 54,
  "rake": 90.0,

  "Scaling": {
    "laws": [
      "Strasser2010_interface",
      "Murotani2013"
    ],
    "magnitude_bins": {
      "mode": "range",
      "min": 8.8,
      "max": 9.2,
      "step": 0.1
    }
  },

  "application": "PTF",

  "Event": {
    "Name": "Tohoku_test",
    "Hypo_LonLat": [
      142.369,
      38.322
    ],
    "Magnitude": 9.0
  },

  "Configure": {
    "shape": "Rectangle",
    "numb_stoch": 5,
    "coupling_shallow_limit": 1.0,
    "coupling_deep_limit": 55.0,
    "minimum_bnd_distance": 0.25,
    "minimum_interdistance": 0.1,
    "Fact_area_scaling": 1.0,
    "variable_mu": 1,
    "Rigidity_file_logic": 0,
    "Stress_drop_var": 0,
    "Fact_rigidity": 0.5,
    "Magnitude_lb": 0.1,
    "Magnitude_ub": 0.1,
    "hypo_baryc_distance": 1.0,
    "mesh_sub_boundary": 0
  }
}
```

---

## 2.3 Notebook configuration

The same input can be configured interactively using:

```text
bin/antifasc_main.ipynb
```

The notebook widget lets the user select:

- mesh source;
- scaling laws;
- magnitude binning;
- application mode;
- event information for Event-based PTF;
- ensemble parameters;
- optional rigidity and sub-boundary files.

The selected options are written to `input.json` before the `Slab` object is created.

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
./antifasc --input input.yaml
```

This searches for:

```text
config_files/Parameters/input.yaml
```

converts it to:

```text
config_files/Parameters/input.json
```

or to the corresponding JSON name, and then runs the code.

Using a custom JSON input:

```bash
./antifasc --input input_tohoku.json
```

This searches for:

```text
config_files/Parameters/input_tohoku.json
```

**Look at [README.md](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/README.md) to see how to run the same pipeline of Windows Powershell**

## Docker notebook run

```bash
./antifasc notebook
```

Then open the JupyterLab URL printed in the terminal and run:

```text
bin/antifasc_main.ipynb
```

**Look at [README.md](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/README.md) to see how to run the same pipeline of Windows Powershell**

## Manual installation

If running without Docker:

```bash
conda activate antifasc
cd bin
python antifasc_main.py
```

or with a custom input:

```bash
python antifasc_main.py --input input.yaml
```

---

# 4 — Expected screen output

During the run, pyANTI-FASc prints the main processing steps.

First, the input file and scaling laws are read:

```text
reading input.json file
reading scaling laws from input.json
Great! You already have nodes and cells, I'm just writing
reading mesh
```

Then the barycenters are selected. Since this is an Event-based PTF example, only magnitude bins around the event magnitude are used:

```text
Barycenter selection
Magnitude bin # ... - Mw=...
Barycenter selection (PTF)
Magnitude bin # ... - Mw=...
```

The selected magnitudes depend on the magnitude binning and on:

```text
Magnitude_lb
Magnitude_ub
```

For this example, the event magnitude is Mw 9.0 and the selected bins are those within the configured tolerance.

After that, rupture areas are computed for each selected magnitude bin and scaling law:

```text
Rupturing area computation
Magnitude bin # ... - Mw=...
Rupturing areas computed!
```

The code then writes the rupture-area inputs and computes the stochastic slip distributions:

```text
Writing Output
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
output/Tohoku_test_M90_E14237_N3832_slip_KuJ/
```

with a structure like:

```text
output/
└── Tohoku_test_M90_E14237_N3832_slip_KuJ
    ├── homogeneous_mu
    │   ├── 8_9000
    │   │   ├── Strasser2010_interface
    │   │   │   ├── Slip4HySea00004_001.dat
    │   │   │   ├── Slip4HySea00004_002.dat
    │   │   │   └── ...
    │   │   └── Murotani2013
    │   │       └── ...
    │   └── ...
    └── variable_mu
        ├── 8_9000
        │   ├── Strasser2010_interface
        │   └── Murotani2013
        └── ...
```

The number of generated slip distributions is controlled by:

```text
Configure.numb_stoch
```

For each rupture area, `numb_stoch` stochastic slip distributions are produced.

---

# 6 — Slip file format

The `Slip4HySea*.dat` files are written in the standard format used as input by [Tsunami-HySea](https://edanya.uma.es/hysea/), one of the most widely used tsunami simulators in the community.

Each row describes one triangular subfault element and contains:

```text
LON1 LAT1 DEPTH1(km) LON2 LAT2 DEPTH2(km) LON3 LAT3 DEPTH3(km) RAKE SLIP(m)
```

Example:

```text
142.209106 36.309288 7.494581 142.289078 36.401661 7.454205 142.143066 36.442127 10.402809 90.000000 7.256678
142.209106 36.309288 7.494581 142.336136 36.256641 5.620936 142.289078 36.401661 7.454205 90.000000 7.095560
142.289078 36.401661 7.454205 142.223038 36.534504 10.362430 142.143066 36.442127 10.402809 90.000000 6.840256
```

For each `Slip4HySea*.dat` file, a corresponding ***GeoJSON*** file containing the same distribution is also produced. It can be promptly used for visualization. See next section

---

# 7 — Post-process

The slip distributions can be plotted using personal scripts, GIS tools, or web services supporting GeoJSON.

An interactive notebook is available at:

```text
bin/interactive_slip_maps.ipynb
```

It allows the user to:

- select one of the output folders;
- visualize slip and rake maps from the GeoJSON files;
- optionally generate and export interactive HTML maps.

![screenshot](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/Screenshot_interactive_plot_JN.png)

*Figure 2 — Screenshot showing the interactive slip plotter Jupyter Notebook.*

