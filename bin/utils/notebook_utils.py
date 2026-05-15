import os
import json
import re
import numpy as np

import ipywidgets as widgets
from IPython.display import display, clear_output

from .scaling_laws import SCALING_LAWS


def get_default_input_config():
    return {
        "zone_name": "",
        "Merc_zone": 33,
        "acronym": "DUM",
        "mesh_gen": 0,
        "rake": 90.0,
        "Scaling": {
            "laws": [list(SCALING_LAWS.keys())[0]],
            "magnitude_bins": {
                "mode": "range",
                "min": 8.8,
                "max": 9.2,
                "step": 0.1
            }
        },
        "application": "PTF",
        "Event": {
            "Name": "event_name",
            "Hypo_LonLat": [0.0, 0.0],
            "Magnitude": 9.0
        },
        "Configure": {
            "shape": "Rectangle",
            "numb_stoch": 1,
            "variable_mu": 0,
            "coupling_shallow_limit": 1.0,
            "coupling_deep_limit": 55.0,
            "mesh_sub_boundary": 0,
            "Magnitude_lb": 0.3,
            "Magnitude_ub": 0.3,
            "hypo_baryc_distance": 1.0,
            "minimum_bnd_distance": 0.25,
            "minimum_interdistance": 0.2,
            "Fact_area_scaling": 1.0,
            "Rigidity_file_logic": 0,
            "Stress_drop_var": 0,
            "Fact_rigidity": 0.5
        }
    }


def create_input_widgets(slab_module, input_file):
    """
    Create widgets to select scaling laws and magnitude bins,
    update input.json, and initialize Slab_obj.

    Parameters
    ----------
    slab_module : module
        The imported utils.slab module.
    input_file : str
        Path to input.json.

    Returns
    -------
    widgets.Output
        Output widget where messages are printed.
    """

    # -----------------------------
    # Mesh widgets
    # -----------------------------
    
    if os.path.isfile(input_file):
        with open(input_file) as f:
            current_cfg = json.load(f)
    else:
        current_cfg = get_default_input_config()

    existing_meshes = list_existing_slab_meshes()
    geojson_meshes = list_geojson_meshes()

    current_zone_name = current_cfg.get("zone_name", "")
    current_mesh_gen = current_cfg.get("mesh_gen", 0)
    

    mesh_gen_widget = widgets.ToggleButtons(
        options=[
            ("Existing mesh", 0),
            ("GeoJSON mesh", 1),
            ("Rectangular fault", 2),
        ],
        value=current_mesh_gen if current_mesh_gen in [0, 1, 2] else 0,
        description="Mesh",
        layout=widgets.Layout(width="650px")
    )

    existing_mesh_widget = widgets.Dropdown(
        options=existing_meshes if existing_meshes else ["NO_AVAILABLE_MESH"],
        value=current_zone_name if current_zone_name in existing_meshes else (
            existing_meshes[0] if existing_meshes else "NO_AVAILABLE_MESH"
        ),
        description="Slab",
        layout=widgets.Layout(width="450px")
    )

    geojson_mesh_widget = widgets.Dropdown(
        options=geojson_meshes if geojson_meshes else ["NO_AVAILABLE_GEOJSON"],
        value=current_zone_name if current_zone_name in geojson_meshes else (
            geojson_meshes[0] if geojson_meshes else "NO_AVAILABLE_GEOJSON"
        ),
        description="GeoJSON",
        layout=widgets.Layout(width="450px")
    )

    lon_c = widgets.FloatText(
        value=float(current_cfg.get("lon_c", 0.0)),
        description="lon_c"
    )

    lat_c = widgets.FloatText(
        value=float(current_cfg.get("lat_c", 0.0)),
        description="lat_c"
    )

    depth_km = widgets.FloatText(
        value=float(current_cfg.get("depth_km", 10.0)),
        description="depth_km"
    )

    length_km = widgets.FloatText(
        value=float(current_cfg.get("length_km", 50.0)),
        description="length_km"
    )

    width_km = widgets.FloatText(
        value=float(current_cfg.get("width_km", 25.0)),
        description="width_km"
    )

    strike = widgets.FloatText(
        value=float(current_cfg.get("strike", 0.0)),
        description="strike"
    )

    dip = widgets.FloatText(
        value=float(current_cfg.get("dip", 30.0)),
        description="dip"
    )

    rake = widgets.FloatText(
        value=float(current_cfg.get("rake", 90.0)),
        description="rake"
    )

    elem_size_km2 = widgets.FloatText(
        value=float(current_cfg.get("elem_size_km2", 5.0)),
        description="elem_size_km2"
    )

    merc_zone = widgets.IntText(
    value=int(current_cfg.get("Merc_zone", 33)),
    description="Mercator projection zone",
    layout=widgets.Layout(width="320px"),
    style={"description_width": "180px"}
    )

    

    acronym = widgets.Text(
        value=str(current_cfg.get("acronym", "DUM")),
        description="acronym"
    )

    zone_name_custom = widgets.Text(
        value=str(current_zone_name) if current_zone_name else "custom_fault",
        description="zone_name"
    )

    existing_mesh_box = widgets.VBox([
        widgets.HTML(
            "<b>Existing mesh:</b> select one slab already available in "
            "<code>../utils/sz_slabs/&lt;name&gt;/subfaults</code>."
        ),
        existing_mesh_widget
    ])

    geojson_mesh_box = widgets.VBox([
        widgets.HTML(
            "<b>GeoJSON mesh:</b> select one <code>*_mesh.geojson</code> "
            "file available in <code>../utils/sz_slabs</code>."
        ),
        geojson_mesh_widget
    ])

    rectangular_fault_box = widgets.VBox([
        widgets.HTML(
            "<b>Rectangular fault:</b> define the center, geometry and mesh resolution."
        ),
        zone_name_custom,
        lon_c,
        lat_c,
        depth_km,
        length_km,
        width_km,
        strike,
        dip,
        elem_size_km2
    ])

    mesh_common_box = widgets.VBox([
        widgets.HTML(
            "<b>Projection and output identifiers</b>"
        ),
        merc_zone,
        acronym,
        rake
    ])
    
    def update_mercator_zone_from_selection(*args):
        mode = mesh_gen_widget.value
    
        inferred_zone = None
    
        if mode == 0:
            zone_name = existing_mesh_widget.value
    
            if zone_name not in ["NO_AVAILABLE_MESH", None]:
                inferred_zone = infer_mercator_zone_from_existing_mesh(zone_name)
    
        elif mode == 1:
            zone_name = geojson_mesh_widget.value
    
            if zone_name not in ["NO_AVAILABLE_GEOJSON", None]:
                inferred_zone = infer_mercator_zone_from_geojson(zone_name)
    
        elif mode == 2:
            inferred_zone = infer_mercator_zone(lon_c.value)
    
        if inferred_zone is not None:
            merc_zone.value = int(inferred_zone) 

            
    def update_mesh_widgets(*args):
        mode = mesh_gen_widget.value

        existing_mesh_box.layout.display = "none"
        geojson_mesh_box.layout.display = "none"
        rectangular_fault_box.layout.display = "none"

        if mode == 0:
            existing_mesh_box.layout.display = "flex"

            selected_name = existing_mesh_widget.value
            acronym.value = make_default_acronym(selected_name)

        elif mode == 1:
            geojson_mesh_box.layout.display = "flex"

            selected_name = geojson_mesh_widget.value
            acronym.value = make_default_acronym(selected_name)

        elif mode == 2:
            rectangular_fault_box.layout.display = "flex"

            if acronym.value in ["", "DUM"]:
                acronym.value = make_default_acronym(zone_name_custom.value)

            #merc_zone.value = infer_mercator_zone(lon_c.value)
        
        update_mercator_zone_from_selection()

    mesh_gen_widget.observe(update_mesh_widgets, names="value")

    def update_existing_mesh_defaults(*args):
        if mesh_gen_widget.value == 0:
            acronym.value = make_default_acronym(existing_mesh_widget.value)

    def update_geojson_mesh_defaults(*args):
        if mesh_gen_widget.value == 1:
            acronym.value = make_default_acronym(geojson_mesh_widget.value)

    def update_rectangular_fault_defaults(*args):
        if mesh_gen_widget.value == 2:
            merc_zone.value = infer_mercator_zone(lon_c.value)
            if acronym.value in ["", "DUM"]:
                acronym.value = make_default_acronym(zone_name_custom.value)

    existing_mesh_widget.observe(update_mercator_zone_from_selection, names="value")
    geojson_mesh_widget.observe(update_mercator_zone_from_selection, names="value")
    lon_c.observe(update_mercator_zone_from_selection, names="value")
    zone_name_custom.observe(update_rectangular_fault_defaults, names="value")

    update_mesh_widgets()

    # -----------------------------
    # Scaling and Mw widgets
    # -----------------------------

    scaling_cfg = current_cfg.get("Scaling", {})

    current_laws = scaling_cfg.get(
        "laws",
        [list(SCALING_LAWS.keys())[0]]
    )
    
    magnitude_bins_cfg = scaling_cfg.get("magnitude_bins", {})
    
    current_mw_mode = magnitude_bins_cfg.get("mode", "range")
    
    law_checkboxes = [
        widgets.Checkbox(
            value=(name in current_laws),
            description=name,
            indent=False
        )
        for name in SCALING_LAWS.keys()
    ]

    law_box = widgets.GridBox(
    law_checkboxes,
    layout=widgets.Layout(
        grid_template_columns="repeat(4, 250px)",
        grid_gap="5px 20px"
    )
)

    mode_widget = widgets.ToggleButtons(
        options=["range", "values"],
        value=current_mw_mode if current_mw_mode in ["range", "values"] else "range",
        #description="Mw bins"
    )

    mw_min = widgets.FloatText(
        value=float(magnitude_bins_cfg.get("min", 8.8)),
        description="Mw min"
    )
    
    mw_max = widgets.FloatText(
        value=float(magnitude_bins_cfg.get("max", 9.2)),
        description="Mw max"
    )
    
    mw_step = widgets.FloatText(
        value=float(magnitude_bins_cfg.get("step", 0.1)),
        description="Step"
    )

    current_mw_values = magnitude_bins_cfg.get("values", [8.8, 8.9, 9.0])

    mw_values = widgets.Text(
        value=", ".join(str(x) for x in current_mw_values),
        description="Mw list",
        placeholder="es: 8.8, 8.9, 9.0"
    )

    range_box = widgets.VBox([
        widgets.HTML(
            "<b>range:</b> <i>select</i> <b>Mw min</b>, "
            "<b>Mw max</b> <i>and</i> <b>Step</b>"
        ),
        mw_min,
        mw_max,
        mw_step
    ])

    values_box = widgets.VBox([
        widgets.HTML(
            "<b>values:</b> <i>use comma-separated values in</i> "
            "<b>Mw list</b> <i>(e.g., 8.8, 8.9, 9.0)</i>"
        ),
        mw_values
    ])


    def update_magnitude_widgets(*args):
        if mode_widget.value == "range":
            range_box.layout.display = "flex"
            values_box.layout.display = "none"
        elif mode_widget.value == "values":
            range_box.layout.display = "none"
            values_box.layout.display = "flex"

    mode_widget.observe(update_magnitude_widgets, names="value")
    update_magnitude_widgets()

    state = {}

    def run_from_widgets(_):
        with out:
            clear_output()

            mesh_mode = mesh_gen_widget.value

            if mesh_mode == 0:
                zone_name = existing_mesh_widget.value

                mesh_block = {
                    "zone_name": zone_name,
                    "mesh_gen": 0,
                    "Merc_zone": int(merc_zone.value),
                    "acronym": acronym.value or make_default_acronym(zone_name),
                    "rake": float(rake.value),
                }

            elif mesh_mode == 1:
                zone_name = geojson_mesh_widget.value

                mesh_block = {
                    "zone_name": zone_name,
                    "mesh_gen": 1,
                    "Merc_zone": int(merc_zone.value),
                    "acronym": acronym.value or make_default_acronym(zone_name),
                    "rake": float(rake.value),
                }

            elif mesh_mode == 2:
                zone_name = zone_name_custom.value or "custom_fault"

                mesh_block = {
                    "zone_name": zone_name,
                    "mesh_gen": 2,
                    "Merc_zone": int(merc_zone.value),
                    "acronym": acronym.value or make_default_acronym(zone_name),
                    "lon_c": float(lon_c.value),
                    "lat_c": float(lat_c.value),
                    "depth_km": float(depth_km.value),
                    "length_km": float(length_km.value),
                    "width_km": float(width_km.value),
                    "strike": float(strike.value),
                    "dip": float(dip.value),
                    "rake": float(rake.value),
                    "elem_size_km2": float(elem_size_km2.value),
                }

            update_input_json_mesh(input_file, mesh_block)
            
            selected_laws = [
                cb.description
                for cb in law_checkboxes
                if cb.value
            ]

            if len(selected_laws) == 0:
                print("ERROR: select at least one scaling law.")
                return

            scaling_block = slab_module.build_scaling_block(
                laws=selected_laws,
                mode=mode_widget.value,
                mw_min=mw_min.value,
                mw_max=mw_max.value,
                mw_step=mw_step.value,
                mw_values=mw_values.value
            )

            update_input_json_scaling(input_file, scaling_block)


            application = application_widget.value

            if application == "PTF":
                event_block = {
                    "Name": event_name.value,
                    "Hypo_LonLat": [
                        float(hypo_lon.value),
                        float(hypo_lat.value)
                    ],
                    "Magnitude": float(event_mw.value)
                }
            else:
                event_block = None
            
            update_input_json_application(
                input_file,
                application,
                event_block=event_block
            )

            print("Scaling saved in input.json")
            print(scaling_block)
            print("Application:", application)

            if application == "PTF":
                print("Event block:")
                print(event_block)

            configure_block = {
                "shape": shape_widget.value,
                "numb_stoch": int(numb_stoch_widget.value),
                "variable_mu": int(variable_mu_widget.value),

                "coupling_shallow_limit": float(coupling_shallow_widget.value),
                "coupling_deep_limit": float(coupling_deep_widget.value),

                "mesh_sub_boundary": int(sub_boundary_logic.value),

                "Magnitude_lb": float(magnitude_lb_widget.value),
                "Magnitude_ub": float(magnitude_ub_widget.value),
                "hypo_baryc_distance": float(hypo_baryc_distance_widget.value),

                "minimum_bnd_distance": float(minimum_bnd_distance_widget.value),
                "minimum_interdistance": float(minimum_interdistance_widget.value),
                "Fact_area_scaling": float(fact_area_scaling_widget.value),

                "Rigidity_file_logic": int(rigidity_file_logic.value) if variable_mu_widget.value else 0,
                "Stress_drop_var": int(stress_drop_widget.value) if variable_mu_widget.value else 0,
                "Fact_rigidity": float(fact_rigidity.value),
            }

            if sub_boundary_logic.value:
                if sub_boundary_file.value == "NO_SUB_BOUNDARY_CSV_FOUND":
                    print("ERROR: no sub-boundary CSV file found in ../config_files/Mesh")
                    return

                configure_block["sub_boundary_file"] = sub_boundary_file.value

            if variable_mu_widget.value and rigidity_file_logic.value:
                if rigidity_file.value == "NO_RIGIDITY_CSV_FOUND":
                    print("ERROR: no rigidity CSV file found in ../config_files/Rigidity")
                    return

                configure_block["Rigidity_file"] = rigidity_file.value
            
            update_input_json_configure(input_file, configure_block)
            
            print("Configure block updated:")
            print(configure_block)

            state["Slab_obj"] = slab_module.Slab(input_file)

            print("Slab created!")
            print("All magnitude bins:", state["Slab_obj"].Magnitude)
            print("Scaling:", state["Slab_obj"].Name_scaling)
            print("Application:", state["Slab_obj"].application)
            print("Selected magnitude bins:", state["Slab_obj"].get_magnitudes())
            


    # -----------------------------
    # Application widgets
    # -----------------------------
    
    current_application = current_cfg.get(
        "application",
        current_cfg.get("Configure", {}).get("application", "PTF")
    )
    
    event_cfg = current_cfg.get("Event", {})
    
    application_widget = widgets.ToggleButtons(
        options=[
            ("Hazard", "Hazard"),
            ("Event-based PTF", "PTF"),
        ],
        value=current_application if current_application in ["PTF", "Hazard"] else "PTF",
        description="Application",
        layout=widgets.Layout(width="650px"),
        style={"description_width": "300px"}
    )
    
    event_name = widgets.Text(
        value=str(event_cfg.get("Name", "event_name")),
        description="Event name",
        layout=widgets.Layout(width="400px"),
        style={"description_width": "120px"}
    )
    
    hypo_lon = widgets.FloatText(
        value=float(event_cfg.get("Hypo_LonLat", [0.0, 0.0])[0]),
        description="Hypo lon",
        layout=widgets.Layout(width="300px"),
        style={"description_width": "120px"}
    )
    
    hypo_lat = widgets.FloatText(
        value=float(event_cfg.get("Hypo_LonLat", [0.0, 0.0])[1]),
        description="Hypo lat",
        layout=widgets.Layout(width="300px"),
        style={"description_width": "120px"}
    )
    
    event_mw = widgets.FloatText(
        value=float(event_cfg.get("Magnitude", 9.0)),
        description="Mw",
        layout=widgets.Layout(width="300px"),
        style={"description_width": "120px"}
    )
    
    event_box = widgets.VBox([
        widgets.HTML(
            "<b>Event parameters</b> "
            "<i>required only for Event-Baesd Probabilistic Tsunami Forecasting</i>"
        ),
        event_name,
        hypo_lon,
        hypo_lat,
        event_mw
    ])

    def update_application_widgets(*args):
        if application_widget.value == "PTF":
            event_box.layout.display = "flex"
        else:
            event_box.layout.display = "none"
    
    
    # -----------------------------
    # Ensemble configuration widgets
    # -----------------------------

    config_cfg = current_cfg.get("Configure", {})

    shape_widget = widgets.Dropdown(
        options=["Rectangle", "Circle"],
        value=config_cfg.get("shape", "Rectangle"),
        description="Shape",
        layout=widgets.Layout(width="320px"),
        style={"description_width": "150px"}
    )

    numb_stoch_widget = widgets.IntText(
        value=int(config_cfg.get("numb_stoch", 1)),
        description="N stochastic",
        layout=widgets.Layout(width="320px"),
        style={"description_width": "150px"}
    )

    variable_mu_widget = widgets.Checkbox(
        value=bool(config_cfg.get("variable_mu", 0)),
        description="Variable rigidity",
        indent=False
    )

    coupling_shallow_widget = widgets.FloatText(
        value=float(config_cfg.get("coupling_shallow_limit", 1.0)),
        description="Coupling shallow",
        layout=widgets.Layout(width="340px"),
        style={"description_width": "150px"}
    )

    coupling_deep_widget = widgets.FloatText(
        value=float(config_cfg.get("coupling_deep_limit", 55.0)),
        description="Coupling deep",
        layout=widgets.Layout(width="340px"),
        style={"description_width": "150px"}
    )

    minimum_bnd_distance_widget = widgets.FloatText(
        value=float(config_cfg.get("minimum_bnd_distance", 0.25)),
        description="Min bnd distance",
        layout=widgets.Layout(width="340px"),
        style={"description_width": "150px"}
    )

    minimum_interdistance_widget = widgets.FloatText(
        value=float(config_cfg.get("minimum_interdistance", 0.2)),
        description="Min interdistance",
        layout=widgets.Layout(width="340px"),
        style={"description_width": "150px"}
    )

    fact_area_scaling_widget = widgets.FloatText(
        value=float(config_cfg.get("Fact_area_scaling", 1.0)),
        description="Area factor",
        layout=widgets.Layout(width="340px"),
        style={"description_width": "150px"}
    )

    common_config_box = widgets.VBox([
        widgets.HTML("<b>Common ensemble configuration</b>"),
        shape_widget,
        numb_stoch_widget,
        coupling_shallow_widget,
        coupling_deep_widget,
        minimum_bnd_distance_widget,
        minimum_interdistance_widget,
        fact_area_scaling_widget
    ])

    magnitude_lb_widget = widgets.FloatText(
        value=float(config_cfg.get("Magnitude_lb", 0.3)),
        description="Magnitude lb",
        layout=widgets.Layout(width="340px"),
        style={"description_width": "150px"}
    )

    magnitude_ub_widget = widgets.FloatText(
        value=float(config_cfg.get("Magnitude_ub", 0.3)),
        description="Magnitude ub",
        layout=widgets.Layout(width="340px"),
        style={"description_width": "150px"}
    )

    hypo_baryc_distance_widget = widgets.FloatText(
        value=float(config_cfg.get("hypo_baryc_distance", 1.0)),
        description="Hypo-baryc dist",
        layout=widgets.Layout(width="340px"),
        style={"description_width": "150px"}
    )

    ptf_config_box = widgets.VBox([
        widgets.HTML("<b>PTF-only configuration</b>"),
        magnitude_lb_widget,
        magnitude_ub_widget,
        hypo_baryc_distance_widget
    ])

    sub_boundary_files = list_csv_files_for_widget("../config_files/Mesh")

    sub_boundary_logic = widgets.Checkbox(
        value=bool(config_cfg.get("mesh_sub_boundary", 0)),
        description="Use sub-boundary file",
        indent=False
    )

    sub_boundary_file = widgets.Dropdown(
        options=sub_boundary_files if sub_boundary_files else ["NO_SUB_BOUNDARY_CSV_FOUND"],
        value=(
            config_cfg.get("sub_boundary_file")
            if config_cfg.get("sub_boundary_file") in sub_boundary_files
            else (sub_boundary_files[0] if sub_boundary_files else "NO_SUB_BOUNDARY_CSV_FOUND")
        ),
        description="Boundary file",
        layout=widgets.Layout(width="750px"),
        style={"description_width": "150px"}
    )

    sub_boundary_box = widgets.VBox([
        widgets.HTML("<b>Sub-boundary</b>"),
        sub_boundary_logic,
        sub_boundary_file
    ])

    rigidity_files = list_csv_files_for_widget("../config_files/Rigidity")

    stress_drop_widget = widgets.Checkbox(
        value=bool(config_cfg.get("Stress_drop_var", 0)),
        description="Stress drop variation",
        indent=False
    )

    rigidity_file_logic = widgets.Checkbox(
        value=bool(config_cfg.get("Rigidity_file_logic", 0)),
        description="Use rigidity file",
        indent=False
    )

    rigidity_file = widgets.Dropdown(
        options=rigidity_files if rigidity_files else ["NO_RIGIDITY_CSV_FOUND"],
        value=(
            config_cfg.get("Rigidity_file")
            if config_cfg.get("Rigidity_file") in rigidity_files
            else (rigidity_files[0] if rigidity_files else "NO_RIGIDITY_CSV_FOUND")
        ),
        description="Rigidity file",
        layout=widgets.Layout(width="750px"),
        style={"description_width": "150px"}
    )

    fact_rigidity = widgets.FloatText(
        value=float(config_cfg.get("Fact_rigidity", 0.5)),
        description="Fact rigidity",
        layout=widgets.Layout(width="340px"),
        style={"description_width": "150px"}
    )

    rigidity_box = widgets.VBox([
        widgets.HTML("<b>Rigidity / stress-drop configuration</b>"),
        variable_mu_widget,
        stress_drop_widget,
        rigidity_file_logic,
        rigidity_file,
        fact_rigidity
    ])

    ensemble_box = widgets.VBox([
        common_config_box,
        rigidity_box,
        ptf_config_box,
        sub_boundary_box
    ])


    
    def update_ensemble_widgets(*args):

        # Event and PTF-only configuration are visible only for Event-based PTF
        if application_widget.value == "PTF":
            event_box.layout.display = "flex"
            ptf_config_box.layout.display = "flex"
        else:
            event_box.layout.display = "none"
            ptf_config_box.layout.display = "none"
    
        # Sub-boundary file appears only when requested
        if sub_boundary_logic.value:
            sub_boundary_file.layout.display = "flex"
        else:
            sub_boundary_file.layout.display = "none"
    
        # Variable-rigidity block is always visible,
        # but detailed options appear only if variable rigidity is enabled.
        rigidity_box.layout.display = "flex"
    
        if variable_mu_widget.value:
            stress_drop_widget.layout.display = "flex"
            rigidity_file_logic.layout.display = "flex"
    
            if rigidity_file_logic.value:
                rigidity_file.layout.display = "flex"
                fact_rigidity.layout.display = "none"
            else:
                rigidity_file.layout.display = "none"
                fact_rigidity.layout.display = "flex"
    
        else:
            stress_drop_widget.layout.display = "none"
            rigidity_file_logic.layout.display = "none"
            rigidity_file.layout.display = "none"
            fact_rigidity.layout.display = "none"

    
    application_widget.observe(update_ensemble_widgets, names="value")
    sub_boundary_logic.observe(update_ensemble_widgets, names="value")
    variable_mu_widget.observe(update_ensemble_widgets, names="value")
    rigidity_file_logic.observe(update_ensemble_widgets, names="value")

    update_ensemble_widgets()

    run_button = widgets.Button(
        description="Save input",
        button_style="success"
    )

    out = widgets.Output()

    run_button.on_click(run_from_widgets)

    display(
        widgets.HTML("<h2>Mesh input</h2>"),
        widgets.HTML(
            "Select how the slab mesh should be defined. "
            "The required fields will appear depending on the selected option."
        ),
        mesh_gen_widget,
        existing_mesh_box,
        geojson_mesh_box,
        rectangular_fault_box,
        mesh_common_box,
        
        widgets.HTML("<h2>Scaling law(s) and magnitude bins selection</h2>"),
        widgets.HTML("<b>Select one or more scaling laws:</b>"),
        law_box,
        widgets.HTML("<b>Magnitude bins</b>"),
        mode_widget,
        range_box,
        values_box,
        widgets.HTML("<h3>Application</h3>"),
        widgets.HTML(
            "Select the workflow. For "
            "<b>Event-based PTF</b> insert event parameters"
        ),
        application_widget,
        event_box,
        widgets.HTML("<h3>Ensemble configuration</h3>"),
        widgets.HTML(
            "Configure ensemble generation, coupling, PTF-only filters, "
            "sub-boundary options, and variable-rigidity settings."
        ),
        ensemble_box,
        run_button,
        out
    )

    return state

def list_csv_files_for_widget(folder):
    """
    List CSV files in a folder for dropdown widgets.
    Returns relative paths as strings.
    """

    if not os.path.isdir(folder):
        return []

    files = []

    for fname in sorted(os.listdir(folder)):
        if fname.lower().endswith(".csv"):
            files.append(os.path.join(folder, fname))

    return files

def reorder_input_config(cfg):
    """
    Reorder input.json keys for readability.

    Desired order:
      1. mesh-related fields
      2. Scaling
      3. application
      4. Event
      5. Configure
      6. any remaining keys
    """

    ordered_keys = [
        # Mesh / geometry
        "zone_name",
        "Merc_zone",
        "acronym",
        "mesh_gen",
        "rake",

        # Rectangular fault parameters, used when mesh_gen = 2
        "lon_c",
        "lat_c",
        "depth_km",
        "length_km",
        "width_km",
        "strike",
        "dip",
        "elem_size_km2",

        # Scaling and magnitude binning
        "Scaling",

        # Workflow
        "application",

        # Event, mainly for PTF
        "Event",

        # Remaining configuration
        "Configure",
    ]

    reordered = {}

    for key in ordered_keys:
        if key in cfg:
            reordered[key] = cfg[key]

    # Preserve any extra keys not listed above
    for key, value in cfg.items():
        if key not in reordered:
            reordered[key] = value

    return reordered

def create_barycenter_plot_widgets(Slab_obj):
    """
    Widget panel to plot rupture barycenters for selected Mw and scaling law.

    Requires:
        Slab_obj.active_barycenters()
        Slab_obj.select_barycenter2()
    """

    import ipywidgets as widgets
    from IPython.display import display, clear_output
    import matplotlib.pyplot as plt

    available_magnitudes = list(Slab_obj.get_magnitudes())
    available_scalings = list(Slab_obj.Name_scaling)

    mw_widget = widgets.Dropdown(
        options=[(f"{mw:.4f}", mw) for mw in available_magnitudes],
        value=available_magnitudes[0],
        description="Mw",
        layout=widgets.Layout(width="250px")
    )

    scaling_widget = widgets.Dropdown(
        options=available_scalings,
        value=available_scalings[0],
        description="Scaling",
        layout=widgets.Layout(width="350px")
    )

    plot_button = widgets.Button(
        description="Plot barycenters",
        button_style="success",
        layout=widgets.Layout(width="220px", height="40px")
    )

    out = widgets.Output()

    def on_plot_click(_):
        with out:
            clear_output(wait=True)

            Mw = mw_widget.value
            scaling = scaling_widget.value

            print(f"PLOT BARYCENTERS FOR Mw={Mw} and SCALING NAME: {scaling}")

            fig = plt.figure(figsize=(15, 15))
            subplots = 111

            fig, ax = Slab_obj.plot_basemap(fig, subplots)
            Slab_obj.plot_slab(ax, fig, colorbar=True)
            Slab_obj.plot_barycenters_mag(Mw, scaling, ax, fig)

            plt.show()

    plot_button.on_click(on_plot_click)

    display(
        widgets.HTML("<h3>Plot barycenters</h3>"),
        widgets.HTML("Select the magnitude and scaling law to plot the selected rupture barycenters."),
        mw_widget,
        scaling_widget,
        plot_button,
        out
    )

def create_rupture_area_plot_widgets(Slab_obj):
    """
    Widget panel to plot a rupture area for selected Mw, scaling law and area ID.

    Requires:
        Slab_obj.rupture_areas()
    """

    import ipywidgets as widgets
    from IPython.display import display, clear_output
    import matplotlib.pyplot as plt

    available_magnitudes = list(Slab_obj.get_magnitudes())
    available_scalings = list(Slab_obj.Name_scaling)

    mw_widget = widgets.Dropdown(
        options=[(f"{mw:.4f}", mw) for mw in available_magnitudes],
        value=available_magnitudes[0],
        description="Mw",
        layout=widgets.Layout(width="250px")
    )

    scaling_widget = widgets.Dropdown(
        options=available_scalings,
        value=available_scalings[0],
        description="Scaling",
        layout=widgets.Layout(width="350px")
    )

    area_id_widget = widgets.Dropdown(
        options=[0],
        value=0,
        description="Area ID",
        layout=widgets.Layout(width="250px")
    )

    colorbar_widget = widgets.Checkbox(
        value=True,
        description="Colorbar",
        indent=False
    )

    plot_button = widgets.Button(
        description="Plot rupture area",
        button_style="success",
        layout=widgets.Layout(width="220px", height="40px")
    )

    out = widgets.Output()

    def update_area_ids(*args):
        if not hasattr(Slab_obj, "Rupturing_areas"):
            area_id_widget.options = [0]
            area_id_widget.value = 0
            return

        Mw = mw_widget.value
        scaling = scaling_widget.value

        mag_idx = Slab_obj.get_magnitude_index(Mw)
        scaling_idx = Slab_obj.get_scaling_index(scaling)

        if mag_idx is None or scaling_idx is None:
            area_id_widget.options = [0]
            area_id_widget.value = 0
            return

        event = Slab_obj.Rupturing_areas[mag_idx][scaling_idx]
        n_areas = len(event)

        if n_areas == 0:
            area_id_widget.options = [0]
            area_id_widget.value = 0
        else:
            area_id_widget.options = list(range(n_areas))
            area_id_widget.value = 0

    def on_plot_click(_):
        with out:
            clear_output(wait=True)

            Mw = mw_widget.value
            scaling = scaling_widget.value
            area_id = area_id_widget.value

            print(f"PLOT RUPTURE AREA FOR Mw={Mw}, SCALING={scaling}, AREA ID={area_id}")

            fig = plt.figure(figsize=(15, 15))
            subplots = 111

            fig, ax = Slab_obj.plot_basemap(fig, subplots)

            Slab_obj.plot_rupture_area(
                Mw,
                scaling,
                area_id,
                ax,
                fig,
                colorbar=colorbar_widget.value
            )

            plt.show()

    mw_widget.observe(update_area_ids, names="value")
    scaling_widget.observe(update_area_ids, names="value")

    plot_button.on_click(on_plot_click)

    update_area_ids()

    display(
        widgets.HTML("<h3>Plot rupture area</h3>"),
        widgets.HTML("Select the magnitude, scaling law and rupture area ID."),
        mw_widget,
        scaling_widget,
        area_id_widget,
        colorbar_widget,
        plot_button,
        out
    )


def create_slip_distribution_plot_widgets(Slab_obj):
    """
    Widget panel to plot slip distribution for selected Mw, scaling law and area ID.

    Requires:
        Slab_obj.rupture_areas()
        Slab_obj.slip_distribution()
    """

    import ipywidgets as widgets
    from IPython.display import display, clear_output
    import matplotlib.pyplot as plt

    available_magnitudes = list(Slab_obj.get_magnitudes())
    available_scalings = list(Slab_obj.Name_scaling)

    mw_widget = widgets.Dropdown(
        options=[(f"{mw:.4f}", mw) for mw in available_magnitudes],
        value=available_magnitudes[0],
        description="Mw",
        layout=widgets.Layout(width="250px")
    )

    scaling_widget = widgets.Dropdown(
        options=available_scalings,
        value=available_scalings[0],
        description="Scaling",
        layout=widgets.Layout(width="350px")
    )

    area_id_widget = widgets.Dropdown(
        options=[0],
        value=0,
        description="Area ID",
        layout=widgets.Layout(width="250px")
    )

    var_widget = widgets.Checkbox(
        value=False,
        description="Variable rigidity",
        indent=False
    )

    colorbar_widget = widgets.Checkbox(
        value=True,
        description="Colorbar",
        indent=False
    )

    plot_button = widgets.Button(
        description="Plot slip distribution",
        button_style="success",
        layout=widgets.Layout(width="260px", height="40px")
    )

    out = widgets.Output()

    def update_area_ids(*args):
        if not hasattr(Slab_obj, "Rupturing_areas"):
            area_id_widget.options = [0]
            area_id_widget.value = 0
            return

        Mw = mw_widget.value
        scaling = scaling_widget.value

        mag_idx = Slab_obj.get_magnitude_index(Mw)
        scaling_idx = Slab_obj.get_scaling_index(scaling)

        if mag_idx is None or scaling_idx is None:
            area_id_widget.options = [0]
            area_id_widget.value = 0
            return

        event = Slab_obj.Rupturing_areas[mag_idx][scaling_idx]
        n_areas = len(event)

        if n_areas == 0:
            area_id_widget.options = [0]
            area_id_widget.value = 0
        else:
            area_id_widget.options = list(range(n_areas))
            area_id_widget.value = 0

    def on_plot_click(_):
        with out:
            clear_output(wait=True)

            Mw = mw_widget.value
            scaling = scaling_widget.value
            area_id = area_id_widget.value
            var = var_widget.value

            print(
                f"PLOT SLIP DISTRIBUTION FOR Mw={Mw}, "
                f"SCALING={scaling}, AREA ID={area_id}, VARIABLE MU={var}"
            )

            fig = plt.figure(figsize=(15, 15))
            subplots = 111

            fig, ax = Slab_obj.plot_basemap(fig, subplots)

            Slab_obj.plot_slip_dist(
                Mw,
                scaling,
                area_id,
                ax,
                fig,
                var=var,
                colorbar=colorbar_widget.value
            )

            plt.show()

    mw_widget.observe(update_area_ids, names="value")
    scaling_widget.observe(update_area_ids, names="value")

    plot_button.on_click(on_plot_click)

    update_area_ids()

    display(
        widgets.HTML("<h3>Plot slip distribution</h3>"),
        widgets.HTML("Select the magnitude, scaling law, rupture area ID and rigidity model."),
        mw_widget,
        scaling_widget,
        area_id_widget,
        var_widget,
        colorbar_widget,
        plot_button,
        out
    )

def infer_mercator_zone(lon):
    """
    Infer Mercator/UTM projection zone from longitude.
    Longitude is expected in degrees.
    """

    lon = float(lon)

    if lon > 180:
        lon -= 360

    return int((lon + 180) // 6) + 1

def infer_mercator_zone_from_existing_mesh(zone_name, base_dir="../utils/sz_slabs"):
    """
    Infer Mercator/UTM projection zone from an existing slab mesh.

    It reads:
        ../utils/sz_slabs/<zone_name>/subfaults/<zone_name>_mesh_nodes.dat

    Expected node file format:
        idx lon lat depth
    """

    nodes_file = os.path.join(
        base_dir,
        zone_name,
        "subfaults",
        f"{zone_name}_mesh_nodes.dat"
    )

    if not os.path.isfile(nodes_file):
        return None

    try:
        data = np.loadtxt(nodes_file)

        # Expected columns: idx lon lat depth
        lons = data[:, 1]

        # Convert longitudes from [0, 360] to [-180, 180] if needed
        lons = np.where(lons > 180, lons - 360, lons)

        lon_ref = np.nanmean(lons)

        return infer_mercator_zone(lon_ref)

    except Exception:
        return None

def infer_mercator_zone_from_geojson(zone_name, base_dir="../utils/sz_slabs"):
    """
    Infer Mercator/UTM projection zone from a GeoJSON mesh file.

    It reads:
        ../utils/sz_slabs/<zone_name>_mesh.geojson
    """

    geojson_file = os.path.join(
        base_dir,
        f"{zone_name}_mesh.geojson"
    )

    if not os.path.isfile(geojson_file):
        return None

    try:
        with open(geojson_file) as f:
            geojson_data = json.load(f)

        lons = []

        for feature in geojson_data["features"]:
            props = feature.get("properties", {})

            for key in ["lon1", "lon2", "lon3"]:
                if key in props:
                    lon = float(props[key])

                    if lon > 180:
                        lon -= 360

                    lons.append(lon)

        if len(lons) == 0:
            return None

        lon_ref = np.nanmean(lons)

        return infer_mercator_zone(lon_ref)

    except Exception:
        return None




def make_default_acronym(zone_name):
    """
    Create a 3-letter dummy acronym from zone_name.
    """
    clean = re.sub(r"[^A-Za-z0-9]", "", str(zone_name))

    if len(clean) >= 3:
        return clean[:3]

    if len(clean) > 0:
        return clean.ljust(3, "X")

    return "DUM"


def list_existing_slab_meshes(base_dir="../utils/sz_slabs"):
    """
    List slab names that already have subfault mesh nodes/faces.
    Used for mesh_gen = 0.
    """
    if not os.path.isdir(base_dir):
        return []

    names = []

    for name in sorted(os.listdir(base_dir)):
        slab_dir = os.path.join(base_dir, name)
        subfault_dir = os.path.join(slab_dir, "subfaults")

        nodes_file = os.path.join(subfault_dir, f"{name}_mesh_nodes.dat")
        faces_file = os.path.join(subfault_dir, f"{name}_mesh_faces.dat")

        if os.path.isfile(nodes_file) and os.path.isfile(faces_file):
            names.append(name)

    return names


def list_geojson_meshes(base_dir="../utils/sz_slabs"):
    """
    List GeoJSON mesh names from files named <zone_name>_mesh.geojson.
    Used for mesh_gen = 1.
    """
    if not os.path.isdir(base_dir):
        return []

    names = []

    for fname in sorted(os.listdir(base_dir)):
        if fname.endswith("_mesh.geojson"):
            names.append(fname.replace("_mesh.geojson", ""))

    return names


def update_input_json_mesh(input_file, mesh_block):
    """
    Update top-level mesh-related fields in input.json.
    """
    if os.path.isfile(input_file):
        with open(input_file) as f:
            cfg = json.load(f)
    else:
        cfg = get_default_input_config()
    
    cfg.update(mesh_block)

    cfg = reorder_input_config(cfg)

    with open(input_file, "w") as f:
        json.dump(cfg, f, indent=2)

    return cfg

def update_input_json_scaling(input_file, scaling_block):
    """
    Update Scaling block in input.json.
    """

    with open(input_file) as f:
        cfg = json.load(f)

    cfg["Scaling"] = scaling_block

    cfg = reorder_input_config(cfg)

    with open(input_file, "w") as f:
        json.dump(cfg, f, indent=2)

    return cfg

def update_input_json_application(input_file, application, event_block=None):
    """
    Update application and, when needed, Event block in input.json.

    application must be:
        - "PTF"
        - "Hazard"

    For backward compatibility, remove Configure["application"] if present.
    """

    with open(input_file) as f:
        cfg = json.load(f)

    cfg["application"] = application

    if "Configure" in cfg and "application" in cfg["Configure"]:
        del cfg["Configure"]["application"]

    if application == "PTF":
        if event_block is None:
            raise ValueError("Event block is required for PTF application")

        cfg["Event"] = event_block

    elif application == "Hazard":
        # Keep Event if already present, but it will not be used.
        # If you prefer to remove it, uncomment the next line:
        # cfg.pop("Event", None)
        pass

    else:
        raise ValueError(f"Unknown application: {application}")

    cfg = reorder_input_config(cfg)

    with open(input_file, "w") as f:
        json.dump(cfg, f, indent=2)

    return cfg

def update_input_json_configure(input_file, configure_block):
    """
    Update Configure block in input.json.
    """

    with open(input_file) as f:
        cfg = json.load(f)

        if "Configure" not in cfg:
            cfg["Configure"] = {}
    
        deprecated_keys = [
            "preprocess",
            "file_baryc",
            "file_baryc_name",
        ]
    
        for key in deprecated_keys:
            cfg["Configure"].pop(key, None)
    
        cfg["Configure"].update(configure_block)
    
    cfg = reorder_input_config(cfg)

    with open(input_file, "w") as f:
        json.dump(cfg, f, indent=2)

    return cfg

