import streamlit as st
import pandas as pd
import numpy as np
import math

# Function to calculate Nγ, Nq, and Nc based on the friction angle using interpolation
def interpolate_n_values(friction_angle):
    table = {
        0: (5.14, 1.00, 0.00),
        5: (6.49, 1.57, 0.45),
        10: (8.35, 2.47, 1.22),
        15: (10.98, 3.94, 2.65),
        20: (14.83, 6.40, 5.39),
        25: (20.72, 10.66, 10.88),
        30: (30.14, 18.40, 22.40),
        35: (46.12, 33.30, 48.03),
        40: (75.31, 64.20, 109.41),
        45: (138.88, 134.88, 271.76),
        50: (266.89, 319.07, 762.89),
    }
    angles = sorted(table.keys())
    for i in range(len(angles) - 1):
        if angles[i] <= friction_angle <= angles[i + 1]:
            angle1, angle2 = angles[i], angles[i + 1]
            n1, nq1, nc1 = table[angle1]
            n2, nq2, nc2 = table[angle2]
            # Linear interpolation
            nc = n1 + (n2 - n1) * (friction_angle - angle1) / (angle2 - angle1)
            nq = nq1 + (nq2 - nq1) * (friction_angle - angle1) / (angle2 - angle1)
            nγ = nc1 + (nc2 - nc1) * (friction_angle - angle1) / (angle2 - angle1)
            return round(nγ, 4), round(nq, 4), round(nc, 4)
    return table[angles[-1]]

# Streamlit App
st.title("Pile Capacity Calculator")

# Input: Number of Boreholes
num_boreholes = st.number_input("Enter the number of boreholes:", min_value=1, step=1)

boreholes_data = []

# Loop through each borehole for user inputs
for i in range(num_boreholes):
    st.header(f"Borehole {i+1}")
    borehole_name = st.text_input(f"Enter the name of Borehole {i+1}:")
    Lpile = st.number_input(f"Length of pile in soil for Borehole {i+1} (in m):", min_value=0.0, step=0.1)
    Dpile = st.number_input(f"Pile Diameter for Borehole {i+1} (in m):", min_value=0.0, step=0.01)
    num_layers = st.number_input(f"Enter the number of soil layers for Borehole {i+1}:", min_value=1, step=1)
    
    layers = []
    PD = 0  # Initialize PD for each borehole
    for j in range(num_layers):
        st.subheader(f"Layer {j+1}")
        Sdepth = st.number_input(f"Depth of Soil Layer {j+1} (in m):", min_value=0.0, step=0.1)
        gamma = st.number_input(f"Bulk Density of Soil Layer {j+1} (in kg/m³):", min_value=0.0, step=0.1) * 9.81  # Convert to kN/m³
        cohesion = st.number_input(f"Cohesion of Soil Layer {j+1} (in kPa):", min_value=0.0, step=0.1)
        friction_angle = st.number_input(f"Angle of Friction for Soil Layer {j+1} (in degrees):", min_value=0, max_value=50, step=1)

        # Get Nγ, Nq, Nc values based on friction angle
        Nγ, Nq, Nc = interpolate_n_values(friction_angle)

        # Calculate overburden pressure (PD)
        PD += gamma * Sdepth

        # Calculate Asi for friction component
        Asi = round(math.pi * Dpile * Sdepth, 9)  # Surface area of pile

        layers.append({
            "Layer Name": f"Layer {j+1}",
            "Depth (m)": Sdepth,
            "γ (kN/m³)": gamma,
            "C (kPa)": cohesion,
            "φ (degrees)": friction_angle,
            "Nγ": Nγ,
            "Nq": Nq,
            "Nc": Nc,
            "PD (kPa)": PD,
            "Asi (m²)": Asi,
        })

    # Calculate Borehole-wide parameters
    Ap = (math.pi / 4) * Dpile ** 2  # Cross-sectional area of pile
    Aps_total = round(math.pi * Dpile * Lpile, 9)  # Total surface area of pile
    I = (math.pi * Dpile ** 4) / 64  # Moment of inertia

    # Factor of Safety input for each borehole with a unique key
    FOS = st.selectbox(f"Select Factor of Safety (FOS) for Borehole {i+1}:", options=[2.5, 3.0], key=f"fos_{i}")

    # Calculate end bearing component (Qub) in kN
    if Lpile <= layers[0]['Depth (m)']:
        d1 = Sdepth
        PD = layers[0]['γ (kN/m³)'] * d1
        Qub = Ap * ((0.5 * Dpile * layers[0]['γ (kN/m³)'] * layers[0]['Nγ']) + (PD * layers[0]['Nq']))
    elif Lpile <= layers[0]['Depth (m)'] + layers[1]['Depth (m)']:
        d2 = Lpile - layers[0]['Depth (m)']
        PD = (layers[0]['γ (kN/m³)'] * layers[0]['Depth (m)']) + (layers[1]['γ (kN/m³)'] * d2)
        Qub = Ap * (0.5 * Dpile * layers[1]['γ (kN/m³)'] * layers[1]['Nγ'] + PD * layers[1]['Nq'])
    else:
        d3 = Lpile - (layers[0]['Depth (m)'] + layers[1]['Depth (m)'])
        PD = (layers[0]['γ (kN/m³)'] * layers[0]['Depth (m)']) + (layers[1]['γ (kN/m³)'] * layers[1]['Depth (m)']) + (layers[2]['γ (kN/m³)'] * d3)
        Qub = Ap * (0.5 * Dpile * layers[2]['γ (kN/m³)'] * layers[2]['Nγ'] + PD * layers[2]['Nq'])

    # Calculate friction component (Quf) in kN
    Quf = 0
    K = 1  # Assuming K is a constant; you may want to make this an input
    PDf1 = 0
    PDf2 = 0
    PDf3 = 0
    if Lpile <= layers[0]['Depth (m)']:
        PDf1 = 0.5 * layers[0]['γ (kN/m³)'] * d1
        Asi1 = math.pi * Dpile * d1
        Quf = K * PDf1 * math.tan(math.radians(layers[0]['φ (degrees)'])) * Asi1
    elif Lpile <= layers[0]['Depth (m)'] + layers[1]['Depth (m)']:
        d2 = Lpile - layers[0]['Depth (m)']
        PDf1 = 0.5 * layers[0]['γ (kN/m³)'] * layers[0]["Depth (m)"]
        PDf2 = layers[0]['γ (kN/m³)'] * layers[0]["Depth (m)"] + 0.5 * layers[1]['γ (kN/m³)'] * d2
        Asi1 = math.pi * Dpile * layers[0]['Depth (m)']
        Asi2 = math.pi * Dpile * d2
        Quf = (K * PDf1 * math.tan(math.radians(layers[0]['φ (degrees)'])) * Asi1) + (K * PDf2 * math.tan(math.radians(layers[1]['φ (degrees)'])) * Asi2)
    else:
        d3 = Lpile - (layers[0]['Depth (m)'] + layers[1]['Depth (m)'])
        PDf1 = 0.5 * layers[0]['γ (kN/m³)'] * layers[0]["Depth (m)"]
        PDf2 = layers[0]['γ (kN/m³)'] * layers[0]["Depth (m)"] + 0.5 * layers[1]['γ (kN/m³)'] * layers[1]["Depth (m)"]
        PDf3 = layers[0]['γ (kN/m³)'] * layers[0]["Depth (m)"] + layers[1]['γ (kN/m³)'] * layers[1]["Depth (m)"] + 0.5 * layers[2]['γ (kN/m³)'] * d3
        Asi1 = math.pi * Dpile * layers[0]['Depth (m)']
        Asi2 = math.pi * Dpile * layers[1]['Depth (m)']
        Asi3 = math.pi * Dpile * d3
        Quf = (K * PDf1 * math.tan(math.radians(layers[0]['φ (degrees)'])) * Asi1) + (K * PDf2 * math.tan(math.radians(layers[1]['φ (degrees)'])) * Asi2) + (K * PDf3 * math.tan(math.radians(layers[2]['φ (degrees)'])) * Asi3)

    # Calculate total vertical load (Qu)
    Qu = Qub + Quf

    # Calculate allowable vertical load (Qall) in Ton
    Qall = round(Qu / FOS / 9.81, 4)  # Convert to Ton by dividing by 9.81 (gravity constant)

    # Compression Capacity (in tons)
    compression_capacity = Qu / 9.81 / FOS  # in tons

    # Tension Capacity (in tons)
    tension_capacity = Qu / 9.81 / FOS  # in tons

    # Tension including self-weight
    cross_section_area = Ap  # Cross-sectional area of pile in m²
    tension_sw = ((tension_capacity) + ((Ap * Lpile * 25) / 9.81))

    boreholes_data.append({
        "Borehole Name": borehole_name,
        "Lpile (m)": Lpile,
        "Dpile (m)": Dpile,
        "Ap (m²)": Ap,
        "Aps Total (m²)": Aps_total,
        "I (m⁴)": I,
        "Qub (kN)": Qub,
        "Quf (kN)": Quf,
        "Qu (kN)": Qu,
        "Qall (Ton)": Qall,
        "Compression Capacity (Ton)": compression_capacity,
        "Tension Capacity (Ton)": tension_capacity,
        "Tension with Self Weight (Ton)": tension_sw,
        "Layers": layers
    })

# Display Results
st.header("Results")

# Create DataFrame for all boreholes
results = []
for borehole in boreholes_data:
    for layer in borehole["Layers"]:
        results.append({
            "Borehole Name": borehole["Borehole Name"],
            "Length of Pile (m)": borehole["Lpile (m)"],
            "Pile Diameter (m)": borehole["Dpile (m)"],
            "Layer Name": layer["Layer Name"],
            "Depth (m)": layer["Depth (m)"],
            "Aps Total (m²)": borehole["Aps Total (m²)"],
            "Ap (m²)": borehole["Ap (m²)"],
            "I (m⁴)": borehole["I (m⁴)"],
            "γ (kN/m³)": layer["γ (kN/m³)"],
            "C (kPa)": layer["C (kPa)"],
            "φ (degrees)": layer["φ (degrees)"],
            "Nγ": layer["Nγ"],
            "Nq": layer["Nq"],
            "Nc": layer["Nc"],
            "PD (kPa)": layer["PD (kPa)"],
            "Asi (m²)": layer["Asi (m²)"],
            "Qub (kN)": borehole["Qub (kN)"],
            "Quf (kN)": borehole["Quf (kN)"],
            "Qu (kN)": borehole["Qu (kN)"],
            "Qall (Ton)": borehole["Qall (Ton)"],
            "Compression Capacity (Ton)": borehole["Compression Capacity (Ton)"],
            "Tension Capacity (Ton)": borehole["Tension Capacity (Ton)"],
            "Tension with Self Weight (Ton)": borehole["Tension with Self Weight (Ton)"]
        })

df = pd.DataFrame(results)
st.write(df)

# Provide an option to download results as a CSV file
csv = df.to_csv(index=False).encode("utf-8")
st.download_button(
    label="Download Results as CSV",
    data=csv,
    file_name="pile_capacity_results.csv",
    mime="text/csv",
)