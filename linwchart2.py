import streamlit as st
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from sympy import symbols, Eq, solve


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
num_boreholes = st.number_input("Enter the number of boreholes:", min_value=1,step=1)
Total_area = st.number_input("Enter the Area:", min_value=1, step=1)
boreholes_data = []

# Loop through each borehole for user inputs
for i in range(num_boreholes):
    st.header(f"Borehole {i+1}")
    borehole_name = st.text_input(f"Enter the name of Borehole {i+1}:")
    Lpile = st.number_input(
        f"Length of pile in soil for Borehole {i+1} (in m):",
        min_value=0.0,
        step=0.1)
    Dpile = st.number_input(f"Pile Diameter for Borehole {i+1} (in m):", min_value=0.0, step=0.01)
    num_layers = st.number_input(f"Enter the number of soil layers for Borehole {i+1}:", min_value=1, step=1)
    layers = []
    PD = 0  # Initialize PD for each borehole
    for j in range(num_layers):
        st.subheader(f"Layer {j+1}")
        Sdepth = st.number_input(f"Depth of Soil Layer {j+1} (in m):", min_value=0.0, step=0.1)
        gamma = st.number_input(f"Bulk Density of Soil Layer {j+1} (in kg/m³):", min_value=0.0, step=0.1) * 9.81  # Convert to kN/m³
        cohesion = st.number_input(f"Cohesion of Soil Layer {j+1} (in kPa):", min_value=0.0, step=0.1) * 98.1
        friction_angle = st.number_input(f"Angle of Friction for Soil Layer {j+1} (in degrees):",min_value=0, max_value=50, step=1)

        # Get Nγ, Nq, Nc values based on friction angle
        Nγ, Nq, Nc = interpolate_n_values(friction_angle)

        # Calculate overburden pressure (PD)
        PD += gamma * Lpile

        # Calculate Asi for friction component
        Asi = round(math.pi * Dpile * Lpile, 9)  # Surface area of pile

        layers.append({
            "Area": Total_area,
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
    Ap = (math.pi / 4) * Dpile**2  # Cross-sectional area of pile
    Aps_total = round(math.pi * Dpile * Lpile, 9)  # Total surface area of pile
    I = (math.pi * Dpile**4) / 64  # Moment of inertia

    # Factor of Safety input for each borehole with a unique key
    FOS = st.selectbox(f"Select Factor of Safety (FOS) for Borehole {i+1}:", options=[2.5, 3.0], key=f"fos_{i}")
    # Calculate end bearing component (Qub) in kN
    if Lpile <= layers[0]['Depth (m)']:
        d1 = Lpile
        PD = layers[0]['γ (kN/m³)'] * d1
        Qub_g = Ap * ((0.5 * Dpile * layers[0]['γ (kN/m³)'] * layers[0]['Nγ']) + (PD * layers[0]['Nq']))
        Qub_c = Ap * layers[0]['C (kPa)'] * layers[0]['Nc']
    elif Lpile <= layers[0]['Depth (m)'] + layers[1]['Depth (m)']:
        d2 = Lpile - layers[0]['Depth (m)']
        PD = (layers[0]['γ (kN/m³)'] *
              layers[0]['Depth (m)']) + (layers[1]['γ (kN/m³)'] * d2)
        Qub_g = Ap * (0.5 * Dpile * layers[1]['γ (kN/m³)'] * layers[1]['Nγ'] +PD * layers[1]['Nq'])
        Qub_c = Ap * layers[1]['C (kPa)'] * layers[1]['Nc']
    else:
        d3 = Lpile - (layers[0]['Depth (m)'] + layers[1]['Depth (m)'])
        PD = (layers[0]['γ (kN/m³)'] * layers[0]['Depth (m)']) + (layers[1]['γ (kN/m³)'] * layers[1]['Depth (m)']) + (layers[2]['γ (kN/m³)'] * d3)
        Qub_g = Ap * (0.5 * Dpile * layers[2]['γ (kN/m³)'] * layers[2]['Nγ'] + PD * layers[2]['Nq'])
        Qub_c = Ap * layers[2]['C (kPa)'] * layers[2]['Nc']

    # Calculate friction component (Quf) in kN
    Quf_g = 0
    Quf_c = 0
    K = 1  # Assuming K is a constant; you may want to make this an input
    PDf1 = 0
    PDf2 = 0
    PDf3 = 0
    if Lpile <= layers[0]['Depth (m)']:
        PDf1 = 0.5 * layers[0]['γ (kN/m³)'] * d1
        Asi1 = math.pi * Dpile * d1
        Quf_g = K * PDf1 * math.tan(math.radians(layers[0]['φ (degrees)'])) * Asi1
        Quf_c = cohesion * math.pi * Dpile * d1
    elif Lpile <= layers[0]['Depth (m)'] + layers[1]['Depth (m)']:
        d2 = Lpile - layers[0]['Depth (m)']
        PDf1 = 0.5 * layers[0]['γ (kN/m³)'] * layers[0]["Depth (m)"]
        PDf2 = layers[0]['γ (kN/m³)'] * layers[0]["Depth (m)"] + 0.5 * layers[1]['γ (kN/m³)'] * d2
        Asi1 = math.pi * Dpile * layers[0]['Depth (m)']
        Asi2 = math.pi * Dpile * d2
        Quf_g = (K * PDf1 * math.tan(math.radians(layers[0]['φ (degrees)'])) * Asi1) + (K * PDf2 * math.tan(math.radians(layers[1]['φ (degrees)'])) * Asi2)
        Quf_c = (layers[0]['C (kPa)'] * math.pi * Dpile * layers[0]['Depth (m)']) + (layers[1]['C (kPa)'] * math.pi * Dpile * d2)
    else:
        d3 = Lpile - (layers[0]['Depth (m)'] + layers[1]['Depth (m)'])
        PDf1 = 0.5 * layers[0]['γ (kN/m³)'] * layers[0]["Depth (m)"]
        PDf2 = layers[0]['γ (kN/m³)'] * layers[0]["Depth (m)"] + 0.5 * layers[1]['γ (kN/m³)'] * layers[1]["Depth (m)"]
        PDf3 = layers[0]['γ (kN/m³)'] * layers[0]["Depth (m)"] + layers[1]['γ (kN/m³)'] * layers[1]["Depth (m)"] + 0.5 * layers[2]['γ (kN/m³)'] * d3
        Asi1 = math.pi * Dpile * layers[0]['Depth (m)']
        Asi2 = math.pi * Dpile * layers[1]['Depth (m)']
        Asi3 = math.pi * Dpile * d3
        Quf_g = (K * PDf1 * math.tan(math.radians(layers[0]['φ (degrees)'])) * Asi1) + (K * PDf2 * math.tan(math.radians(layers[1]['φ (degrees)'])) * Asi2) + (K * PDf3 * math.tan(math.radians(layers[2]['φ (degrees)'])) * Asi3)
        Quf_c = (layers[0]['C (kPa)'] * math.pi * Dpile * layers[0]['Depth (m)']) + (layers[1]['C (kPa)'] * math.pi * Dpile * layers[1]['Depth (m)']) + (layers[2]['C (kPa)'] * math.pi * Dpile * d3)

    # Calculate total vertical load (Qu)
    Qug = Qub_g + Quf_g
    Quc = Qub_c + Quf_c
    Qu = Qug + Quc

    # Calculate allowable vertical load (Qall) in Ton
    Qall = round(Qug / FOS / 9.81, 4)  # Convert to Ton by dividing by 9.81 (gravity constant)
    # Compression Capacity (in tons)
    compression_capacity = Qu / 9.81 / FOS  # in tons

    # Tension Capacity (in tons)
    tension_capacity = (Quf_c + Quf_g) / 9.81 / FOS  # in tons

    # Tension including self-weight
    cross_section_area = Ap  # Cross-sectional area of pile in m²
    tension_sw = ((tension_capacity) + (Ap / 9.81) * (Lpile + 0.15) * 25)

    # Constants
    Dsteel = 78.5  # kN/m³
    Dconc = 25.0  # kN/m³
    Bcov = 0.05  # m
    e = 0.15  # m

    # Get the first two layers
    layer1 = layers[0]
    layer2 = layers[1]

    Sdepth1 = layer1['Depth (m)']
    Sdepth2 = layer2['Depth (m)']
    Y1 = layer1['γ (kN/m³)']
    Y2 = layer2['γ (kN/m³)']
    C1 = layer1['C (kPa)'] / 98.1
    C2 = layer2['C (kPa)'] / 98.1
    φ1 = layer1['φ (degrees)']
    φ2 = layer2['φ (degrees)']

    # Calculate unit weight of soil
    Y = ((Sdepth1 * Y1) + (Sdepth2 * Y2)) / (Sdepth1 + Sdepth2)

    # Calculate cohesion of soil
    Cu = C1 + C2

    # Calculate angle of internal friction
    Cphi = ((Sdepth1 * φ1) + (Sdepth2 * φ2)) / (Sdepth1 + Sdepth2)

    # Calculate passive pressure factor
    kp = (1 + math.sin(math.radians(Cphi))) / (1 - math.sin(math.radians(Cphi)))

    # Calculate area of pile
    P_area = 0.785 * (Dpile**2)

    L = Sdepth1 + Sdepth2
    L1 = L - Bcov

    #define H 
    H = symbols('H')
    # Define H 
    H = symbols('H')

    # Check if cohesion (Cu) is greater than 0
    if Cu > 0:
        # Calculate H1
        equation = Eq(H * (e + (1.5 * Dpile) + ((0.5 * H)/(9 * Cu * Dpile))), 
                        2.25 * Cu * Dpile * (L1 - (1.5 * Dpile) - (H/((9*Cu*Dpile)*(9*Cu*Dpile)))))
        
        # Solve equations for H
        solution = solve(equation, H)
        H1 = solution[0]
    else:
        # If Cu is 0, set H1 to 0 or some default value
        H1 = 0

    # Calculate H2
    H2 = (0.5 * Dpile * L1**3 * kp * Y) / (e + L1)

    # Calculate total H
    H = H1 + H2

    # Calculate lateral load carrying capacity of pile with FOS
    Hmax = H / FOS  # in kN
    lateral_capacity = Hmax / 9.81  # Convert to tons

    #constants
    configuration = "2PX54"
    module_wp = 600 #module power in watts
    dc_capacity = 64.8 # DC capacity in kW
    num_tables = 4977 #Number of Tables
    structure_type = "Single Pole Structure"
    num_piles_per_table = 23 # number of piles per table
    per_m3_cost = 7145   # per m3 cost of M25 concrete
    Capacity = Total_area / 4   #calculate capacity
    Total_Pile_no = num_tables * num_piles_per_table  # Total Number of pile
    
    Lpile1=Lpile + 0
    volume_1 = (3.14 * Dpile * Dpile) * Lpile1
    cost_per_pile_1 = volume_1 * per_m3_cost
    net_cost_1 = cost_per_pile_1 * Total_Pile_no

    Lpile2=Lpile + 0.05
    volume_2 = (3.14 * Dpile * Dpile) * (Lpile2)
    cost_per_pile_2 = volume_2 * per_m3_cost
    net_cost_2 = cost_per_pile_2 * Total_Pile_no
    
    cost_difference = net_cost_2 - net_cost_1
    
    boreholes_data.append({
    "Borehole Name": borehole_name,
    "Area": Total_area,
    "Lpile (m)": Lpile,
    "Dpile (m)": Dpile,
    "Ap (m²)": Ap,
    "Aps Total (m²)": Aps_total,
    "I (m⁴)": I,
    "Qub_g (kN)": Qub_g,
    "Qub_c (kN)": Qub_c,
    "Quf_g (kN)": Quf_g,
    "Quf_c (kN)": Quf_c,
    "Qu (kN)": Qu,
    "Qall (Ton)": Qall,
    "Compression Capacity (Ton)": compression_capacity,
    "Tension Capacity (Ton)": tension_capacity,
    "Tension with Self Weight (Ton)": tension_sw,
    "Unit Weight of Soil (kN/m3)": Y,
    "Cohesion of Soil (kN/m3)": Cu,
    "Angle of Friction (Cphi)": Cphi,
    "Passive Pressure Factor (kp)": kp,
    "Area of Pile (m2)": P_area,
    "Total length of pile below NGL (m)": L,
    "length of pile below NGL (m)": L1,
    "H1 (kN)": H1,
    "H2 (kN)": H2,
    "Total H (kN)": H,
    "Hmax (kN)": Hmax,
    "Lateral Load Capacity (Ton)": lateral_capacity,
    "Layers": layers,
    "Capacity": Capacity,
    "Total Pile No": Total_Pile_no,
    "LPile 1":Lpile1,
    "Volume1": volume_1,
    "Cost Per Pile (1)": cost_per_pile_1,
    "Net Cost 1": net_cost_1,
    "LPile 2":Lpile2,
    "Volume2": volume_2,
    "Cost Per Pile (2)": cost_per_pile_2,
    "Net Cost 2": net_cost_2,
    "Cost Difference Per 0.05": cost_difference
})
    
# Display Results
st.header("Results")

# Create DataFrame for all boreholes
results = []
for borehole in boreholes_data:
    for layer in borehole["Layers"]:
        results.append({
            "Borehole Name":
            borehole["Borehole Name"],
            "Area": 
            borehole["Area"],
            "Length of Pile (m)":
            borehole["Lpile (m)"],
            "Pile Diameter (m)":
            borehole["Dpile (m)"],
            "Layer Name":
            layer["Layer Name"],
            "Depth (m)":
            layer["Depth (m)"],
            "Aps Total (m²)":
            borehole["Aps Total (m²)"],
            "Ap (m²)":
            borehole["Ap (m²)"],
            "I (m⁴)":
            borehole["I (m⁴)"],
            "γ (kN/m³)":
            layer["γ (kN/m³)"],
            "C (kPa)":
            layer["C (kPa)"],
            "φ (degrees)":
            layer["φ (degrees)"],
            "Nγ":
            layer["Nγ"],
            "Nq":
            layer["Nq"],
            "Nc":
            layer["Nc"],
            "PD (kPa)":
            layer["PD (kPa)"],
            "Asi (m²)":
            layer["Asi (m²)"],
            "Qub_g (kN)":
            borehole["Qub_g (kN)"],
            "Qub_c (kN)":
            borehole["Qub_c (kN)"],
            "Quf_g (kN)":
            borehole["Quf_g (kN)"],
            "Quf_c (kN)":
            borehole["Quf_c (kN)"],
            "Qu (kN)":
            borehole["Qu (kN)"],  # Corrected from "Qub" to "Qu"
            "Qall (Ton)":
            borehole["Qall (Ton)"],
            "PDf1":
            PDf1,
            "PDf2":
            PDf2,
            "PDf3":
            PDf3,
            "Compression Capacity (Ton)":
            borehole["Compression Capacity (Ton)"],
            "Tension Capacity (Ton)":
            borehole["Tension Capacity (Ton)"],
            "Tension with Self Weight (Ton)":
            borehole["Tension with Self Weight (Ton)"],
            "Unit Weight of Soil (kN/m3)":
            borehole["Unit Weight of Soil (kN/m3)"],
            "Cohesion of Soil (kN/m3)":
            borehole["Cohesion of Soil (kN/m3)"],
            "Angle of Friction (Cphi)":
            borehole["Angle of Friction (Cphi)"],
            "Passive Pressure Factor (kp)":
            borehole["Passive Pressure Factor (kp)"],
            "Area of Pile (m2)":
            borehole["Area of Pile (m2)"],
            "Total length of pile below NGL (m)":
            borehole["Total length of pile below NGL (m)"],
            "length of pile below NGL (m)":
            borehole["length of pile below NGL (m)"],
            "H1 (kN)":
            borehole["H1 (kN)"],
            "H2 (kN)":
            borehole["H2 (kN)"],
            "Total H (kN)":
            borehole["Total H (kN)"],
            "Hmax (kN)":
            borehole["Hmax (kN)"],
            "Lateral Load Capacity (Ton)":
            borehole["Lateral Load Capacity (Ton)"],
            "Capacity": 
            borehole['Capacity'],
            "Total No. of pile": 
            borehole['Total Pile No'],  
            "LPile-1":
            borehole['LPile 1'],
            "Volume1": 
            borehole['Volume1'],
            "Cost Per Pile (1)": 
            borehole['Cost Per Pile (1)'],
            "Net Cost 1": 
            borehole['Net Cost 1'],
            "LPile-2":
            borehole['LPile 2'],
            "Volume2": 
            borehole['Volume2'],
            "Cost Per Pile (2)": 
            borehole['Cost Per Pile (2)'],
            "Net Cost 2": 
            borehole['Net Cost 2'],
            "Cost Difference Per 0.05": 
            borehole['Cost Difference Per 0.05']
        })

df = pd.DataFrame(results)
st.write(df)

# Create a selection box for the user to choose a borehole
selected_borehole = st.selectbox("Select a Borehole to View Load Carrying Capacity Chart:", df["Borehole Name"].unique())

# Filter the DataFrame for the selected borehole
borehole_data = df[df["Borehole Name"] == selected_borehole]

# User-defined pile depth
pile_depth_start = Lpile  # Example starting pile depth
pile_depth_end = pile_depth_start + 1.0  # Range till next 1m
pile_depth = np.arange(pile_depth_start, pile_depth_end + 0.05, 0.05)  # Interval of 0.05m

# Assuming compression_03, lateral_03, and tension_03 are based on the same pile depth range
compression_capacity = borehole_data["Compression Capacity (Ton)"].values[0]  # Get the first value
lateral_capacity = borehole_data["Lateral Load Capacity (Ton)"].values[0]  # Get the first value
tension_capacity = borehole_data["Tension Capacity (Ton)"].values[0]  # Get the first value

# Create arrays for plotting
compression_03 = np.full_like(pile_depth, compression_capacity)  # Fill with the same compression capacity
lateral_03 = np.full_like(pile_depth, lateral_capacity)  # Fill with the same lateral capacity
tension_03 = np.full_like(pile_depth, tension_capacity)  # Fill with the same tension capacity

# Plot
plt.figure(figsize=(8, 5))
plt.plot(pile_depth, compression_03, marker='o', linestyle='-', color='deepskyblue', label='COMPRESSION (0.3DIA)')
plt.plot(pile_depth, lateral_03, marker='o', linestyle='-', color='purple', label='LATERAL (0.3DIA)')
plt.plot(pile_depth, tension_03, marker='o', linestyle='-', color='green', label='TENSION (0.3DIA)')

# Labels and title
plt.xlabel("Pile Depth (m)")
plt.ylabel("Safe Pile Load Carrying Capacity (T)")
plt.title(f"Borehole: {selected_borehole}")
plt.legend()
plt.grid()

# Show plot in Streamlit
st.pyplot(plt)

# Provide an option to download results as a CSV file
csv = df.to_csv(index=False).encode("utf-8")
st.download_button(
    label="Download Results as CSV",
    data=csv,
    file_name="pile_capacity_results.csv",
    mime="text/csv",
)
