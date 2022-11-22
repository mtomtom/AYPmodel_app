import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import plotly.figure_factory as ff
import plotly.express as px
from plotly.subplots import make_subplots
from p2_functions.simulation import Simulation
from p2_functions.leafhoppers_class import Leafhoppers
from p2_functions.plants_class import Plants

st.set_page_config(layout="wide")
st.title('Leafhopper pathogen model')

## Add in the parameter options - organise into sections
with st.expander("Main parameters"):
    cols = st.columns(4)
    timestep = float(cols[0].text_input("Timestep (days)", value = 0.1))
    sim_length = int(cols[0].text_input("Simulation length (days)", 365))
    world_edge_x = int(cols[0].text_input("Field x (m)", value = 400))
    world_edge_y = int(cols[0].text_input("Field y (m)", value = 300))
    infected_plant_initial_conc1 = float(cols[1].text_input("Infected plant", value = 0.0))
   # infected_plant_initial_conc2 = float(cols[1].text_input("Infected plant 2", value = 0.0))
    mort_PI = float(cols[1].text_input("Virulence",value = 0.02))
    kGrow = float(cols[1].text_input("kGrow", value = 9))
    grow_PU = float(cols[1].text_input("grow_PU", value = 0.05))
    infect = float(cols[2].text_input("Transmission", value = 0.1))
    leafhopper_initial_conc = float(cols[2].text_input("Lu initial", value = 0.01))
    infected_leafhopper_initial_conc1 = float(cols[2].text_input("Li initial",value = 0.001))
    #infected_leafhopper_initial_conc2 = float(cols[2].text_input("Li2 initial",value = 0.0))
    mort = float(cols[2].text_input("Mort",value = 0.025))
    ddmort = float(cols[2].text_input("DDmort", value = 0.9))
    rep = float(cols[3].text_input("rep", value = 0.1))
    step = float(cols[3].text_input("Step", value = 3.0))
    transmission = cols[3].radio("Infection transmission",("density","frequency"))
    shuffle_patches = cols[3].checkbox("Shuffle patches", value = False)
    
with st.expander("Additional parameters"):
    resolution = int(st.text_input("Resolution", value = 2))
    lh_direct = int(st.text_input("lh_direct", value = 10))
    kl = float(st.text_input("kl", value = 0.015))
    n = float(st.text_input("n", value = 2.0))
    sigma = float(st.text_input("sigma", value = 10.0))
    

sim_steps = int(sim_length / timestep)

## Run the model
cols = st.columns(2)

## There is no option for two pathogen strains
infected_plant_initial_conc2 = 0.0
infected_leafhopper_initial_conc2 = 0.0

if st.checkbox("Run simulation"):
    sim = Simulation()	
    sim.setUpPlants(world_edge_x,world_edge_y,infect,resolution,infected_plant_initial_conc1,infected_plant_initial_conc2,1,1,1,kGrow, mort_PI, grow_PU, timestep)
    sim.setUpLeafhoppers(world_edge_x,world_edge_y,step,infect, leafhopper_initial_conc, infected_leafhopper_initial_conc1,infected_leafhopper_initial_conc2,resolution,
    1,1,1, rep, mort, ddmort, lh_direct, kl, n, timestep)

    ## Create the figure output
    time_update = st.empty()
    #cols = st.columns([1.5, 1])
    plot_spot1 = st.empty()
    plot_spot2 = st.empty()

    def show_snapshot(Hi1, Hi2, Hu, Li1, Li2, Lu, df):
        fig = px.imshow(Hi1.T,color_continuous_scale='Aggrnyl',width=650, height=450, origin="lower")
        fig.add_scatter(x=df["Li1_x"] * sim.leaf.lhop_to_patch_x, y=df["Li1_y"] * sim.leaf.lhop_to_patch_y,mode='markers', line=dict(
                color='Red'), name="Li")
        fig.add_scatter(x=df["Lu_x"] * sim.leaf.lhop_to_patch_x, y=df["Lu_y"] * sim.leaf.lhop_to_patch_y,mode='markers',line=dict(
                color='Blue'), name="Lu")
        # set colorbar position for X-axis
        #fig.update_layout(coloraxis_colorbar_x=-1.5)
        
        # set colorbar position for Y-axis
        #fig.update_layout(coloraxis_colorbar_y=0.5)
        fig.update_layout(legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=-0.2
    ))
        # Plot
        st.plotly_chart(fig, use_container_width=False)
    
    def plot_data(df_main):
        subfig = make_subplots(specs=[[{"secondary_y": True}]])

        # create two independent figures with px.line each containing data from multiple columns
        fig = px.line(df_main, x="time_point", y=["Li1","Lu"], width=800, height=600)
        fig2 = px.line(df_main, x="time_point", y=["Hi1","Hu"], width=800, height=600)

        fig2.update_traces(yaxis="y2")

        subfig.add_traces(fig.data + fig2.data)
        subfig.layout.xaxis.title="Days"
        subfig.layout.yaxis.title="Leafhopper density"
        subfig.layout.yaxis2.title="Plant density"
        # recoloring is necessary otherwise lines from fig und fig2 would share each color
        # e.g. Linear-, Log- = blue; Linear+, Log+ = red... we don't want this
        subfig.for_each_trace(lambda t: t.update(line=dict(color=t.marker.color)))
        st.plotly_chart(subfig)

    ## Create a dataframe to hold the main output
    df_columns=["Li1","Li2","Lu","Hi1","Hi2","Hu","time_point"]
    df_main = pd.DataFrame(0.0, index=np.arange(sim_steps), columns=df_columns)
    
    ## Run the model
    for i in range(sim_steps):
        sim.moveLeafhoppers()
        sim.updateLeafhoppersPlants()
        if transmission == 'density':
            sim.infectLeafhoppers_density()
        if transmission == 'frequency':
            sim.infectLeafhoppers_freq()
        sim.reproduction()
        sim.mortLeafhoppers()
        sim.updatePlants()
        sim.addNewLeafs()
        if shuffle_patches: 
            sim.shufflePlants()
            sim.updateLeafhoppersPlants()
        ## Dynamic plot of state of infection
        Hi1 = sim.plants.matrix[:,:,sim.plants.p_dict["Hi1"]]
        Hi2 = sim.plants.matrix[:,:,sim.plants.p_dict["Hi2"]]
        Hu = sim.plants.matrix[:,:,sim.plants.p_dict["Hu"]]
        Li1 = sim.plants.matrix[:,:,sim.plants.p_dict["Li1"]]
        Li2 = sim.plants.matrix[:,:,sim.plants.p_dict["Li2"]]
        Lu = sim.plants.matrix[:,:,sim.plants.p_dict["Lu"]]

        ## Make leafhopper datapoints
        df = pd.DataFrame([sim.leaf.inf_xcoord1,sim.leaf.inf_ycoord1, sim.leaf.inf_xcoord2,sim.leaf.inf_ycoord2,sim.leaf.uninf_xcoord, sim.leaf.uninf_ycoord])
        df = df.T
        df.columns=["Li1_x", "Li1_y","Li2_x","Li2_y","Lu_x","Lu_y"]

        ## Add to output
        df_main.loc[i]["Li1"] = np.mean(Li1)
        df_main.loc[i]["Li2"] = np.mean(Li2)
        df_main.loc[i]["Lu"] = np.mean(Lu)
        df_main.loc[i]["Hi1"] = np.mean(Hi1)
        df_main.loc[i]["Hi2"] = np.mean(Hi2)
        df_main.loc[i]["Hu"] = np.mean(Hu)
        df_main.loc[i]["time_point"] = i * timestep

        if i%10 == 0:
            with plot_spot1:
                show_snapshot(Hi1, Hi2, Hu, Li1, Li2, Lu, df)
            with plot_spot2:
                plot_data(df_main.loc[0:i])
            time_text = str(round(i * timestep,3)) + " days"
            time_update.text(time_text)
            #time.sleep(0.1)
    @st.cache
    def convert_df(df):
        # IMPORTANT: Cache the conversion to prevent computation on every rerun
        return df.to_csv()

    csv = convert_df(df_main)

    st.download_button(
        label="Download data as CSV",
        data=csv,
        file_name='results.csv',
        mime='text/csv',
    )
