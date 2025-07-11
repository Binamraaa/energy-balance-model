import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# --- Constants ---
SIGMA = 5.67e-8  # W/m² K^4
WATER_DENSITY = 1000  # kg/m³
WATER_SPECIFIC_HEAT = 4184  # J/(kg·K)
SECS_PER_YEAR = 365 * 24 * 3600
EARTH_RADIUS = 6.371e6 # meters (for distance calculation in diffusion)

# --- Temperature Conversion & Bounds ---
K_TO_C_OFFSET = 273.15

def k_to_c(temp_k):
    """Converts Kelvin to Celsius"""
    return temp_k - K_TO_C_OFFSET

def c_to_k(temp_c):
    """Converts Celsius to Kelvin"""
    return temp_c + K_TO_C_OFFSET

# Define plausible temperature bounds (internal use in K)
TEMP_MIN_PLAUSIBLE_K = 150.0 # K
TEMP_MAX_PLAUSIBLE_K = 450.0 # K
# Corresponding bounds in C for plot scaling
TEMP_MIN_PLAUSIBLE_C = k_to_c(TEMP_MIN_PLAUSIBLE_K)
TEMP_MAX_PLAUSIBLE_C = k_to_c(TEMP_MAX_PLAUSIBLE_K)


# --- Helper Functions (No changes needed below) ---
def get_latitude_bands(num_bands):
    """Generates center latitudes (deg, rad), edges (deg), and band width (rad)."""
    band_width_deg = 180.0 / num_bands
    edges_deg = np.linspace(-90.0, 90.0, num_bands + 1)
    centers_deg = edges_deg[:-1] + band_width_deg / 2.0
    centers_rad = np.deg2rad(centers_deg)
    band_width_rad = np.deg2rad(band_width_deg)
    delta_y = EARTH_RADIUS * band_width_rad # meters
    return centers_deg, centers_rad, edges_deg, delta_y

def annual_avg_insolation(latitude_rad, L=1361.0):
    """ Simple polynomial fit for annual average insolation distribution. """
    x = np.sin(latitude_rad)
    p2 = 0.5 * (3 * x**2 - 1)
    s2 = -0.48 # Empirical coefficient
    insolation_factor = np.maximum(0, 1.0 + s2 * p2)
    return (L / 4.0) * insolation_factor

def latitude_albedo(latitude_deg, base_albedo, ice_latitude, polar_albedo):
    """ Albedo increases sharply poleward of the ice_latitude. """
    albedo = np.full_like(latitude_deg, base_albedo)
    ice_mask = np.abs(latitude_deg) >= ice_latitude
    albedo[ice_mask] = polar_albedo
    return albedo

def latitude_heat_capacity(latitude_deg, ocean_frac_tropics, ocean_frac_poles, depth_ocean, depth_land):
    """ Models heat capacity based on a simple ocean fraction profile. """
    abs_lat_norm = np.abs(latitude_deg) / 90.0
    ocean_fraction = ocean_frac_tropics * (1 - abs_lat_norm) + ocean_frac_poles * abs_lat_norm
    ocean_fraction = np.clip(ocean_fraction, 0, 1)
    depth_eff = depth_ocean * ocean_fraction + depth_land * (1 - ocean_fraction)
    heat_capacity_local = depth_eff * WATER_DENSITY * WATER_SPECIFIC_HEAT
    return heat_capacity_local

def latitude_clear_sky_epsilon(latitude_deg, epsilon_tropics, epsilon_poles):
    """ Models lower epsilon (more trapping) in warmer/wetter tropics. """
    abs_lat_norm = np.abs(latitude_deg) / 90.0
    epsilon_local = epsilon_tropics * (1 - abs_lat_norm) + epsilon_poles * abs_lat_norm
    return np.clip(epsilon_local, 0.01, 1.0)

def latitude_aerosol_forcing(latitude_deg, base_forcing, nh_midlat_peak_forcing, peak_latitude, peak_width):
    """ Models aerosol forcing with a peak, e.g., in NH mid-latitudes. """
    forcing = base_forcing + (nh_midlat_peak_forcing - base_forcing) * \
              np.exp(-((latitude_deg - peak_latitude)**2) / (2 * peak_width**2))
    return forcing

def get_latitude_zone(latitude_deg):
    """Returns a descriptive string for the latitude zone."""
    abs_lat = abs(latitude_deg)
    if abs_lat <= 10:
        return "Equatorial"
    elif abs_lat <= 30:
        return "Tropics"
    elif abs_lat <= 60:
        return "Mid-latitudes"
    else:
        return "Polar"

# --- Simulation Function (Input/Output remains KELVIN internally) ---
def run_latitude_simulation_enhanced(
    L, initial_temp_profile_k, # *** Expects initial temp in Kelvin ***
    dt_years, n_steps, num_latitude_bands,
    ocean_frac_tropics, ocean_frac_poles, depth_ocean, depth_land,
    heat_diffusivity_coeff,
    base_albedo, ice_latitude, polar_albedo,
    epsilon_tropics, epsilon_poles,
    base_aerosol_forcing, nh_midlat_peak_forcing, peak_latitude, peak_width,
    additional_ghg_forcing
):

    lat_centers_deg, lat_centers_rad, lat_edges_deg, delta_y = get_latitude_bands(num_latitude_bands)
    heat_capacity_local = latitude_heat_capacity(lat_centers_deg, ocean_frac_tropics, ocean_frac_poles, depth_ocean, depth_land)
    insolation_local = annual_avg_insolation(lat_centers_rad, L)
    albedo_local = latitude_albedo(lat_centers_deg, base_albedo, ice_latitude, polar_albedo)
    clear_sky_epsilon_local = latitude_clear_sky_epsilon(lat_centers_deg, epsilon_tropics, epsilon_poles)
    aerosol_forcing_local = latitude_aerosol_forcing(lat_centers_deg, base_aerosol_forcing, nh_midlat_peak_forcing, peak_latitude, peak_width)
    effective_albedo_local = albedo_local
    effective_epsilon_local = clear_sky_epsilon_local
    absorbed_solar_local = insolation_local * (1.0 - effective_albedo_local)
    total_direct_forcing_local = aerosol_forcing_local + additional_ghg_forcing

    if np.any(heat_capacity_local <= 0):
        st.error("Error: Calculated Heat Capacity is non-positive in some bands.")
        # *** Return signature must match expected output ***
        return None, None, None, "Input Error", None, None, True # Indicate error

    # --- Initialize with Kelvin Temperature ---
    if isinstance(initial_temp_profile_k, (int, float)):
        # Input MUST be Kelvin here
        T_current = np.full(num_latitude_bands, max(0.0, float(initial_temp_profile_k)))
    elif isinstance(initial_temp_profile_k, (np.ndarray, list)) and len(initial_temp_profile_k) == num_latitude_bands:
         # Input MUST be Kelvin here
        T_current = np.maximum(0.0, np.array(initial_temp_profile_k))
    else:
        st.error(f"Internal error: Initial temperature profile has wrong length or type.")
        return None, None, None, "Input Error", None, None, True

    heatContent = T_current * heat_capacity_local
    dt_seconds = dt_years * SECS_PER_YEAR
    time_list = np.arange(n_steps + 1) * dt_years
    temp_history_k = np.zeros((n_steps + 1, num_latitude_bands)) # Store history in Kelvin
    temp_history_k[0, :] = T_current
    stability_status = "Stable"
    blew_up = False
    unphysical_temps_detected = False
    T_plus = np.zeros_like(T_current)
    T_minus = np.zeros_like(T_current)
    d2T_dy2_approx = np.zeros_like(T_current)
    startup_steps = min(10, n_steps // 10)

    for i in range(n_steps):
        # --- Calculations remain in KELVIN ---
        outgoing_flux_local = effective_epsilon_local * SIGMA * (T_current**4)
        d2T_dy2_approx[1:-1] = T_current[2:] - 2 * T_current[1:-1] + T_current[:-2]
        d2T_dy2_approx[0]  = (T_current[1] - T_current[0]) * 2
        d2T_dy2_approx[-1] = (T_current[-2] - T_current[-1]) * 2
        transport_heating = heat_diffusivity_coeff * d2T_dy2_approx
        dH_dt = absorbed_solar_local - outgoing_flux_local + transport_heating + total_direct_forcing_local
        heatContent = heatContent + dH_dt * dt_seconds
        T_current = np.maximum(0.0, heatContent / heat_capacity_local) # T_current is KELVIN

        # --- Stability checks remain in KELVIN ---
        if not np.all(np.isfinite(T_current)) or np.any(T_current > 10000):
            stability_status = "Unstable/Blowup - Reduce Time Step or Diffusivity!"
            blew_up = True
            temp_history_k[i+1:, :] = np.nan
            break

        if i > startup_steps and not unphysical_temps_detected:
             # Check against Kelvin bounds
             if np.any(T_current < TEMP_MIN_PLAUSIBLE_K) or np.any(T_current > TEMP_MAX_PLAUSIBLE_K):
                 stability_status = "Oscillatory/Unstable - Unphysical Temps Occurred"
                 unphysical_temps_detected = True

        temp_history_k[i+1, :] = T_current # Store Kelvin history

    final_T_profile_k = temp_history_k[-1, :] if not blew_up else np.full(num_latitude_bands, np.nan)

    # --- Return results (still in KELVIN) ---
    return (time_list, temp_history_k, final_T_profile_k, stability_status,
            lat_centers_deg, lat_edges_deg, blew_up)


# --- Streamlit App ---
st.set_page_config(page_title="1D Climate Sim (°C)", layout="wide") # Changed title
st.title("🌍 Interactive 1D Climate Model") # Changed title
st.markdown("""
Simulate temperature (displayed in **Degrees Celsius**) across **latitude bands** with more realistic physics:
- **Latitude-dependent** heat capacity, insolation, albedo, emissivity, and aerosols.
- **Heat Transport via Diffusion** driven by temperature gradients.
Compare regions, view profiles, and explore climate dynamics! Parameters are grouped in expandable sections below.
""")

# --- Sidebar Inputs ---

with st.sidebar.expander("⚙️ Core & Grid Settings", expanded=True):
    L_input = st.number_input("Solar Constant (L, W/m²)", min_value=0.0, value=1361.0, step=1.0)
    num_bands_input = st.slider("Number of Latitude Bands", 4, 36, 18, 2, help="More bands = higher resolution, slower calculation.")
    # *** Initial Temp Input in Celsius ***
    initial_temp_input_c = st.number_input("Initial Temperature (°C)", value=7.0, min_value=-273.15, max_value=200.0, step=1.0, help="Uniform starting temperature for all bands.")

# (Other expanders remain largely the same, just update help text if needed)
with st.sidebar.expander("🏝️ Surface Heat Capacity", expanded=False):
    ocean_frac_tropics_input = st.slider("Ocean Fraction @ Equator", 0.0, 1.0, 0.75)
    ocean_frac_poles_input = st.slider("Ocean Fraction @ Poles", 0.0, 1.0, 0.5)
    depth_ocean_input = st.number_input("Ocean Effective Depth (m)", value=100.0, min_value=1.0, step=10.0)
    depth_land_input = st.number_input("Land Effective Depth (m)", value=2.0, min_value=0.1, step=0.5)

with st.sidebar.expander("💨 Heat Transport (Diffusion)", expanded=False):
    diffusion_input = st.number_input("Heat Diffusivity Coefficient (D)", value=0.55, min_value=0.0, max_value=5.0, step=0.05, format="%.3f", help="Controls efficiency of heat transport. Higher D = flatter T profile.")

with st.sidebar.expander("☀️ Albedo Profile", expanded=False):
    base_albedo_input = st.slider("Base Albedo (Non-Ice)", 0.0, 1.0, 0.20)
    ice_latitude_input = st.slider("Ice Line Latitude (°N/S)", 0, 90, 65)
    polar_albedo_input = st.slider("Polar Albedo (Ice/Snow)", 0.0, 1.0, 0.65)

with st.sidebar.expander("♨️ Clear Sky Emissivity Profile", expanded=False):
    epsilon_tropics_input = st.slider("Emissivity (ε) @ Equator", 0.01, 1.0, 0.55)
    epsilon_poles_input = st.slider("Emissivity (ε) @ Poles", 0.01, 1.0, 0.70)

with st.sidebar.expander("🌫️ Aerosol Forcing Profile", expanded=False):
    base_aerosol_input = st.number_input("Base Aerosol Forcing (W/m²)", value=-0.3, min_value=-5.0, max_value=5.0, step=0.1)
    nh_peak_aerosol_input = st.number_input("NH Mid-Lat Peak Forcing (W/m²)", value=-1.0, min_value=-5.0, max_value=5.0, step=0.1)
    peak_aerosol_lat_input = st.slider("Peak Forcing Latitude (°N)", 0, 90, 45)
    peak_aerosol_width_input = st.slider("Peak Forcing Width (°)", 5, 45, 15)

with st.sidebar.expander("🌍 Additional GHG Forcing", expanded=False):
    ghg_forcing_input = st.number_input("Uniform GHG Forcing (W/m²)", value=1.0, min_value=-5.0, max_value=10.0, step=0.1)

with st.sidebar.expander("⏱️ Simulation Control", expanded=True):
    dt_years_input = st.number_input("Time Step (years)", value=0.1, min_value=0.001, step=0.01, format="%.3f", help="Crucial for stability! Reduce if unstable.")
    n_steps_input = st.number_input("Number of Steps", value=500, min_value=1, max_value=20000, step=10)

total_sim_time = dt_years_input * n_steps_input
st.sidebar.info(f"Total Simulation Time: {total_sim_time:.1f} years")

with st.sidebar.expander("📈 Comparison Bands", expanded=True):
    # (Latitude zone selection logic remains the same)
    dummy_lats, _, _, _ = get_latitude_bands(num_bands_input)
    lat1_value = 0.0
    lat2_value = 0.0
    zone1 = "N/A"
    zone2 = "N/A"
    if len(dummy_lats) > 0:
        lat_options = {f"{lat:.1f}°": i for i, lat in enumerate(dummy_lats)}
        default_idx1 = len(lat_options) // 2
        default_idx2 = len(lat_options) - 1
        if default_idx1 == default_idx2 and default_idx1 > 0: default_idx1 -=1
        elif len(lat_options) == 1: default_idx1 = default_idx2 = 0
        band1_key = st.selectbox("Select Latitude Band 1", options=list(lat_options.keys()), index=default_idx1)
        band2_key = st.selectbox("Select Latitude Band 2", options=list(lat_options.keys()), index=default_idx2)
        idx1 = lat_options.get(band1_key, 0)
        idx2 = lat_options.get(band2_key, 0)
        try:
            lat1_value = float(band1_key.split('°')[0])
            zone1 = get_latitude_zone(lat1_value)
            st.caption(f"Band 1 Zone: {zone1}")
        except: zone1 = "Error"
        try:
            lat2_value = float(band2_key.split('°')[0])
            zone2 = get_latitude_zone(lat2_value)
            st.caption(f"Band 2 Zone: {zone2}")
        except: zone2 = "Error"
    else:
        st.sidebar.warning("Number of bands too small for selection.")
        idx1, idx2 = 0, 0


# --- Main Area ---
run_button = st.button("🚀 Run 1D Simulation")

results_container = st.container()
plot_container = st.container()
explanation_container = st.container()

if run_button:
    # *** Convert Celsius input temp to Kelvin before passing to simulation ***
    initial_temp_k = c_to_k(initial_temp_input_c)

    sim_results = run_latitude_simulation_enhanced(
        L=L_input, initial_temp_profile_k=initial_temp_k, # Pass K value
        dt_years=dt_years_input, n_steps=n_steps_input, num_latitude_bands=num_bands_input,
        ocean_frac_tropics=ocean_frac_tropics_input, ocean_frac_poles=ocean_frac_poles_input,
        depth_ocean=depth_ocean_input, depth_land=depth_land_input,
        heat_diffusivity_coeff=diffusion_input,
        base_albedo=base_albedo_input, ice_latitude=ice_latitude_input, polar_albedo=polar_albedo_input,
        epsilon_tropics=epsilon_tropics_input, epsilon_poles=epsilon_poles_input,
        base_aerosol_forcing=base_aerosol_input, nh_midlat_peak_forcing=nh_peak_aerosol_input,
        peak_latitude=peak_aerosol_lat_input, peak_width=peak_aerosol_width_input,
        additional_ghg_forcing=ghg_forcing_input
    )

    if sim_results is not None:
        # *** Results arrive in KELVIN ***
        (time_list, temp_history_k, final_T_profile_k, stability,
         lat_centers_deg, lat_edges_deg, blew_up) = sim_results

        # *** Convert outputs to CELSIUS for display ***
        final_T_profile_c = k_to_c(final_T_profile_k)
        temp_history_c = k_to_c(temp_history_k)

        with results_container:
            st.subheader("Simulation Results")
            # (Status display remains the same)
            if "Unstable/Blowup" in stability: st.error(f"Status: {stability}")
            elif "Unphysical Temps" in stability: st.warning(f"Status: {stability} - Check parameters or reduce time step.")
            elif "Input Error" in stability: st.error(f"Status: {stability} - Check input parameters.")
            else: st.success(f"Status: {stability}")

            # *** Calculate and display metrics in CELSIUS ***
            avg_final_temp_c = np.nanmean(final_T_profile_c)
            T1_final_c = final_T_profile_c[idx1] if not blew_up and final_T_profile_c.shape[0]>idx1 else np.nan
            T2_final_c = final_T_profile_c[idx2] if not blew_up and final_T_profile_c.shape[0]>idx2 else np.nan

            col1, col2, col3 = st.columns(3)
            col1.metric("Global Avg Final Temp", f"{avg_final_temp_c:.2f} °C" if np.isfinite(avg_final_temp_c) else "N/A")
            col2.metric(f"Final Temp @ {band1_key} ({zone1})", f"{T1_final_c:.2f} °C" if np.isfinite(T1_final_c) else "N/A")
            col3.metric(f"Final Temp @ {band2_key} ({zone2})", f"{T2_final_c:.2f} °C" if np.isfinite(T2_final_c) else "N/A")

        with plot_container:
            st.subheader("Visualizations")
            if "Input Error" not in stability:
                tab1, tab2, tab3 = st.tabs(["Comparison Plot", "Final Latitude Profile", "Temperature Heatmap"])
                can_plot = not blew_up # True if simulation didn't return NaN/inf

                # Define plot range in Celsius using converted bounds
                plot_ymin = TEMP_MIN_PLAUSIBLE_C - 20
                plot_ymax = TEMP_MAX_PLAUSIBLE_C + 50

                with tab1:
                    if can_plot:
                        fig_comp = go.Figure()
                        # *** Plot using Celsius history ***
                        fig_comp.add_trace(go.Scatter(x=time_list, y=temp_history_c[:, idx1], mode='lines', name=f'Temp @ {band1_key} ({zone1})'))
                        fig_comp.add_trace(go.Scatter(x=time_list, y=temp_history_c[:, idx2], mode='lines', name=f'Temp @ {band2_key} ({zone2})'))
                        fig_comp.update_layout(title=f'Temperature Evolution: {band1_key} ({zone1}) vs {band2_key} ({zone2})',
                                               xaxis_title='Time (years)', yaxis_title='Temperature (°C)') # Y-axis label change
                        fig_comp.update_yaxes(range=[plot_ymin, plot_ymax])
                        st.plotly_chart(fig_comp, use_container_width=True)
                        st.markdown("Shows how temperature changes over time for the two selected latitude bands.")
                    else:
                        st.warning("Cannot plot comparison - simulation unstable.")

                with tab2:
                    if can_plot:
                        fig_prof = go.Figure()
                         # *** Plot using Celsius profile ***
                        fig_prof.add_trace(go.Scatter(x=lat_centers_deg, y=final_T_profile_c, mode='lines+markers', name='Final Temperature'))
                        fig_prof.update_layout(title='Final Temperature vs. Latitude', xaxis_title='Latitude (°)', yaxis_title='Temperature (°C)') # Y-axis label change
                        fig_prof.update_yaxes(range=[plot_ymin, plot_ymax])
                        st.plotly_chart(fig_prof, use_container_width=True)
                        st.markdown("Displays the temperature profile across all latitudes at the *end* of the simulation.")
                    else:
                        st.warning("Cannot plot profile - simulation unstable.")

                with tab3:
                    if can_plot:
                        if temp_history_c.ndim == 2 and temp_history_c.shape[0] > 0 and temp_history_c.shape[1] > 0:
                             fig_heat = go.Figure(data=go.Heatmap(
                                z=temp_history_c.T, # Use Celsius data
                                x=time_list,
                                y=lat_centers_deg,
                                colorscale='Viridis',
                                colorbar_title='Temp (°C)', # Color bar label change
                                zmin=plot_ymin, zmax=plot_ymax # Color bar range in Celsius
                            ))
                             fig_heat.update_layout(title='Temperature Heatmap (Latitude vs. Time)',
                                                      xaxis_title='Time (years)', yaxis_title='Latitude (°)')
                             st.plotly_chart(fig_heat, use_container_width=True)
                             st.markdown("Visualizes the temperature in all latitude bands over the entire simulation duration.")
                        else:
                            st.warning("Cannot plot heatmap - insufficient data generated.")
                    else:
                         st.warning("Cannot plot heatmap - simulation unstable.")
            else:
                st.warning("Cannot display plots due to input error.")
    else:
         with results_container:
            # Error message should have been displayed by st.error within simulation function
            pass


# --- Footer / Explanations ---
with explanation_container:
    st.markdown("---")
    st.header("Model Details & Caveats")
    # *** Update text to mention Celsius display ***
    st.markdown("""
    - **1D Energy Balance Model (EBM):** Averages properties longitudinally within each latitude band. It captures the fundamental equator-to-pole temperature gradient. Temperatures are displayed in **Degrees Celsius (°C)**.
    - **Heat Capacity:** Varies with latitude based on a simplified ocean/land fraction profile specified by the user. Oceans store much more heat than land over seasonal/annual cycles.
    - **Heat Transport (Diffusion):** Heat flows from warmer to colder regions driven by the temperature *curvature*. The **Diffusivity Coefficient (D)** controls how effectively heat is mixed; higher values lead to a smaller equator-pole temperature difference. *Note: The scaling/units of D here are simplified for the interface; tuning may be needed depending on the exact physics desired.*
    - **Emissivity (`ε`):** Varies with latitude to represent higher infrared trapping (lower `ε`) by water vapor in the tropics, based on user-defined equator and pole values.
    - **Aerosol Forcing:** Can be set to vary with latitude using a Gaussian profile centered at a specified latitude, representing, e.g., industrial pollution. Negative values typically represent cooling.
    - **GHG Forcing:** Applied uniformly as it represents well-mixed gases like CO₂.
    - **Stability:** Diffusion models often require **smaller time steps (`dt_years`)** than simpler models, especially with higher diffusivity or fewer latitude bands. If the simulation becomes unstable (`NaN` values, extreme temperatures > 10000K) or shows unphysical oscillations (temps outside approx. -123°C to 177°C), the status will indicate this. **Reduce the time step** or the diffusivity coefficient if instability occurs. (Note: Internal physics calculations use Kelvin).
    - **Simplifications:** Clouds are currently highly simplified (effectively ignored). No seasons, ocean currents, complex atmospheric dynamics, or dynamic ice-albedo feedback (ice line is fixed) are included.
    """)
