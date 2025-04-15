# 🌍 Interactive 1D Climate Model Explorer

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://energy-balance-model.streamlit.app/) 

Welcome! This project provides an interactive web application built with Streamlit to explore a **1-Dimensional Energy Balance Model (EBM)** of planetary climate.

Climate science is complex, involving intricate interactions across space and time. This tool simplifies some of this complexity by focusing on the fundamental **energy balance** of the planet across different **latitude bands**, from the equator to the poles. It allows you to visually explore how factors like incoming sunlight, surface properties (albedo, heat capacity), atmospheric effects (emissivity, forcing), and crucially, **heat transport**, shape the planet's temperature profile.

Whether you're a student learning about climate dynamics, an educator looking for a demonstration tool, or just curious about climate modeling, this app provides a hands-on way to build intuition about the Earth's climate system (or even hypothetical planets!).

## 🚀 Live Demo

**[Access the interactive application here!](https://energy-balance-model.streamlit.app/)** <!-- *** Replace with your actual deployment URL! *** -->

*(Please allow a few moments for the application to load initially.)*

## ✨ Features

*   **1D Latitude Resolution:** Simulates temperature in distinct bands from -90°S to +90°N.
*   **Energy Balance Core:** Models the fundamental balance between absorbed solar radiation, outgoing infrared radiation, and energy storage.
*   **Latitude-Dependent Physics:** Includes realistic variations in:
    *   **Insolation:** Sunlight intensity naturally decreases towards the poles.
    *   **Albedo:** Higher reflectivity over icy polar regions.
    *   **Heat Capacity:** Accounts for differences between land and ocean dominance across latitudes.
    *   **Emissivity:** Represents variations in atmospheric trapping (e.g., due to water vapor).
    *   **Aerosol Forcing:** Allows specifying geographically varying aerosol effects.
*   **Heat Transport:** Implements meridional (equator-to-pole) heat transport using a **diffusion** mechanism, crucial for moderating temperature gradients.
*   **Interactive Parameter Control:** Use sidebar widgets to adjust numerous planetary and simulation parameters in expandable sections.
*   **Multiple Visualizations:**
    *   **Comparison Plot:** Track temperature evolution over time for two user-selected latitude bands.
    *   **Final Latitude Profile:** See the equilibrium temperature distribution across all latitudes.
    *   **Temperature Heatmap:** Visualize the full temperature field `T(latitude, time)` throughout the simulation.
*   **Stability Monitoring:** Provides feedback on the numerical stability of the simulation based on the chosen parameters.

## 🤔 How Does It Work? (Simplified Science Background)

At its heart, this is an **Energy Balance Model (EBM)**. The core idea is that the temperature of a region changes based on the net energy it receives:

`Change in Heat = Absorbed Sunlight - Outgoing Infrared + Heat Transport In/Out + Direct Forcing`

*   **Absorbed Sunlight:** Depends on the incoming Solar Constant and how much is reflected (Albedo). Varies strongly with latitude (sun angle).
*   **Outgoing Infrared:** How the planet radiates heat back to space. Depends on temperature (Stefan-Boltzmann Law: `σT⁴`) and atmospheric trapping (Emissivity, `ε`). Lower emissivity means more trapping (stronger greenhouse effect).
*   **Heat Transport (Diffusion):** Without this, the equator would be scorching and the poles frozen solid! Heat naturally flows from warmer regions to colder regions. This model uses diffusion, where the flow is proportional to the *curvature* of the temperature profile (`d²T/dy²`). The `Diffusivity Coefficient (D)` controls how efficiently this happens.
*   **Direct Forcing:** Represents additional energy inputs/outputs, like the warming from extra Greenhouse Gases (GHGs) or the cooling/warming effect of Aerosols.

The model calculates these energy fluxes for each latitude band and steps forward in time, updating the temperature based on the net energy change and the region's **Heat Capacity** (how much energy it takes to change its temperature, influenced by land vs. ocean).

## 🎮 How to Use the App

1.  **Access the Live Demo:** Click the link above.
2.  **Explore the Sidebar:** Parameters are grouped into expandable sections:
    *   **⚙️ Core & Grid Settings:** Adjust the overall solar input, number of latitude bands (resolution), and starting temperature.
    *   **🏝️ Surface Heat Capacity:** Define how much ocean vs. land influences heat storage at different latitudes.
    *   **💨 Heat Transport (Diffusion):** Control the efficiency of heat mixing between latitudes with the `D` coefficient.
    *   **☀️ Albedo Profile:** Set the reflectivity of ice-free areas and define where the high-reflectivity ice/snow line begins.
    *   **♨️ Clear Sky Emissivity Profile:** Define the baseline greenhouse effect, allowing it to vary between the potentially wetter tropics (more trapping) and drier poles.
    *   **🌫️ Aerosol Forcing Profile:** Simulate the effect of aerosols, potentially peaking them in certain latitude bands (e.g., NH mid-latitudes).
    *   **🌍 Additional GHG Forcing:** Add a uniform warming effect representing well-mixed greenhouse gases.
    *   **⏱️ Simulation Control:** Set the **Time Step** (crucial for stability!) and the **Number of Steps** (total simulation duration).
    *   **📈 Comparison Bands:** Choose two specific latitude bands (and see their climate zone) to compare in detail.
3.  **Run Simulation:** Click the "🚀 Run Enhanced 1D Simulation" button.
4.  **Analyze Results:**
    *   Check the **Status** message (Stable, Oscillatory, Unstable).
    *   Observe the **Final Temperatures** (Global Average and selected bands).
    *   Explore the **Visualizations** using the tabs:
        *   See how the selected bands evolve over time.
        *   Examine the final temperature gradient across all latitudes.
        *   Analyze the heatmap to see the spatial and temporal patterns.
5.  **Experiment!** Change parameters (e.g., increase GHG forcing, decrease Diffusivity, move the ice line) and re-run the simulation to see how the climate responds.

## 🛠️ Running Locally (Optional)

If you want to run or modify the code yourself:

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/your-username/your-repo-name.git 
    # Replace with your actual repo URL
    cd your-repo-name
    ```
2.  **Install Requirements:**
    ```bash
    pip install -r requirements.txt
    ```
3.  **Run the App:**
    ```bash
    streamlit run app.py 
    # Or the specific name of your main Python script
    ```

## ⚠️ Limitations

This model, while illustrative, has several simplifications:

*   **Zero-Dimensional Atmosphere:** Atmospheric effects (emissivity, forcing) are parameterized simply; there's no explicit atmospheric model.
*   **No Seasons:** Uses annual average insolation.
*   **Simplified Heat Transport:** Diffusion is a basic representation; real transport involves complex atmospheric and oceanic currents.
*   **No Ocean Dynamics:** Heat capacity uses an "effective depth"; no explicit ocean circulation.
*   **Fixed Ice Line:** The albedo change happens at a fixed latitude; real ice extent responds dynamically to temperature (ice-albedo feedback).
*   **Simplified Clouds:** Cloud effects on albedo and emissivity are currently minimal or ignored (a potential area for future work).
*   **Numerical Method:** Uses a basic Forward Euler time-stepping scheme, requiring careful choice of time step (`dt`) for stability, especially with diffusion.

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🙏 Acknowledgments

*   Global Warming I & II: A Course by University of Chicago; taught by [David Archer](https://www.coursera.org/instructor/davidarcher)
*   Streamlit, Google AI studio

## 📧 Contact

Binamra Bhusal | M.Sc. Climate Change & Development | https://binamra4bhusal@gmail.com

---
