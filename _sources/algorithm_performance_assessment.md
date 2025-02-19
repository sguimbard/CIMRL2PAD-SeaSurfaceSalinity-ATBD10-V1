# Algorithm Performance Assessment

To assess the performance of the previously described OWV retrieval algorithm, we used the SCEPS prototype simulator to generate a very large test card (~6000 km × 4000 km), with 1 km resolution including strong SSS gradients due to the Amazon River plume, and strong surface wind speed gradients with the presence of two hurricanes (see {numref}`Figure25`).

As described in RD-2, the simulator consists of a single statically-linked Linux executable program, **cimrProject**, which produces a Level 1b antenna temperature product from input scene brightness information. The program takes as input an EO-CFI orbit definition file (Team, 2022), a file containing the scene brightness temperature in the form of a *test card*, complete GRASP-generated antenna patterns provided by industry, and a file that specifies the sample times at which to produce antenna temperatures. 

The scene brightness includes both the isotropic and anisotropic complete modified Stokes vector ($T_h$, $T_p$, $U$, $V$) in the surface polarization basis. The isotropic part can vary with incidence angle, while, for practical reasons in our analyses, the anisotropic part is assumed to be independent of incidence angle.

```{figure} Figure24.png
---
name: Figure24
---
Schematic of the retrieval algorithm assessment using SCEPS
```

A schematic of the overall retrieval algorithm assessment using SCEPS is provided in {numref}`Figure24`. The assessment is based upon these successive steps:

1. The scene is first described by geophysical fields, including: SSS, SST, OWV at ~1 km spatial resolution and then interpolated (in space and time) atmospheric vertical parameters from ERA5.

2. The TOA brightness temperatures are then evaluated for all CIMR frequencies and polarization using the previously described RTM for an ensemble of incidence angles. The azimuthal harmonic contributions to total TOA temperature are also evaluated separately.

3. The computation of (noisy) simulated CIMR antenna temperature from the scene brightness temperature field involves several aspects:
   - 3.1 Calculation of orbit and viewing geometry,
   - 3.2 Calculation of the scene brightness Stokes vector,
   - 3.3 Integration of the scene brightness over the antenna patterns (both near and far from the feed boresights) for each horn.
   - 3.4 Addition of NEDT noise to the antenna temperatures,
   - 3.5 Calculation and application of the Antenna Pattern Correction (APC) matrices for all feeds, to reduce/remove cross-pol contamination. For these tests, this APC has the polarization basis rotation *built-in*;
   - 3.6 Production of a preliminary version of the Level 1b product.

4. Resampled brightness temperature are then computed by applying the Backus-Gilbert method to express in the surface polarization basis:
   - SCEPS uses elliptical power patterns (on Earth) to compute the BG inner products.
   - Fore and aft sampled can be combined or kept separated in fore and aft images.
   - The Target grid is the original test card at 1 km resolution and the target patterns are all circular Gaussians with FWHM for L, C, X, Ku, and Ka band defined as follows: 60, 30, 30, 8, and 8 km.
   - A single application of BG yields weights that can be applied to all antenna pattern contributions.
   - Maps of error variance amplification factors are computed from the BG weights.

5. The noisy simulated resampled CIMR L1b antenna temperature from the scene brightness temperature field are then used for simulation of the future CIMR remapped measured L1B Tbs.


## Geophysical Data Used to Describe the Ocean Scene

The synthetic wind field used to model the scene is a merge between ERA5 re-analyses used as the background wind field, and more local and high-resolution surface winds from Sentinel-1/2 SAR data around two hurricanes. The SSS is from the MERCATOR model. The SST is from ERA5 re-analyses, which were also used to simulate the atmospheric vertical profiles. The test scene includes several coasts to analyze the land-contamination.

```{figure} Figure25.png
---
name: Figure25
---
Input geophysical data used to simulate the test card: (top left): SSS from Mercator; (top right): ERA5 SST and (bottom) merged ERA5 & SAR images (over the two hurricanes) surface wind vectors.
```

## Simulated TOA Tbs over the Scene

Using these geophysical data as input, we generated TOA L- to Ka-band TOA Tbs using the previously described RTM at an ensemble of EIA with values of 0°, 20°, 40°, 50°, 55°, 60°, 70°, 80° including azimuthal harmonics for the full Stokes vector. Note: In the Forward Model simulations, we assumed the anisotropic roughness effects to be similar at all EIA values. ERA5 vertical profiles for Cloud Liquid, Ice Water Content, air temperature, humidity were used to evaluate the atmospheric contributions at each frequency. Land Tbs > 300K are from the SCEPS project (courtesy Carlos Jimenez). Example of such Tb fields are shown in {numref}`Figure26`.

```{figure} Figure26.png
---
name: Figure26
---
(Left) L-band H-pol; (right) L-band V-pol based on the Radiative forward model (EIA = 52°).
```

As described in previous sections, the roughness-induced emission $\Delta T_p$ can be decomposed into isotropic and anisotropic emission components using a second-order azimuthal harmonic expansion as a function of $\tilde{\varphi}$, the relative azimuth between the radiometer look direction and the wind direction. 

The Stokes vector of the total brightness temperature emitted by the sea surface at location $(x_1, x_2)$ and at an EIA of $\theta_s$ thus reads:

$$
T_h(x_1, x_2, \theta_s) = T_{hso} + T_{ho}(U_{10}) + T_{hc1}(U_{10}) \cos(\tilde{\varphi}) + T_{hc2}(U_{10}) \cos(2\tilde{\varphi})
$$

$$
T_v(x_1, x_2, \theta_s) = T_{vso} + T_{vo}(U_{10}) + T_{vc1}(U_{10}) \cos(\tilde{\varphi}) + T_{vc2}(U_{10}) \cos(2\tilde{\varphi})
$$

$$
U(x_1, x_2, \theta_s) = U_{hs1}(U_{10}) \sin(\tilde{\varphi}) + U_{hs2}(U_{10}) \sin(2\tilde{\varphi})
$$

$$
V(x_1, x_2, \theta_s) = V_{hs1}(U_{10}) \sin(\tilde{\varphi}) + V_{hs2}(U_{10}) \sin(2\tilde{\varphi})
$$


As found, only the horizontal and vertical polarization exhibit isotropic contributions (0 order azimuthal harmonics), namely: $T_{ho}(U_{10})$ and $T_{vo}(U_{10})$. The amplitude of these isotropic components is significantly larger than the first-order azimuthal harmonics: $(T_{hc1}(U_{10}), T_{vc1}(U_{10}), U_{hs1}(U_{10}), V_{hs1}(U_{10}))$ and the second-order azimuthal harmonics: $(T_{hc2}(U_{10}), T_{vc2}(U_{10}), U_{hs2}(U_{10}), V_{hs2}(U_{10}))$. Note that the third $U$ and fourth $V$ Stokes vector components are non-zero only because of anisotropic features of the sea surface roughness.

For each frequency, SCEPS integrates the harmonic contributions (see {numref}`Figure27`): 

$$
T_{hc1}(U_{10}), T_{vc1}(U_{10}), U_{hs1}(U_{10}), V_{hs1}(U_{10}), T_{hc2}(U_{10}), T_{vc2}(U_{10}), U_{hs2}(U_{10}), V_{hs2}(U_{10})
$$

to total antenna temperature separately and stores them in an extended L1B netCDF file.

```{figure} Figure27.png
---
name: Figure27
---
Example of the TOA scene azimuthal harmonics ($T_h$, $T_p$, $U$, $V$) simulated from the surface wind field in {numref}`Figure25` for L-band.
```

Here, we specified wind direction-relative harmonics together with the wind direction at each grid point. Note that the harmonics other than 0 cannot vary with incidence angle and are provided as a single value at each horizontal grid point.

## Simulated Instrument Tbs Over the Scene

The computation of (noisy) simulated CIMR antenna temperature from the scene brightness temperature field then involves several aspects:

1. Calculation of orbit and viewing geometry,
2. Calculation of the scene brightness Stokes vector for an ensemble of incidence angles, and of the wind-relative azimuthal harmonics components,
3. Integration of the scene brightness over the antenna patterns (both near and far from the feed boresights),
4. Addition of NEDT noise to the antenna temperatures,
5. Calculation and application of the antenna pattern correction matrices for all feeds,
6. Production of a preliminary version of the Level 1B product,
7. Resampling of the TOA Tbs using Backus-Gilbert for all bands except at L-band for which a weighting average is applied,
8. An OZA adjustment is then simply computed by interpolating the test card Tbs (without harmonics) to the reference and actual OZA for each sample/measurement using those values directly without integrating over the antenna pattern. Then, the resampler resamples the interpolated test card Tbs at the OZA values as well as the unadjusted antenna Tbs separately so as to keep track of the adjustment. The OZA correction is then applied after resampling as a simple addition of the resampled adjustment.

Details about each step of the SCEPS simulator can be found in [RD-2]. The simulator includes an orbit propagator to simulate the orbiting satellite positions in time and an antenna geometry routine to simulate the antenna scanning geometry and integrate the brightness fields over the antenna power patterns. 

Noting that most of the power originates near the feed boresight, the integration domain is divided into an inner and an outer domain. These domains are defined for each feed (see {numref}`Figure28`): The inner domain includes the boresight and extends out to an ellipse corresponding to an azimuthally averaged power 40 dB below the peak (at the feed boresight), where the azimuthal average is computed in the tilted direction cosine coordinate system. The outer domain covers the rest of the Front Half-Space (FHS) of antenna (containing boresight) in either the tilted or un-tilted coordinate system.

```{figure} Figure28.png
---
name: Figure28
---
Inner (left) and outer (right) integration grids in untilted director cosine (u, v) coordinates for L. Only the directivities down to -40 dB and -60 dB below the peak are shown for the inner and outer grids, respectively. The color indicates the antenna pattern power relative to peak.
```


The CIMR instrument consists of 25 feeds, with one at L-band, 4 at C/X bands, and 8 at Ku/Ka bands. {numref}`Figure29` shows the measurement (full integration time) footprints for a portion of three successive scans for the L-, C-band, X-band, and Ku-band horns. The figure illustrates the extent of overlap between successive footprints both along and across the scan near the swath center. Although overlap can be increased (at the expense of increased noise) along the scan by using individual samples, the across-scan overlap cannot be changed. To obtain a resolution consistent with the antenna pattern footprint size, distance between neighboring measurements should not exceed half the Half-Power Beam Width (by the Whittaker–Nyquist–Shannon sampling theorem).

As discussed in [RD-2], this condition is not satisfied across scan in L and Ku bands, but is satisfied in C and X bands. For Ku-band, there is no overlap in the half-power footprints across scan, which may be a limiting factor in the resolution of images obtained from resampled antenna temperatures.

The figure illustrates the extent of overlap between successive footprints both along and across scan near the swath center. Although overlap can be increased (at the expense of increased noise) along scan by using individual samples, the across-scan overlap cannot be changed.

```{figure} Figure29.png
---
name: Figure29
---
(a): Projection of the L-band feed half-power footprints and boresight locations onto the ground at two closely separated times for three successive scans.
```
In SCEPS, the total antenna temperature integral over all space is separated into several parts, with different integration methods for each:

- Inner domain: Extending from each feed boresight out to a circle (in tilted director cosine coordinates) corresponding to an isotropic directivity of -40 dB relative to the maximum at feed center (HBS).
- Outer domain: Extending from the boundary of the inner domain out to the Earth limb in the untilted frame (OUT).
- Contribution from distributed foreign sources: Cold sky and galactic radiation in nominal instrument pointing (sky).
- Contribution from localized sources: Bright spots in the test card scene (bs).
- Contribution from the direct sun and moon: Without reflection from the Earth.

In SCEPS, the default reflector rotation rate is set to the fixed value of 7.8 rpm (clockwise looking towards the instrument along the spin axis from behind the reflector), and the radiometer integration times are set to the values shown in {numref}`tablerecparam`. 

SCEPS can be configured to produce antenna temperatures for only the full measurements or for all samples. In the present simulation, there are five samples per measurement.

```{table} Receiver characteristics. Both the antenna temperature and the effective receiver noise temperature ($T_{rn}$) corresponding to the NEDT values are taken to be 150 K in SCEPS. There are five samples per measurement for all bands.
:name: tablerecparam
| Band Name | Number of Feeds | Center Frequency (GHz) | Maximum Bandwidth (MHz) | Sample Integration Time (ms) | Measurement NEDT at $T_a$ = $T_{rn}$ = 150 K |
|----|----------|-------|---------|-----------|----------|
| L   | 1     | 1.4135  | 50.0 | 11.12  | 0.3   |
```

The reflector rotation rate is sufficiently fast that the antenna patterns move substantially during the measurement integration times, and this must be accounted for when SCEPS is configured to compute antenna temperatures for the full measurements. To do so, SCEPS extends the antenna pattern integration to include integration over the integration time of the receiver.

For L-band, one uses a weighted average of the original scene $T_Bs$ based on a Gaussian window of FWHM = 30 km. Examples of simulated TOA brightness temperatures integrated over CIMR antenna patterns for Aft view of two successive ascending passes intercepting the scene are shown in {numref}`Figure30`.

Note that instrument measurements are simulated from scene $T_p$ using SCEPS with several options:

$\Rightarrow$ L1B TOA + NEDT (constant over the antenna pattern),

$\Rightarrow$ L1B TOA without NEDT.

$\Rightarrow$ L1B TOA integrated over inner antenna patterns (HBS): Boresight to -40 dB,

$\Rightarrow$ L1B TOA integrated over full antenna patterns (TOT): 0 dB to -60 dB,

$\Rightarrow$ L1B TOA $T_p$ are then remapped into L1C TOA $T_p$ using Backus-Gilbert interpolation (all frequencies except L-band) or using « gaussian weighted-average » for L-band (do not respect Nyquist sampling)

$\Rightarrow$ L1C TOA $T_p$ are then generated for both Aft and Fore views and for 2 ascending passes intercepting the scene.


```{figure} Figure30.png
---
name: Figure30
---
Simulated TOA L-band brightness temperatures integrated over CIMR antenna patterns for Aft view of two successive ascending passes intercepting the scene (left) V-pol; (right) H-pol.
```

## SSS Retrieval simulations

Examples of SSS retrieved over the scene is shown in {numref}`Figure31`. In these retrieval simulations, the instrument Tbs were derived by integrating the scene over the full antenna patterns (HBS+outer patterns). SCEPS provides the option to add NEDT noise to the antenna temperatures: we perform retrieval simulations with and without NEDT, and in the following, we present only results with NEDT. As shown in {numref}`Figure31`, the $\Delta SSS=SSS_{retrieved}-SSS_{input}$ fields exhibit high resolution features (North Brazilian current induced SSS fronts and Tropical Cyclone center) which represent area where the retrieved SSS are not well resolved by the ~60 km footprint of CIMR L-band data. The ∆SSS distribution is also skewed with large negative values along the coasts. We analyzed the $\Delta$SSS as a function of distance to nearest coasts in bins of 10 km. As illustrated in {numref}`Figure32`, a significant fresh bias and increased standard deviation of the ∆SSS is observed when the data are acquired at distance to nearest coasts less than ~70 km.

```{figure} Figure31.png
---
name: Figure31
---
 (Top Left) Input 4km resolution SSS field (Top Right) CIMR retrived SSS from SCEPS generated L1C data. (Bottom Left) Differences between input and retrieved SSS (Bottom Right) Histrogram of the differences between input and retrieved SSS.
```

```{figure} Figure32.png
---
name: Figure32
---
Median (Left) and standard deviation (Right) of the differences between input and retrieved SSS binned as a function of the distance to nearest coasts.
```

The closer to the coast, the higher the fresh bias and standard deviation of the $\Delta$SSS.

```{figure} Figure33.png
---
name: Figure33
---
Filtered $\Delta$SSS when data with distance to nearest coasts less than 70 km are removed. The histogram of ∆SSS is shown in the right panel.
```

While a land-contamination correction will be potentially derived from the future CIMR data (e.g., using the approach of {cite:p}`Meissner2017,meissner2018salinity`, flagging/removing the data within a distance of 70 km from the nearest coasts, will solve most of the land-contamination issues. When removing data within 70 km from the coasts, the statistics of the ∆SSS becomes more gaussian with a significantly reduced STD of ~0.2 pss (see {numref}`Figure33`), compared to ~1.7 pss when considering all data including the SSS retrieved within a band of 70 km from nearest coasts (see {numref}`Figure31`).

## SSS Retrieval simulations: comparison with MRD

As stated in the CIMR MRD, the mission objectives for SSS retrieval are the following:
“SEC-OBJ-9. Measure Sea Surface Salinity (SSS) over the global ocean from space [RD-3],[RD-4] with a target gridded spatial resolution of 40 km and uncertainty $\leq$0.3 pss over monthly time-scales [RD-3]”

The accuracy of weekly or monthly SSS products from orbiting satellites depends on the number of satellite passes over a given grid cell within the period of time. To evaluate such number for the CIMR L-band data collected over one month, in the MACRAD project [RD-3], we have conducted a sampling analysis with the satellite instrument retrieval simulator for the mission: the CIMR SCEPS. We estimated the number of L1B samples per day in 0.36°x 0.36° rectangular grid (~40km resolution) boxes for the complete month of January 2029. We reproduced some of the results hereafter.


### Error Estimates for individual CIMR L2 SSS

In the following, we first estimated the accuracy of the CIMR satellite L2 SSS measurements from the effects of the various error sources at a nominal incidence angle of 52° that were described in detail in [RD-3]. Table 1 summarizes the influence of various geophysical error sources for four water temperature ranges:

1. Very low SSTs with -4°C < $T_s$ < 0°C
2. Low SSTs: 0°C < $T_s$ < 10°C
3. Moderate SSTs: 10°C < $T_s$ < 20°C
4. Warm SSTs: 20°C < $T_s$ < 30°C

And for each SST range, we investigated the errors for the following three wind speed ranges:

5. Low wind speed: $U_{10} < 5$ m/s
6. Moderate wind speed: $5 \leq U_{10} < 15$ m/s
7. High wind speed: $U_{10} > 15$ m/s

In [RD-3], we performed sensitivity analyses from the RTM presented in this ATBD. For each temperature and wind speed range, the estimated brightness temperature errors were listed for $T_v$, $T_h$, and $(T_v + T_h)/2$. The estimates for the average brightness temperature $(T_v + T_h)/2$ were included because it is insensitive to the Faraday rotation. 

We provided the sensitivity of the brightness temperatures to SSS and SST, indicating a varying degree of sensitivity versus polarization and water temperature, as well as the sensitivity to surface wind speed and direction. We also estimated the sensitivity to air surface temperature, surface pressure, columnar water vapor content, and cloud liquid water content. Additionally, we provided the residual brightness temperature errors after data corrections or under an assumed threshold for data flags related to external sources (e.g., galactic and solar glints). Land-sea contamination was filtered (e.g., all data within a 70 km band from coasts were not considered).

The net effect of geophysical errors is described in {numref}`Figure34` under the assumption that all geophysical error sources are uncorrelated.

```{figure} Figure34.png
---
name: Figure34
---
Distribution of uncertainties in the modeling of the V-polarization L-band Tbs at incidence angle of 52° in a) cold sea condition and b) warm sea conditions. The wind speed is moderate with 5 m/s $\leq$ $U_{10}$ < 15 m/s.
```

For individual L-band measurements, the square Root of the Sum of Squares (RSS) of quantity is varying from about 0.5 K to 1.2 K, depending on SST and Wind Speed and on polarization. One of the dominant error sources is the surface roughness (e.g., wind speed). The second most impacting effect is NEDT and radiometer stability. Contrarily to V-polarization (see {numref}`Figure35`), the RSS do not increase monotonically with increasing wind speed for H-polarization because the sensitivity to wind speed is higher for low wind speed than for moderate and higher winds (see RD-3). As a reference, 0.1 K corresponds to a wind speed uncertainty of 0.5 m/s- for vertical polarization and 0.25 m/s- for horizontal polarization.

```{figure} Figure35.png
---
name: Figure35
---
RSS of total Tb errors for a) V-polarization; b) H-polarization and c) First Stokes parameter for the 3 wind speed and 4 SST regimes.
```

The sensor errors include the calibration errors and the radiometer noise-equivalent-delta-T (NEDT). The instrument calibration stability does not include the bias of the instrument calibration errors, and only accounts for the temporal variability of the instrument calibration. It is assumed that the instrument calibration bias is a constant and hence can be removed by a comparative analysis of in-situ measurements. The dominant error sources are the antenna pointing angle, antenna emissivities, and the stability of the calibration device, such as the noise source, for the electronics. The radiometer NEDT is assumed to be 0.3 K for each polarization.

As found, the vertical polarization outperforms the other polarization and combinations because of its superior sensitivity to salinity. In moderate SST and wind conditions, it is about ~1 pss. It increases with increasing wind speed and decreasing SST. The average brightness temperature $(T_v + T_h)/2$ has the advantage of a negligible Faraday rotation error, but its performance is slightly worse than the vertical polarization because of a slight degradation of SSS sensitivity.

```{figure} Figure36.png
---
name: Figure36
---
RSS of total SSS errors for a) V-polarization; b) H-polarization and c) First Stokes parameter for the 3 wind speed and 4 SST regimes.
```

In V-polarization, the accuracy ranges from ~0.9 pss to ~3.5 pss depending on wind speed and SST conditions. In H-pol, it is significantly higher, ranging from ~2.5 to 9 pss. The accuracy indicated here is for one satellite observation and approaches 2 pss (vertical polarization) for cold waters and improves with increasing water temperatures to reach about 0.9 psu for warm waters.

### Error Estimates for Level 3 CIMR SSS

To achieve better accuracy, multiple independent observations are needed. The independent measurements can possibly result from temporal or spatial averaging. The temporal averaging involves the average of data from multiple satellite passes over a given surface grid cell. The revisit time is a few hours for high latitudes and 1 to 3 days for equatorial regions; therefore, it is reasonable to assume that the data from different satellite passes are uncorrelated. In contrast, the spatial averaging over the measurements from adjacent antenna footprints may not reduce the error because most of the geophysical error sources and instrument calibration uncertainty are likely to have a spatial correlation.

For the CIMR conical scanner, a correlation time of a few minutes for the instrument calibration error will make the measurements across the entire swath and along the track of about 1800 km correlated. The only error source that definitely can be improved by spatial averaging is the sensor NEDT, which is random. To be conservative, we assume only the measurements from different satellite passes as independent estimates. Under this assumption, the number of satellite passes required to reduce the error to 0.3 pss has been evaluated in RD-3. It is about 50 for very cold waters with vertical polarization and diminishes to ~5 for equatorial regions.

It is indicated in {numref}`Figure37` and {numref}`Figure38` that one month of sampling with the CIMR polar-orbiting L-band instrument results in about 100 passes. This suggests that even during strong sunglint periods (when more than 25$\%$ of the data might be contaminated in the Northern Hemisphere), a monthly averaged SSS of 0.3 pss is achievable with the vertical polarization data under the assumptions stated above. An average of the SSS estimates from both polarizations should yield better accuracy.

```{figure} Figure37.png
---
name: Figure37
---
Total number of L-band L1B samples in 0.36°x 0.36° for 1 day (top panel), 1 week (middle panel) and 1 month (bottom panel).
```

```{figure} Figure38.png
---
name: Figure38
---
Median number of L-band L1B samples per 0.36°x0.36° as a function of latitude for 1 week (orange) and 1 month (blue curve).
```