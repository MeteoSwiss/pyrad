"""
Configuration file for the Python ARM Radar Toolkit (Py-ART)

The values for a number of Py-ART parameters and the default metadata created
when reading files, correcting fields, etc. is controlled by this single
Python configuration file.

Py-ART's configuration can be modified by setting the environment variable
PYART_CONFIG to point to a configuration file with formatting similar to this
file.

The recommended method for changing these defaults is for users to copy this
file into their home directory, rename it to pyart_config.py, make any
desired changes, and adjust their login scripts to set the PYART_CONFIG
environment variable to point to pyart_config.py in their home directory.

Py-ART's configuration can also be modified within a script or shell session
using the load_config functions, such modification will last until the end of
the script/session or until a new configuration is loaded.

"""

##############################################################################
##############################################################################
# Simple configuration
#
# Adjust the values of the variable (right hand side of the equal sign) in
# this section for an easy method of customizing Py-ART.  Do not change the
# variable names (the left hand side of the equal sign). More advanced
# settings are based upon these variables.  Most users will find that
# adjusting this section is all that is needed.
##############################################################################
##############################################################################

# The default fill value for masked arrays and _FillValue keys
fill_value = -9999.0

# The default equivalent earth radius factor used in atm. refraction
# computations
KE = 1.25

# Field names used when reading in radar files and in the various correction
# and retrieval algorithms. The comments in this section provide additional
# information about the fields in that section.

# Radar reflectivity fields
total_power = "total_power"

unfiltered_reflectivity = "unfiltered_reflectivity"
unfiltered_reflectivity_vv = "unfiltered_reflectivity_vv"

reflectivity = "reflectivity"
reflectivity_vv = "reflectivity_vv"

corrected_unfiltered_reflectivity = "corrected_unfiltered_reflectivity"

corrected_reflectivity = "corrected_reflectivity"
corrected_reflectivity_vv = "corrected_reflectivity_vv"

clutter_correction_ratio_hh = "clutter_correction_ratio_hh"
clutter_correction_ratio_vv = "clutter_correction_ratio_vv"

reflectivity_bias = "reflectivity_bias"

volumetric_reflectivity = "volumetric_reflectivity"
volumetric_reflectivity_vv = "volumetric_reflectivity_vv"

radar_cross_section_hh = "radar_cross_section_hh"
radar_cross_section_vv = "radar_cross_section_vv"

signal_power_hh = "signal_power_hh"
signal_power_vv = "signal_power_vv"

signal_to_noise_ratio = "signal_to_noise_ratio"
signal_to_noise_ratio_hh = "signal_to_noise_ratio_hh"
signal_to_noise_ratio_vv = "signal_to_noise_ratio_vv"

noisedBZ_hh = "noisedBZ_hh"
noisedBZ_vv = "noisedBZ_vv"

noisedBm_hh = "noisedBm_hh"
noisedBm_vv = "noisedBm_vv"

noisedBADU_hh = "noisedBADU_hh"
noisedBADU_vv = "noisedBADU_vv"

noiseADU_hh = "noiseADU_hh"
noiseADU_vv = "noiseADU_vv"

noise_pos_h = "noise_pos_h"
noise_pos_v = "noise_pos_v"

stat_test_lag1 = "stat_test_lag1"
stat_test_lag2 = "stat_test_lag2"
wide_band_noise = "wide_band_noise"

fields_difference = "fields_difference"
field_mask = "field_mask"
field_texture = "field_texture"

transmitted_signal_power_h = "transmitted_signal_power_h"
transmitted_signal_power_v = "transmitted_signal_power_v"

# spectral data
complex_spectra_hh_ADU = "complex_spectra_hh_ADU"
complex_spectra_vv_ADU = "complex_spectra_vv_ADU"

spectral_power_hh_ADU = "spectral_power_hh_ADU"
spectral_power_vv_ADU = "spectral_power_vv_ADU"

spectral_power_hh_dBADU = "spectral_power_hh_dBADU"
spectral_power_vv_dBADU = "spectral_power_vv_dBADU"

spectral_power_hh_dBm = "spectral_power_hh_dBm"
spectral_power_vv_dBm = "spectral_power_vv_dBm"

spectral_noise_power_hh_dBZ = "spectral_noise_power_hh_dBZ"
spectral_noise_power_vv_dBZ = "spectral_noise_power_vv_dBZ"
spectral_noise_power_hh_dBm = "spectral_noise_power_hh_dBm"
spectral_noise_power_vv_dBm = "spectral_noise_power_vv_dBm"
spectral_noise_power_hh_dBADU = "spectral_noise_power_hh_dBADU"
spectral_noise_power_vv_dBADU = "spectral_noise_power_vv_dBADU"
spectral_noise_power_hh_ADU = "spectral_noise_power_hh_ADU"
spectral_noise_power_vv_ADU = "spectral_noise_power_vv_ADU"

spectral_phase_hh = "spectral_phase_hh"
spectral_phase_vv = "spectral_phase_vv"

spectral_reflectivity_hh = "spectral_reflectivity_hh"
spectral_reflectivity_vv = "spectral_reflectivity_vv"
spectral_differential_reflectivity = "spectral_differential_reflectivity"
spectral_differential_phase = "spectral_differential_phase"
spectral_copolar_correlation_coefficient = "spectral_copolar_correlation_coefficient"

unfiltered_complex_spectra_hh_ADU = "unfiltered_complex_spectra_hh_ADU"
unfiltered_complex_spectra_vv_ADU = "unfiltered_complex_spectra_vv_ADU"

unfiltered_spectral_power_hh_ADU = "unfiltered_spectral_power_hh_ADU"
unfiltered_spectral_power_vv_ADU = "unfiltered_spectral_power_vv_ADU"

unfiltered_spectral_power_hh_dBADU = "unfiltered_spectral_power_hh_dBADU"
unfiltered_spectral_power_vv_dBADU = "unfiltered_spectral_power_vv_dBADU"

unfiltered_spectral_power_hh_dBm = "unfiltered_spectral_power_hh_dBm"
unfiltered_spectral_power_vv_dBm = "unfiltered_spectral_power_vv_dBm"

unfiltered_spectral_phase_hh = "unfiltered_spectral_phase_hh"
unfiltered_spectral_phase_vv = "unfiltered_spectral_phase_vv"

unfiltered_spectral_reflectivity_hh = "unfiltered_spectral_reflectivity_hh"
unfiltered_spectral_reflectivity_vv = "unfiltered_spectral_reflectivity_vv"
unfiltered_spectral_differential_reflectivity = (
    "unfiltered_spectral_differential_reflectivity"
)
unfiltered_spectral_differential_phase = "unfiltered_spectral_differential_phase"
unfiltered_spectral_copolar_correlation_coefficient = (
    "unfiltered_spectral_copolar_correlation_coefficient"
)

# IQ data
IQ_hh_ADU = "IQ_hh_ADU"
IQ_vv_ADU = "IQ_vv_ADU"
IQ_noise_power_hh_dBZ = "IQ_noise_power_hh_dBZ"
IQ_noise_power_vv_dBZ = "IQ_noise_power_vv_dBZ"
IQ_noise_power_hh_dBm = "IQ_noise_power_hh_dBm"
IQ_noise_power_vv_dBm = "IQ_noise_power_vv_dBm"
IQ_noise_power_hh_dBADU = "IQ_noise_power_hh_dBADU"
IQ_noise_power_vv_dBADU = "IQ_noise_power_vv_dBADU"
IQ_noise_power_hh_ADU = "IQ_noise_power_hh_ADU"
IQ_noise_power_vv_ADU = "IQ_noise_power_vv_ADU"

# Mean Doppler velocity fields, VEL
velocity = "velocity"
corrected_velocity = "corrected_velocity"
unfiltered_velocity = "unfiltered_velocity"
simulated_velocity = "simulated_velocity"

velocity_vv = "velocity_vv"
unfiltered_velocity_vv = "unfiltered_velocity_vv"

dealiased_corrected_velocity = "dealiased_corrected_velocity"
dealiased_velocity = "dealiased_velocity"

# retrieved Doppler velocity (for VAD)
retrieved_velocity = "retrieved_velocity"
retrieved_velocity_std = "retrieved_velocity_std"
velocity_difference = "velocity_difference"

# Spectral width fields, SW
spectrum_width = "spectrum_width"
unfiltered_spectrum_width = "unfiltered_spectrum_width"
corrected_spectrum_width = "corrected_spectrum_width"

unfiltered_spectrum_width_vv = "unfiltered_spectrum_width_vv"
spectrum_width_vv = "spectrum_width_vv"

turbulence = "turbulence"

# Differential reflectivity fields, ZDR
differential_reflectivity = "differential_reflectivity"
corrected_differential_reflectivity = "corrected_differential_reflectivity"

unfiltered_differential_reflectivity = "unfiltered_differential_reflectivity"

differential_reflectivity_in_precipitation = (
    "differential_reflectivity_in_precipitation"
)

differential_reflectivity_in_snow = "differential_reflectivity_in_snow"

differential_reflectivity_column_height = "differential_reflectivity_column_height"

# Cross correlation ratio, correlation coefficient, RhoHV
cross_correlation_ratio = "cross_correlation_ratio"

unfiltered_cross_correlation_ratio = "unfiltered_cross_correlation_ratio"
uncorrected_cross_correlation_ratio = "uncorrected_cross_correlation_ratio"
corrected_cross_correlation_ratio = "corrected_cross_correlation_ratio"
logarithmic_cross_correlation_ratio = "logarithmic_cross_correlation_ratio"

cross_correlation_ratio_in_rain = "cross_correlation_ratio_in_rain"

# Normalized coherent power, signal quality index, SQI, NCP
normalized_coherent_power = "normalized_coherent_power"

signal_quality_index = "signal_quality_index"
signal_quality_index_vv = "signal_quality_index_vv"
unfiltered_signal_quality_index = "unfiltered_signal_quality_index"
unfiltered_signal_quality_index_vv = "unfiltered_signal_quality_index_vv"

# visibility processing fields
terrain_altitude = "terrain_altitude"
bent_terrain_altitude = "bent_terrain_altitude"
terrain_slope = "terrain_slope"
terrain_aspect = "terrain_aspect"
elevation_angle = "elevation_angle"
visibility = "visibility"
min_vis_altitude = "min_vis_altitude"
min_vis_height_above_ground = "min_vis_height_above_ground"
min_rad_vis_height_above_ground = "min_rad_vis_height_above_ground"
min_vis_elevation = "min_vis_elevation"
incident_angle = "incident_angle"
effective_area = "effective_area"
sigma_0 = "sigma_0"
rcs_clutter = "rcs_clutter"
dBm_clutter = "dBm_clutter"
dBZ_clutter = "dBZ_clutter"
visibility_polar = "visibility_polar"


# Differential phase shift, PhiDP
mean_phase = "mean_phase"
differential_phase = "differential_phase"
unfolded_differential_phase = "unfolded_differential_phase"
corrected_differential_phase = "corrected_differential_phase"

uncorrected_differential_phase = "uncorrected_differential_phase"
uncorrected_unfiltered_differential_phase = "uncorrected_unfiltered_differential_phase"

system_differential_phase = "system_differential_phase"
first_gate_differential_phase = "first_gate_differential_phase"

# Specific differential phase shift, KDP
specific_differential_phase = "specific_differential_phase"
corrected_specific_differential_phase = "corrected_specific_differential_phase"

uncorrected_specific_differential_phase = "uncorrected_specific_differential_phase"
uncorrected_unfiltered_specific_differential_phase = (
    "uncorrected_unfiltered_specific_differential_phase"
)

# Linear depolarization ration (h - horizontal, v - vertical), LDR
linear_depolarization_ratio = "linear_polarization_ratio"
linear_depolarization_ratio_h = "linear_depolarization_ratio_h"
linear_depolarization_ratio_v = "linear_depolarization_ratio_v"

circular_depolarization_ratio = "circular_depolarization_ratio"

# attenuation
specific_attenuation = "specific_attenuation"
corrected_specific_attenuation = "corrected_specific_attenuation"
path_integrated_attenuation = "path_integrated_attenuation"
corrected_path_integrated_attenuation = "corrected_path_integrated_attenuation"
specific_differential_attenuation = "specific_differential_attenuation"
corrected_specific_differential_attenuation = (
    "corrected_specific_differential_attenuation"
)
path_integrated_differential_attenuation = "path_integrated_differential_attenuation"
corrected_path_integrated_differential_attenuation = (
    "corrected_path_integrated_differential_attenuation"
)

# Textures
differential_phase_texture = "differential_phase_texture"
cross_correlation_ratio_texture = "cross_correlation_ratio_texture"
differential_reflectivity_texture = "differential_reflectivity_texture"
reflectivity_texture = "reflectivity_texture"

# Misc fields
radar_echo_id = "radar_echo_id"
clutter_exit_code = "clutter_exit_code"
number_of_samples = "number_of_samples"
sum = "sum"
sum_squared = "sum_squared"
standard_deviation = "standard_deviation"
time_avg_flag = "time_avg_flag"
colocated_gates = "colocated_gates"
occurrence = "occurrence"
frequency_of_occurrence = "frequency_of_occurrence"

# icon data fields
temperature = "temperature"
iso0 = "iso0"
height_over_iso0 = "height_over_iso0"
icon_index = "icon_index"
hzt_index = "hzt_index"
iso0_height = "iso0_height"

# DEM fields
visibility = "visibility"
mininum_visible_elevation = "mininum_visible_elevation"

# precipitation
rain_rate = "rain_rate"
radar_estimated_rain_rate = "radar_estimated_rain_rate"
corrected_radar_estimated_rain_rate = "corrected_radar_estimated_rain_rate"
rainfall_accumulation = "rainfall_accumulation"

# melting layer
melting_layer = "melting_layer"
melting_layer_height = "melting_layer_height"

# hydroclass
radar_echo_classification = "radar_echo_classification"
corrected_radar_echo_classification = "corrected_radar_echo_classification"
hydroclass_entropy = "hydroclass_entropy"
proportion_AG = "proportion_AG"
proportion_CR = "proportion_CR"
proportion_LR = "proportion_LR"
proportion_RP = "proportion_RP"
proportion_RN = "proportion_RN"
proportion_VI = "proportion_VI"
proportion_WS = "proportion_WS"
proportion_MH = "proportion_MH"
proportion_IH = "proportion_IH"

# rad4alp products
probability_of_hail = "probability_of_hail"
thunderstorm_forecast = "thunderstorm_forecast"
maximum_expected_severe_hail_size = "maximum_expected_severe_hail_size"
maximum_echo = "maximum_echo"
maximum_echo_height = "maximum_echo_height"
echo_top_15dBZ = "echo_top_15dBZ"
echo_top_20dBZ = "echo_top_20dBZ"
echo_top_45dBZ = "echo_top_45dBZ"
echo_top_50dBZ = "echo_top_50dBZ"
vertically_integrated_liquid = "vertically_integrated_liquid"


# Wind retrieval fields
eastward_wind_component = "eastward_wind_component"
northward_wind_component = "northward_wind_component"
vertical_wind_component = "vertical_wind_component"
azimuthal_horizontal_wind_component = "azimuthal_horizontal_wind_component"
vertical_wind_shear = "vertical_wind_shear"
wind_speed = "wind_speed"
wind_direction = "wind_direction"

# Sun signal fields
sun_hit_power_h = "sun_hit_power_h"
sun_hit_power_v = "sun_hit_power_v"
sun_hit_differential_reflectivity = "sun_hit_differential_reflectivity"

sun_est_power_h = "sun_est_power_h"
sun_est_power_v = "sun_est_power_v"
sun_est_differential_reflectivity = "sun_est_differential_reflectivity"

sun_hit_h = "sun_hit_h"
sun_hit_v = "sun_hit_v"
sun_hit_zdr = "sun_hit_zdr"

# birds signal
bird_density = "bird_density"
bird_reflectivity = "bird_reflectivity"

# profile variables
height = "height"
interpolated_profile = "interpolated_profile"

# wind lidar variables
corrected_radial_wind_speed = "corrected_radial_wind_speed"
radial_wind_speed = "radial_wind_speed"
radial_wind_speed_ci = "radial_wind_speed_ci"
radial_wind_speed_status = "radial_wind_speed_status"

doppler_spectrum_width = "doppler_spectrum_width"
doppler_spectrum_mean_error = "doppler_spectrum_mean_error"

atmospherical_structures_type = "atmospherical_structures_type"
relative_beta = "relative_beta"
absolute_beta = "absolute_beta"

cnr = "cnr"

avg_reflectivity = "avg_reflectivity"
npoints_reflectivity = "npoints_reflectivity"
quant05_reflectivity = "quant05_reflectivity"
quant10_reflectivity = "quant10_reflectivity"
quant20_reflectivity = "quant20_reflectivity"
quant50_reflectivity = "quant50_reflectivity"
quant80_reflectivity = "quant80_reflectivity"
quant90_reflectivity = "quant90_reflectivity"
quant95_reflectivity = "quant95_reflectivity"

avg_radar_estimated_rain_rate = "avg_radar_estimated_rain_rate"
npoints_radar_estimated_rain_rate = "npoints_radar_estimated_rain_rate"
quant05_radar_estimated_rain_rate = "quant05_radar_estimated_rain_rate"
quant10_radar_estimated_rain_rate = "quant10_radar_estimated_rain_rate"
quant20_radar_estimated_rain_rate = "quant20_radar_estimated_rain_rate"
quant50_radar_estimated_rain_rate = "quant50_radar_estimated_rain_rate"
quant80_radar_estimated_rain_rate = "quant80_radar_estimated_rain_rate"
quant90_radar_estimated_rain_rate = "quant90_radar_estimated_rain_rate"
quant95_radar_estimated_rain_rate = "quant95_radar_estimated_rain_rate"

avg_velocity = "avg_velocity"
npoints_velocity = "npoints_velocity"
quant05_velocity = "quant05_velocity"
quant10_velocity = "quant10_velocity"
quant20_velocity = "quant20_velocity"
quant50_velocity = "quant50_velocity"
quant80_velocity = "quant80_velocity"
quant90_velocity = "quant90_velocity"
quant95_velocity = "quant95_velocity"

avg_corrected_velocity = "avg_corrected_velocity"
npoints_corrected_velocity = "npoints_corrected_velocity"
quant05_corrected_velocity = "quant05_corrected_velocity"
quant10_corrected_velocity = "quant10_corrected_velocity"
quant20_corrected_velocity = "quant20_corrected_velocity"
quant50_corrected_velocity = "quant50_corrected_velocity"
quant80_corrected_velocity = "quant80_corrected_velocity"
quant90_corrected_velocity = "quant90_corrected_velocity"
quant95_corrected_velocity = "quant95_corrected_velocity"


# Satellite products:
# all 10*um except CTH
# all 3 km resolution except HRV

# infra-red
IR_039 = "IR_039"
WV_062 = "WV_062"  # Water vapour channel
WV_073 = "WV_073"  # Water vapour channel
IR_087 = "IR_087"
IR_097 = "IR_097"
IR_108 = "IR_108"
IR_120 = "IR_120"
IR_134 = "IR_134"

# visible
HRV = "HRV"
VIS006 = "VIS006"
VIS008 = "VIS008"
IR_016 = "IR_016"  # This is near-infrared

# Cloud top height
CTH = "CTH"  # m

# visible normalized by sun elevation
HRV_norm = "HRV_norm"
VIS006_norm = "VIS006_norm"
VIS008_norm = "VIS008_norm"
IR_016_norm = "IR_016_norm"  # This is near-infrared


# End of Simple Configuration section

##############################################################################
##############################################################################
# Advanced Configuration
#
# Most users will not want to make any changes in this section.  For users
# who want a more fine-grained control over Py-ART's configuration this
# section provides access to these controls.  The layout of this section can
# be changed, the only requirement for a valid configuration file is that
# the ALL CAPITALIZED variable must must be present with the formatting
# present in this file.  These required variables are:
#
# FILL_VALUE, DEFAULT_METADATA, FILE_SPECIFIC_METADATA, FIELD_MAPPINGS,
# DEFAULT_FIELD_NAMES
#
# This section makes generous use of the variables in the Simple Configuration
# section, this is not required, but simplifies and enforces uniformity on
# the configuration.
##############################################################################
##############################################################################


##############################################################################
# Parameters
#
# Various parameters used in Py-ART.
##############################################################################

# the default fill value for masked arrays and the _FillValue key
FILL_VALUE = fill_value

# The DEFAULT_FIELD_NAMES controls the field names which are used in the
# correction and retrieval algorithms in Py-ART. The keys of the dictionary
# are "internal" names which cannot change, the values are the field names
# which will be used in the algorithms by default. For best results use the
# names defined by the variables in simple configuration section which are
# also used in the DEFAULT_METADATA and FIELD_MAPPINGS variable. If you
# choose to change a field name the names should also be changed in the
# DEFAULT_METADATA and FIELD_MAPPINGS variable. This is not required but
# highly suggested.

DEFAULT_FIELD_NAMES = {
    # Internal field name (do not change): field name used (can change)
    "reflectivity": reflectivity,
    "bird_reflectivity": bird_reflectivity,
    "corrected_reflectivity": corrected_reflectivity,
    "total_power": total_power,
    "unfiltered_reflectivity": unfiltered_reflectivity,
    "corrected_unfiltered_reflectivity": corrected_unfiltered_reflectivity,
    "reflectivity_vv": reflectivity_vv,
    "corrected_reflectivity_vv": corrected_reflectivity_vv,
    "unfiltered_reflectivity_vv": unfiltered_reflectivity_vv,
    "clutter_correction_ratio_hh": clutter_correction_ratio_hh,
    "clutter_correction_ratio_vv": clutter_correction_ratio_vv,
    "reflectivity_bias": reflectivity_bias,
    "signal_power_hh": signal_power_hh,
    "signal_power_vv": signal_power_vv,
    "volumetric_reflectivity": volumetric_reflectivity,
    "volumetric_reflectivity_vv": volumetric_reflectivity_vv,
    "radar_cross_section_hh": radar_cross_section_hh,
    "radar_cross_section_vv": radar_cross_section_vv,
    "sun_hit_power_h": sun_hit_power_h,
    "sun_hit_power_v": sun_hit_power_v,
    "sun_hit_differential_reflectivity": sun_hit_differential_reflectivity,
    "sun_est_power_h": sun_est_power_h,
    "sun_est_power_v": sun_est_power_v,
    "sun_est_differential_reflectivity": sun_est_differential_reflectivity,
    "velocity": velocity,
    "simulated_velocity": simulated_velocity,
    "corrected_velocity": corrected_velocity,
    "unfiltered_velocity": unfiltered_velocity,
    "velocity_vv": velocity_vv,
    "unfiltered_velocity_vv": unfiltered_velocity_vv,
    "dealiased_corrected_velocity": dealiased_corrected_velocity,
    "dealiased_velocity": dealiased_velocity,
    "retrieved_velocity": retrieved_velocity,
    "retrieved_velocity_std": retrieved_velocity_std,
    "velocity_difference": velocity_difference,
    "spectrum_width": spectrum_width,
    "corrected_spectrum_width": corrected_spectrum_width,
    "unfiltered_spectrum_width": unfiltered_spectrum_width,
    "spectrum_width_vv": spectrum_width_vv,
    "unfiltered_spectrum_width_vv": unfiltered_spectrum_width_vv,
    "turbulence": turbulence,
    "terrain_altitude": terrain_altitude,
    "bent_terrain_altitude": bent_terrain_altitude,
    "terrain_slope": terrain_slope,
    "terrain_aspect": terrain_aspect,
    "elevation_angle": elevation_angle,
    "visibility": visibility,
    "min_vis_altitude": min_vis_altitude,
    "min_vis_height_above_ground": min_vis_height_above_ground,
    "min_rad_vis_height_above_ground": min_rad_vis_height_above_ground,
    "min_vis_elevation": min_vis_elevation,
    "incident_angle": incident_angle,
    "effective_area": effective_area,
    "sigma_0": sigma_0,
    "rcs_clutter": rcs_clutter,
    "dBm_clutter": dBm_clutter,
    "dBZ_clutter": dBZ_clutter,
    "visibility_polar": visibility_polar,
    "differential_reflectivity": differential_reflectivity,
    "corrected_differential_reflectivity": corrected_differential_reflectivity,
    "unfiltered_differential_reflectivity": (unfiltered_differential_reflectivity),
    "differential_reflectivity_in_precipitation": (
        differential_reflectivity_in_precipitation
    ),
    "differential_reflectivity_in_snow": differential_reflectivity_in_snow,
    "differential_reflectivity_column_height": (
        differential_reflectivity_column_height
    ),
    "cross_correlation_ratio": cross_correlation_ratio,
    "corrected_cross_correlation_ratio": corrected_cross_correlation_ratio,
    "unfiltered_cross_correlation_ratio": unfiltered_cross_correlation_ratio,
    "uncorrected_cross_correlation_ratio": uncorrected_cross_correlation_ratio,
    "logarithmic_cross_correlation_ratio": logarithmic_cross_correlation_ratio,
    "cross_correlation_ratio_in_rain": cross_correlation_ratio_in_rain,
    "normalized_coherent_power": normalized_coherent_power,
    "wide_band_noise": wide_band_noise,
    "fields_difference": fields_difference,
    "field_mask": field_mask,
    "field_texture": field_texture,
    "stat_test_lag1": stat_test_lag1,
    "stat_test_lag2": stat_test_lag2,
    "mean_phase": mean_phase,
    "differential_phase": differential_phase,
    "unfolded_differential_phase": unfolded_differential_phase,
    "corrected_differential_phase": corrected_differential_phase,
    "uncorrected_differential_phase": uncorrected_differential_phase,
    "uncorrected_unfiltered_differential_phase": (
        uncorrected_unfiltered_differential_phase
    ),
    "system_differential_phase": system_differential_phase,
    "specific_differential_phase": specific_differential_phase,
    "corrected_specific_differential_phase": (corrected_specific_differential_phase),
    "uncorrected_specific_differential_phase": (
        uncorrected_specific_differential_phase
    ),
    "uncorrected_unfiltered_specific_differential_phase": (
        uncorrected_unfiltered_specific_differential_phase
    ),
    "linear_depolarization_ratio": linear_depolarization_ratio,
    "linear_depolarization_ratio_h": linear_depolarization_ratio_h,
    "linear_depolarization_ratio_v": linear_depolarization_ratio_v,
    "circular_depolarization_ratio": circular_depolarization_ratio,
    "signal_quality_index": signal_quality_index,
    "signal_quality_index_vv": signal_quality_index_vv,
    "unfiltered_signal_quality_index": unfiltered_signal_quality_index,
    "unfiltered_signal_quality_index_vv": unfiltered_signal_quality_index_vv,
    "signal_to_noise_ratio": signal_to_noise_ratio,
    "signal_to_noise_ratio_hh": signal_to_noise_ratio_hh,
    "signal_to_noise_ratio_vv": signal_to_noise_ratio_vv,
    "noisedBZ_hh": noisedBZ_hh,
    "noisedBZ_vv": noisedBZ_vv,
    "noisedBm_hh": noisedBm_hh,
    "noisedBm_vv": noisedBm_vv,
    "noisedBADU_hh": noisedBADU_hh,
    "noisedBADU_vv": noisedBADU_vv,
    "noiseADU_hh": noiseADU_hh,
    "noiseADU_vv": noiseADU_vv,
    "noise_pos_h": noise_pos_h,
    "noise_pos_v": noise_pos_v,
    "transmitted_signal_power_h": transmitted_signal_power_h,
    "transmitted_signal_power_v": transmitted_signal_power_v,
    "complex_spectra_hh_ADU": complex_spectra_hh_ADU,
    "complex_spectra_vv_ADU": complex_spectra_vv_ADU,
    "spectral_power_hh_ADU": spectral_power_hh_ADU,
    "spectral_power_vv_ADU": spectral_power_vv_ADU,
    "spectral_power_hh_dBADU": spectral_power_hh_dBADU,
    "spectral_power_vv_dBADU": spectral_power_vv_dBADU,
    "spectral_power_hh_dBm": spectral_power_hh_dBm,
    "spectral_power_vv_dBm": spectral_power_vv_dBm,
    "spectral_noise_power_hh_dBZ": spectral_noise_power_hh_dBZ,
    "spectral_noise_power_vv_dBZ": spectral_noise_power_vv_dBZ,
    "spectral_noise_power_hh_dBm": spectral_noise_power_hh_dBm,
    "spectral_noise_power_vv_dBm": spectral_noise_power_vv_dBm,
    "spectral_noise_power_hh_dBADU": spectral_noise_power_hh_dBADU,
    "spectral_noise_power_vv_dBADU": spectral_noise_power_vv_dBADU,
    "spectral_noise_power_hh_ADU": spectral_noise_power_hh_ADU,
    "spectral_noise_power_vv_ADU": spectral_noise_power_vv_ADU,
    "spectral_phase_hh": spectral_phase_hh,
    "spectral_phase_vv": spectral_phase_vv,
    "spectral_reflectivity_hh": spectral_reflectivity_hh,
    "spectral_reflectivity_vv": spectral_reflectivity_vv,
    "spectral_differential_reflectivity": spectral_differential_reflectivity,
    "spectral_differential_phase": spectral_differential_phase,
    "spectral_copolar_correlation_coefficient": (
        spectral_copolar_correlation_coefficient
    ),
    "unfiltered_complex_spectra_hh_ADU": unfiltered_complex_spectra_hh_ADU,
    "unfiltered_complex_spectra_vv_ADU": unfiltered_complex_spectra_vv_ADU,
    "unfiltered_spectral_power_hh_ADU": unfiltered_spectral_power_hh_ADU,
    "unfiltered_spectral_power_vv_ADU": unfiltered_spectral_power_vv_ADU,
    "unfiltered_spectral_power_hh_dBADU": unfiltered_spectral_power_hh_dBADU,
    "unfiltered_spectral_power_vv_dBADU": unfiltered_spectral_power_vv_dBADU,
    "unfiltered_spectral_power_hh_dBm": unfiltered_spectral_power_hh_dBm,
    "unfiltered_spectral_power_vv_dBm": unfiltered_spectral_power_vv_dBm,
    "unfiltered_spectral_phase_hh": unfiltered_spectral_phase_hh,
    "unfiltered_spectral_phase_vv": unfiltered_spectral_phase_vv,
    "unfiltered_spectral_reflectivity_hh": unfiltered_spectral_reflectivity_hh,
    "unfiltered_spectral_reflectivity_vv": unfiltered_spectral_reflectivity_vv,
    "unfiltered_spectral_differential_reflectivity": (
        unfiltered_spectral_differential_reflectivity
    ),
    "unfiltered_spectral_differential_phase": (unfiltered_spectral_differential_phase),
    "unfiltered_spectral_copolar_correlation_coefficient": (
        unfiltered_spectral_copolar_correlation_coefficient
    ),
    "IQ_hh_ADU": IQ_hh_ADU,
    "IQ_vv_ADU": IQ_vv_ADU,
    "IQ_noise_power_hh_dBZ": IQ_noise_power_hh_dBZ,
    "IQ_noise_power_vv_dBZ": IQ_noise_power_vv_dBZ,
    "IQ_noise_power_hh_dBm": IQ_noise_power_hh_dBm,
    "IQ_noise_power_vv_dBm": IQ_noise_power_vv_dBm,
    "IQ_noise_power_hh_dBADU": IQ_noise_power_hh_dBADU,
    "IQ_noise_power_vv_dBADU": IQ_noise_power_vv_dBADU,
    "IQ_noise_power_hh_ADU": IQ_noise_power_hh_ADU,
    "IQ_noise_power_vv_ADU": IQ_noise_power_vv_ADU,
    "rain_rate": rain_rate,
    "bird_density": bird_density,
    "sun_hit_h": sun_hit_h,
    "sun_hit_v": sun_hit_v,
    "sun_hit_zdr": sun_hit_zdr,
    "radar_estimated_rain_rate": radar_estimated_rain_rate,
    "corrected_radar_estimated_rain_rate": corrected_radar_estimated_rain_rate,
    "rainfall_accumulation": rainfall_accumulation,
    "radar_echo_classification": radar_echo_classification,
    "corrected_radar_echo_classification": corrected_radar_echo_classification,
    "hydroclass_entropy": hydroclass_entropy,
    "proportion_AG": proportion_AG,
    "proportion_CR": proportion_CR,
    "proportion_LR": proportion_LR,
    "proportion_RP": proportion_RP,
    "proportion_RN": proportion_RN,
    "proportion_VI": proportion_VI,
    "proportion_WS": proportion_WS,
    "proportion_MH": proportion_MH,
    "proportion_IH": proportion_IH,
    "radar_echo_id": radar_echo_id,
    "clutter_exit_code": clutter_exit_code,
    "melting_layer": melting_layer,
    "melting_layer_height": melting_layer_height,
    "probability_of_hail": probability_of_hail,
    "thunderstorm_forecast": thunderstorm_forecast,
    "maximum_expected_severe_hail_size": maximum_expected_severe_hail_size,
    "maximum_echo": maximum_echo,
    "maximum_echo_height": maximum_echo_height,
    "echo_top_15dBZ": echo_top_15dBZ,
    "echo_top_20dBZ": echo_top_20dBZ,
    "echo_top_45dBZ": echo_top_45dBZ,
    "echo_top_50dBZ": echo_top_50dBZ,
    "vertically_integrated_liquid": vertically_integrated_liquid,
    "specific_attenuation": specific_attenuation,
    "path_integrated_attenuation": path_integrated_attenuation,
    "specific_differential_attenuation": specific_differential_attenuation,
    "path_integrated_differential_attenuation": (
        path_integrated_differential_attenuation
    ),
    "corrected_specific_attenuation": corrected_specific_attenuation,
    "corrected_path_integrated_attenuation": (corrected_path_integrated_attenuation),
    "corrected_specific_differential_attenuation": (
        corrected_specific_differential_attenuation
    ),
    "corrected_path_integrated_differential_attenuation": (
        corrected_path_integrated_differential_attenuation
    ),
    "temperature": temperature,
    "iso0": iso0,
    "height_over_iso0": height_over_iso0,
    "iso0_height": iso0_height,
    "icon_index": icon_index,
    "hzt_index": hzt_index,
    "mininum_visible_elevation": mininum_visible_elevation,
    "differential_phase_texture": differential_phase_texture,
    "cross_correlation_ratio_texture": cross_correlation_ratio_texture,
    "differential_reflectivity_texture": differential_reflectivity_texture,
    "reflectivity_texture": reflectivity_texture,
    "eastward_wind_component": eastward_wind_component,
    "northward_wind_component": northward_wind_component,
    "vertical_wind_component": vertical_wind_component,
    "azimuthal_horizontal_wind_component": azimuthal_horizontal_wind_component,
    "vertical_wind_shear": vertical_wind_shear,
    "wind_speed": wind_speed,
    "wind_direction": wind_direction,
    "height": height,
    "number_of_samples": number_of_samples,
    "standard_deviation": standard_deviation,
    "sum": sum,
    "sum_squared": sum_squared,
    "colocated_gates": colocated_gates,
    "time_avg_flag": time_avg_flag,
    "occurrence": occurrence,
    "frequency_of_occurrence": frequency_of_occurrence,
    "interpolated_profile": interpolated_profile,
    "corrected_radial_wind_speed": corrected_radial_wind_speed,
    "radial_wind_speed": radial_wind_speed,
    "radial_wind_speed_ci": radial_wind_speed_ci,
    "radial_wind_speed_status": radial_wind_speed_status,
    "doppler_spectrum_width": doppler_spectrum_width,
    "doppler_spectrum_mean_error": doppler_spectrum_mean_error,
    "atmospherical_structures_type": atmospherical_structures_type,
    "relative_beta": relative_beta,
    "absolute_beta": absolute_beta,
    "cnr": cnr,
    "avg_reflectivity": avg_reflectivity,
    "npoints_reflectivity": npoints_reflectivity,
    "quant05_reflectivity": quant05_reflectivity,
    "quant10_reflectivity": quant10_reflectivity,
    "quant20_reflectivity": quant20_reflectivity,
    "quant50_reflectivity": quant50_reflectivity,
    "quant80_reflectivity": quant80_reflectivity,
    "quant90_reflectivity": quant90_reflectivity,
    "quant95_reflectivity": quant95_reflectivity,
    "avg_radar_estimated_rain_rate": avg_radar_estimated_rain_rate,
    "npoints_radar_estimated_rain_rate": npoints_radar_estimated_rain_rate,
    "quant05_radar_estimated_rain_rate": quant05_radar_estimated_rain_rate,
    "quant10_radar_estimated_rain_rate": quant10_radar_estimated_rain_rate,
    "quant20_radar_estimated_rain_rate": quant20_radar_estimated_rain_rate,
    "quant50_radar_estimated_rain_rate": quant50_radar_estimated_rain_rate,
    "quant80_radar_estimated_rain_rate": quant80_radar_estimated_rain_rate,
    "quant90_radar_estimated_rain_rate": quant90_radar_estimated_rain_rate,
    "quant95_radar_estimated_rain_rate": quant95_radar_estimated_rain_rate,
    "avg_velocity": avg_velocity,
    "npoints_velocity": npoints_velocity,
    "quant05_velocity": quant05_velocity,
    "quant10_velocity": quant10_velocity,
    "quant20_velocity": quant20_velocity,
    "quant50_velocity": quant50_velocity,
    "quant80_velocity": quant80_velocity,
    "quant90_velocity": quant90_velocity,
    "quant95_velocity": quant95_velocity,
    "avg_corrected_velocity": avg_corrected_velocity,
    "npoints_corrected_velocity": npoints_corrected_velocity,
    "quant05_corrected_velocity": quant05_corrected_velocity,
    "quant10_corrected_velocity": quant10_corrected_velocity,
    "quant20_corrected_velocity": quant20_corrected_velocity,
    "quant50_corrected_velocity": quant50_corrected_velocity,
    "quant80_corrected_velocity": quant80_corrected_velocity,
    "quant90_corrected_velocity": quant90_corrected_velocity,
    "quant95_corrected_velocity": quant95_corrected_velocity,
    "CTH": CTH,
    "HRV": HRV,
    "VIS006": VIS006,
    "VIS008": VIS008,
    "IR_016": IR_016,
    "IR_039": IR_039,
    "WV_062": WV_062,
    "WV_073": WV_073,
    "IR_087": IR_087,
    "IR_097": IR_097,
    "IR_108": IR_108,
    "IR_120": IR_120,
    "IR_134": IR_134,
    "HRV_norm": HRV_norm,
    "VIS006_norm": VIS006_norm,
    "VIS008_norm": VIS008_norm,
    "IR_016_norm": IR_016_norm,
}


##############################################################################
# Default metadata
#
# The DEFAULT_METADATA dictionary contains dictionaries which provide the
# default radar attribute and field metadata. When reading in a file with
# Py-ART the FILE_SPECIFIC_METADATA variable is first queued for a metadata
# dictionary, if it is not found then the metadata in DEFAULT_METADATA is
# utilized.
##############################################################################

DEFAULT_METADATA = {
    # Metadata for radar attributes. These closely follow the CF/Radial
    # standard
    "azimuth": {
        "units": "degrees",
        "standard_name": "beam_azimuth_angle",
        "long_name": "azimuth_angle_from_true_north",
        "axis": "radial_azimuth_coordinate",
        "comment": "Azimuth of antenna relative to true north",
    },
    "elevation": {
        "units": "degrees",
        "standard_name": "beam_elevation_angle",
        "long_name": "elevation_angle_from_horizontal_plane",
        "axis": "radial_elevation_coordinate",
        "comment": "Elevation of antenna relative to the horizontal plane",
    },
    "number_of_pulses": {
        "units": "-",
        "standard_name": "number_of_pulses",
        "long_name": "number of pulses per ray",
        "axis": "radial_pulses_coordinate",
    },
    "scan_rate": {
        "units": "degrees_per_second",
        "long_name": "Antenna angle scan rate",
    },
    "range": {
        "units": "meters",
        "standard_name": "projection_range_coordinate",
        "long_name": "range_to_measurement_volume",
        "axis": "radial_range_coordinate",
        "spacing_is_constant": "true",
        "comment": ("Coordinate variable for range. Range to center of each bin."),
    },
    "time": {
        "units": "seconds",
        "standard_name": "time",
        "long_name": "time_in_seconds_since_volume_start",
        "calendar": "gregorian",
        "comment": (
            "Coordinate variable for time. "
            "Time at the center of each ray, in fractional seconds "
            "since the global variable time_coverage_start"
        ),
    },
    "metadata": {
        "Conventions": "CF/Radial instrument_parameters",
        "version": "1.3",
        "title": "",
        "institution": "",
        "references": "",
        "source": "",
        "history": "",
        "comment": "",
        "instrument_name": "",
    },
    # Metadata for radar sweep information dictionaries
    "sweep_number": {
        "units": "count",
        "standard_name": "sweep_number",
        "long_name": "Sweep number",
    },
    "sweep_mode": {
        "units": "unitless",
        "standard_name": "sweep_mode",
        "long_name": "Sweep mode",
        "comment": (
            'Options are: "sector", "coplane", "rhi", '
            '"vertical_pointing", "idle", "azimuth_surveillance", '
            '"elevation_surveillance", "sunscan", "pointing", '
            '"manual_ppi", "manual_rhi"'
        ),
    },
    "fixed_angle": {
        "long_name": "Target angle for sweep",
        "units": "degrees",
        "standard_name": "target_fixed_angle",
    },
    "sweep_start_ray_index": {
        "long_name": "Index of first ray in sweep, 0-based",
        "units": "count",
    },
    "sweep_end_ray_index": {
        "long_name": "Index of last ray in sweep, 0-based",
        "units": "count",
    },
    "rays_per_sweep": {"long_name": "Number of rays in each sweep", "units": "count"},
    "target_scan_rate": {
        "long_name": "Target scan rate for sweep",
        "units": "degrees_per_second",
    },
    "rays_are_indexed": {
        "long_name": "Flag for indexed rays",
        "units": "unitless",
        "options": (
            "true: rays are indexed to a regular grid, "
            + "false: rays are not indexed to a regular grid"
        ),
    },
    "ray_angle_res": {
        "long_name": "Angular resolution between rays",
        "units": "degrees",
        "comment": "Only applicable when rays_are_indexed variable is true",
    },
    # Metadata for radar location attributes
    "latitude": {
        "long_name": "Latitude",
        "standard_name": "Latitude",
        "units": "degrees_north",
    },
    "longitude": {
        "long_name": "Longitude",
        "standard_name": "Longitude",
        "units": "degrees_east",
    },
    "altitude": {
        "long_name": "Altitude",
        "standard_name": "Altitude",
        "units": "meters",
        "positive": "up",
    },
    "gate_x": {
        "long_name": "Cartesian x location of gate with origin at the radar",
        "units": "meters",
    },
    "gate_y": {
        "long_name": "Cartesian y location of gate with origin at the radar",
        "units": "meters",
    },
    "gate_z": {
        "long_name": "Cartesian z location of gate with origin at the radar",
        "units": "meters",
    },
    "gate_edge_x": {
        "long_name": "Cartesian x location of the edges of each gate",
        "units": "meters",
    },
    "gate_edge_y": {
        "long_name": "Cartesian y location of the edges of each gate",
        "units": "meters",
    },
    "gate_edge_z": {
        "long_name": "Cartesian z location of the edges of each gate",
        "units": "meters",
    },
    "gate_longitude": {
        "long_name": "Longitude of radar gate.",
        "units": "degrees_north",
    },
    "gate_latitude": {"long_name": "Latitude of radar gate", "units": "degrees_east"},
    "gate_altitude": {"long_name": "Altitude of radar gate", "units": "meters"},
    # Metadata for instrument_parameter dictionary
    "prt_mode": {
        "comments": (
            'Pulsing mode Options are: "fixed", "staggered", '
            '"dual". Assumed "fixed" if missing.'
        ),
        "meta_group": "instrument_parameters",
        "long_name": "Pulsing mode",
        "units": "unitless",
    },
    "nyquist_velocity": {
        "units": "m/s",
        "comments": "Unambiguous velocity",
        "meta_group": "instrument_parameters",
        "long_name": "Nyquist velocity",
    },
    "prt": {
        "units": "seconds",
        "comments": (
            "Pulse repetition time. For staggered prt, " "also see prt_ratio."
        ),
        "meta_group": "instrument_parameters",
        "long_name": "Pulse repetition time",
    },
    "unambiguous_range": {
        "units": "meters",
        "comments": "Unambiguous range",
        "meta_group": "instrument_parameters",
        "long_name": "Unambiguous range",
    },
    "pulse_width": {
        "units": "seconds",
        "comments": "Pulse width",
        "meta_group": "instrument_parameters",
        "long_name": "Pulse width",
    },
    "prt_ratio": {
        "units": "unitless",
        "meta_group": "instrument_parameters",
        "long_name": "Pulse repetition frequency ratio",
    },
    "frequency": {
        "units": "s-1",
        "meta_group": "instrument_parameters",
        "long_name": "Radiation frequency",
    },
    "n_samples": {
        "units": "unitless",
        "meta_group": "instrument_parameters",
        "long_name": "Number of samples used to compute moments",
    },
    "radar_antenna_gain_h": {
        "units": "dB",
        "meta_group": "instrument_parameters",
        "long_name": "Antenna gain H polarization",
    },
    "radar_antenna_gain_v": {
        "units": "dB",
        "meta_group": "instrument_parameters",
        "long_name": "Antenna gain V polarization",
    },
    # metadata for radar calibration constant
    "calibration_constant_hh": {
        "units": "dB",
        "meta_group": "radar_calibration",
        "long_name": " radar calibration constant H polarization",
    },
    "calibration_constant_vv": {
        "units": "dB",
        "meta_group": "radar_calibration",
        "long_name": " radar calibration constant V polarization",
    },
    "transmit_power_h": {
        "units": "dBm",
        "meta_group": "radar_calibration",
        "long_name": " transmit power H channel",
    },
    "transmit_power_v": {
        "units": "dBm",
        "meta_group": "radar_calibration",
        "long_name": " transmit power V channel",
    },
    "dBADU_to_dBm_hh": {
        "units": "dBm",
        "meta_group": "radar_calibration",
        "long_name": "dBADU to dBm H polarization",
    },
    "dBADU_to_dBm_vv": {
        "units": "dBm",
        "meta_group": "radar_calibration",
        "long_name": "dBADU to dBm V polarization",
    },
    "matched_filter_loss_h": {
        "units": "dB",
        "meta_group": "radar_calibration",
        "long_name": "matched filter loss H polarization",
    },
    "matched_filter_loss_v": {
        "units": "dB",
        "meta_group": "radar_calibration",
        "long_name": "matched filter loss V polarization",
    },
    "path_attenuation": {
        "units": "dB/km",
        "meta_group": "radar_calibration",
        "long_name": "matched filter loss",
    },
    # non-standard parameter for specifying the PRF high/low for each ray
    "prf_flag": {
        "units": "unitless",
        "comments": "PRF used to collect ray. 0 for high PRF, 1 for low PRF.",
        "meta_group": "instrument_parameters",
        "long_name": "PRF flag",
    },
    # Metadata for radar_parameter sub-convention
    "radar_beam_width_h": {
        "units": "degrees",
        "meta_group": "radar_parameters",
        "long_name": "Antenna beam width H polarization",
    },
    "radar_beam_width_v": {
        "units": "degrees",
        "meta_group": "radar_parameters",
        "long_name": "Antenna beam width V polarization",
    },
    # Metadata for airborne radar parameters
    "rotation": {
        "units": "degrees",
        "standard_name": "ray_rotation_angle_relative_to_platform",
        "long_name": "Ray rotation angle relative to platform",
    },
    "tilt": {
        "units": "degrees",
        "standard_name": "ray_tilt_angle_relative_to_platform",
        "long_name": "Ray tilt angle relative to platform",
    },
    "roll": {
        "units": "degrees",
        "standard_name": "platform_roll_angle",
        "long_name": "Platform roll angle",
    },
    "drift": {
        "units": "degrees",
        "standard_name": "platform_drift_angle",
        "long_name": "Platform drift angle",
    },
    "heading": {
        "units": "degrees",
        "standard_name": "platform_heading_angle",
        "long_name": "Platform heading angle",
    },
    "pitch": {
        "units": "degrees",
        "standard_name": "platform_pitch_angle",
        "long_name": "Platform pitch angle",
    },
    "georefs_applied": {
        "units": "unitless",
        "standard_name": "georefs_have_been_applied_to_ray",
        "long_name": "Geoferences have been applied to ray",
        "comment": "1 if georefs have been applied, 0 otherwise",
    },
    # Metadata for Doppler spectra
    "Doppler_frequency": {
        "units": "Hz",
        "standard_name": "Doppler_frequency",
        "long_name": "Doppler frequency",
    },
    "Doppler_velocity": {
        "units": "m/s",
        "standard_name": "Doppler_velocity",
        "long_name": "Doppler velocity",
    },
    # Reflectivity fields
    reflectivity: {
        "units": "dBZ",
        "standard_name": "horizontal_reflectivity",
        "long_name": "Horizontal Reflectivity",
        "coordinates": "elevation azimuth range",
        "scale_factor": 0.5,
        "add_offset": -32.0,
    },
    avg_reflectivity: {
        "units": "dBZ",
        "standard_name": "horizontal_reflectivity",
        "long_name": "Horizontal Reflectivity",
        "coordinates": "elevation azimuth range",
        "scale_factor": 0.5,
        "add_offset": -32.0,
        "_Write_as_dtype": "uint8",
    },
    quant05_reflectivity: {
        "units": "dBZ",
        "standard_name": "horizontal_reflectivity",
        "long_name": "Horizontal Reflectivity",
        "coordinates": "elevation azimuth range",
        "scale_factor": 0.5,
        "add_offset": -32.0,
        "_Write_as_dtype": "uint8",
    },
    quant10_reflectivity: {
        "units": "dBZ",
        "standard_name": "horizontal_reflectivity",
        "long_name": "Horizontal Reflectivity",
        "coordinates": "elevation azimuth range",
        "scale_factor": 0.5,
        "add_offset": -32.0,
        "_Write_as_dtype": "uint8",
    },
    quant20_reflectivity: {
        "units": "dBZ",
        "standard_name": "horizontal_reflectivity",
        "long_name": "Horizontal Reflectivity",
        "coordinates": "elevation azimuth range",
        "scale_factor": 0.5,
        "add_offset": -32.0,
        "_Write_as_dtype": "uint8",
    },
    quant50_reflectivity: {
        "units": "dBZ",
        "standard_name": "horizontal_reflectivity",
        "long_name": "Horizontal Reflectivity",
        "coordinates": "elevation azimuth range",
        "scale_factor": 0.5,
        "add_offset": -32.0,
        "_Write_as_dtype": "uint8",
    },
    quant80_reflectivity: {
        "units": "dBZ",
        "standard_name": "horizontal_reflectivity",
        "long_name": "Horizontal Reflectivity",
        "coordinates": "elevation azimuth range",
        "scale_factor": 0.5,
        "add_offset": -32.0,
        "_Write_as_dtype": "uint8",
    },
    quant90_reflectivity: {
        "units": "dBZ",
        "standard_name": "horizontal_reflectivity",
        "long_name": "Horizontal Reflectivity",
        "coordinates": "elevation azimuth range",
        "scale_factor": 0.5,
        "add_offset": -32.0,
        "_Write_as_dtype": "uint8",
    },
    quant95_reflectivity: {
        "units": "dBZ",
        "standard_name": "horizontal_reflectivity",
        "long_name": "Horizontal Reflectivity",
        "coordinates": "elevation azimuth range",
        "scale_factor": 0.5,
        "add_offset": -32.0,
        "_Write_as_dtype": "uint8",
    },
    bird_reflectivity: {
        "units": "dBZ",
        "standard_name": "bird_reflectivity",
        "long_name": "Bird Reflectivity",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_reflectivity: {
        "units": "dBZ",
        "standard_name": "horizontal_reflectivity",
        "long_name": "Unfiltered Horizontal Reflectivity",
        "coordinates": "elevation azimuth range",
    },
    corrected_reflectivity: {
        "units": "dBZ",
        "standard_name": "horizontal_reflectivity",
        "long_name": "Corrected Horizontal Reflectivity",
        "coordinates": "elevation azimuth range",
    },
    reflectivity_vv: {
        "units": "dBZ",
        "standard_name": "vertical_reflectivity",
        "long_name": "Vertical Reflectivity",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_reflectivity_vv: {
        "units": "dBZ",
        "standard_name": "vertical_reflectivity",
        "long_name": "Unfiltered Vertical Reflectivity",
        "coordinates": "elevation azimuth range",
    },
    corrected_reflectivity_vv: {
        "units": "dBZ",
        "standard_name": "vertical_reflectivity",
        "long_name": "Corrected Vertical Reflectivity",
        "coordinates": "elevation azimuth range",
    },
    clutter_correction_ratio_hh: {
        "units": "dB",
        "standard_name": "clutter_correction_ratio_hh",
        "long_name": "Horizontal Clutter Correction Ratio",
        "valid_min": 0.0,
        "coordinates": "elevation azimuth range",
    },
    clutter_correction_ratio_vv: {
        "units": "dB",
        "standard_name": "clutter_correction_ratio_vv",
        "long_name": "Vertical Clutter Correction Ratio",
        "valid_min": 0.0,
        "coordinates": "elevation azimuth range",
    },
    volumetric_reflectivity: {
        "units": "10log10(cm2/km3)",
        "standard_name": "volumetric_reflectivity",
        "long_name": "Volumetric Reflectivity",
        "coordinates": "elevation azimuth range",
    },
    volumetric_reflectivity_vv: {
        "units": "10log10(cm2/km3)",
        "standard_name": "volumetric_reflectivity_vv",
        "long_name": "Vertical Volumetric Reflectivity",
        "coordinates": "elevation azimuth range",
    },
    radar_cross_section_hh: {
        "units": "dBsm",
        "standard_name": "radar_cross_section_hh",
        "long_name": "Horizontal Radar Cross-Section",
        "coordinates": "elevation azimuth range",
    },
    radar_cross_section_vv: {
        "units": "dBsm",
        "standard_name": "radar_cross_section_vv",
        "long_name": "Vertical Radar Cross-Section",
        "coordinates": "elevation azimuth range",
    },
    total_power: {
        "units": "dBZ",
        "standard_name": "equivalent_reflectivity_factor",
        "long_name": "Total power",
        "coordinates": "elevation azimuth range",
    },
    reflectivity_bias: {
        "units": "dB",
        "standard_name": "reflectivity_bias",
        "long_name": "Reflectivity bias",
        "coordinates": "elevation azimuth range",
    },
    signal_power_hh: {
        "units": "dBm",
        "standard_name": "signal_power_hh",
        "long_name": "Signal power horizontal",
        "coordinates": "elevation azimuth range",
    },
    signal_power_vv: {
        "units": "dBm",
        "standard_name": "signal_power_vv",
        "long_name": "Signal power vertical",
        "coordinates": "elevation azimuth range",
    },
    sun_hit_power_h: {
        "units": "dBm",
        "standard_name": "sun_hit_power_h",
        "long_name": "sun hit power horizontal",
        "coordinates": "elevation azimuth",
    },
    sun_hit_power_v: {
        "units": "dBm",
        "standard_name": "sun_hit_power_v",
        "long_name": "sun hit power vertical",
        "coordinates": "elevation azimuth",
    },
    sun_hit_differential_reflectivity: {
        "units": "dB",
        "standard_name": "sun_hit_differential_reflectivity",
        "long_name": "sun hit differential reflectivity",
        "coordinates": "elevation azimuth",
    },
    sun_est_power_h: {
        "units": "dBm",
        "standard_name": "sun_est_power_h",
        "long_name": "estimated sun power horizontal",
        "coordinates": "elevation azimuth",
    },
    sun_est_power_v: {
        "units": "dBm",
        "standard_name": "sun_est_power_v",
        "long_name": "estimated sun power vertical",
        "coordinates": "elevation azimuth",
    },
    sun_est_differential_reflectivity: {
        "units": "dB",
        "standard_name": "sun_est_differential_reflectivity",
        "long_name": "estimated sun differential reflectivity",
        "coordinates": "elevation azimuth",
    },
    stat_test_lag1: {
        "units": "dB",
        "standard_name": "stat_test_lag1",
        "long_name": "Statistical test one lag fluctuation",
        "coordinates": "elevation azimuth range",
    },
    stat_test_lag2: {
        "units": "dB",
        "standard_name": "stat_test_lag2",
        "long_name": "Statistical test two lag fluctuation",
        "coordinates": "elevation azimuth range",
    },
    wide_band_noise: {
        "units": "dB",
        "standard_name": "wide_band_noise",
        "long_name": "Wide band noise",
        "coordinates": "elevation azimuth range",
    },
    fields_difference: {
        "units": "-",
        "standard_name": "fields_difference",
        "long_name": "Fields difference",
        "coordinates": "elevation azimuth range",
    },
    field_texture: {
        "units": "-",
        "standard_name": "field_texture",
        "long_name": "Field texture",
        "coordinates": "elevation azimuth range",
    },
    field_mask: {
        "units": "-",
        "standard_name": "field_mask",
        "long_name": "Field mask",
        "labels": ["BELOW", "ABOVE"],
        "ticks": [1, 2],
        "boundaries": [-0.5, 1.5, 2.5],
        "coordinates": "elevation azimuth range",
    },
    mean_phase: {
        "units": "deg",
        "standard_name": "mean_phase",
        "long_name": "Mean phase",
        "coordinates": "elevation azimuth range",
    },
    # Velocity fields
    velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    simulated_velocity: {
        "units": "meters_per_second",
        "standard_name": ("radial_velocity_of_scatterers_away_from_instrument"),
        "long_name": "Simulated mean doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    corrected_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Corrected mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Unfiltered mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    avg_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    quant05_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    quant10_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    quant20_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    quant50_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    quant80_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    quant90_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    quant95_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    avg_corrected_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    quant05_corrected_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    quant10_corrected_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    quant20_corrected_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    quant50_corrected_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    quant80_corrected_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    quant90_corrected_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    quant95_corrected_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    dealiased_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Dealiased mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    dealiased_corrected_velocity: {
        "units": "m/s",
        "standard_name": "mean_Doppler_velocity",
        "long_name": "Dealiased corrected mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    retrieved_velocity: {
        "units": "m/s",
        "standard_name": "retrieved_mean_Doppler_velocity",
        "long_name": "Retrieved mean Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    retrieved_velocity_std: {
        "units": "m/s",
        "standard_name": "retrieved_mean_Doppler_velocity_std",
        "long_name": "Retrieved mean Doppler velocity standard deviation",
        "coordinates": "elevation azimuth range",
    },
    velocity_difference: {
        "units": "m/s",
        "standard_name": "retrieved_mean_Doppler_velocity_difference",
        "long_name": "Difference between retrieved and measured Doppler velocity",
        "coordinates": "elevation azimuth range",
    },
    # Spectrum width fields
    spectrum_width: {
        "units": "m/s",
        "standard_name": "Doppler_spectrum_width",
        "long_name": "Doppler spectrum width",
        "coordinates": "elevation azimuth range",
    },
    corrected_spectrum_width: {
        "units": "m/s",
        "standard_name": "Doppler_spectrum_width",
        "long_name": "Corrected Doppler spectrum width",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_spectrum_width: {
        "units": "m/s",
        "standard_name": "Doppler_spectrum_width",
        "long_name": "Unfiltered Doppler spectrum width",
        "coordinates": "elevation azimuth range",
    },
    turbulence: {
        "units": "m^2/3 s^-1",
        "standard_name": "EDR^1/3",
        "long_name": "Cubic Root of Eddy Dissipation Rate",
        "coordinates": "elevation azimuth range",
    },
    # Visibility processing fields
    terrain_altitude: {
        "units": "meters",
        "standard_name": "terrain_altitude",
        "long_name": "Altitude",
        "coordinates": "x y",
    },
    bent_terrain_altitude: {
        "units": "meters",
        "standard_name": "bent_terrain_altitude",
        "long_name": "Altitude after refraction correction",
        "coordinates": "x y",
    },
    terrain_slope: {
        "units": "degrees",
        "standard_name": "terrain_slope",
        "long_name": "Slope of the terrain",
        "coordinates": "x y",
    },
    terrain_aspect: {
        "units": "degrees",
        "standard_name": "terrain_aspect",
        "long_name": "Orientation of the terrain",
        "coordinates": "x y",
    },
    elevation_angle: {
        "units": "degrees",
        "standard_name": "elevation_angle",
        "long_name": "Elevation angle at terrain [deg]",
        "coordinates": "x y",
    },
    visibility: {
        "units": "-",
        "standard_name": "visibility",
        "long_name": "Radar visibility",
        "coordinates": "x y",
    },
    min_vis_altitude: {
        "units": "meters",
        "standard_name": "min_vis_altitude",
        "long_name": "Minimum visible altitude",
        "coordinates": "x y",
    },
    min_vis_height_above_ground: {
        "units": "meters",
        "standard_name": "min_vis_height_above_ground",
        "long_name": "Minimum visible height above ground",
        "coordinates": "x y",
    },
    min_rad_vis_height_above_ground: {
        "units": "meters",
        "standard_name": "min_rad_vis_height_above_ground",
        "long_name": "Minimum radar visible height above ground",
        "coordinates": "x y",
    },
    min_vis_elevation: {
        "units": "degrees",
        "standard_name": "min_vis_elevation",
        "long_name": "Minimum visible elevation angle",
        "coordinates": "x y",
    },
    incident_angle: {
        "units": "degrees",
        "standard_name": "incident_angle",
        "long_name": "Angle of incidence at terrain",
        "coordinates": "x y",
    },
    effective_area: {
        "units": "m",
        "standard_name": "effective_area",
        "long_name": "Effective backscattering area",
        "coordinates": "x y",
    },
    sigma_0: {
        "units": "-",
        "standard_name": "sigma_0",
        "long_name": "ratio RCS/backscattering area",
        "coordinates": "x y",
    },
    rcs_clutter: {
        "units": "square meters",
        "standard_name": "rcs_clutter",
        "long_name": "RCS of the ground clutter",
        "coordinates": "elevation azimuth elevation",
    },
    dBm_clutter: {
        "units": "dBm",
        "standard_name": "dBm_clutter",
        "long_name": "Ground clutter power signal",
        "coordinates": "elevation azimuth elevation",
    },
    dBZ_clutter: {
        "units": "dBZ",
        "standard_name": "dBZ_clutter",
        "long_name": "Ground clutter reflectivity",
        "coordinates": "elevation azimuth elevation",
    },
    visibility_polar: {
        "units": "%",
        "standard_name": "visibility_polar",
        "long_name": "Visibility",
        "coordinates": "elevation azimuth range",
    },
    # Dual-polarization fields
    differential_reflectivity: {
        "units": "dB",
        "standard_name": "differential_reflectivity",
        "long_name": "Differential reflectivity",
        "coordinates": "elevation azimuth range",
        "scale_factor": 0.062254902,
        "add_offset": -7.9375,
    },
    corrected_differential_reflectivity: {
        "units": "dB",
        "standard_name": "differential_reflectivity",
        "long_name": "Corrected differential reflectivity",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_differential_reflectivity: {
        "units": "dB",
        "standard_name": "differential_reflectivity",
        "long_name": "Unfiltered differential reflectivity",
        "coordinates": "elevation azimuth range",
    },
    differential_reflectivity_in_precipitation: {
        "units": "dB",
        "standard_name": "differential_reflectivity_in_precip",
        "long_name": "Differential reflectivity in precipitation",
        "coordinates": "elevation azimuth range",
    },
    differential_reflectivity_in_snow: {
        "units": "dB",
        "standard_name": "differential_reflectivity_in_snow",
        "long_name": "Differential reflectivity in snow",
        "coordinates": "elevation azimuth range",
    },
    differential_reflectivity_column_height: {
        "units": "Km above freezing level",
        "standard_name": "differential_reflectivity_column_height",
        "long_name": "Differential reflectivity column height above freezing level",
        "coordinates": "elevation azimuth range",
    },
    cross_correlation_ratio: {
        "units": "-",
        "standard_name": "copolar_correlation_coefficient",
        "long_name": "Copolar correlation coefficient (RHOHV)",
        "coordinates": "elevation azimuth range",
    },
    corrected_cross_correlation_ratio: {
        "units": "-",
        "standard_name": "copolar_correlation_coefficient",
        "long_name": "Corrected copolar correlation coefficient (RHOHV)",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_cross_correlation_ratio: {
        "units": "-",
        "standard_name": "copolar_correlation_coefficient",
        "long_name": "Unfiltered copolar correlation coefficient (RHOHV)",
        "coordinates": "elevation azimuth range",
    },
    uncorrected_cross_correlation_ratio: {
        "units": "-",
        "standard_name": "copolar_correlation_coefficient",
        "long_name": "Uncorrected copolar correlation coefficient (RHOHV)",
        "coordinates": "elevation azimuth range",
    },
    logarithmic_cross_correlation_ratio: {
        "units": "dB",
        "standard_name": "logarithmic_copolar_correlation_coefficient",
        "long_name": "Logarithmic copolar correlation coefficient (L)",
        "coordinates": "elevation azimuth range",
    },
    cross_correlation_ratio_in_rain: {
        "units": "-",
        "standard_name": "copolar_correlation_coefficient_in_rain",
        "long_name": "copolar correlation coefficient in rain",
        "coordinates": "elevation azimuth range",
    },
    normalized_coherent_power: {
        "units": "-",
        "standard_name": "normalized_coherent_power",
        "long_name": "Normalized coherent power",
        "valid_max": 1.0,
        "valid_min": 0.0,
        "comment": "Also know as signal quality index (SQI)",
        "coordinates": "elevation azimuth range",
    },
    differential_phase: {
        "units": "deg",
        "standard_name": "differential_phase",
        "long_name": "Differential propagation phase (PhiDP)",
        "valid_max": 180.0,
        "valid_min": -180.0,
        "coordinates": "elevation azimuth range",
    },
    unfolded_differential_phase: {
        "units": "deg",
        "standard_name": "differential_phase",
        "long_name": "Unfolded differential propagation phase",
        "coordinates": "elevation azimuth range",
    },
    corrected_differential_phase: {
        "units": "deg",
        "standard_name": "differential_phase",
        "long_name": "Corrected differential propagation phase",
        "coordinates": "elevation azimuth range",
    },
    uncorrected_differential_phase: {
        "units": "deg",
        "standard_name": "differential_phase",
        "long_name": "Uncorrected differential propagation phase",
        "coordinates": "elevation azimuth range",
    },
    uncorrected_unfiltered_differential_phase: {
        "units": "deg",
        "standard_name": "differential_phase",
        "long_name": "Uncorrected unfiltered differential propagation phase",
        "coordinates": "elevation azimuth range",
    },
    system_differential_phase: {
        "units": "deg",
        "standard_name": "system_differential_phase",
        "long_name": "System differential phase (PhiDP0)",
        "coordinates": "elevation azimuth range",
    },
    first_gate_differential_phase: {
        "units": "gate index",
        "standard_name": "first_gate_differential_phase",
        "long_name": "First valid differential phase gate",
        "coordinates": "elevation azimuth",
    },
    specific_differential_phase: {
        "units": "deg/km",
        "standard_name": "specific_differential_phase",
        "long_name": "Specific differential phase (KDP)",
        "coordinates": "elevation azimuth range",
    },
    corrected_specific_differential_phase: {
        "units": "deg/km",
        "standard_name": "specific_differential_phase",
        "long_name": "Corrected specific differential phase (KDP)",
        "coordinates": "elevation azimuth range",
    },
    # Depolarization ratio fields
    linear_depolarization_ratio: {
        "units": "dB",
        "standard_name": "linear_depolarization_ratio",
        "long_name": "Linear depolarization ratio",
        "coordinates": "elevation azimuth range",
    },
    linear_depolarization_ratio_h: {
        "units": "dB",
        "standard_name": "linear_depolarization_ratio_h",
        "long_name": "Linear depolarization ratio horizontal",
        "coordinates": "elevation azimuth range",
    },
    linear_depolarization_ratio_v: {
        "units": "dB",
        "standard_name": "linear_depolarization_ratio_v",
        "long_name": "Linear depolarization ratio vertical",
        "coordinates": "elevation azimuth range",
    },
    circular_depolarization_ratio: {
        "units": "dB",
        "standard_name": "circular_depolarization_ratio",
        "long_name": "Circular depolarization ratio",
        "coordinates": "elevation azimuth range",
    },
    # Misc fields
    signal_to_noise_ratio: {
        "units": "dB",
        "standard_name": "signal_to_noise_ratio",
        "long_name": "Signal to noise ratio",
        "coordinates": "elevation azimuth range",
    },
    signal_to_noise_ratio_hh: {
        "units": "dB",
        "standard_name": "signal_to_noise_ratio_hh",
        "long_name": "Signal to noise ratio horizontal",
        "coordinates": "elevation azimuth range",
    },
    signal_to_noise_ratio_vv: {
        "units": "dB",
        "standard_name": "signal_to_noise_ratio_vv",
        "long_name": "Signal to noise ratio vertical",
        "coordinates": "elevation azimuth range",
    },
    noisedBZ_hh: {
        "units": "dBZ",
        "standard_name": "noisedBZ_hh",
        "long_name": "noise in dBZ horizontal",
        "coordinates": "elevation azimuth range",
    },
    noisedBZ_vv: {
        "units": "dBZ",
        "standard_name": "noisedBZ_vv",
        "long_name": "noise in dBZ vertical",
        "coordinates": "elevation azimuth range",
    },
    noisedBm_hh: {
        "units": "dBm",
        "standard_name": "noisedBm_hh",
        "long_name": "noise in dBm horizontal",
        "coordinates": "elevation azimuth range",
    },
    noisedBm_vv: {
        "units": "dBm",
        "standard_name": "noisedBm_vv",
        "long_name": "noise in dBm vertical",
        "coordinates": "elevation azimuth range",
    },
    noisedBADU_hh: {
        "units": "dBADU",
        "standard_name": "noisedBADU_hh",
        "long_name": "noise in dBADU horizontal",
        "coordinates": "elevation azimuth range",
    },
    noisedBADU_vv: {
        "units": "dBADU",
        "standard_name": "noisedBADU_vv",
        "long_name": "noise in dBADU vertical",
        "coordinates": "elevation azimuth range",
    },
    noiseADU_hh: {
        "units": "ADU",
        "standard_name": "noiseADU_hh",
        "long_name": "noise in ADU horizontal",
        "coordinates": "elevation azimuth range",
    },
    noiseADU_vv: {
        "units": "dBADU",
        "standard_name": "noiseADU_vv",
        "long_name": "noise in ADU vertical",
        "coordinates": "elevation azimuth range",
    },
    noise_pos_h: {
        "units": "-",
        "standard_name": "noise_pos_h",
        "long_name": "noisy radar bins horizontal polarization",
        "labels": ["OTHER", "NOISE"],
        "ticks": [1, 2],
        "boundaries": [0.5, 1.5, 2.5],
        "coordinates": "elevation azimuth range",
        "scale_factor": 1,
        "add_offset": 0,
        "_FillValue": 0,
        "_Write_as_dtype": "uint8",
    },
    noise_pos_v: {
        "units": "-",
        "standard_name": "noise_pos_v",
        "long_name": "noisy radar bins vertical polarization",
        "labels": ["OTHER", "NOISE"],
        "ticks": [1, 2],
        "boundaries": [0.5, 1.5, 2.5],
        "coordinates": "elevation azimuth range",
        "scale_factor": 1,
        "add_offset": 0,
        "_FillValue": 0,
        "_Write_as_dtype": "uint8",
    },
    transmitted_signal_power_h: {
        "units": "kW",
        "standard_name": "transmitted_signal_power_h",
        "long_name": "Transmitted signal power horizontal",
        "coordinates": "elevation azimuth range",
    },
    transmitted_signal_power_v: {
        "units": "kW",
        "standard_name": "transmitted_signal_power_v",
        "long_name": "Transmitted signal power vertical",
        "coordinates": "elevation azimuth range",
    },
    complex_spectra_hh_ADU: {
        "units": "ADU",
        "standard_name": "complex_spectra_hh_ADU",
        "long_name": "Complex spectra horizontal",
        "coordinates": "elevation azimuth range",
    },
    complex_spectra_vv_ADU: {
        "units": "ADU",
        "standard_name": "complex_spectra_vv_ADU",
        "long_name": "Complex spectra vertical",
        "coordinates": "elevation azimuth range",
    },
    spectral_power_hh_ADU: {
        "units": "ADU",
        "standard_name": "spectral_power_hh_ADU",
        "long_name": "spectral power horizontal",
        "coordinates": "elevation azimuth range",
    },
    spectral_power_vv_ADU: {
        "units": "ADU",
        "standard_name": "spectral_power_vv_ADU",
        "long_name": "spectral power vertical",
        "coordinates": "elevation azimuth range",
    },
    spectral_power_hh_dBADU: {
        "units": "dBADU",
        "standard_name": "spectral_power_hh_dBADU",
        "long_name": "spectral power horizontal",
        "coordinates": "elevation azimuth range",
    },
    spectral_power_vv_dBADU: {
        "units": "dBADU",
        "standard_name": "spectral_power_vv_dBADU",
        "long_name": "spectral power vertical",
        "coordinates": "elevation azimuth range",
    },
    spectral_power_hh_dBm: {
        "units": "dBm",
        "standard_name": "spectral_power_hh_dBm",
        "long_name": "spectral power horizontal",
        "coordinates": "elevation azimuth range",
    },
    spectral_power_vv_dBm: {
        "units": "dBm",
        "standard_name": "spectral_power_vv_dBm",
        "long_name": "spectral power vertical",
        "coordinates": "elevation azimuth range",
    },
    spectral_noise_power_hh_dBZ: {
        "units": "dBZ",
        "standard_name": "spectral_noise_power_hh_dBZ",
        "long_name": "spectral noise power horizontal",
        "coordinates": "elevation azimuth range",
    },
    spectral_noise_power_vv_dBZ: {
        "units": "dBZ",
        "standard_name": "spectral_noise_power_vv_dBZ",
        "long_name": "spectral noise power vertical",
        "coordinates": "elevation azimuth range",
    },
    spectral_noise_power_hh_dBm: {
        "units": "dBm",
        "standard_name": "spectral_noise_power_hh_dBm",
        "long_name": "spectral noise power horizontal",
        "coordinates": "elevation azimuth range",
    },
    spectral_noise_power_vv_dBm: {
        "units": "dBm",
        "standard_name": "spectral_noise_power_hh_dBm",
        "long_name": "spectral noise power vertical",
        "coordinates": "elevation azimuth range",
    },
    spectral_noise_power_hh_dBADU: {
        "units": "dBADU",
        "standard_name": "spectral_noise_power_hh_dBADU",
        "long_name": "spectral noise power horizontal",
        "coordinates": "elevation azimuth range",
    },
    spectral_noise_power_vv_dBADU: {
        "units": "dBADU",
        "standard_name": "spectral_noise_power_hh_dBADU",
        "long_name": "spectral noise power vertical",
        "coordinates": "elevation azimuth range",
    },
    spectral_noise_power_hh_ADU: {
        "units": "ADU",
        "standard_name": "spectral_noise_power_hh_ADU",
        "long_name": "spectral noise power horizontal",
        "coordinates": "elevation azimuth range",
    },
    spectral_noise_power_vv_ADU: {
        "units": "ADU",
        "standard_name": "spectral_noise_power_hh_ADU",
        "long_name": "spectral noise power vertical",
        "coordinates": "elevation azimuth range",
    },
    spectral_phase_hh: {
        "units": "deg",
        "standard_name": "spectral_phase_hh",
        "long_name": "spectral phase horizontal",
        "coordinates": "elevation azimuth range",
    },
    spectral_phase_vv: {
        "units": "deg",
        "standard_name": "spectral_phase_vv",
        "long_name": "spectral phase vertical",
        "coordinates": "elevation azimuth range",
    },
    spectral_reflectivity_hh: {
        "units": "dBZ",
        "standard_name": "spectral_reflectivity_hh",
        "long_name": "Spectral Horizontal Reflectivity",
        "coordinates": "elevation azimuth range",
    },
    spectral_reflectivity_vv: {
        "units": "dBZ",
        "standard_name": "spectral_reflectivity_vv",
        "long_name": "Spectral Vertical Reflectivity",
        "coordinates": "elevation azimuth range",
    },
    spectral_differential_reflectivity: {
        "units": "dBZ",
        "standard_name": "spectral_differential_reflectivity",
        "long_name": "Spectral Differential Reflectivity",
        "coordinates": "elevation azimuth range",
    },
    spectral_differential_phase: {
        "units": "deg",
        "standard_name": "spectral_differential_phase",
        "long_name": "Spectral Differential Phase",
        "coordinates": "elevation azimuth range",
    },
    spectral_copolar_correlation_coefficient: {
        "units": "-",
        "standard_name": "spectral_copolar_correlation_coefficient",
        "long_name": "Spectral copolar correlation coefficient (RHOHV)",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_complex_spectra_hh_ADU: {
        "units": "ADU",
        "standard_name": "unfiltered_complex_spectra_hh_ADU",
        "long_name": "Unfiltered complex spectra horizontal",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_complex_spectra_vv_ADU: {
        "units": "ADU",
        "standard_name": "unfiltered_complex_spectra_vv_ADU",
        "long_name": "Unfiltered complex spectra vertical",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_spectral_power_hh_ADU: {
        "units": "ADU",
        "standard_name": "unfiltered_spectral_power_hh_ADU",
        "long_name": "Unfiltered spectral power horizontal",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_spectral_power_vv_ADU: {
        "units": "ADU",
        "standard_name": "unfiltered_spectral_power_vv_ADU",
        "long_name": "Unfiltered spectral power vertical",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_spectral_power_hh_dBADU: {
        "units": "dBADU",
        "standard_name": "unfiltered_spectral_power_hh_dBADU",
        "long_name": "Unfiltered spectral power horizontal",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_spectral_power_vv_dBADU: {
        "units": "dBADU",
        "standard_name": "unfiltered_spectral_power_vv_dBADU",
        "long_name": "Unfiltered spectral power vertical",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_spectral_power_hh_dBm: {
        "units": "dBm",
        "standard_name": "unfiltered_spectral_power_hh_dBm",
        "long_name": "Unfiltered spectral power horizontal",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_spectral_power_vv_dBm: {
        "units": "dBm",
        "standard_name": "unfiltered_spectral_power_vv_dBm",
        "long_name": "Unfiltered spectral power vertical",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_spectral_phase_hh: {
        "units": "deg",
        "standard_name": "unfiltered_spectral_phase_hh",
        "long_name": "Unfiltered spectral phase horizontal",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_spectral_phase_vv: {
        "units": "deg",
        "standard_name": "unfiltered_spectral_phase_vv",
        "long_name": "Unfiltered spectral phase vertical",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_spectral_reflectivity_hh: {
        "units": "dBZ",
        "standard_name": "unfiltered_spectral_reflectivity_hh",
        "long_name": "Unfiltered Spectral Horizontal Reflectivity",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_spectral_reflectivity_vv: {
        "units": "dBZ",
        "standard_name": "unfiltered_spectral_reflectivity_vv",
        "long_name": "Unfiltered Spectral Vertical Reflectivity",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_spectral_differential_reflectivity: {
        "units": "dBZ",
        "standard_name": "unfiltered_spectral_differential_reflectivity",
        "long_name": "Unfiltered Spectral Differential Reflectivity",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_spectral_differential_phase: {
        "units": "deg",
        "standard_name": "unfiltered_spectral_differential_phase",
        "long_name": "Unfiltered Spectral Differential Phase",
        "coordinates": "elevation azimuth range",
    },
    unfiltered_spectral_copolar_correlation_coefficient: {
        "units": "-",
        "standard_name": "unfiltered_spectral_copolar_correlation_coefficient",
        "long_name": "Unfiltered Spectral copolar correlation coefficient (RHOHV)",
        "coordinates": "elevation azimuth range",
    },
    IQ_hh_ADU: {
        "units": "ADU",
        "standard_name": "IQ_hh_ADU",
        "long_name": "IQ signal horizontal",
        "coordinates": "elevation azimuth range",
    },
    IQ_vv_ADU: {
        "units": "ADU",
        "standard_name": "IQ_vv_ADU",
        "long_name": "IQ signal vertical",
        "coordinates": "elevation azimuth range",
    },
    IQ_noise_power_hh_dBZ: {
        "units": "dBZ",
        "standard_name": "IQ_noise_power_hh_dBZ",
        "long_name": "IQ noise power horizontal",
        "coordinates": "elevation azimuth range",
    },
    IQ_noise_power_vv_dBZ: {
        "units": "dBZ",
        "standard_name": "IQ_noise_power_vv_dBZ",
        "long_name": "IQ noise power vertical",
        "coordinates": "elevation azimuth range",
    },
    IQ_noise_power_hh_dBm: {
        "units": "dBm",
        "standard_name": "IQ_noise_power_hh_dBm",
        "long_name": "IQ noise power horizontal",
        "coordinates": "elevation azimuth range",
    },
    IQ_noise_power_vv_dBm: {
        "units": "dBm",
        "standard_name": "IQ_noise_power_vv_dBm",
        "long_name": "IQ noise power vertical",
        "coordinates": "elevation azimuth range",
    },
    IQ_noise_power_hh_dBADU: {
        "units": "dBADU",
        "standard_name": "IQ_noise_power_hh_dBADU",
        "long_name": "IQ noise power horizontal",
        "coordinates": "elevation azimuth range",
    },
    IQ_noise_power_vv_dBADU: {
        "units": "dBADU",
        "standard_name": "IQ_noise_power_vv_dBADU",
        "long_name": "IQ noise power vertical",
        "coordinates": "elevation azimuth range",
    },
    IQ_noise_power_hh_ADU: {
        "units": "ADU",
        "standard_name": "IQ_noise_power_hh_ADU",
        "long_name": "IQ noise power horizontal",
        "coordinates": "elevation azimuth range",
    },
    IQ_noise_power_vv_ADU: {
        "units": "ADU",
        "standard_name": "IQ_noise_power_vv_ADU",
        "long_name": "IQ noise power vertical",
        "coordinates": "elevation azimuth range",
    },
    rain_rate: {
        "units": "mm/h",
        "standard_name": "rain_rate",
        "long_name": "Rain rate",
        "coordinates": "elevation azimuth range",
    },
    bird_density: {
        "units": "birds/km3",
        "standard_name": "bird_density",
        "long_name": "Birds Density",
        "coordinates": "elevation azimuth range",
    },
    radar_estimated_rain_rate: {
        "units": "mm/h",
        "standard_name": "radar_estimated_rain_rate",
        "long_name": "Radar estimated rain rate",
        "labels": [
            "0.",
            "0.4",
            "0.63",
            "1.",
            "1.6",
            "2.5",
            "4.0",
            "6.3",
            "10.",
            "16.",
            "25.",
            "40.",
            "63.",
            "100.",
            "160.",
            "250.",
        ],
        "ticks": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
        ],
        "boundaries": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
            500.0,
        ],
        "coordinates": "elevation azimuth range",
    },
    corrected_radar_estimated_rain_rate: {
        "units": "mm/h",
        "standard_name": "corrected_radar_estimated_rain_rate",
        "long_name": "Radar estimated rain rate",
        "labels": [
            "0.",
            "0.4",
            "0.63",
            "1.",
            "1.6",
            "2.5",
            "4.0",
            "6.3",
            "10.",
            "16.",
            "25.",
            "40.",
            "63.",
            "100.",
            "160.",
            "250.",
        ],
        "ticks": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
        ],
        "boundaries": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
            500.0,
        ],
        "coordinates": "elevation azimuth range",
    },
    avg_radar_estimated_rain_rate: {
        "units": "mm/h",
        "standard_name": "radar_estimated_rain_rate",
        "long_name": "Radar estimated rain rate",
        "labels": [
            "0.",
            "0.4",
            "0.63",
            "1.",
            "1.6",
            "2.5",
            "4.0",
            "6.3",
            "10.",
            "16.",
            "25.",
            "40.",
            "63.",
            "100.",
            "160.",
            "250.",
        ],
        "ticks": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
        ],
        "boundaries": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
            500.0,
        ],
        "coordinates": "elevation azimuth range",
    },
    quant05_radar_estimated_rain_rate: {
        "units": "mm/h",
        "standard_name": "radar_estimated_rain_rate",
        "long_name": "Radar estimated rain rate",
        "labels": [
            "0.",
            "0.4",
            "0.63",
            "1.",
            "1.6",
            "2.5",
            "4.0",
            "6.3",
            "10.",
            "16.",
            "25.",
            "40.",
            "63.",
            "100.",
            "160.",
            "250.",
        ],
        "ticks": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
        ],
        "boundaries": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
            500.0,
        ],
        "coordinates": "elevation azimuth range",
    },
    quant10_radar_estimated_rain_rate: {
        "units": "mm/h",
        "standard_name": "radar_estimated_rain_rate",
        "long_name": "Radar estimated rain rate",
        "labels": [
            "0.",
            "0.4",
            "0.63",
            "1.",
            "1.6",
            "2.5",
            "4.0",
            "6.3",
            "10.",
            "16.",
            "25.",
            "40.",
            "63.",
            "100.",
            "160.",
            "250.",
        ],
        "ticks": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
        ],
        "boundaries": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
            500.0,
        ],
        "coordinates": "elevation azimuth range",
    },
    quant20_radar_estimated_rain_rate: {
        "units": "mm/h",
        "standard_name": "radar_estimated_rain_rate",
        "long_name": "Radar estimated rain rate",
        "labels": [
            "0.",
            "0.4",
            "0.63",
            "1.",
            "1.6",
            "2.5",
            "4.0",
            "6.3",
            "10.",
            "16.",
            "25.",
            "40.",
            "63.",
            "100.",
            "160.",
            "250.",
        ],
        "ticks": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
        ],
        "boundaries": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
            500.0,
        ],
        "coordinates": "elevation azimuth range",
    },
    quant50_radar_estimated_rain_rate: {
        "units": "mm/h",
        "standard_name": "radar_estimated_rain_rate",
        "long_name": "Radar estimated rain rate",
        "labels": [
            "0.",
            "0.4",
            "0.63",
            "1.",
            "1.6",
            "2.5",
            "4.0",
            "6.3",
            "10.",
            "16.",
            "25.",
            "40.",
            "63.",
            "100.",
            "160.",
            "250.",
        ],
        "ticks": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
        ],
        "boundaries": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
            500.0,
        ],
        "coordinates": "elevation azimuth range",
    },
    quant80_radar_estimated_rain_rate: {
        "units": "mm/h",
        "standard_name": "radar_estimated_rain_rate",
        "long_name": "Radar estimated rain rate",
        "labels": [
            "0.",
            "0.4",
            "0.63",
            "1.",
            "1.6",
            "2.5",
            "4.0",
            "6.3",
            "10.",
            "16.",
            "25.",
            "40.",
            "63.",
            "100.",
            "160.",
            "250.",
        ],
        "ticks": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
        ],
        "boundaries": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
            500.0,
        ],
        "coordinates": "elevation azimuth range",
    },
    quant90_radar_estimated_rain_rate: {
        "units": "mm/h",
        "standard_name": "radar_estimated_rain_rate",
        "long_name": "Radar estimated rain rate",
        "labels": [
            "0.",
            "0.4",
            "0.63",
            "1.",
            "1.6",
            "2.5",
            "4.0",
            "6.3",
            "10.",
            "16.",
            "25.",
            "40.",
            "63.",
            "100.",
            "160.",
            "250.",
        ],
        "ticks": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
        ],
        "boundaries": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
            500.0,
        ],
        "coordinates": "elevation azimuth range",
    },
    quant95_radar_estimated_rain_rate: {
        "units": "mm/h",
        "standard_name": "radar_estimated_rain_rate",
        "long_name": "Radar estimated rain rate",
        "labels": [
            "0.",
            "0.4",
            "0.63",
            "1.",
            "1.6",
            "2.5",
            "4.0",
            "6.3",
            "10.",
            "16.",
            "25.",
            "40.",
            "63.",
            "100.",
            "160.",
            "250.",
        ],
        "ticks": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
        ],
        "boundaries": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
            500.0,
        ],
        "coordinates": "elevation azimuth range",
    },
    rainfall_accumulation: {
        "units": "mm",
        "standard_name": "rainfall_accumulation",
        "long_name": "Rainfall accumulation",
        "labels": [
            "0.",
            "0.4",
            "0.63",
            "1.",
            "1.6",
            "2.5",
            "4.0",
            "6.3",
            "10.",
            "16.",
            "25.",
            "40.",
            "63.",
            "100.",
            "160.",
            "250.",
        ],
        "ticks": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
        ],
        "boundaries": [
            0.0,
            0.4,
            0.63,
            1.0,
            1.6,
            2.5,
            4.0,
            6.3,
            10.0,
            16.0,
            25.0,
            40.0,
            63.0,
            100.0,
            160.0,
            250.0,
            500.0,
        ],
        "coordinates": "elevation azimuth range",
    },
    sun_hit_h: {
        "units": "-",
        "standard_name": "sun_hit_h",
        "long_name": "sun hit radar bins horizontal polarization",
        "labels": ["OTHER", "SUN"],
        "ticks": [1, 2],
        "boundaries": [0.5, 1.5, 2.5],
        "coordinates": "elevation azimuth range",
        "scale_factor": 1,
        "add_offset": 0,
        "_FillValue": 0,
        "_Write_as_dtype": "uint8",
    },
    sun_hit_v: {
        "units": "-",
        "standard_name": "sun_hit_v",
        "long_name": "sun hit radar bins vertical polarization",
        "labels": ["OTHER", "SUN"],
        "ticks": [1, 2],
        "boundaries": [0.5, 1.5, 2.5],
        "coordinates": "elevation azimuth range",
        "scale_factor": 1,
        "add_offset": 0,
        "_FillValue": 0,
        "_Write_as_dtype": "uint8",
    },
    sun_hit_zdr: {
        "units": "-",
        "standard_name": "sun_hit_zdr",
        "long_name": "sun hit radar bins differential reflectivity",
        "labels": ["OTHER", "SUN"],
        "ticks": [1, 2],
        "boundaries": [0.5, 1.5, 2.5],
        "coordinates": "elevation azimuth range",
        "scale_factor": 1,
        "add_offset": 0,
        "_FillValue": 0,
        "_Write_as_dtype": "uint8",
    },
    radar_echo_classification: {
        "units": "-",
        "standard_name": "radar_echo_classification",
        "long_name": "Radar echo classification",
        "labels": ["NC", "AG", "CR", "LR", "RP", "RN", "VI", "WS", "MH", "IH/HDG"],
        "ticks": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        "boundaries": [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5],
        "coordinates": "elevation azimuth range",
        "scale_factor": 1,
        "add_offset": 0,
        "_FillValue": 0,
        "_Write_as_dtype": "uint8",
    },
    corrected_radar_echo_classification: {
        "units": "-",
        "standard_name": "corrected_radar_echo_classification",
        "long_name": "Radar echo classification",
        "labels": ["NC", "AG", "CR", "LR", "RP", "RN", "VI", "WS", "MH", "IH/HDG"],
        "ticks": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        "boundaries": [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5],
        "coordinates": "elevation azimuth range",
        "scale_factor": 1,
        "add_offset": 0,
        "_FillValue": 0,
        "_Write_as_dtype": "uint8",
    },
    hydroclass_entropy: {
        "units": "-",
        "standard_name": "hydroclass_entropy",
        "long_name": "Semi-supervised hydrometeor classification entropy",
        "coordinates": "elevation azimuth range",
    },
    proportion_AG: {
        "units": "percent",
        "standard_name": "proportion_AG",
        "long_name": "Aggregates proportion",
        "coordinates": "elevation azimuth range",
    },
    proportion_CR: {
        "units": "percent",
        "standard_name": "proportion_CR",
        "long_name": "Crystals proportion",
        "coordinates": "elevation azimuth range",
    },
    proportion_LR: {
        "units": "percent",
        "standard_name": "proportion_LR",
        "long_name": "Light Rain proportion",
        "coordinates": "elevation azimuth range",
    },
    proportion_RP: {
        "units": "percent",
        "standard_name": "proportion_RP",
        "long_name": "Rimed particles proportion",
        "coordinates": "elevation azimuth range",
    },
    proportion_RN: {
        "units": "percent",
        "standard_name": "proportion_RN",
        "long_name": "Rain proportion",
        "coordinates": "elevation azimuth range",
    },
    proportion_VI: {
        "units": "percent",
        "standard_name": "proportion_VI",
        "long_name": "Vertical Ice Crystals proportion",
        "coordinates": "elevation azimuth range",
    },
    proportion_WS: {
        "units": "percent",
        "standard_name": "proportion_WS",
        "long_name": "Wet Snow proportion",
        "coordinates": "elevation azimuth range",
    },
    proportion_MH: {
        "units": "percent",
        "standard_name": "proportion_MH",
        "long_name": "Melting Hail proportion",
        "coordinates": "elevation azimuth range",
    },
    proportion_IH: {
        "units": "percent",
        "standard_name": "proportion_IH",
        "long_name": "Ice Hail proportion",
        "coordinates": "elevation azimuth range",
    },
    radar_echo_id: {
        "units": "-",
        "standard_name": "radar_echo_id",
        "long_name": "Radar Echo Identification",
        "labels": ["NOISE", "CLT", "PREC"],
        "ticks": [1, 2, 3],
        "boundaries": [0.5, 1.5, 2.5, 3.5],
        "coordinates": "elevation azimuth range",
        "scale_factor": 1,
        "add_offset": 0,
        "_FillValue": 0,
        "_Write_as_dtype": "uint8",
    },
    clutter_exit_code: {
        "units": "-",
        "standard_name": "clutter_exit_code",
        "long_name": "Clutter Exit Code",
        "labels": ["NOISE", "PREC", "CLT"],
        "ticks": [1, 50, 150],
        "boundaries": [0.5, 1.5, 99.5, 200.0],
        "coordinates": "elevation azimuth range",
    },
    melting_layer: {
        "units": "-",
        "standard_name": "melting_layer",
        "long_name": "Position of the range bin respect to the melting layer",
        "labels": ["BELOW", "ENTERING", "INSIDE", "EXITING", "ABOVE"],
        "ticks": [1, 2, 3, 4, 5],
        "boundaries": [0.5, 1.5, 2.5, 3.5, 4.5, 5.5],
        "coordinates": "elevation azimuth range",
        "scale_factor": 1,
        "add_offset": 0,
        "_FillValue": 0,
        "_Write_as_dtype": "uint8",
    },
    melting_layer_height: {
        "units": "m MSL",
        "standard_name": "melting_layer_height",
        "long_name": "Top and bottom melting layer height",
        "coordinates": "elevation azimuth",
    },
    probability_of_hail: {
        "units": "percent",
        "standard_name": "probability_of_hail",
        "long_name": "Probability of hail",
        "labels": ["0", "10", "20", "30", "40", "50", "60", "70", "80", "90"],
        "ticks": [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0],
        "boundaries": [
            0.0,
            10.0,
            20.0,
            30.0,
            40.0,
            50.0,
            60.0,
            70.0,
            80.0,
            90.0,
            100.0,
        ],
    },
    thunderstorm_forecast: {
        "units": "-",
        "standard_name": "thunderstorm_forecast",
        "long_name": "TRT Thunderstorm Forecast Ranks",
        "labels": [str(val) for val in range(41)],
        "ticks": list(range(41)),
        "boundaries": list(range(41)),
    },
    maximum_expected_severe_hail_size: {
        "units": "cm",
        "standard_name": "maximum_expected_severe_hail_size",
        "long_name": "Maximum expected severe hail size",
    },
    maximum_echo: {
        "units": "dBZ",
        "standard_name": "maximum_echo",
        "long_name": "Maximum echo",
    },
    maximum_echo_height: {
        "units": "km ASL",
        "standard_name": "maximum_echo_height",
        "long_name": "Maximum echo height",
    },
    echo_top_15dBZ: {
        "units": "km ASL",
        "standard_name": "echo_top_15_dBZ",
        "long_name": "Echo top 15 dBZ",
    },
    echo_top_20dBZ: {
        "units": "km ASL",
        "standard_name": "echo_top_20_dBZ",
        "long_name": "Echo top 20 dBZ",
    },
    echo_top_45dBZ: {
        "units": "km ASL",
        "standard_name": "echo_top_45_dBZ",
        "long_name": "Echo top 45 dBZ",
    },
    echo_top_50dBZ: {
        "units": "km ASL",
        "standard_name": "echo_top_50_dBZ",
        "long_name": "Echo top 50 dBZ",
    },
    vertically_integrated_liquid: {
        "units": "kg/m2",
        "standard_name": "vertically_integrated_liquid",
        "long_name": "Vertically integrated liquid",
    },
    specific_attenuation: {
        "units": "dB/km",
        "standard_name": "specific_attenuation",
        "long_name": "Specific attenuation",
        "coordinates": "elevation azimuth range",
    },
    path_integrated_attenuation: {
        "units": "dB",
        "standard_name": "path_integrated_attenuation",
        "long_name": "Path integrated attenuation",
        "coordinates": "elevation azimuth range",
    },
    specific_differential_attenuation: {
        "units": "dB/km",
        "standard_name": "specific_differential_attenuation",
        "long_name": "Specific differential attenuation",
        "coordinates": "elevation azimuth range",
    },
    path_integrated_differential_attenuation: {
        "units": "dB",
        "standard_name": "path_integrated_differential_attenuation",
        "long_name": "Path integrated differential attenuation",
        "coordinates": "elevation azimuth range",
    },
    corrected_specific_attenuation: {
        "units": "dB/km",
        "standard_name": "specific_attenuation",
        "long_name": "Corrected specific attenuation",
        "coordinates": "elevation azimuth range",
    },
    corrected_path_integrated_attenuation: {
        "units": "dB",
        "standard_name": "path_integrated_attenuation",
        "long_name": "Corrected path integrated attenuation",
        "coordinates": "elevation azimuth range",
    },
    corrected_specific_differential_attenuation: {
        "units": "dB/km",
        "standard_name": "specific_differential_attenuation",
        "long_name": "Corrected specific differential attenuation",
        "coordinates": "elevation azimuth range",
    },
    corrected_path_integrated_differential_attenuation: {
        "units": "dB",
        "standard_name": "path_integrated_differential_attenuation",
        "long_name": "Corrected path integrated differential attenuation",
        "coordinates": "elevation azimuth range",
    },
    number_of_samples: {
        "units": "count",
        "standard_name": "number_of_samples",
        "long_name": "Number of samples",
        "valid_min": 0,
        "coordinates": "elevation azimuth range",
    },
    standard_deviation: {
        "units": "-",
        "standard_name": "standard_deviation",
        "long_name": "Standard deviation",
        "valid_min": 0,
        "coordinates": "elevation azimuth range",
    },
    sum: {
        "units": "-",
        "standard_name": "sum",
        "long_name": "Sum",
        "coordinates": "elevation azimuth range",
    },
    sum_squared: {
        "units": "-",
        "standard_name": "sum_squared",
        "long_name": "Sum squared",
        "coordinates": "elevation azimuth range",
    },
    occurrence: {
        "units": "count",
        "standard_name": "occurrence",
        "long_name": "occurrence",
        "coordinates": "elevation azimuth range",
    },
    frequency_of_occurrence: {
        "units": "percent",
        "standard_name": "frequency_of_occurrence",
        "long_name": "Frequency of occurrence",
        "coordinates": "elevation azimuth range",
    },
    time_avg_flag: {
        "units": "count",
        "standard_name": "time_avg_flag",
        "long_name": "Time average flag",
        "coordinates": "elevation azimuth range",
    },
    colocated_gates: {
        "units": "flag",
        "standard_name": "colocated_gates",
        "long_name": "Colocated gates",
        "labels": ["FALSE", "TRUE"],
        "ticks": [1, 2],
        "boundaries": [0.5, 1.5, 2.5],
        "coordinates": "elevation azimuth range",
        "scale_factor": 1,
        "add_offset": 0,
        "_FillValue": 0,
        "_Write_as_dtype": "uint8",
    },
    # icon model data
    temperature: {
        "units": "deg Celsius",
        "standard_name": "temperature",
        "long_name": "Temperature",
        "coordinates": "elevation azimuth range",
    },
    iso0: {
        "units": "-",
        "standard_name": "iso0",
        "long_name": "Position of the range bin respect to the iso0 level",
        "labels": ["BELOW", "INSIDE", "ABOVE"],
        "ticks": [1, 2, 3],
        "boundaries": [0.5, 1.5, 2.5, 3.5],
        "coordinates": "elevation azimuth range",
        "scale_factor": 1,
        "add_offset": 0,
        "_FillValue": 0,
        "_Write_as_dtype": "uint8",
    },
    height_over_iso0: {
        "units": "m",
        "standard_name": "height_over_iso0",
        "long_name": "Height of the range bin respect to the iso0 level",
        "coordinates": "elevation azimuth range",
    },
    iso0_height: {
        "units": "m MSL",
        "standard_name": "iso0_height",
        "long_name": "iso0 height",
    },
    icon_index: {
        "units": "bin index",
        "standard_name": "icon_index",
        "long_name": ("indices of the icon model corresponding to each radar gate"),
        "coordinates": "elevation azimuth range",
    },
    hzt_index: {
        "units": "bin index",
        "standard_name": "hzt_index",
        "long_name": ("indices of the HZT corresponding to each radar gate"),
        "coordinates": "elevation azimuth range",
    },
    visibility: {
        "units": "percent",
        "standard_name": "visibility",
        "long_name": "visibility",
        "coordinates": "elevation azimuth range",
    },
    mininum_visible_elevation: {
        "units": "deg",
        "standard_name": "mininum_visible_elevation",
        "long_name": "Minimum visible elevation",
        "coordinates": "elevation azimuth range",
    },
    # Textures
    differential_phase_texture: {
        "units": "deg",
        "standard_name": "differential_phase_texture",
        "long_name": "Differential phase texture (PhiDP)",
        "coordinates": "elevation azimuth range",
    },
    differential_reflectivity_texture: {
        "units": "dB",
        "standard_name": "differential_reflectivity_texture",
        "long_name": "differential reflectivity texture",
        "coordinates": "elevation azimuth range",
    },
    reflectivity_texture: {
        "units": "dB",
        "standard_name": "reflectivity_texture",
        "long_name": "reflectivity texture",
        "coordinates": "elevation azimuth range",
    },
    cross_correlation_ratio_texture: {
        "units": "-",
        "standard_name": "copolar_correlation_coefficient_texture",
        "long_name": "Copolar correlation coefficient texture",
        "coordinates": "elevation azimuth range",
    },
    # Wind retrieval fields
    eastward_wind_component: {
        "units": "m/s",
        "standard_name": "eastward_wind_component",
        "long_name": "Eastward wind component",
    },
    northward_wind_component: {
        "units": "m/s",
        "standard_name": "northward_wind_component",
        "long_name": "Northward wind component",
    },
    vertical_wind_component: {
        "units": "m/s",
        "standard_name": "vertical_wind_component",
        "long_name": "Vertical wind component",
    },
    azimuthal_horizontal_wind_component: {
        "units": "m/s",
        "standard_name": "azimuthal_horizontal_wind_component",
        "long_name": "Azimuthal horizontal wind component",
    },
    vertical_wind_shear: {
        "units": "m/s/km",
        "standard_name": "vertical_wind_shear",
        "long_name": "Vertical wind shear",
    },
    wind_speed: {
        "units": "m/s",
        "standard_name": "horizontal_wind_speed",
        "long_name": "Horizontal wind speed",
    },
    wind_direction: {
        "units": "deg_from_north",
        "standard_name": "horizontal_wind_direction",
        "long_name": "Horizontal wind direction",
    },
    # profile variables
    height: {
        "long_name": "Height of radar beam",
        "standard_name": "height",
        "units": "meters",
    },
    interpolated_profile: {
        "long_name": "Interpolated profile",
        "standard_name": "interpolated_profile",
        "units": "unknown",
    },
    # wind lidar variables
    atmospherical_structures_type: {
        "units": "-",
        "standard_name": "atmospherical_structures_type",
        "long_name": (
            "Atmospherical structures detected out of the planetary boundary" " layer."
        ),
        # ND: No data, RL: residual layer, ML: mixed layer,
        # UC: unclassified cloud, IC: ice cloud, WC: water cloud,
        # UA: unclassified aerosol, SA: spherical aerosol,
        # AA: aspherical aerosol
        "labels": ["ND", "RL", "ML", "UC", "IC", "WC", "UA", "SA", "AA"],
        "ticks": [0, 20, 30, 142, 300, 400, 1475, 3000, 4000],
        "boundaries": [-15, 15, 25, 35, 250, 350, 450, 2500, 3500, 4500],
        "coordinates": "elevation azimuth range",
        "scale_factor": 1,
        "add_offset": 0,
        "_FillValue": -999,
        "_Write_as_dtype": "uint8",
    },
    radial_wind_speed_status: {
        "units": "flag",
        "standard_name": "radial_wind_speed_status",
        "long_name": "radial_wind_speed_status",
        "labels": ["Rejected", "Accepted"],
        "ticks": [0, 1],
        "boundaries": [-0.5, 0.5, 1.5],
        "coordinates": "elevation azimuth range",
        "scale_factor": 1,
        "add_offset": 0,
        "_FillValue": 0,
        "_Write_as_dtype": "uint8",
    },
    HRV: {
        "units": "%",
        "standard_name": "HRV",
        "long_name": "SEVIRI High Resolution Visible Reflectance",
    },
    VIS006: {
        "units": "%",
        "standard_name": "VIS006",
        "long_name": "SEVIRI Visible 0.6 um Reflectance",
    },
    VIS008: {
        "units": "%",
        "standard_name": "VIS008",
        "long_name": "SEVIRI Visible 0.8 um Reflectance",
    },
    IR_016: {
        "units": "%",
        "standard_name": "IR_016",
        "long_name": "SEVIRI Near-Infrared 1.6 um Reflectance",
    },
    IR_039: {
        "units": "K",
        "standard_name": "IR_039",
        "long_name": "SEVIRI Infrared 3.9 um",
    },
    WV_062: {
        "units": "K",
        "standard_name": "WV_062",
        "long_name": "SEVIRI Water Vapour 6.2 um",
    },
    WV_073: {
        "units": "K",
        "standard_name": "WV_073",
        "long_name": "SEVIRI Water Vapour 7.3 um",
    },
    IR_087: {
        "units": "K",
        "standard_name": "IR_087",
        "long_name": "SEVIRI Infrared 8.7 um",
    },
    IR_097: {
        "units": "K",
        "standard_name": "IR_097",
        "long_name": "SEVIRI Infrared 9.7 um",
    },
    IR_108: {
        "units": "K",
        "standard_name": "IR_108",
        "long_name": "SEVIRI Infrared 10.8 um",
    },
    IR_120: {
        "units": "K",
        "standard_name": "IR_120",
        "long_name": "SEVIRI Infrared 12 um",
    },
    IR_134: {
        "units": "K",
        "standard_name": "IR_134",
        "long_name": "SEVIRI Infrared 13.4 um",
    },
    CTH: {"units": "m", "standard_name": "CTH", "long_name": "Cloud Top Height"},
    HRV_norm: {
        "units": "%",
        "standard_name": "HRV_norm",
        "long_name": "Normalized SEVIRI High Resolution Visible Reflectance",
    },
    VIS006_norm: {
        "units": "%",
        "standard_name": "VIS006_norm",
        "long_name": "Normalized SEVIRI Visible 0.6 um Reflectance",
    },
    VIS008_norm: {
        "units": "%",
        "standard_name": "VIS008_norm",
        "long_name": "Normalized SEVIRI Visible 0.8 um Reflectance",
    },
    IR_016_norm: {
        "units": "%",
        "standard_name": "IR_016_norm",
        "long_name": "Normalized SEVIRI Near-Infrared 1.6 um Reflectance",
    },
    # Grid metadata
    "grid_time": {
        "units": "seconds",
        "standard_name": "time",
        "long_name": "Time of grid",
        "calendar": "gregorian",
    },
    "origin_longitude": {
        "long_name": "Longitude at grid origin",
        "units": "degrees_east",
        "standard_name": "longitude",
        "valid_min": -180.0,
        "valid_max": 180.0,
    },
    "origin_latitude": {
        "long_name": "Latitude at grid origin",
        "units": "degrees_north",
        "standard_name": "latitude",
        "valid_min": -90.0,
        "valid_max": 90.0,
    },
    "origin_altitude": {
        "long_name": "Altitude at grid origin",
        "units": "m",
        "standard_name": "altitude",
    },
    "x": {
        "standard_name": "projection_x_coordinate",
        "long_name": "X distance on the projection plane from the origin",
        "axis": "X",
        "units": "m",
    },
    "y": {
        "standard_name": "projection_y_coordinate",
        "long_name": "Y distance on the projection plane from the origin",
        "axis": "Y",
        "units": "m",
    },
    "z": {
        "standard_name": "projection_z_coordinate",
        "long_name": "Z distance on the projection plane from the origin",
        "axis": "Z",
        "units": "m",
        "positive": "up",
    },
    "point_x": {
        "long_name": "Cartesian x distance of each grid point from the origin",
        "units": "meters",
    },
    "point_y": {
        "long_name": "Cartesian y distance of each grid point from the origin",
        "units": "meters",
    },
    "point_z": {
        "long_name": "Cartesian z distance of each grid point from the origin",
        "positive": "up",
        "units": "meters",
    },
    "point_longitude": {
        "long_name": "Longitude of each grid point",
        "units": "degrees_north",
    },
    "point_latitude": {
        "long_name": "Latitude of each grid point",
        "units": "degrees_east",
    },
    "point_altitude": {"long_name": "Altitude of each grid point", "units": "meters"},
    "radar_latitude": {
        "long_name": "Latitude of radars used to make the grid.",
        "units": "degrees_north",
    },
    "radar_longitude": {
        "long_name": "Longitude of radars used to make the grid.",
        "units": "degrees_east",
    },
    "radar_altitude": {
        "long_name": "Altitude of radars used to make the grid.",
        "units": "m",
    },
    "radar_time": {
        "calendar": "gregorian",
        "long_name": "Time in seconds of the volume start for each radar",
    },
    "radar_name": {
        "long_name": "Name of radar used to make the grid",
    },
}


##############################################################################
# File specific metadata
#
# These dictionaries define metadata that is to be used only when reading in
# a given type of file.  This metadata is used in place of the
# DEFAULT_METADATA when it is avialable.  The main use of these variable
# is to define field specific data, it is safe to leave some/all of these
# empty if the default metadata is acceptable.
##############################################################################

# Metadata for Sigmet/IRIS files
sigmet_metadata = {}

# Metadata for NEXRAD Level II files (Archive and CDM files)
nexrad_metadata = {
    reflectivity: {
        "units": "dBZ",
        "standard_name": "equivalent_reflectivity_factor",
        "long_name": "Reflectivity",
        "valid_max": 94.5,
        "valid_min": -32.0,
        "coordinates": "elevation azimuth range",
    },
    velocity: {
        "units": "meters_per_second",
        "standard_name": "radial_velocity_of_scatterers_away_from_instrument",
        "long_name": "Mean doppler Velocity",
        "valid_max": 95.0,
        "valid_min": -95.0,
        "coordinates": "elevation azimuth range",
    },
    spectrum_width: {
        "units": "meters_per_second",
        "standard_name": "doppler_spectrum_width",
        "long_name": "Spectrum Width",
        "valid_max": 63.0,
        "valid_min": -63.5,
        "coordinates": "elevation azimuth range",
    },
    differential_reflectivity: {
        "units": "dB",
        "standard_name": "log_differential_reflectivity_hv",
        "long_name": "log_differential_reflectivity_hv",
        "valid_max": 7.9375,
        "valid_min": -7.8750,
        "coordinates": "elevation azimuth range",
    },
    differential_phase: {
        "units": "degrees",
        "standard_name": "differential_phase_hv",
        "long_name": "differential_phase_hv",
        "valid_max": 360.0,
        "valid_min": 0.0,
        "coordinates": "elevation azimuth range",
    },
    cross_correlation_ratio: {
        "units": "ratio",
        "standard_name": "cross_correlation_ratio_hv",
        "long_name": "Cross correlation_ratio (RHOHV)",
        "valid_max": 1.0,
        "valid_min": 0.0,
        "coordinates": "elevation azimuth range",
    },
}

# Metadata for NEXRAD Level 3 Products
nexrad_level3_metadata = {
    radar_estimated_rain_rate: {
        "units": "inches",
        "standard_name": "radar_estimated_rain_rate",
        "long_name": "Radar estimated rain rate",
        "coordinates": "elevation azimuth range",
    },
    radar_echo_classification: {
        "units": "legend",
        "standard_name": "radar_echo_classification",
        "long_name": "Radar echo classification",
        "options": (
            "0: Below Threshold (ND), "
            "10: Biological (BI), "
            "20: Anomalous Propagation/Ground Clutter (GC), "
            "30: Ice Crystals (IC), "
            "40: Dry Snow (DS), "
            "50: Wet Snow (WS), "
            "60: Light and/or Moderate Rain (RA), "
            "70: Heavy Rain (HR), "
            "80: Big Drops (rain) (BD), "
            "90: Graupel (GR), "
            "100: Hail, possibly with rain (HA), "
            "140: Unknown Classification (UK), "
            "150: Range Folded (RH)"
        ),
        "coordinates": "elevation azimuth range",
    },
}

# Metadata for CF/Radial files
cfradial_metadata = {}

# Metadata for MDV files
mdv_metadata = {}

# Metadata for RSL files
rsl_metadata = {}

# Metadata for CSU-CHILL, CHL files
chl_metadata = {}

FILE_SPECIFIC_METADATA = {  # Required
    "sigmet": sigmet_metadata,
    "nexrad_archive": nexrad_metadata,
    "nexrad_cdm": nexrad_metadata,
    "nexrad_level3": nexrad_level3_metadata,
    "cfradial": cfradial_metadata,
    "mdv": mdv_metadata,
    "rsl": rsl_metadata,
    "chl": chl_metadata,
}

##############################################################################
# Field name mapping
#
# These dictionaries map file field names or data types to a radar field
# name.  These are used to populate the radar.fields dictionary during a read
# in Py-ART.  A value of None will not include that field in the radar object.
# These can be over-ridden on a per-read basis using the field_mapping
# parameter, or using setting the file_field_names parameter to True.
##############################################################################

# Sigmet/IRIS file field mapping
# Note that multiple sigmet fields map to the same radar field, if
# more than one of these fields are present the radar field will be
# overwritten with the last sigmet field.
sigmet_field_mapping = {
    # Sigmet data type :field name              # (Data_type) Description
    "XHDR": None,  # (0) Extended Header
    "DBT": total_power,  # (1) Total Power
    "DBZ": reflectivity,  # (2) Reflectivity
    "VEL": velocity,  # (3) Velocity
    "WIDTH": spectrum_width,  # (4) Width
    "ZDR": differential_reflectivity,  # (5) Diff. reflectivity
    "DBZC": corrected_reflectivity,  # (7) Corrected reflectivity
    "DBT2": total_power,  # (8) Total Power
    "DBZ2": reflectivity,  # (9) Reflectivity
    "VEL2": velocity,  # (10) Velocity
    "WIDTH2": spectrum_width,  # (11) Width
    "ZDR2": differential_reflectivity,  # (12) Diff. reflectivity
    "RAINRATE2": radar_estimated_rain_rate,  # (13) Rainfall rate
    "KDP": specific_differential_phase,  # (14) KDP (diff. phase)
    "KDP2": specific_differential_phase,  # (15) KDP (diff. phase)
    "PHIDP": differential_phase,  # (16) PhiDP (diff. phase)
    "VELC": corrected_velocity,  # (17) Corrected velocity
    "SQI": normalized_coherent_power,  # (18) SQI
    "RHOHV": cross_correlation_ratio,  # (19) RhoHV
    "RHOHV2": cross_correlation_ratio,  # (20) RhoHV
    "DBZC2": corrected_reflectivity,  # (21) Corrected Reflec.
    "VELC2": corrected_velocity,  # (21) Corrected Velocity
    "SQI2": normalized_coherent_power,  # (23) SQI
    "PHIDP2": differential_phase,  # (24) PhiDP (diff. phase)
    "LDRH": linear_depolarization_ratio_h,  # (25) LDR xmt H, rcv V
    "LDRH2": linear_depolarization_ratio_h,  # (26) LDR xmt H, rcv V
    "LDRV": linear_depolarization_ratio_v,  # (27) LDR xmt V, rcv H
    "LDRV2": linear_depolarization_ratio_v,  # (28) LDR xmt V, rcv H
    "HEIGHT": None,  # (32) Height (1/10 km)
    "VIL2": None,  # (33) Linear Liquid
    "RAW": None,  # (34) Raw Data
    "SHEAR": None,  # (35) Wind Shear
    "DIVERGE2": None,  # (36) Divergence
    "FLIQUID2": None,  # (37) Floated liquid
    "USER": None,  # (38) User type
    "OTHER": None,  # (39) Unspecified
    "DEFORM2": None,  # (40) Deformation
    "VVEL2": None,  # (41) Vertical velocity
    "HVEL2": None,  # (42) Horizontal velocity
    "HDIR2": None,  # (43) Horiz. wind direction
    "AXDIL2": None,  # (44) Axis of dilation
    "TIME2": None,  # (45) Time in seconds
    "RHOH": None,  # (46) Rho, xmt H, rcv V
    "RHOH2": None,  # (47) Rho, xmt H, rcv V
    "RHOV": None,  # (48) Rho, xmt V, rcv H
    "RHOV2": None,  # (49) Rho, xmt V, rcv H
    "PHIH": None,  # (50) Phi, xmt H, rcv V
    "PHIH2": None,  # (51) Phi, xmt H, rcv V
    "PHIV": None,  # (52) Phi, xmt V, rcv H
    "PHIV2": None,  # (53) Phi, xmt V, rcv H
    "USER2": None,  # (54) User type
    "HCLASS": radar_echo_classification,  # (55) Hydrometeor class
    "HCLASS2": radar_echo_classification,  # (56) Hydrometeor class
    "ZDRC": corrected_differential_reflectivity,
    # (57) Corrected diff. refl.
    "ZDRC2": corrected_differential_reflectivity,
    # (58) Corrected diff. refl.
    "UNKNOWN_59": None,  # Unknown field
    "UNKNOWN_60": None,  # Unknown field
    "UNKNOWN_61": None,  # Unknown field
    "UNKNOWN_62": None,  # Unknown field
    "UNKNOWN_63": None,  # Unknown field
    "UNKNOWN_64": None,  # Unknown field
    "UNKNOWN_65": None,  # Unknown field
    "UNKNOWN_66": None,  # Unknown field
    # there may be more field, add as needed
}


# NEXRAD Level II Archive files
nexrad_archive_field_mapping = {
    # NEXRAD field: radar field name
    "REF": reflectivity,
    "VEL": velocity,
    "SW": spectrum_width,
    "ZDR": differential_reflectivity,
    "PHI": differential_phase,
    "RHO": cross_correlation_ratio,
}

# NEXRAD Level II CDM files
nexrad_cdm_field_mapping = {
    # CDM variable name (without _HI): radar field name
    "Reflectivity": reflectivity,
    "RadialVelocity": velocity,
    "SpectrumWidth": spectrum_width,
    "DifferentialReflectivity": differential_reflectivity,
    "DifferentialPhase": differential_phase,
    "CorrelationCoefficient": cross_correlation_ratio,
}

# NEXRAD Level 3 Product files.
nexrad_level3_mapping = {
    # Message code : field name         # Product name
    19: reflectivity,  # Base Reflectivity
    20: reflectivity,  # Base Reflectivity
    25: velocity,  # Base Velocity
    27: velocity,  # Base Velocity
    28: spectrum_width,  # Base Spectrum Width
    30: spectrum_width,  # Base Spectrum Width
    32: reflectivity,  # Digital Hybrid Scan Reflectivity
    34: None,  # Clutter Filter Control
    56: velocity,  # Storm Relative Mean Radial Velocity
    78: radar_estimated_rain_rate,  # Surface Rainfall Accum. (1 hr)
    79: radar_estimated_rain_rate,  # Surface Rainfall Accum. (3 hr)
    80: radar_estimated_rain_rate,  # Storm Total Rainfall Accumulation
    94: reflectivity,  # Base Reflectivity Data Array
    99: velocity,  # Base Velocity Data Array
    134: None,  # High Resolution VIL
    135: None,  # Enhanced Echo Tops
    138: radar_estimated_rain_rate,  # Digital Storm Total Precipitation
    159: differential_reflectivity,  # Digital Differential Reflectivity
    161: cross_correlation_ratio,  # Digital Correlation Coefficient
    163: specific_differential_phase,  # Digital Specific Differential Phase
    165: radar_echo_classification,  # Digital Hydrometeor Classification
    169: radar_estimated_rain_rate,  # One Hour Accumulation
    170: radar_estimated_rain_rate,  # Digital Accumulation Array
    171: radar_estimated_rain_rate,  # Storm Total Accumulation
    172: radar_estimated_rain_rate,  # Digital Storm Total Accumulation
    173: radar_estimated_rain_rate,  # Digital User-Selectable Accum.
    174: radar_estimated_rain_rate,  # Digital 1 hr Diff. Accum.
    175: radar_estimated_rain_rate,  # Digital Storm Total Diff. Accum.
    176: radar_estimated_rain_rate,  # Digital Inst. Precipitation Rate
    177: radar_echo_classification,  # Hybrid Hydrometeor Classification
    181: reflectivity,  # Base Reflectivity
    182: velocity,  # Base Velocity
    186: reflectivity,  # Base Reflectivity
}

# MDV files
mdv_field_mapping = {
    # MDV moment: radar field name
    "DBZ_F": reflectivity,
    "VEL_F": velocity,
    "WIDTH_F": spectrum_width,
    "ZDR_F": differential_reflectivity,
    "RHOHV_F": cross_correlation_ratio,
    "NCP_F": normalized_coherent_power,
    "KDP_F": specific_differential_phase,
    "PHIDP_F": differential_phase,
    "VEL_COR": corrected_velocity,
    "PHIDP_UNF": unfolded_differential_phase,
    "KDP_SOB": corrected_specific_differential_phase,
    "DBZ_AC": corrected_reflectivity,
    # repeated integer moments
    "DBZ": reflectivity,
    "VEL": velocity,
    "WIDTH": spectrum_width,
    "ZDR": differential_reflectivity,
    "RHOHV": cross_correlation_ratio,
    "NCP": normalized_coherent_power,
    "KDP": specific_differential_phase,
    "PHIDP": differential_phase,
}

# CF/Radial files
cfradial_field_mapping = {}

# RSL files
# Note that multiple RSL field map to the same radar field, if
# more than one of these fields are present in the RSL data structure
# the radar field will be overwritten with the last field.
rsl_field_mapping = {
    # RSL 2 letter field: radar field           # RSL description
    "DZ": reflectivity,  # reflectivity
    "VR": velocity,  # velocity
    "SW": spectrum_width,  # spectrum width
    "CZ": corrected_reflectivity,  # corrected reflectivity
    "ZT": reflectivity,  # uncorrected reflectivity
    "DR": differential_reflectivity,  # differential reflectivity
    "LR": differential_reflectivity,  # another diff. reflectivity
    "ZD": differential_reflectivity,  # another diff. reflectivity
    "DM": None,  # received power
    "RH": cross_correlation_ratio,  # RhoHV
    "PH": differential_phase,  # PhiDP
    "XZ": None,  # X-band reflectivity
    "CD": corrected_differential_reflectivity,  # Corrected DR.
    "MZ": None,  # DZ mask
    "MD": None,  # DR Mask
    "ZE": corrected_reflectivity,  # edited reflectivity
    "VE": corrected_velocity,  # edited velocity
    "KD": specific_differential_phase,  # specific diff. phase
    "TI": None,  # TIME (unknown)
    "DX": None,  # ???
    "CH": None,  # ???
    "AH": None,  # ???
    "CV": None,  # ???
    "AV": None,  # ???
    "SQ": normalized_coherent_power,  # Signal Quality Index (sigmet)
    "VS": None,  # Radial Vel. combined
    "VL": None,  # Radial Vel. combined
    "VG": None,  # Radial Vel. combined
    "VT": None,  # Radial Vel. combined
    "NP": normalized_coherent_power,  # Normalized Coherent Power
    "HC": radar_echo_classification,  # Hydroclass
    "VC": None,  # Radial Vel. Corrected.
    "V2": None,  # Radial Vel cut 2
    "S2": None,  # Spectrum width cut 2
    "V3": None,  # Radial Vel cut 3
    "S3": None,  # Spectrum width cut 3
}

chl_field_mapping = {
    # Chill field name : radar field name
    "Z": reflectivity,
    "V": velocity,
    "W": spectrum_width,
    "ZDR": differential_reflectivity,
    "LDRH": linear_depolarization_ratio_h,
    "LDRV": linear_depolarization_ratio_v,
    b"\xce\xa8 DP".decode("utf-8"): differential_phase,
    "KDP": specific_differential_phase,
    b"\xcf\x81 HV".decode("utf-8"): cross_correlation_ratio,
    "NCP": normalized_coherent_power,
    # These fields are not mapped by default
    "H Re(lag 1)": None,  # Real part of lag-1 correlation, H Channel
    "V Re(lag 2)": None,  # Real part of lag-2 correlation, V Channel
    "VAvgQ": None,  # Average Q, V Channel
    "V Im(lag 1)": None,  # Imaginary part of lag-1 correlation, V Channel
    "HAvgQ": None,  # Average Q, H Channel
    "H Im(lag 2)": None,  # Imaginary part of lag-2 correlation, H Channel
    "V lag 0": None,  # Absolute value of lag-0 correlation, V Channel
    "H lag 0": None,  # Absolute value of lag-0 correlation, H Channel
    "H lag 0 cx": None,  # Absolute value of lag-0 cross correlation,
    # H Channel
    "H Im(lag 1)": None,  # Imaginary part of lag-1 correlation, H Channel
    "H Re(lag 2)": None,  # Real part of lag-2 correlation, H Channel
    "V lag 0 cx": None,  # Absolute value of lag-0 cross correlation,
    # V Channel
    "V Re(lag 1)": None,  # Real part of lag-1 correlation, V Channel
    "V Im(lag 2)": None,  # Imaginary part of lag-2 correlation, V Channel
    "HV lag 0 I": None,  # Real part of cross channel correlation at lag-0
    "HV lag 0 Q": None,  # Imaginary part of cross channel correlation at
    # lag-0
    "VAvgI": None,  # Average I, V Channel
    "HAvgI": None,  # Average I, H Channel
    b"\xcf\x81 HCX".decode("utf-8"): None,  # H Co to Cross Correlation
    b"\xcf\x81 VCX".decode("utf-8"): None,  # V Co to Cross Correlation
}

# GAMIC HDF5 files
gamic_field_mapping = {
    # Description of GAMIC fields taken from Radx source code file:
    # GamicHdf5RadxFile.cc
    # General notes on field names
    # ----------------------------
    # * An 'U' prefix indicated the moment was corrected from a timeseries
    #   that was not clutter corrected.  Fields without a 'U' prefix are
    #   calculated from a timeseries that has been clutter corrected.
    # * An 'A' prefix indicates that the moment has been corrected
    #   for rainfall attenuation.
    # * 'h' and 'v' endings indicate a the horizontal and vertical channels
    # GAMIC field name: radar field name
    "I": None,  # In phase signal
    "Q": None,  # Quadrature signal
    "Z": corrected_reflectivity,
    # Reflectivity, corrected for 2nd trip & clutter
    "UZ": reflectivity,
    # Uncorrected reflectivity
    "Zh": corrected_reflectivity,
    # Corrected reflectivity, horizontal channel
    "Zv": None,  # Corrected reflectivity, vertical channel
    "UZh": reflectivity,
    # Uncorrected reflectivity, horizontal channel
    "UZv": None,  # Uncorrected reflectivity, vertical channel
    "AZh": None,  # Refl., rainfall atten. & clutter corrected, h chan.
    # 'F' in velocity fields indicate folded velocities (no de-aliasing),
    # velocities fields without an F have been unfolded.
    "V": corrected_velocity,
    # Unfolded velocity from corrected timeseries
    "VF": velocity,  # Folded velocity from corr. t.s.
    "UV": None,  # Unfolded velocity from uncorrected timeseries
    "UVF": None,  # Folded velcity from uncorr. t.s.
    "Vh": corrected_velocity,
    # Velocity from corr. t.s., horizontal channel
    "Vv": None,  # Velocity from corr. t.s., vertical channel
    "UVh": None,  # Velocity from uncorr. t.s., horizontal channel
    "UVv": None,  # Velocity from uncorr. t.s., vertical channel
    "VFh": velocity,
    # Folded velocity, corr. t.s., horizontal channel
    "VFv": None,  # Folded velocity, corr. t.s., vertical channel
    "UnV": None,  # Folded velocity from uncorrected timeseries
    "UnVFh": None,  # Folded velocity, uncorr. t.s., horizontal channel
    "UnVFv": None,  # Folded velocity, uncorr. t.s., vertical channel
    # 'C' in spectral width fields indicate that the field has been
    # corrected for decorrelation causes by antenna rotation
    "W": spectrum_width,
    # Spectral width from clutter corrected time series.
    "UW": None,  # Spectral width from uncorrected timeseries.
    "CW": None,  # Spec. width, antenna rotation corrected, corr. t.s.
    "UCW": None,  # Spec. width, antenna rotation corrected, uncorr t.s.
    "Wh": spectrum_width,
    # Spectral width, corr t.s., horizontal channel
    "Wv": None,  # Spectral width, corr t.s., vertical channel
    "UWh": None,  # Spectral width, uncorr t.s., horizontal channel
    "UWv": None,  # Spectral width, uncorr t.s., vertical channel
    "CWh": None,  # Spec. width, antenna rot. corr., horizontal channel
    "CWv": None,  # Spec. width, antenna rot. corr., vertical channel
    "SQI": normalized_coherent_power,
    # Signal quality index
    "SQIh": normalized_coherent_power,
    # Signal quality index, horizontal channel
    "SQIv": None,  # Signal quality index, vertical channel
    "CCOR": None,  # Clutter power correction
    "CCORh": None,  # Clutter power correction, horizontal channel
    "CCORv": None,  # Clutter power correction, vertical channel
    "SIGPOW": None,  # Singal Power
    "SNR": None,  # Raw signal to noise ratio
    "SNRh": None,  # Raw signal to noise ratio, horizontal channel
    "SNRv": None,  # Raw signal to noise ration, vertical channel
    "DFT": None,  # Signal spectrum amplitude
    "DFTh": None,  # Signal spectrum amplitude, horizontal channel
    "DFTv": None,  # Signal spectrum amplitude, vertical channel
    "LOG": None,  # Logarithmic amplitude 10*log Isq + Qsq
    "LOGh": None,  # Logarithmic amplitude, horizontal channel
    "LOGv": None,  # Logarithmic amplitude, vertical channel
    "CMAP": None,  # Censor map
    # A '1' ending on differential reflectivity fields indicates that the
    # moment has calculated using a 1st LAG algorithm.
    "ZDR": corrected_differential_reflectivity,
    # Differential reflectivity from corrected timeseries
    "UZDR": differential_reflectivity,
    # Diff. refl. from uncorrected timeseries
    "AZDR": None,  # Diff. refl., rainfall atten. corr., corr t.s.
    "ZDR1": None,  # Diff. refl., corr. t.s., 1st LAG algo.
    "UZDR1": None,  # Diff. refl., uncorr. t.s., 1st LAG algo.
    "AZDR1": None,  # Diff. refl., rain. atten. corr., corr. t.s., 1st LAG
    "PHI": corrected_differential_phase,
    # Differential phase from corrected timeseries
    "PHIDP": corrected_differential_phase,
    # Differential phase from corrected timeseries
    "UPHIDP": differential_phase,
    # Diff. phase from uncorrected timeseries.
    "PHIH": None,  # Diff. phase, corr. t.s., horizontal channel
    "UPHIH": None,  # Diff. phase, uncorr t.s., hortizontal channel
    "KDP": specific_differential_phase,
    # Specific differential phase
    "RHO": cross_correlation_ratio,
    # Cross correlation coefficient from corrected t.s.
    "RHOHV": cross_correlation_ratio,
    # Cross correlation coefficient from corrected t.s.
    "URHOHV": None,  # Cross correlation coefficient from uncorr. t.s.
    "RHOH": None,  # Cross corr., corr. t.s., horizontal transmit only
    "URHOH": None,  # Cross corr., uncorr t.s., horizontal transmit only
    "LDR": linear_depolarization_ratio,
    # Linear depolarization ratio from corr. t.s.
    "ULDR": None,  # Linear depolarization ratio  from uncorr. t.s.
}

# UF field mappings
uf_field_mapping = {
    # UF 2 letter field: radar field
    "DZ": reflectivity,
    "CZ": corrected_reflectivity,
    "ZE": corrected_reflectivity,
    "ZT": total_power,
    "UZ": total_power,
    "VR": velocity,
    "VE": corrected_velocity,
    "SW": spectrum_width,
    "ZD": differential_reflectivity,
    "DR": corrected_differential_reflectivity,
    "CD": corrected_differential_reflectivity,
    "LR": linear_depolarization_ratio,
    "PH": differential_phase,
    "KD": specific_differential_phase,
    "RH": cross_correlation_ratio,
    "SQ": normalized_coherent_power,
    "NP": normalized_coherent_power,
    "HC": radar_echo_classification,
}

# Mapping used when writing UF files
write_uf_mapping = {
    # UF 2 letter field: radar field
    reflectivity: "DZ",
    corrected_reflectivity: "CZ",
    total_power: "ZT",
    velocity: "VR",
    corrected_velocity: "VE",
    spectrum_width: "SW",
    differential_reflectivity: "ZD",
    corrected_differential_reflectivity: "DR",
    linear_depolarization_ratio: "LR",
    differential_phase: "PH",
    specific_differential_phase: "KD",
    cross_correlation_ratio: "RH",
    normalized_coherent_power: "SQ",
    radar_echo_classification: "HC",
}

FIELD_MAPPINGS = {  # Required variable
    "sigmet": sigmet_field_mapping,
    "nexrad_archive": nexrad_archive_field_mapping,
    "nexrad_cdm": nexrad_cdm_field_mapping,
    "nexrad_level3": nexrad_level3_mapping,
    "cfradial": cfradial_field_mapping,
    "mdv": mdv_field_mapping,
    "rsl": rsl_field_mapping,
    "chl": chl_field_mapping,
    "gamic": gamic_field_mapping,
    "uf": uf_field_mapping,
    "write_uf": write_uf_mapping,
}


def velocity_limit(container=None, selection=0):
    import pyart

    if isinstance(container, pyart.core.Radar):
        try:
            if selection >= 0 and selection < container.nsweeps:
                vel = container.get_nyquist_vel(selection, check_uniform=False)
            else:
                vel = container.get_nyquist_vel(0, check_uniform=False)
            return (-vel, vel)
        except LookupError:
            # return (-42., 42.)
            return (-15.8, 15.8)
    else:
        # return (-42., 42.)
        return (-15.8, 15.8)


def spectrum_width_limit(container=None, selection=0):
    import pyart

    if isinstance(container, pyart.core.Radar):
        try:
            if selection >= 0 and selection < container.nsweeps:
                vel = container.get_nyquist_vel(selection, check_uniform=False)
            else:
                vel = container.get_nyquist_vel(0, check_uniform=False)
            return (0, vel)
        except LookupError:
            # return (0., 4.)
            return (0.0, 15.8)
    else:
        # return (0., 4.)
        return (0.0, 15.8)


DEFAULT_FIELD_COLORMAP = {
    # field name : colormap
    reflectivity: "ChaseSpectral",
    corrected_reflectivity: "ChaseSpectral",
    total_power: "ChaseSpectral",
    unfiltered_reflectivity: "ChaseSpectral",
    corrected_unfiltered_reflectivity: "ChaseSpectral",
    reflectivity_vv: "ChaseSpectral",
    corrected_reflectivity_vv: "ChaseSpectral",
    unfiltered_reflectivity_vv: "ChaseSpectral",
    reflectivity_bias: "ChaseSpectral",
    signal_power_hh: "ChaseSpectral",
    signal_power_vv: "ChaseSpectral",
    volumetric_reflectivity: "ChaseSpectral",
    volumetric_reflectivity_vv: "ChaseSpectral",
    bird_density: "ChaseSpectral",
    bird_reflectivity: "ChaseSpectral",
    radar_cross_section_hh: "ChaseSpectral",
    radar_cross_section_vv: "ChaseSpectral",
    spectral_reflectivity_hh: "ChaseSpectral",
    spectral_reflectivity_vv: "ChaseSpectral",
    unfiltered_spectral_reflectivity_hh: "ChaseSpectral",
    unfiltered_spectral_reflectivity_vv: "ChaseSpectral",
    avg_reflectivity: "ChaseSpectral",
    quant05_reflectivity: "ChaseSpectral",
    quant10_reflectivity: "ChaseSpectral",
    quant20_reflectivity: "ChaseSpectral",
    quant50_reflectivity: "ChaseSpectral",
    quant80_reflectivity: "ChaseSpectral",
    quant90_reflectivity: "ChaseSpectral",
    quant95_reflectivity: "ChaseSpectral",
    noisedBm_hh: "viridis",
    noisedBm_vv: "viridis",
    noisedBZ_hh: "ChaseSpectral",
    noisedBZ_vv: "ChaseSpectral",
    noisedBADU_hh: "ChaseSpectral",
    noisedBADU_vv: "ChaseSpectral",
    noiseADU_hh: "ChaseSpectral",
    noiseADU_vv: "ChaseSpectral",
    noise_pos_h: "LangRainbow12",
    noise_pos_v: "LangRainbow12",
    signal_to_noise_ratio: "Carbone17",
    signal_to_noise_ratio_hh: "Carbone17",
    signal_to_noise_ratio_vv: "Carbone17",
    clutter_correction_ratio_hh: "Carbone17",
    clutter_correction_ratio_vv: "Carbone17",
    visibility: "Carbone17",
    frequency_of_occurrence: "Carbone17",
    occurrence: "Carbone17",
    stat_test_lag1: "ChaseSpectral",
    stat_test_lag2: "ChaseSpectral",
    wide_band_noise: "ChaseSpectral",
    fields_difference: "BuDRd18",
    field_mask: "Greys_r",
    field_texture: "Greys_r",
    sun_hit_power_h: "ChaseSpectral",
    sun_hit_power_v: "ChaseSpectral",
    sun_hit_differential_reflectivity: "RefDiff",
    sun_est_power_h: "ChaseSpectral",
    sun_est_power_v: "ChaseSpectral",
    sun_est_differential_reflectivity: "RefDiff",
    simulated_velocity: "BuDRd18",
    velocity: "BuDRd18",
    corrected_velocity: "BuDRd18",
    unfiltered_velocity: "BuDRd18",
    dealiased_corrected_velocity: "BuDRd18",
    dealiased_velocity: "BuDRd18",
    velocity_vv: "BuDRd18",
    unfiltered_velocity_vv: "BuDRd18",
    eastward_wind_component: "BuDRd18",
    northward_wind_component: "BuDRd18",
    vertical_wind_component: "BuDRd18",
    azimuthal_horizontal_wind_component: "BuDRd18",
    vertical_wind_shear: "BuDRd18",
    retrieved_velocity: "BuDRd18",
    retrieved_velocity_std: "ChaseSpectral",
    velocity_difference: "BuDRd18",
    wind_speed: "ChaseSpectral",
    wind_direction: "Wild25",
    avg_velocity: "BuDRd18",
    quant05_velocity: "BuDRd18",
    quant10_velocity: "BuDRd18",
    quant20_velocity: "BuDRd18",
    quant50_velocity: "BuDRd18",
    quant80_velocity: "BuDRd18",
    quant90_velocity: "BuDRd18",
    quant95_velocity: "BuDRd18",
    avg_corrected_velocity: "BuDRd18",
    quant05_corrected_velocity: "BuDRd18",
    quant10_corrected_velocity: "BuDRd18",
    quant20_corrected_velocity: "BuDRd18",
    quant50_corrected_velocity: "BuDRd18",
    quant80_corrected_velocity: "BuDRd18",
    quant90_corrected_velocity: "BuDRd18",
    quant95_corrected_velocity: "BuDRd18",
    spectrum_width: "NWS_SPW",
    corrected_spectrum_width: "NWS_SPW",
    unfiltered_spectrum_width: "NWS_SPW",
    spectrum_width_vv: "NWS_SPW",
    unfiltered_spectrum_width_vv: "NWS_SPW",
    turbulence: "NWS_SPW",
    normalized_coherent_power: "Carbone17",
    differential_reflectivity: "RefDiff",
    corrected_differential_reflectivity: "RefDiff",
    unfiltered_differential_reflectivity: "RefDiff",
    differential_reflectivity_in_precipitation: "RefDiff",
    differential_reflectivity_in_snow: "RefDiff",
    differential_reflectivity_column_height: "RefDiff",
    spectral_differential_reflectivity: "RefDiff",
    unfiltered_spectral_differential_reflectivity: "RefDiff",
    cross_correlation_ratio: "RefDiff",
    corrected_cross_correlation_ratio: "RefDiff",
    unfiltered_cross_correlation_ratio: "RefDiff",
    uncorrected_cross_correlation_ratio: "RefDiff",
    logarithmic_cross_correlation_ratio: "RefDiff",
    cross_correlation_ratio_in_rain: "RefDiff",
    spectral_copolar_correlation_coefficient: "RefDiff",
    unfiltered_spectral_copolar_correlation_coefficient: "RefDiff",
    mean_phase: "Wild25",
    differential_phase: "Wild25",
    unfolded_differential_phase: "Wild25",
    corrected_differential_phase: "Wild25",
    uncorrected_differential_phase: "Wild25",
    uncorrected_unfiltered_differential_phase: "Wild25",
    system_differential_phase: "Wild25",
    spectral_phase_hh: "Wild25",
    spectral_phase_vv: "Wild25",
    spectral_differential_phase: "Wild25",
    unfiltered_spectral_phase_hh: "Wild25",
    unfiltered_spectral_phase_vv: "Wild25",
    unfiltered_spectral_differential_phase: "Wild25",
    specific_differential_phase: "Theodore16",
    corrected_specific_differential_phase: "Theodore16",
    uncorrected_specific_differential_phase: "Theodore16",
    uncorrected_unfiltered_specific_differential_phase: "Theodore16",
    linear_depolarization_ratio: "SCook18",
    linear_depolarization_ratio_h: "SCook18",
    linear_depolarization_ratio_v: "SCook18",
    circular_depolarization_ratio: "SCook18",
    rain_rate: "RRate11",
    radar_estimated_rain_rate: "RRate11",
    corrected_radar_estimated_rain_rate: "RRate11",
    avg_radar_estimated_rain_rate: "RRate11",
    quant05_radar_estimated_rain_rate: "RRate11",
    quant10_radar_estimated_rain_rate: "RRate11",
    quant20_radar_estimated_rain_rate: "RRate11",
    quant50_radar_estimated_rain_rate: "RRate11",
    quant80_radar_estimated_rain_rate: "RRate11",
    quant90_radar_estimated_rain_rate: "RRate11",
    quant95_radar_estimated_rain_rate: "RRate11",
    rainfall_accumulation: "RRate11",
    sun_hit_h: "LangRainbow12",
    sun_hit_v: "LangRainbow12",
    sun_hit_zdr: "LangRainbow12",
    radar_echo_classification: "LangRainbow12",
    corrected_radar_echo_classification: "LangRainbow12",
    hydroclass_entropy: "LangRainbow12",
    proportion_AG: "LangRainbow12",
    proportion_CR: "LangRainbow12",
    proportion_LR: "LangRainbow12",
    proportion_RP: "LangRainbow12",
    proportion_RN: "LangRainbow12",
    proportion_VI: "LangRainbow12",
    proportion_WS: "LangRainbow12",
    proportion_MH: "LangRainbow12",
    proportion_IH: "LangRainbow12",
    radar_echo_id: "LangRainbow12",
    clutter_exit_code: "LangRainbow12",
    melting_layer: "LangRainbow12",
    height_over_iso0: "BuDRd18",
    specific_attenuation: "Carbone17",
    path_integrated_attenuation: "Carbone17",
    specific_differential_attenuation: "Carbone17",
    path_integrated_differential_attenuation: "Carbone17",
    corrected_specific_attenuation: "Carbone17",
    corrected_path_integrated_attenuation: "Carbone17",
    corrected_specific_differential_attenuation: "Carbone17",
    corrected_path_integrated_differential_attenuation: "Carbone17",
    differential_phase_texture: "BlueBrown11",
    differential_reflectivity_texture: "BlueBrown11",
    reflectivity_texture: "BlueBrown11",
    cross_correlation_ratio_texture: "BlueBrown11",
    height: "SCook18",
    interpolated_profile: "SCook18",
    corrected_radial_wind_speed: "BuDRd18",
    radial_wind_speed: "BuDRd18",
    # radial_wind_speed_ci:
    # radial_wind_speed_status:
    doppler_spectrum_width: "NWS_SPW",
    # doppler_spectrum_mean_error:
    atmospherical_structures_type: "LangRainbow12",
    relative_beta: "ChaseSpectral",
    absolute_beta: "ChaseSpectral",
    cnr: "Carbone17",
    # appropriate colors are matplotlib Greys, gray, etc.
    HRV: "Greys_r",  # 'Gray9',
    VIS006: "Greys_r",  # 'Gray9',
    VIS008: "Greys_r",  # 'Gray9',
    IR_016: "Greys_r",  # 'Gray9_r',
    IR_039: "ChaseSpectral",
    WV_062: "ChaseSpectral",
    WV_073: "ChaseSpectral",
    IR_087: "ChaseSpectral",
    IR_097: "ChaseSpectral",
    IR_108: "ChaseSpectral",
    IR_120: "ChaseSpectral",
    IR_134: "ChaseSpectral",
    CTH: "ChaseSpectral",
    HRV_norm: "Greys_r",  # 'Gray9',
    VIS006_norm: "Greys_r",  # 'Gray9',
    VIS008_norm: "Greys_r",  # 'Gray9',
    IR_016_norm: "Greys_r",  # 'Gray9_r',
    # Additional reflectivity like fields
    "CZ": "ChaseSpectral",
    "DZ": "ChaseSpectral",
    "AZ": "ChaseSpectral",
    "Z": "ChaseSpectral",
    "dbz": "ChaseSpectral",
    "DBZ": "ChaseSpectral",
    "dBZ": "ChaseSpectral",
    "DBZH": "ChaseSpectral",
    "DBZ_S": "ChaseSpectral",
    "DBZ_K": "ChaseSpectral",
    "reflectivity_horizontal": "ChaseSpectral",
    "corr_reflectivity": "ChaseSpectral",
}

# map each field to a limit or a limit function

DEFAULT_FIELD_LIMITS = {
    # field name : limits
    # reflectivity: (-30., 75.),
    reflectivity: (-30.0, 85.0),
    avg_reflectivity: (-30.0, 75.0),
    quant05_reflectivity: (-30.0, 75.0),
    quant10_reflectivity: (-30.0, 75.0),
    quant20_reflectivity: (-30.0, 75.0),
    quant50_reflectivity: (-30.0, 75.0),
    quant80_reflectivity: (-30.0, 75.0),
    quant90_reflectivity: (-30.0, 75.0),
    quant95_reflectivity: (-30.0, 75.0),
    bird_reflectivity: (-30.0, 75.0),
    corrected_reflectivity: (-30.0, 75.0),
    total_power: (-30.0, 75.0),
    # unfiltered_reflectivity: (-30., 75.),
    unfiltered_reflectivity: (-30.0, 85.0),
    corrected_unfiltered_reflectivity: (-30.0, 75.0),
    # reflectivity_vv: (-30., 75.),
    reflectivity_vv: (-30.0, 85.0),
    corrected_reflectivity_vv: (-30.0, 75.0),
    # unfiltered_reflectivity_vv: (-30., 75.),
    unfiltered_reflectivity_vv: (-30.0, 85.0),
    signal_to_noise_ratio: (-5.0, 30.0),
    signal_to_noise_ratio_hh: (-5.0, 30.0),
    signal_to_noise_ratio_vv: (-5.0, 30.0),
    noisedBZ_hh: (-40.0, 10.0),
    noisedBZ_vv: (-40.0, 10.0),
    reflectivity_bias: (-30.0, 30.0),
    volumetric_reflectivity: (20.0, 60.0),
    volumetric_reflectivity_vv: (20.0, 60.0),
    bird_density: (0.0, 400.0),
    radar_cross_section_hh: (-50.0, 55.0),
    radar_cross_section_vv: (-50.0, 55.0),
    stat_test_lag1: (0.0, 10.0),
    stat_test_lag2: (0.0, 10.0),
    wide_band_noise: (0.0, 25.0),
    # spectral_reflectivity_hh: (-60., 45.),
    # spectral_reflectivity_vv: (-60., 45.),
    # unfiltered_spectral_reflectivity_hh: (-60., 45.),
    # unfiltered_spectral_reflectivity_vv: (-60., 45.),
    signal_power_hh: (-130.0, 0.0),
    signal_power_vv: (-130.0, 0.0),
    sun_hit_power_h: (-120.0, -90.0),
    sun_hit_power_v: (-120.0, -90.0),
    sun_hit_differential_reflectivity: (-2.0, 2.0),
    sun_est_power_h: (-120.0, -90.0),
    sun_est_power_v: (-120.0, -90.0),
    sun_est_differential_reflectivity: (-2.0, 2.0),
    velocity: velocity_limit,
    avg_velocity: velocity_limit,
    quant05_velocity: velocity_limit,
    quant10_velocity: velocity_limit,
    quant20_velocity: velocity_limit,
    quant50_velocity: velocity_limit,
    quant80_velocity: velocity_limit,
    quant90_velocity: velocity_limit,
    quant95_velocity: velocity_limit,
    avg_corrected_velocity: velocity_limit,
    quant05_corrected_velocity: velocity_limit,
    quant10_corrected_velocity: velocity_limit,
    quant20_corrected_velocity: velocity_limit,
    quant50_corrected_velocity: velocity_limit,
    quant80_corrected_velocity: velocity_limit,
    quant90_corrected_velocity: velocity_limit,
    quant95_corrected_velocity: velocity_limit,
    corrected_velocity: velocity_limit,
    unfiltered_velocity: velocity_limit,
    velocity_vv: velocity_limit,
    unfiltered_velocity_vv: velocity_limit,
    dealiased_corrected_velocity: velocity_limit,
    dealiased_velocity: velocity_limit,
    eastward_wind_component: velocity_limit,
    northward_wind_component: velocity_limit,
    vertical_wind_component: velocity_limit,
    azimuthal_horizontal_wind_component: velocity_limit,
    radial_wind_speed: velocity_limit,
    corrected_radial_wind_speed: velocity_limit,
    vertical_wind_shear: velocity_limit,
    retrieved_velocity: velocity_limit,
    retrieved_velocity_std: (0.0, 15.0),
    velocity_difference: (-15.0, 15.0),
    wind_speed: (0.0, 50.0),
    wind_direction: (0.0, 360.0),
    spectrum_width: spectrum_width_limit,
    corrected_spectrum_width: spectrum_width_limit,
    unfiltered_spectrum_width: spectrum_width_limit,
    spectrum_width_vv: spectrum_width_limit,
    unfiltered_spectrum_width_vv: spectrum_width_limit,
    turbulence: (0.0, 0.2),
    normalized_coherent_power: (0.0, 1.0),
    # differential_reflectivity: (-1., 8.),
    differential_reflectivity: (-8.0, 12.0),
    corrected_differential_reflectivity: (-1.0, 8.0),
    # unfiltered_differential_reflectivity: (-1., 8.),
    unfiltered_differential_reflectivity: (-8.0, 12.0),
    differential_reflectivity_in_precipitation: (-10.0, 10.0),
    differential_reflectivity_in_snow: (-10.0, 10.0),
    differential_reflectivity_column_height: (0.0, 6.0),
    spectral_differential_reflectivity: (-1.0, 8.0),
    unfiltered_spectral_differential_reflectivity: (-1.0, 8.0),
    # cross_correlation_ratio: (0.7, 1.),
    cross_correlation_ratio: (0.0, 1.0),
    corrected_cross_correlation_ratio: (0.7, 1.0),
    # unfiltered_cross_correlation_ratio: (0.7, 1.),
    unfiltered_cross_correlation_ratio: (0.0, 1.0),
    uncorrected_cross_correlation_ratio: (0.7, 1.0),
    logarithmic_cross_correlation_ratio: (0, 4),
    cross_correlation_ratio_in_rain: (0.9, 1.0),
    # spectral_copolar_correlation_coefficient: (0.7, 1),
    # unfiltered_spectral_copolar_correlation_coefficient: (0.7, 1),
    mean_phase: (-180.0, 180.0),
    differential_phase: (-180.0, 180.0),
    unfolded_differential_phase: (-180.0, 180.0),
    corrected_differential_phase: (-180.0, 180.0),
    uncorrected_differential_phase: (-180, 180.0),
    uncorrected_unfiltered_differential_phase: (-180, 180.0),
    system_differential_phase: (-180.0, 180.0),
    spectral_phase_hh: (-180.0, 180.0),
    spectral_phase_vv: (-180.0, 180.0),
    spectral_differential_phase: (-180.0, 180.0),
    unfiltered_spectral_phase_hh: (-180.0, 180.0),
    unfiltered_spectral_phase_vv: (-180.0, 180.0),
    unfiltered_spectral_differential_phase: (-180.0, 180.0),
    specific_differential_phase: (-1.0, 2.0),
    corrected_specific_differential_phase: (-1.0, 2.0),
    uncorrected_specific_differential_phase: (-1.0, 2.0),
    uncorrected_unfiltered_specific_differential_phase: (-1.0, 2.0),
    linear_depolarization_ratio: (-40.0, 0.0),
    linear_depolarization_ratio_h: (-40.0, 0.0),
    linear_depolarization_ratio_v: (-40.0, 0.0),
    circular_depolarization_ratio: (-40.0, 0.0),
    rain_rate: (0.0, 10.0),
    radar_estimated_rain_rate: (0.02, 10.0),
    rainfall_accumulation: (0.0, 100.0),
    radar_echo_classification: (0.0, 9.0),
    hydroclass_entropy: (0.0, 1.0),
    proportion_AG: (0.0, 100.0),
    proportion_CR: (0.0, 100.0),
    proportion_LR: (0.0, 100.0),
    proportion_RP: (0.0, 100.0),
    proportion_RN: (0.0, 100.0),
    proportion_VI: (0.0, 100.0),
    proportion_WS: (0.0, 100.0),
    proportion_MH: (0.0, 100.0),
    proportion_IH: (0.0, 100.0),
    radar_echo_id: (0, 3),
    melting_layer: (0, 5),
    clutter_exit_code: (0, 200),
    echo_top_15dBZ: (0, 10),
    echo_top_20dBZ: (0, 10),
    echo_top_45dBZ: (0, 10),
    echo_top_50dBZ: (0, 10),
    probability_of_hail: (0.0, 100.0),
    vertically_integrated_liquid: (0.0, 30.0),
    sun_hit_h: (0, 1),
    sun_hit_v: (0, 1),
    sun_hit_zdr: (0, 1),
    specific_attenuation: (0.0, 1.0),
    path_integrated_attenuation: (0.0, 20.0),
    specific_differential_attenuation: (0.0, 0.3),
    path_integrated_differential_attenuation: (0.0, 3.0),
    corrected_specific_attenuation: (0.0, 1.0),
    corrected_path_integrated_attenuation: (0.0, 20.0),
    corrected_specific_differential_attenuation: (0.0, 0.3),
    corrected_path_integrated_differential_attenuation: (0.0, 3.0),
    differential_phase_texture: (0, 180.0),
    height: (0, 20000),
    interpolated_profile: (0, 10000),
    visibility: (0, 100),
    frequency_of_occurrence: (0, 100),
    temperature: (-60, 30),
    height_over_iso0: (-6000.0, 6000.0),
    iso0_height: (0.0, 5000.0),
    # Additional reflectivity like fields
    "CZ": (-10.0, 65.0),
    "DZ": (-10.0, 65.0),
    "AZ": (-10.0, 65.0),
    "Z": (-10.0, 65.0),
    "dbz": (-10.0, 65.0),
    "DBZ": (-10.0, 65.0),
    "dBZ": (-10.0, 65.0),
    "DBZH": (-10.0, 65.0),
    "DBZ_S": (-10.0, 65.0),
    "DBZ_K": (-10.0, 65.0),
    "reflectivity_horizontal": (-10.0, 65.0),
    "corr_reflectivity": (-10.0, 65.0),
    HRV: (0.0, 100.0),
    VIS006: (0.0, 85.0),
    VIS008: (0.0, 90.0),
    IR_016: (0.0, 30.0),
    IR_039: (210.0, 340.0),
    WV_062: (210.0, 260.0),
    WV_073: (190.0, 280.0),
    IR_087: (205.0, 320.0),
    IR_097: (210.0, 285.0),
    IR_108: (205.0, 320.0),
    IR_120: (205.0, 320.0),
    IR_134: (205.0, 290.0),
    CTH: (0.0, 12000.0),
    HRV_norm: (0.0, 100.0),
    VIS006_norm: (0.0, 85.0),
    VIS008_norm: (0.0, 90.0),
    IR_016_norm: (0.0, 30.0),
    # radial_wind_speed: 'BuDRd18',
    # radial_wind_speed_ci:
    # radial_wind_speed_status:
    # doppler_spectrum_width: 'NWS_SPW',
    # doppler_spectrum_mean_error:
    # atmospherical_structures_type: 'LangRainbow12',
    # relative_beta: (0., 4e-6),
    # absolute_beta: (0., 4e-6),
    # cnr: 'Carbone17',
}
