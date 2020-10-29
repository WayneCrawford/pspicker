# ps_picker

Seismological P- and S- wave picker using the modified Kurtosis method

Python port of the picker described in Baillard et al., 2014 

## Methodology
The picker is based around the Kurtosis, but also uses energy levels, polarity,
clustering and phase association in a 3-step process:

### Step 1: define a global pick window

The *Kurtosis* is calculated for all stations.  The global window
surrounds the most densely clustered region of triggers.

### Step 2: pick P and S arrivals on each station individually

For each station:
    - calculate the *Kurtosi*s over coarse to fine scales.
    - Identify candidates on the coarse scale and refine their times using
      the finier scales
    - Choose P- and S- candidates based on the *signal-to-noise level* of
      each pick
    - Verify the candidates using the waveform *polarity*, if possible
       * polarity is only used if one of the picks has a dip of > 30 degrees

### Step 3: associate picks
    - Calculate origin times for each trace, based on the P-S delay and
      a simple velocity model (could I use a single Vp/Vs value?)
    - If at least 3 origin times are clustered, use their average origin time
      to validate all candidates, possibly dipping into the pool of unused
      candidates for replacemene P and S picks
    - If less than 3 origin times are clustered, reject bad P- and S- picks
      based on clustering of P-pick times, S-pick times and P-S delays

## Examples

To pick one event from a database in `/SEISAN/MAYOBS`:
```python
from ps_picker import PSPicker

picker = PSPicker('/SEISAN/MAYOBS/WAV/MAYOB', 'parameters_C.yaml')
picker.run_one('/SEISAN/MAYOBS/REA/MAYOB/2019/05/19-0607-59L.S201905')
```
To pick events from May to September 2019 in the same database:
```python
from ps_picker import PSPicker

picker = PSPicker('/SEISAN/MAYOBS/WAV/MAYOB/', 'my_params.yaml')
picker.run_many('/SEISAN/MAYOBS/REA/MAYOB/', '201905', '201909')
```
## Parameters
Picker parameters are passed in a
[YAML](https://tools.ietf.org/id/draft-pbryan-zyp-json-ref-03.html) file with
the following fields (fields with values shown have defaults and are not
required in the file):
```yaml
---
global_window: # Parameters affecting the initial selection of a global pick window across all stations using the distribution of kurtosis extrema)
    kurt_frequency_band:       # Kurtosis cutoff frequencies [low, high] for kurtosis calculation
    kurt_window_length:        # Kurtosis sliding window length in seconds for kurtosis calculation
    kurt_extrema_smoothing: 40 # Kurtosis number of samples to smooth extrema by when looking for pick
    distri_secs:        # size of window in seconds in which to look for the maximum # of picks
    offsets:            # final window offset in seconds [left, right] from peak distribution
    end_cutoff: 0.9     # don't look for extrema beyond this fraction of the overall time
    n_extrema: 5        # number of kurtosis extrema to pick for each trace
SNR: # Parameters affecting the signal-to-noise level calculation and use
    noise_window:              # seconds to use for noise window
    signal_window:             # seconds to use for signal_window
    max_threshold_crossings: 2 # Maximum allowed crossings of SNR threshold within global window
    quality_thresholds:        # [4-list] of SNR levels associated with quality levels '3', '2', '1' and '0'
    threshold_parameter: 0.2   # Controls the SNR_threshold for SNR-based quality evaluation
                               # if between 0 and 1, then SNR_threshold = max(SNR)*threshold_parameter
                               # if < 0, then SNR_threshold = -threshold_parameter
channel_parameters: # Parameters affecting the choice of channels to pick on and save to
    compZ: 'Z3'               # Component names that will be interpreted as 'Z'
    compN: 'N1Y'              # Component names that will be interpreted as 'N'
    compE: 'E2X'              # Component names that will be interpreted as 'E'
    compH: 'HF'               # Component names that will be interpreted as 'H'
    S_write_cmp: 'N'          # Assign S picks to this component (or equivalent as defined above)
    P_write_cmp: 'Z'          # Assign P picks to this component (or equivalent as defined above)
    P_write_phase: 'Pg'       # Give this phase hint to P picks
    S_write_phase: 'Sg'       # Give this phase hint to S picks
    band_order: 'GFDCEHSBMLV' # If multiple traces have the same component, chose the one with the earliest listed band code
                              # 'GFDCEHSBMLV' prioritizes high sampling rates over low, and short period over broadband
polarity: # polarity analyses parameters (mostly related to dip_rect, or DR, see Baillard et al 2014)
    DR_threshold_P: 0.4   # minimum DR to assign 'P'
    DR_threshold_S: -0.4  # maximum DR to assign 'S'
    calculate_window: 2.  # number of seconds after a pick over which to calculate dip_rect
    analyze_window: 4.    # number of seconds around a calc point to calculate polarity
    smooth_length: 1.     # smoothing window to apply to dip and rectilinearity when calculating DR
association: # Parameters affecting the association between different stations
    cluster_windows_P:     # Window length in seconds for cluster-based rejection of P arrivals
    cluster_windows_S:     # Window length in seconds for cluster-based rejection of S arrivals
    distri_min_values: 4   # minimum number of values (P picks, S picks, or PS-times) needed for distribution-based rejection
    distri_nstd_picks: 3.2 # reject picks outside of this number of standard deviations
    distri_nstd_delays: 4  # reject delays outside of this number of standard deviations
response_filetype: '' # 'GSE' or '': the latter means a Baillard PoleZeros-type format
station_parameters:  # List of objects with key = station_type
    - station_type1
        P_comp:    # string of all components (one letter each, selected from 'ZNEH') used for P-picks
        S_comp:    # string of all components (one letter each, selected from 'ZNEH') used for S-picks
        f_energy:  # frequency band [low, high] used for SNR and energy calculations
        kurt_frequency bands:    # Kurtosis list of frequency bands over which to run Kurtosis, e.g.[[3, 15], [8, 30]]
        kurt_window_lengths:     # Kurtosis list of window lengths in seconds, e.g. [0.3, 0.5, 1, 2, 4, 8]
        kurt_extrema_smoothings: # Kurtosis list of smoothing sequences in samples, e.g. [2, 4, 6, 8, 10, 20, 30, 40, 50]
        use_polarity:    # Use polarities (dip_rect thresholds) to assign P and S picks
        nrg_win:  # only look at data from t-nrg_win to t when evaluating energy, where t is the time of the peak waveform energy.
                  # If == 0, don't use energy criteria.
        n_extrema: 5 # number of extrema to follow
    - station2_name
      ...
    - station3_name
      ...
    ...
stations:  # List of stations with their station_parameters and responsefiles
    station1_name: {parameters: "station_typeN", response: "responsefilename"}
    station2_name: {parameters: "station_typeM", response: "responsefilename"}
    station2_name: {parameters: "station_typeM", response: "responsefilename"}
    ...    
```

## To Do

    - Make associator switch to cluster-based if origin-time based doesn't work
        - and make sure origin-time based parameters are in paramter file
    - Fix magnitude calculations
        - Allow SAC PZ files for response files
    - Figure out why S-picks aren't being saved
    - Put pick uncertainties into Nordic files
        - short-term solution: by creating an associated arrival and setting
          its time_weight to the appropriate number
        - long-term solution: integrating my recommended time-weight procedure
          into obspy nordic:
            - add a parameter uncertainty-mapping to write_select. This could
              be a list of the maximum time_errors for a given nordic_weight,
              for example [0.1, 0.2, 0.4, 0.8] would give a weight of "0" for
              uncertainties less than 0.1s, "1" for 0.1-0.2s, "2" for 0.2-0.4,
              "3" for 0.4-0.8 and "4" for the rest. No weight is given for
              a pick without a time_error
            - if this parameter is not set, then the pick weight will be based
              on the Arrival.time_weight according to the following formula:
                ```
                if all(time_weight <= 1): 
                    nordic_weight = round(4. * (1-time_weight))  # alternatively floor(4.9999. * (1-time_weight))
                else:
                    nordic_weight = round(4. * (max(time_weights) - 1)  # or the alternate form
            - Should I also make nordic outputs include I lines (all in one
              integrated variable including waveform filenames)?

    - Add event location-based acceptance of solitary P- and S- candidates
    - In P-, S- and P-S clustering stage, allow unused candidates to be
      substituted for rejected picks