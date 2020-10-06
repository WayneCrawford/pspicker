# ps_picker

Seismological P- and S- wave picker using the modified Kurtosis method

Python port of the picker described in Baillard et al., 2014 

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
    kurtosis:
        frequency_band:       # cutoff frequencies [low, high] for kurtosis calculation
        window_length:        # sliding window length in seconds for kurtosis calculation
        extrema_smoothing: 40 # number of samples to smooth extrema by when looking for pick
    distri_secs:        # size of window in seconds in which to look for the maximum # of picks
    offsets:            # final window offset in seconds [left, right] from peak distribution
    end_cutoff: 0.9     # don't look for extrema beyond this fraction of the overall time
    n_picks: 5          # number of picks to use for each trace
SNR: # Parameters affecting the signal-to-noise level calculation and use
    noise_window_length:       # seconds to use for noise window
    signal_window_length:      # seconds to use for signal_window
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
dip_rect_thresholds: # minimum rectilinearity thresholds needed to assign 'P' or 'S' to an onset (P positive, S negative)
    P: 0.4
    S: -0.4
kurtosis: # Parameters affecting Kurtosis calculations (except in inital global window selection)
    frequency_bands:     # object with one or more "keys", each followed by a list of frequency bands over which to run Kurtosis
                         # e.g. {A: [[3, 15], [8, 30]]}
    window_lengths:      # object with one or more "keys", each followed by a list of window lengths in seconds
                         # e.g. {A: [0.3, 0.5, 1, 2, 4, 8]}
    extrema_smoothings:  # object with one or more "keys", each followed by a list of smoothing sequences in samples
                         # e.g. {A: [2, 4, 6, 8, 10, 20, 30, 40, 50]}
association: # Parameters affecting the association between different stations
    cluster_windows_P:     # Window length in seconds for cluster-based rejection of P arrivals
    cluster_windows_S:     # Window length in seconds for cluster-based rejection of S arrivals
    distri_min_values: 4   # minimum number of values (P picks, S picks, or PS-times) needed for distribution-based rejection
    distri_nstd_picks: 3.2 # reject picks outside of this number of standard deviations
    distri_nstd_delays: 4  # reject delays outside of this number of standard deviations
responses:
    filetype: '' # 'GSE' or '': the latter means a Baillard PoleZeros-type format
    filenames:   # object with one or more "keys", each followed by a filename
                 # e.g. {A: 'SPOBS2_response.txt', B: 'micrOBS_G1_response.txt'}
station_parameters:  # List of objects with key = station_name
    - station1_name
        P_comp:   # string of all components (one letter each) used for P-picks
        S_comp:   # string of all components (one letter each) used for S-picks
        f_nrg:    # frequency band [low, high] used for SNR and energy calculations
        k_parms:  # Kurtosis parameters
            freqs:   # key from kurtosis:frequency_bands
            wind:    # key from kurtosis:window_lengths
            smooth:  # key from kurtosis:extrema_smoothings
        polar:    # Use polarities (dip_rect thresholds) to assign P and S picks
        nrg_win:  # only look at data from t-nrg_win to t when evaluating energy, where t is the time of the peak waveform energy.
                  # If == 0, don't use energy criteria.
        n_follow: # number of extrema to follow (1 or 2).  Generally use 2 (S and P) unless data are problematic
        resp:     # key from responses:filename
    - station2_name
      ...
    - station3_name
      ...
    ...
```

## More information

[TO DO](ToDo.rst)
[JSONref](https://tools.ietf.org/id/draft-pbryan-zyp-json-ref-03.html)
