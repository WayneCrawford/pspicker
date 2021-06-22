
The parameter file
========================================

The parameter file is written in 
`YAML <https://tools.ietf.org/id/draft-pbryan-zyp-json-ref-03.html>`_
If the program can't read your's, try using an online YAML validator, like
`this one <https://codebeautify.org/yaml-validator>`_

Example file
-------------------------------

```yaml

    ---
    global_window:
        kurtosis:
            frequency_bands: [[5, 30]]
            window_lengths: [20]
        distri_secs: 5
        offsets: [-10, 10]
        end_cutoff: 0.9
        max_candidates: 5
    SNR:
        noise_window: 2.
        signal_window: 1.
        quality_thresholds: [1.5, 2.5, 4, 6]
        threshold_parameter: -3.
    polarity:
        calculate_window: 1.
        analyze_window: 1.
    association:
        cluster_window_otime: 1.
        otime_vp_vs: 1.70
        cluster_window_P: 3.
        cluster_window_S: 5.
        cluster_window_otime: 1.
    station_parameters:
        SPOBS:
            picking_components:
                P: 'Z'
                S: 'ZNE'
            SNR_energy:
                frequency_band: [3, 30]
                window: 20
            kurtosis:
                frequency_bands: [[3, 15], [8, 30]]
                window_lengths: [0.3, 0.5, 1, 2, 4, 8]
                extrema_smoothings: [2, 4, 6, 8, 10, 20, 30, 40, 50]
            use_polarity: true
        BBLAND:
            picking_components:
                P: 'Z'
                S: 'ZNE'
            SNR_energy:
                frequency_band: [3, 30]
                window: 20
            kurtosis:
                frequency_bands: [[3, 15], [8, 30]]
                window_lengths: [0.3, 0.5, 1, 2, 4, 8]
                extrema_smoothings: [2, 4, 6, 8, 10, 20, 30, 40, 50]
            use_polarity: true
    stations:
        MOCA: {parameters: 'SPOBS', resp_file: 'SPOBS2_response.txt'}
        MOFA: {parameters: 'SPOBS', resp_file: 'SPOBS2_response.txt'}
        MONA: {parameters: 'SPOBS', resp_file: 'SPOBS2_response.txt'}
        MODA: {parameters: 'SPOBS', resp_file: 'SPOBS2_response.txt'}
        MOSA: {parameters: 'SPOBS', resp_file: 'SPOBS2_response.txt'}
        MOVA: {parameters: 'SPOBS', resp_file: 'SPOBS2_response.txt'}
        IF1A: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF2A: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF3A: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF4A: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF5A: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF6A: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF7A: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF8A: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF1B: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF2B: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF3B: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF4B: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF5B: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF6B: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF7B: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        IF8B: {parameters: 'SPOBS', resp_file: 'micrOBS_G1_response.txt'}
        KNKL: {parameters: 'BBLAND', resp_file: 'KNKL_BBOBS1_1.response.txt'}
```

A description of every line
-------------------------------

The values provided on some lines are defaults.  If you don't want
to change them, you don't have to include them in your parameter file.

```yaml

    ---
    global_window: # Parameters affecting the initial selection of a global pick window across all stations using the distribution of kurtosis extrema)
        kurtosis:           # Kurtosis calculation parameters
            frequency_bands:         # list of frequency bands [low, high] to use
            window_lengths:          # list of sliding window lengths to use
            extrema_smoothings: [40] # list of number of samples to smooth extrema by when looking for pick
        distri_secs:        # size of window in seconds in which to look for the maximum # of picks
        offsets:            # final window offset in seconds [left, right] from peak distribution
        end_cutoff: 0.9     # don't look for extrema beyond this fraction of the overall time
        max_candidates: 5   # maxium number of pick candidates for each trace
    SNR: # Parameters affecting the signal-to-noise level calculation and use
        noise_window:              # seconds to use for noise window
        signal_window:             # seconds to use for signal_window
        quality_thresholds:        # [4-list] of SNR levels associated with quality levels '3', '2', '1' and '0'
        threshold_parameter: 0.2   # Controls the SNR_threshold for SNR-based quality evaluation
                                   # if between 0 and 1: SNR_threshold = max(max(SNR)*threshold_parameter, quality_thresholds[0])
                                   # if < 0:  SNR_threshold = -threshold_parameter
        max_threshold_crossings: 5 # Maximum allowed crossings of SNR threshold within global window
    channel_parameters: # Parameters affecting the choice of channels to pick on and save to
        component_orientation_codes: # Component names to map to Z, N, E and H
            Z: 'Z3'
            N: 'N1Y'
            E: 'E2X'
            H: 'HF'
        write_components_phases:  # Channel to write to ('Z', 'N', 'E' or 'H', mapped as above)
                                  # and phase name to give for picks
            S: ['N', 'Sg']  # S-wave picks
            S: ['Z', 'Pg']  # P-wave picks
        band_order: 'GFDCEHSBMLV' # If multiple traces have the same component, chose the one with the earliest listed band code
                                  # 'GFDCEHSBMLV' prioritizes high sampling rates over low, and short period over broadband
    polarity: # polarity analyses parameters (mostly related to dip_rect, or DR, see Baillard et al 2014)
        DR_threshold_P: 0.4   # minimum DR to assign 'P'
        DR_threshold_S: -0.4  # maximum DR to assign 'S'
        DR_smooth_length: 1.  # smoothing window to apply to dip and rectilinearity when calculating DR
        calculate_window: 2.  # number of seconds after a pick over which to calculate dip_rect
        analyze_window: 4.    # number of seconds around a calc point to calculate polarity
    association: # Parameters affecting the association between different stations
        method: 'origin_time'  # Preferred association method: ['origin_time', 'arrival_time']
        cluster_window_otime:  # Window length in seconds for cluster-based rejection of origin times
        otime_vp_vs: 1.75      # Vp/Vs value to use for origin time calculations
        cluster_window_P:      # Window length in seconds for cluster-based rejection of P arrivals
        cluster_window_S:      # Window length in seconds for cluster-based rejection of S arrivals
        distri_min_values: 4   # minimum number of values (P picks, S picks, or PS-times) needed for distribution-based rejection
        distri_nstd_picks: 3.2 # reject picks outside of this number of standard deviations
        distri_nstd_delays: 4  # reject delays outside of this number of standard deviations
    response_file_type: ''  # 'GSE', 'SACPZ', 'JSON_PZ', 'STATIONXML' or '': the latter means Baillard PoleZero format
    station_parameters:  # List of objects with key = station_type
        - station_type1
            picking_components:  # components to use for picks (selected from 'ZNEH')
                P: 'Z'           # P-picks
                S: 'ZNE'         # S-picks
            SNR_energy:
                frequency_band:   # frequency band [low, high] for SNR and energy calculations
                window:           # only look at data from t-nrg_win to t when evaluating energy, where t is the time of the peak waveform energy.
                                     # If == 0, don't use energy criteria.
            kurtosis:                # Kurtosis calculation parameters
                frequency_bands:         # list of frequency bands [low, high] to use
                window_lengths:          # list of sliding window lengths to use
                extrema_smoothings:      # list of number of samples to smooth extrema by when looking for pick
            kurt_frequency bands:    # Kurtosis list of frequency bands over which to run Kurtosis, e.g.[[3, 15], [8, 30]]
            kurt_window_lengths:     # Kurtosis list of window lengths in seconds, e.g. [0.3, 0.5, 1, 2, 4, 8]
            kurt_extrema_smoothings: # Kurtosis list of smoothing sequences in samples, e.g. [2, 4, 6, 8, 10, 20, 30, 40, 50]
            use_polarity:            # Use polarities (dip_rect thresholds) to assign P and S picks
            max_candidates: 5        # number of candidates to pick (a big number allows alternate candidates)
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

Response files
========================================


Response files are needed to calculate local amplitudes.  They can be provided
in 'SACPZ', 'GSE2', "Baillard" or 'JSON' format (the latter is just a cleaner
version of the "Baillard' format).  The code either reads an absolute gain
or calculates it from  a "passband" gain provided at a reference frequency.  
The "A0" normalization constant, needed to calculate the absolute gain from the
passband gain, is directly calculated from the poles and zeros such that A0 times
the pole-zero formula equals 1.0 at the reference frequency. The
parameters used for each format are:

 format   | gain     | passband gain @ ref_freq | poles | zeros | input units |
----------|----------|--------------------------|-------|-------|-------------|
 SAC PZ   | CONSTANT |                          | POLES | ZEROS |  meters     |
 GSE2     |          |  1/sensitivity at f_ref  (values from CAL2 line) | poles | zeros |  nm         |
 Baillard |          | 1/sensitivity (line 1) at f_req(line 2)  | poles | zeros |  nm         |
 JSON     |          | 1/sensitivity at f_ref   | poles | zeros |  nm         |
