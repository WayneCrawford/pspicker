# ps_picker

Seismological P- and S- wave picker using the modified Kurtosis method

Python port of the picker described in Baillard et al., 2014

debugging information is saved to the local file run_{datetime}.log

## Methodology
The picker is based around the Kurtosis, but also uses energy levels, polarity,
clustering and phase association in a 3-step process:

### Step 1: define a global pick window

The *Kurtosis* is calculated for all stations.  The global window
surrounds the most densely clustered region of triggers.

### Step 2: pick P and S arrivals on each station individually

For each station:
    - calculate the *Kurtosis* over coarse to fine scales.
    - Identify candidates on the coarse scale and refine their times using
      the finier scales
    - Choose P- and S- candidates based on the *signal-to-noise level* of
      each pick
    - Verify the candidates using the waveform *polarity*, if possible
       - polarity is only used if one of the picks has a dip of > 30 degrees

### Step 3: associate picks
    - Calculate origin times for each trace, based on the P-S delay and
      a simple velocity model (could I use a single Vp/Vs value?)
    - If at least 3 origin times are clustered, use their average origin time
      to validate all candidates, possibly dipping into the pool of unused
      candidates for replacemene P and S picks
    - If less than 3 origin times are clustered, reject bad P- and S- picks
      based on clustering of P-pick times, S-pick times and P-S delays

## Code and parameter file examples
Are located [here](code_examples.md):


## Example workflow

### Start by autopicking a few events, with all the bells and whistles on:

To pick one event from a database in `/SEISAN/MAYOBS`:

    from ps_picker import PSPicker

    picker = PSPicker('parameters_C.yaml', '/SEISAN/MAYOBS/WAV/MAYOB', 
                      '/SEISAN/MAYOBS/REA/MAYOB')
    picker.run_one('19-0607-59L.S201905', plot_global=True, plot_stations=True,
                   verbose=True)

Look at all of the plots and verify that the picks and association are as
you expect.  If not, change the paramters and run again.

### Next, pick several events with only the global plots on

The bells and whistles text will be saved to a log file named
run_{DATETIME}.log

To pick events from May 5th to 25th in the same database:

    from ps_picker import PSPicker

    picker = PSPicker('parameters_C.yaml', '/SEISAN/MAYOBS/WAV/MAYOB', 
                      '/SEISAN/MAYOBS/REA/MAYOB')
    picker.run_many('20190505', '20190525', plot_global=True)

### Finally, run the whole database without plots

*(run_{DATETIME}.log is always created)*

To pick events from May 26th 2019 May 1st 2020:

    from ps_picker import PSPicker

    picker = PSPicker('parameters_C.yaml', '/SEISAN/MAYOBS/WAV/MAYOB', 
                      '/SEISAN/MAYOBS/REA/MAYOB')
    picker.run_many('20190526', '20200501')

## Parameters
Picker parameters are passed in a
[YAML](https://tools.ietf.org/id/draft-pbryan-zyp-json-ref-03.html) file with
the following fields (fields with values shown have defaults and are not
required in the file):

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
        quality_thresholds:        # [4-list] of SNR levels associated with quality levels '3', '2', '1' and '0'
        threshold_parameter: 0.2   # Controls the SNR_threshold for SNR-based quality evaluation
                                   # if between 0 and 1, then SNR_threshold = max(SNR)*threshold_parameter
                                   # if < 0, then SNR_threshold = -threshold_parameter
        max_threshold_crossings: 2 # Maximum allowed crossings of SNR threshold within global window
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
        DR_smooth_length: 1.  # smoothing window to apply to dip and rectilinearity when calculating DR
        calculate_window: 2.  # number of seconds after a pick over which to calculate dip_rect
        analyze_window: 4.    # number of seconds around a calc point to calculate polarity
    association: # Parameters affecting the association between different stations
        cluster_window_otime:  # Window length in seconds for cluster-based rejection of origin times
        otime_vp_vs: 1.75      # Vp/Vs value to use for origin time calculations
        cluster_window_P:      # Window length in seconds for cluster-based rejection of P arrivals
        cluster_window_S:      # Window length in seconds for cluster-based rejection of S arrivals
        distri_min_values: 4   # minimum number of values (P picks, S picks, or PS-times) needed for distribution-based rejection
        distri_nstd_picks: 3.2 # reject picks outside of this number of standard deviations
        distri_nstd_delays: 4  # reject delays outside of this number of standard deviations
    response_filetype: '' # 'GSE' or '': the latter means a Baillard PoleZeros-type format
    station_parameters:  # List of objects with key = station_type
        - station_type1
            P_comp:                  # components (one letter each, selected from 'ZNEH') to use for P-picks
            S_comp:                  # components (one letter each, selected from 'ZNEH') to use for S-picks
            energy_frequency_band:   # frequency band [low, high] used for SNR and energy calculations
            energy_window:           # only look at data from t-nrg_win to t when evaluating energy, where t is the time of the peak waveform energy.
                                     # If == 0, don't use energy criteria.
            kurt_frequency bands:    # Kurtosis list of frequency bands over which to run Kurtosis, e.g.[[3, 15], [8, 30]]
            kurt_window_lengths:     # Kurtosis list of window lengths in seconds, e.g. [0.3, 0.5, 1, 2, 4, 8]
            kurt_extrema_smoothings: # Kurtosis list of smoothing sequences in samples, e.g. [2, 4, 6, 8, 10, 20, 30, 40, 50]
            use_polarity:            # Use polarities (dip_rect thresholds) to assign P and S picks
            n_extrema: 5             # number of candidates to pick (a big number allows alternate candidates)
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

## To Do

- Add event location-based acceptance of solitary P- and S- candidates
- In P-, S- and P-S clustering stage, allow unused candidates to be
  substituted for rejected picks
- More in [ToDo.md](ToDo.md)
    
Also see the [profiling file](profiling.md)
