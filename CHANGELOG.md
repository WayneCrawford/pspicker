# Version 0.2:

- Changed date representation in run_many()
- Removed asctime and caller from logging.basicConfig format, because matplotlib
  calls to logging.debug do not fill in these values, which gives an error
- added VERBOSE logging level (using verboselogs package)
- Fixed bug that set all amplitude pick components (and associated phase pick)
  to the amplitude level and component of the last picked station
  
## v0.2.2:
- Added `associator:method` attribute (can use to specify 'arrival_time'
  instead of the default 'origin_time'
- Replaced `verbose` argument to run_one() and run_many() by `log_level`
  (values=('info', 'verbose', or 'debug'))
- Fixed a bug that could create a false kurtosis "hit" 50 points into a
  time series

# v0.3

- Parameter file changes:

    - replaced ``kurt_frequency_band(s)``, ``kurt_window_length(s)``
      and ``kurt_extrema_smoothing(s)`` to
      ``kurtosis:{frequency_bands, window_lengths, extrema_smoothings``
    - renamed ``n_extrema`` to ``max_candidates``

- Changed some internal kurtosis code to simplify and better align with
  Matlab code outputs
- Fixed a (probably unimportant) bug that kept kurtosis gradients from 
  reaching zero
- added a prepend and append to second np.diff in kurtosis._loca_ext(),
  to retain array length and index alignment
- simplified smooth_filter, eliminating a bug that put zero in the first index

## v0.3.1:

cleaned up parameter file reading, fixed log_level use in run_many()

# v0.4

- json_pz files now use m/s as default units
- Improved PAZ class
- Sorts input files by filename within each day

## v0.4.1:
 - Updated pyyaml and fixed some pathlib errors
 
# v0.5:

- Updated PAZ to have separate initialization routines for pole-zero gain and
bandpass gain.

- Updated parameter file format:

```yaml
  channel_parameters:
    compZ: 'Z3'
    compN: 'N1Y'
    compE: 'E2X'
    compH: 'HF'
    S_write_cmp: 'N'
    P_write_cmp: 'Z'
    P_write_phase: 'Pg'
    S_write_phase: 'Sg'
```

becomes:

```yaml
  channel_parameters:
    component_orientation_codes:
        Z: 'Z3'
        N: 'N1Y'
        E: 'E2X'
        H: 'HF'
    write_components_phases:
        S: ['N', 'Sg']
        P: ['Z', 'Pg']
```

```yaml

    station_parameters:
        SPOBS:
            P_comp: 'Z'
            S_comp: 'ZNE'
            energy_frequency_band: [3, 30]
            energy_window: 20  # What does this really do?
```

becomes

```yaml

    station_parameters:
        SPOBS:
            picking_components
                P: 'Z'
                S: 'ZNE'
            SNR_energy:
                frequency_band: [3, 30]
                window: 20  # What does this really do?
```


- Got rid of sampling_rate attribute, a relic of C Baillard's code 
  which doesn't belong in a Pole-Zero

- Changed behavior of setting input_unit to an unknown value: previously
  accepted the name, now keeps the old name if it can't convert units (was
  a problem when saving unknown units to json_pz files, which forced units to 
  m/s, even if they were Pascals)
  
- Changed fast_kurtosis to reduce edge effect without throwing away first
  win_samples points, by calculating using a win_secs length buffer of repeating
  first second of data_
## v0.5.1:
 - Fixed bugs(?) in amplitude calculation
 
