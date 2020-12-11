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
