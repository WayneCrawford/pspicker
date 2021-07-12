Critical
------------------------

- Make origin-time association works

  - Even if I increase cluster_window_otime to 5 seconds, doesn't increase
    number of events selected!
    
Non-critical
------------------------

- Allow read from SDS (or other client) rather than SEISAN structure.  Changes
  PSPicker arguments:
    - rename:
        - ``database_path_in`` to ``catalog_in`` (file rather than directory)
        - ``database_path_out`` to ``catalog_out`` (file appended to after each
          event?, or stick with a directory?)
        - ``wav_base_path`` to ``wave_db_base``. Can be a URL or string
          (for FDSN) or a path (for SEISAN or SDS)
    - add:
        - catalog_in_format (allow read_events arguments)
        - catalog_out_format (allow Catalog.write arguments)
        - wave_db_format ('SEISAN', 'SDS', 'FDSN')

- Integrate station name wildcard search into station_parameters list object
  (currently runs search each time, in pspicker.py)

- Read responses before running (instead of each time and at the end)

- Allow catalog reads (automatically link station names to their responses)

- Allow existing manual picks to be used (and kept as basis for other picks)

- Allow manual specification of the global window

- Check that the polarity analysis works like that in Ppol (and make a test case)

- Change smoothing parameters to seconds rather than samples? Or set a
  max_seconds for the smoothing parameters?

- Allow user to set ChannelMappingRules and make them more generic
  (by default, only select based on the component code).  If more
  than one channel is selected, keep the one with the highest sample rate.
  If two with highest sample rate choose least broadband(?)

Code cleanup
------------------------

- Clean up picker_parameters.from_yaml_file()!

-  Get rid of all "isnan" parts in code (remnant of Matlab version that didn't
   cut traces but just set values to NaN
   
   - get rid of nan testing everywhere and stop using np.nanmax() and np.nanargmax())

### ChannelMapping

- Integrate "ChannelMappingRules" into "ChannelMapping"?
- Make select_traces() part of ChannelMapping?