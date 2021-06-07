Critical
------------------------

- Make origin-time association works

  - Even if I increase cluster_window_otime to 5 seconds, doesn't increase
    number of events selected!

Non-critical
------------------------

- Read responses before running (instead of each time and at the end)

- Make converter to different PAZ formats

- Make run_many create a single log file (looks like run_one makes a new one
  each time)

- Allow wildcards in station names (would simplify Mayotte to ``MO*``, ``MOOH``,
  ``IF*``, ``IF*E``, ``KNKL``, ``PMZI``, ``MTSB``)

  - Needs to intelligently handle conflicting station names (taking the most
    specific)

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

- Replace os.path calls with pathlib

-  Get rid of all "isnan" parts in code (remnant of Matlab version that didn't
   cut traces but just set values to NaN
   
   - get rid of nan testing everywhere and stop using np.nanmax() and np.nanargmax())

### ChannelMapping

- Integrate "ChannelMappingRules" into "ChannelMapping"?
- Make select_traces() part of ChannelMapping?
    
### Parameter file

Replace:

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

by:

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

Replace:

```yaml

    station_parameters:
        SPOBS:
            P_comp: 'Z'
            S_comp: 'ZNE'
            energy_frequency_band: [3, 30]
            energy_window: 20  # What does this really do?
```

by:

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
