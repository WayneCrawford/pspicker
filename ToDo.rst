Critical
------------------------

- Make origin-time association works

  - Even if I increase cluster_window_otime to 5 seconds, doesn't increase
    number of events selected!

Non-critical
------------------------

- Allow catalog reads (automatically link station names to their responses)

- Allow existing manual picks to be used (and kept as basis for other picks)

- Allow manual specification of the global window

- Check that the polarity analysis works like that in Ppol (and make a test case)

- Change smoothing parameters to seconds rather than samples? Or set a
  max_seconds for the smoothing parameters?

-  Get rid of all "isnan" parts in code (remnant of Matlab version that didn't
   cut traces but just set values to NaN
   
   - get rid of nan testing everywhere and stop using np.nanmax() and np.nanargmax())

- Allow user to set ChannelMappingRules and make them more generic
  (by default, only select based on the component code).  If more
  than one channel is selected, keep the one with the highest sample rate.
  If two with highest sample rate choose least broadband(?)

ChannelMapping
------------------------

- Integrate "ChannelMappingRules" into "ChannelMapping"?
- Make select_traces() part of ChannelMapping?
    
Parameter file
------------------------

Replace:

.. code:: yaml
  channel_parameters:
    compZ: 'Z3'
    compN: 'N1Y'
    compE: 'E2X'
    compH: 'HF'
    S_write_cmp: 'N'
    P_write_cmp: 'Z'
    P_write_phase: 'Pg'
    S_write_phase: 'Sg'

by:

.. code:: yaml
  channel_parameters:
    component_orientation_codes:
        Z: 'Z3'
        N: 'N1Y'
        E: 'E2X'
        H: 'HF'
    write_components_phases:
        S: ['N', 'Sg']
        P: ['Z', 'Pg']

Replace:

.. code:: yaml
    station_parameters:
        SPOBS:
            P_comp: 'Z'
            S_comp: 'ZNE'
            energy_frequency_band: [3, 30]
            energy_window: 20  # What does this really do?

by:

.. code:: yaml
    station_parameters:
        SPOBS:
            picking_components
                P: 'Z'
                S: 'ZNE'
            SNR_energy:
                frequency_band: [3, 30]
                window: 20  # What does this really do?
