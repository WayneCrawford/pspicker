#Critical
Make amplitude picks print out on a different line from S picks.
    - Initially set Amplitude.pick_id = pick used to find it, but it keeps the
      phase hint instead of changing to IAML and may overwrite the pick time
    - need to make a new "pick" for each amplitude
        - Amplitude can have magnitude_hint: 'ML', category: 'point'
        - Pick should have time:UTCDateTime, phase_hint='IAML'
- Associator should do arrival-time based clustering if origin-time based clustering
  does not work
- Origin-time based clustering should be possible simply with a Vp_over_Vs ratio

#Non-critical

Allow catalog reads (automatically link station names to their responses)

Check that the polarity analysis works like that in Ppol (and make a test case)

Change smoothing parameters to seconds rather than samples? Or set a
max_seconds for the smoothing parameters?

cut traces rather than setting their "out of bounds" values to nan?  This could
reduce memory usage and simplify the code (get rid of nan testing everywhere
and stop using np.nanmax() and np.nanargmax())

Allow user to set ChannelMappingRules and make them more generic
(by default, only select based on the component code).  If more
than one channel is selected, keep the one with the highest sample rate.
If two with highest sample rate choose least broadband(?)

Integrate "ChannelMappingRules" into "ChannelMapping"?  Make select_traces() part
of ChannelMapping?

