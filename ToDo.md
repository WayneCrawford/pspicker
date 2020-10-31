# Critical
- Make Associator  do arrival-time based clustering if origin-time based
  clustering does not work
- Allow origin-time based clustering simply with a Vp_over_Vs ratio

# Non-critical

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

# ChannelMapping
    - Integrate "ChannelMappingRules" into "ChannelMapping"?
    - Make select_traces() part of ChannelMapping?

