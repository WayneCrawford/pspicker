from paz import PAZ

print('Entering simplified STS-2 response, no explicit norm_factor')
paz = PAZ(poles=[0.037+0.037j, -0.037-0.037j], zeros=[0, 0],
          passband_gain=1500, ref_freq=1.)
# calculates norm factor if it was not provided, so that
# norm_factor * pole-zero equation = 1.0 at ref_freq
print(paz)
print('Changing input_units to nm')
paz.input_units = 'nm'
# norm_factor is changed
print(paz)
print()

print('Entering simplified STS-2 response, explicit norm_factor')
paz = PAZ(poles=[0.037+0.037j, -0.037-0.037j], zeros=[0, 0],
          passband_gain=1500, ref_freq=1., norm_factor=1.0)
print(paz)
print('Changing input_units to nm')
paz.input_units = 'nm'
print(paz)
