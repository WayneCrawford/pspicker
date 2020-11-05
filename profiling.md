## method:

Disable plots in runPicker_C.py

    python -m cProfile -o run.prof runPicker_C.py
    python -m pstats run.prof
    [run.prof% strip
    [run.prof% sort cumtime
    [run.prof% stats 20

## First tests: trace.times('utcdatetime') is a HOG!

#### Profiling `run_one()`:
- 24.5 seconds total
- 23.0 in ps_picker.run_one()
- 21.7 in pick_one_station() [13 calls at 1.668s per call]
- 17.4 in trace.times()   [68 calls at 0.256s/call]
- 17.0 in polarity.py/calc_dip_rect()
    - 17.0 in polarity.py.polar_analysis()
        - 14.5 in polarity._calc_indices [12 calls at 1.2s/call]
- 15.0 in trace.<listcomp> [68 calls at 0.231s/call]
- 15.0 in utcdatetime.__add__ [1 276 962 calls at 0.000s/call]

### ACTION: replace times('utcdatetime') by times('timestamp')

Also took trace.times calulation out of loop in polarity._calc_indices()

#### RESULT:
- cumtime is now 7.1s:
  - 1.5s in __init__()???
     -  No significant time in reading the YAML file, another level of __init__?
  - 0.4s to read_waveforms()
  - 0.7s to choose global window
     - ~0.6s in kurtosis:pick_trace (25*0.073)
  - 4.2s in _pick_one_station [13 * 0.325s]
      - ~1.2s in kurtosis.pick_trace (mostly calc_kurtocum)
      - 2.6s in _polarity_analysis() [12 * 0.219s]
         -  2.6s in _calc_dip_rect()
              - 2.5s in polar_analysis()
  - 0.0s to associate
  - 0.2s to save picks
    
# TEST

Ran run_many for all events on 2019-07-21 (5 stations

## Profiling

### ordered by cumtime (time in function and its calls

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.000    0.000  336.084  336.084 runPicker_C.py:1(<module>)
        1    0.000    0.000  333.541  333.541   ps_picker.py:85(run_many)
        1    0.017    0.017  333.532  333.532     ps_picker.py:217(_run_one_day)
      305    0.122    0.000  332.210    1.089       ps_picker.py:131(run_one) 
     1525    0.103    0.000  257.235    0.169         ps_picker.py:279(_pick_one_station)
     1033    0.080    0.000  135.800    0.131           ps_picker.py:544(_polarity_analysis)
     1033    0.173    0.000  127.213    0.123             polarity.py:82(calc_dip_rect)
     1033   27.990    0.027  124.706    0.121               polarity.py:127(polar_analysis)
     2558    0.144    0.000  118.976    0.047           kurtosis.py:84(pick_trace)
  3608934   10.872    0.000  112.345    0.000 {built-in method numpy.core._multiarray_umath.implement_array_function}
     1033    0.032    0.000   85.400    0.083           ps_picker.py:478(_run_Kurtosis)
     2558    0.072    0.000   80.705    0.032             kurtosis.py:113(calc_kurtocum)
     7723    0.228    0.000   65.685    0.009             kurtosis.py:158(calc_cum_kurtoses)
   149991   12.627    0.000   58.785    0.000 copy.py:132(deepcopy)
   151868    0.924    0.000   52.094    0.000 attribdict.py:138(__deepcopy__)
   149991    3.316    0.000   50.536    0.000 copy.py:236(_deepcopy_dict)
   541838    0.561    0.000   48.401    0.000 <__array_function__ internals>:2(cov)
   215913    1.771    0.000   47.292    0.000 copy.py:268(_reconstruct)
   541838   12.726    0.000   46.688    0.000 function_base.py:2270(cov)
    13921    2.894    0.000   45.827    0.003            kurtosis.py:320(_fast_kurtosis)
      305    0.014    0.000   44.215    0.145         ps_picker.py:385(_choose_global_window)
   107243    0.156    0.000   42.991    0.000 trace.py:2238(copy)
   541838   17.904    0.000   40.196    0.000 decomp.py:115(eig)
   100111    1.197    0.000   39.385    0.000 trace.py:250(_add_processing_info)
    60449    1.244    0.000   38.424    0.001 signaltools.py:1719(lfilter)
     2558    0.048    0.000   38.127    0.015 kurtosis.py:198(follow_extrem)
     2558    0.132    0.000   36.604    0.014 kurtosis.py:263(_get_extrema)
      305    0.061    0.000   36.226    0.119 ps_picker.py:412(_gw_get_distri)
    60449    0.078    0.000   35.974    0.001 <__array_function__ internals>:2(apply_along_axis)

### Ordered by tottime (time spent in the function itself (not its calls))
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
   168803   29.457    0.000   29.457    0.000 {built-in method numpy.core._multiarray_umath.correlate}
     1033   27.990    0.027  124.706    0.121 polarity.py:127(polar_analysis)
   541838   17.904    0.000   40.196    0.000 decomp.py:115(eig)
   541838   12.726    0.000   46.688    0.000 function_base.py:2270(cov)
   149991   12.627    0.000   58.785    0.000 copy.py:132(deepcopy)
  3608934   10.872    0.000  112.345    0.000 {built-in method numpy.core._multiarray_umath.implement_array_function}
  2015977   10.739    0.000   10.739    0.000 {method 'reduce' of 'numpy.ufunc' objects}
  5957649    9.067    0.000    9.067    0.000 {built-in method numpy.array}
 34345008    8.048    0.000   18.631    0.000 {built-in method builtins.isinstance}
   113140    7.461    0.000   18.315    0.000 inspect.py:714(getmodule)
  6796345    6.288    0.000   10.586    0.000 abc.py:180(__instancecheck__)
   571182    5.654    0.000   11.906    0.000 _methods.py:143(_mean)
   542871    5.533    0.000    8.302    0.000 stride_tricks.py:114(_broadcast_to)
 32579429    5.196    0.000    5.196    0.000 {method 'get' of 'dict' objects}
 11626891    5.172    0.000    5.206    0.000 {built-in method builtins.hasattr}
  2150855    4.924    0.000   30.342    0.000 trace.py:177(__setitem__)
  3477026    4.369    0.000   18.575    0.000 attribdict.py:84(__setitem__)
  2250998    4.296    0.000    7.164    0.000 utcdatetime.py:1266(__setattr__)
 13362197    4.288    0.000    4.288    0.000 _weakrefset.py:70(__contains__)
   149991    3.316    0.000   50.536    0.000 copy.py:236(_deepcopy_dict)
      305    3.102    0.010    3.937    0.013 ps_picker.py:741(center_distri)
    13921    2.894    0.000   45.827    0.003 kurtosis.py:320(_fast_kurtosis)
   349493    2.872    0.000    2.872    0.000 {built-in method posix.stat}
   541838    2.646    0.000   25.398    0.000 function_base.py:280(average)
 10969445    2.611    0.000    4.037    0.000 inspect.py:64(ismodule)
   541845    2.589    0.000    9.481    0.000 _util.py:226(_asarray_validated)
  1125499    2.584    0.000   10.525    0.000 utcdatetime.py:291(__init__)
   937462    2.584    0.000   12.139    0.000 utcdatetime.py:985(__add__)
   541838    2.381    0.000    3.988    0.000 lapack.py:928(_compute_lwork)
   541838    2.255    0.000    6.071    0.000 function_base.py:422(asarray_chkfinite)
   648959    2.102    0.000   29.711    0.000 attribdict.py:143(update)
   850855    1.975    0.000    6.637    0.000 fromnumeric.py:70(_wrapreduction)
     1525    1.872    0.001    5.784    0.004 energy_snr.py:45(_calc_energy)
   215913    1.771    0.000   47.292    0.000 copy.py:268(_reconstruct)
   566623    1.754    0.000    1.754    0.000 {method 'argsort' of 'numpy.ndarray' objects}
  3573534    1.654    0.000    1.654    0.000 {built-in method builtins.getattr}
  1125499    1.631    0.000    3.238    0.000 utcdatetime.py:521(_set_ns)
   101623    1.507    0.000    4.085    0.000 inspect.py:2102(_signature_from_function)
  1754542    1.452    0.000    1.871    0.000 copy.py:252(_keep_alive)
    60449    1.442    0.000   35.758    0.001 shape_base.py:267(apply_along_axis)
    27301    1.424    0.000   22.861    0.001 smooth_filter.py:8(smooth_filter)
 13832724    1.402    0.000    1.402    0.000 {built-in method builtins.id}
   100276    1.302    0.000    6.499    0.000 inspect.py:1089(getfullargspec)
    60449    1.244    0.000   38.424    0.001 signaltools.py:1719(lfilter)
   541841    1.210    0.000    1.531    0.000 blas.py:370(getter)
   100111    1.197    0.000   39.385    0.000 trace.py:250(_add_processing_info)
   298094    1.176    0.000   10.155    0.000 trace.py:466(__setattr__)
   571237    1.133    0.000    1.602    0.000 _methods.py:59(_count_reduce_items)
  1345221    1.125    0.000    4.669    0.000 fromnumeric.py:52(_wrapfunc)
   541838    1.108    0.000    1.108    0.000 lapack.py:961(_check_work_float)

### Summary

- By code blocks (cumtime):

    - 40% is in polarity.py (called 1033 times)
        - 92% of that is in polarity_analysis() (20%) and its calls (72%)
    - 35% is in kurtosis.py (called 2558 times)
        - 13% is in choose_global window (mostly kurtosis calls)
    -  9% is in logger (mostly inspect)
    -  3% is in read_waveforms()
    -  2% is in local_amplitude
    -  1% is in associator.py
    -  1% is in get_nordic_wavefile_name
    -  0.6% is in save_waveform

- By functions (tottime)

    - 9% in  168803 calls to {built-in method numpy.core._multiarray_umath.correlate}
    - 8% in    1033 calls to polarity.py(polar_analysis)
    - 5% in  541838 calls to decomp.py(eig)
    - 4% in  541838 calls to function_base.py(cov)
    - 4% in  149991 calls to copy.py(deepcopy)
    - 3% in 3608934 calls to {built-in method numpy.core._multiarray_umath.implement_array_function}
    - 3% in 2015977 calls to {method 'reduce' of 'numpy.ufunc' objects}
