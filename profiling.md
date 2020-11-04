# Profiling `run_one()` without plots:
- 24.5 seconds total
- 23.0 in ps_picker.run_one()
- 21.7 in pick_one_station() [13 calls at 1.668s per call]
- 17.4 in trace.times()   [68 calls at 0.256s/call]
- 17.0 in polarity.py/calc_dip_rect()
    - 17.0 in polarity.py.polar_analysis()
        - 14.5 in polarity._calc_indices [12 calls at 1.2s/call]
- 15.0 in trace.<listcomp> [68 calls at 0.231s/call]
- 15.0 in utcdatetime.__add__ [1 276 962 calls at 0.000s/call]

## method:
    python -m cProfile -o run.prof runPicker_C.py
    python -m pstats run.prof
    [run.prof% strip
    [run.prof% sort cumtime
    [run.prof% stats 20

## ACTION: In polarity:_calc_indices(): took calculation of trace.times out of loop

### Result:
    - cumtime reduced to 15.7s
    - time in trace.times() reduced to 7.1s [24 calls at 0.296s/call]
### Profiling `run_one()` without plots:
- **15.7** seconds total
- 13.9 in ps_picker.run_one()
- 12.3 in pick_one_station() [13 calls at 0.924s per call]
- **7.1** in trace.times()   [**24** calls at 0.296s/call]
- 6.9 in polarity.py/calc_dip_rect()
    - 6.9 in polarity.py.polar_analysis()
        - **3.8** in polarity._calc_indices [12 calls at 0.32s/call]
- **6.4** in trace.<listcomp> [**24** calls at 0.231s/call]
- **6.2** in utcdatetime.__add__ [**456 108** calls at 0.000s/call]

### Profiling `run_one()` without polarity analysis:
- 8.8 seconds total
- 7.0 in ps_picker.run_one()
- 5.5 in pick_one_station() [13 calls at 0.924s per call]
- 5.0 in run_Kurtosis [12 calls at 0.420s]
- 3.5 in trace.times [12 calls at 0.296s]
- 3.1 in utcdatetime(__add___ [230802 calls at 0.000s]
- 2.1s in kurtosis.pick_trace() [25 calls at 0.085s]
- 1.6s in kutosis.calc_kurtocum() [25 calls at 0.064s]

## ACTION: replace `times = energy.snr.times('utcdatetime')`
- by `times = energy.snr.times('timestamp')`
- in _run_Kurtosis
### RESULT:
- cumtime reduced from 8.8s to 5.0s!
- cumtime with polarity is now 10.7
    
## ACTION: replace times('utcdatetime') by times('timestamp') in polarity()

### RESULT:
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
    
