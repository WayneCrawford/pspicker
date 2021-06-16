ps_picker
===========

Seismological P- and S- wave picker using the modified Kurtosis method

Python port of the picker described in Baillard et al., 2014

debugging information is saved to the local file ``run_{datetime}.log``

Methodology
------------

The picker is based around the Kurtosis, but also uses energy levels, polarity,
clustering and phase association in a 3-step process:

### Step 1: define a global pick window


The *Kurtosis* is calculated for all stations.  The global window
surrounds the most densely clustered region of triggers.

### Step 2: pick P and S arrivals on each station individually

For each station:
    - calculate the *Kurtosis* over coarse to fine scales.
    - Identify candidates on the coarse scale and refine their times using
      the finier scales
    - Choose P- and S- candidates based on the *signal-to-noise level* of
      each pick
    - Verify the candidates using the waveform *polarity*, if possible
       - polarity is only used if one of the picks has a dip of > 30 degrees

### Step 3: associate picks

- Calculate origin times for each trace, based on the P-S delay and
  a simple velocity model (could I use a single Vp/Vs value?)
- If at least 3 origin times are clustered, use their average origin time
  to validate all candidates, possibly dipping into the pool of unused
  candidates for replacemene P and S picks
- If less than 3 origin times are clustered, reject bad P- and S- picks
  based on clustering of P-pick times, S-pick times and P-S delays



Database and waveform files
---------------------------

Are assumed to be in SEISAN structure:
  - Database files: NORDIC format, in ``database_path_in``/``YEAR``/``MONTH``/
    (except run_one, for which the file may be local)
  - Waveform files: one miniseed file per event.  Filename is read from the
    database file and assumed to start with ``YEAR``-``MONTH``.  File is read
    from ``waveform_path_in``/``YEAR``/``MONTH``/
  
 
Example workflow
----------------

### Start by autopicking a few events, with all bells and whistles on:

To pick one event from a database in ``/SEISAN/MAYOBS``:

```python
from pspicker import PSPicker
picker = PSPicker('parameters_C.yaml', '/SEISAN/MAYOBS/WAV/MAYOB',  '/SEISAN/MAYOBS/REA/MAYOB')
picker.run_one('19-0607-59L.S201905', plot_global=True, plot_stations=True, log_level='verbose')
```

Look at all of the plots and verify that the picks and association are as
you expect.  If not, change the paramters and run again.

### Next, pick several events with only the global plots on

The bells and whistles text will be saved to a log file named
run_{DATETIME}.log

To pick events from May 5th to 25th in the same database:

```python
from pspicker import PSPicker
picker = PSPicker('parameters_C.yaml', '/SEISAN/MAYOBS/WAV/MAYOB',  '/SEISAN/MAYOBS/REA/MAYOB')
picker.run_many('20190505', '20190525', plot_global=True)
```

### Finally, run the whole database without plots

*(run_{DATETIME}.log is always created)*

To pick events from May 26th 2019 May 1st 2020:

```python
from pspicker import PSPicker
picker = PSPicker('parameters_C.yaml', '/SEISAN/MAYOBS/WAV/MAYOB', '/SEISAN/MAYOBS/REA/MAYOB')
picker.run_many('20190526', '20200501')
```

The three main methods:
-----------------------

```python
def __init__(self, parm_file, wav_base_path, database_path_in,
             database_path_out='Sfile_directory', database_format='NORDIC'):
    """
    :param parm_file: path/name of the parameter file
    :param wav_base_path: absolute basepath to the waveform files (just before
                          the YEAR/MONTH subdirectories)
    :param database_path_in: absolute basepath to the database/catalog file(s)
                             (just before the YEAR/MONTH subdirectories)
    :param database_path_out: path to output database files
    :param database_format: 'NORDIC' is the only choice for now
        'NORDIC': Use SEISAN conventions for waveform  and database files
                  (naming, and location in YEAR/MONTH subdirectories)
    """
```
```python
def run_one(self, database_filename, plot_global=True, plot_stations=False,
            assoc=None, log_level="verbose", plot_debug=None):
    """
    Picks P and S arrivals on one waveform, using the Kurtosis

    Information in the database file will be appended with the picks.
    :param database_filename: database file to read
    :param plot_global: show global and overall pick plots
    :param plot_stations: show individual station plots
    :param assoc: Associator object (used by run_many())
    :param log_level: console log level (choices = 'debug', 'verbose',
        'info', 'warning', 'error', 'critical'), default='info'
    :param plot_debug: show some debugging plots
    """
```
```python
def run_many(self, start_date, end_date, plot_global=False,
    plot_stations=False, ignore_fails=False, log_level='info'):
    """
    Loops over events in a date range

    :param start_date: "YYYYMMDD" or "YYYYMMDDHHMM" of first data to process
    :param end_date: "YYYYMMDD" of last data to process
    :param plot_global: show global and overall pick plots
    :param plot_stations: show individual station plots
    :param ignore_fails: keep going if one run fails
    :param log_level: console log level (choices = 'debug', 'verbose',
                      'info', 'warning', 'error', 'critical'), default='info'        
    """
```

Parameter and response files 
-----------------------------

[Are documented here](file_examples.md)

To get the same results as with the old Matlab program, set the following
values:

- set ``association:method`` to **"arrival_time"**
- set ``station_parameters:{type}:max_candidates`` to **2**
- set ``SNR:threshold_parameter`` to **0.2**
- set ``SNR:max_threshold_crossings`` to **5**
- set ``global_window:max_candidates`` to **2**

Event amplitudes 
-----------------

Event amplitudes calculations need accurate instrument responses.  The
instrument response filename(s) are input in the parameter file.  If you have
as stationxml file, you can make a pspicker_compatible json_pz file like this:

```python
paz = PAZ.read_stationxml(filename, channel=xxx[, station=xxxx])
paz.write_json_pz (ps_filename)
```

If you have a response in another format that you can read in using obspy,
you can output it to a pspicker-compatible json_pz file like this:

```python
paz = PAZ.from_obspy_response(resp)
paz.write_json_pz(pz_filename)
```

In both cases, you can look at the response using `paz.plot(min_freq=xxx)`, or
you could compare it to the obspy_response using:

```python
fig = resp.plot(min_freq=xxx, label='obspy', show=False)
paz = PAZ.from_obspy_response(resp)
paz.plot(min_freq=xxx, axes=fig.axes, label='PAZ', sym='g.')
```

To Do
-------

- Add event location-based acceptance of solitary P- and S- candidates
- In P-, S- and P-S clustering stage, allow unused candidates to be
  substituted for rejected picks
- Dedicated [To Do file](ToDo.md)
    
Also see the [profiling file](profiling.md)
