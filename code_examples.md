Had to put these in a separate file from README because PyPI can't handle
enhanced Markdown

## Docstrings of the three main methods:

```python
    def __init__(self, parm_file, wav_base_path, database_path_in,
        database_path_out=\'Sfile_directory\', database_format=\'NORDIC\',
        verbose=True, debug_plots=False):
        """
        *parm_file*: path/name of the parameter file
        *wav_base_path*: absolute basepath to the waveform files (just before
        the YEAR/MONTH subdirectories)
        *database_path_in*: absolute basepath to the database/catalog file(s)
        (just before the YEAR/MONTH subdirectories)
        *database_path_out*: path to output database files
        *database_format*: 'NORDIC': Use SEISAN conventions for waveform 
        and database files (naming, and location in YEAR/MONTH subdirectories)
        *verbose*: output \'verbose\' and \'debug\' logs to console (will be 
        flagged DEBUG because logging module has no VERBOSE level)
        *debug_plots*: show debugging plots
        """
```

```python
    def run_one(self, database_filename, plot_global=True, plot_stations=False,
        assoc=None, verbose=False, debug_plots=None):
        """
        Picks P and S arrivals on one waveform, using the Kurtosis
    
        Information in the database file will be appended with the picks.
        *database_filename*: database file to read
        *plot_global*: show global and overall pick plots
        *plot_stations*: show individual station plots
        *assoc*: Associator object (used by run_many())
        *verbose*: same as in creator
        *debug_plots*: same as in creator
        """
```

```python
    def run_many(self, start_date, end_date, plot_global=False,
        plot_stations=False, verbose=False, ignore_fails=False):
        """
        Loops over events in a date range
    
        :param start_date: "YYYYMMDD" or "YYYYMMDDHHMM" of first data to process
        :param end_date: "YYYYMMDD" of last data to process
        :param plot_global: show global and overall pick plots
        :param plot_stations: show individual station plots
        :param ignore_fails: keep going if one run fails
        
        """
```
