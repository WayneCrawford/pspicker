
## Version 0.2:

- Changed date representation in run_many()
- Removed asctime and caller from logging.basicConfig format, because matplotlib
  calls to logging.debug do not fill in these values, which gives an error
- added VERBOSE logging level (using verboselogs package)
- Fixed bug that set all amplitude pick components (and associated phase pick)
  to the amplitude level and component of the last picked station
  
ToDo: add PyPI verboselogs package to allow logger.verbose()?