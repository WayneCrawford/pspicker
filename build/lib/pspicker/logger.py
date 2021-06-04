import logging
import verboselogs
# import sys
import inspect
from datetime import datetime

# from obspy.core import UTCDateTime
#
logger_name = 'pspicker'


# Set up Filter to not output logs from other programs to file handler
# (because it requires the "caller" variable)
class NoParsingFilter(logging.Filter):
    def filter(self, record):
        return record.name == logger_name


def setup_log(stream_log_level=logging.INFO,
              file_log_level=logging.VERBOSE):
    """
    Set console and log file logging levels

    :param stream_log_level: log level for console output
    :kind stream_log_level: verboselogs level, str or None
    :param file_log_level: log level for log file.  If higher than the
        stream_log_level, will be set to stream_log_level

    For information, the log levels are:
        DEBUG: 10
        VERBOSE: 15
        INFO: 20
        WARNING: 30
        ERROR: 40
        CRITICAL: 50
    """
    if stream_log_level is None:
        stream_log_level = logging.INFO
    elif isinstance(stream_log_level, str):
        stream_log_level = getattr(logging, stream_log_level.upper())
    if stream_log_level not in [logging.DEBUG, logging.VERBOSE, logging.INFO,
                                logging.WARNING, logging.ERROR,
                                logging.CRITICAL]:
        stream_log_level = logging.INFO
    ts = datetime.today().strftime('%Y%m%dT%H%M')

    file_log_level = min(file_log_level, stream_log_level)

    verboselogs.install()   # Sets VerboseLogger as the default Logger

    # Set up Formatters
    lf = logging.Formatter('%(asctime)s %(levelname)-8s - '
                           '%(message)s (in %(caller)s)')
    cf = logging.Formatter('%(levelname)-8s %(message)s')

    # Set up Handlers
    lh = logging.FileHandler(f'run_{ts}.log', 'w')
    lh.setLevel(file_log_level)
    lh.setFormatter(lf)
    lh.addFilter(NoParsingFilter())
    ch = logging.StreamHandler()
    ch.setLevel(stream_log_level)
    ch.setFormatter(cf)

    global logger
    logger = logging.getLogger(logger_name)
    logger.setLevel(file_log_level)  # THIS IS NECESSARY!!!!
    logger.addHandler(lh)
    logger.addHandler(ch)
    logging.captureWarnings(True)
    print(f'{logger=}')
    print(f'{logger.handlers=}')


def log(string, level="info"):
    """
    Prints a string and logs to file

    :param level: the log level
    """
    caller = inspect.stack()[1]
    caller_str = caller[3] + '()'
    level = level.upper()
    extra = {'caller': caller_str}

    level = getattr(logging, level.upper())

    global logger
    logger.log(level, string, extra=extra)
