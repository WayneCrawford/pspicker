import logging, verboselogs
import sys
import inspect
from datetime import datetime

from obspy.core import UTCDateTime

def setup_log(stream_log_level=logging.INFO):
    global logger
    # verboselogs.install()
    ts = datetime.today().strftime('%Y%m%dT%H%M')
    logging.basicConfig(filename=f'run_{ts}.log',
                        # format='%(asctime)s %(caller)-25s %(levelname)-8s %(message)s',
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filemode='w',
                        level=logging.VERBOSE)
    console = logging.StreamHandler()
    console.setLevel(stream_log_level)
    f = logging.Formatter('%(levelname)-8s %(message)s')
    console.setFormatter(f)
    logger = verboselogs.VerboseLogger('')
    logger.addHandler(logging.getLogger(''))
    # logger = logging.getLogger('')
    if len(logger.handlers) == 1:
        logger.addHandler(console)
    else:
        logger.handlers[1] = console

def log(string, level="info"):
    """
    Prints a string and logs to file

    :param level: the log level
    """
    global logger
    caller = inspect.stack()[1]
    caller_str = caller[3] + '()'
    level = level.upper()
    extra = {'caller':caller_str}
    if level == "INFO":
        #print(bc.BRIGHTGREEN + full_str + bc.RESET)
        logger.info(string, extra=extra)
    elif level == "VERBOSE":
        logger.verbose(string, extra=extra)
    elif level == "DEBUG":
        logger.debug(string, extra=extra)
    elif level == "ERROR":
        logger.error(string, extra=extra)
    elif level == 'WARNING':
        logger.warning(string, extra=extra)
    elif level == 'CRITICAL':
        logger.critical(string, extra=extra)
    else:
        logger.warning(string, extra=extra)
    # sys.stdout.flush()


class bc:
    RED = '\033[30m'
    BRIGHTRED = '\033[1;31m'
    BRIGHTGREEN = '\033[0;32m'
    BRIGHTYELLOW = '\033[1;33m'
    BRIGHTBLUE = '\033[1;34m'
    BRIGHTMAGENTA = '\033[1;35m'
    BRIGHTCYAN = '\033[1;36m'
    BRIGHTWHITE = '\033[1;37m'
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    RESET = '\033[0m'
    ENDC = '\033[1;m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
