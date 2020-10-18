import logging
import sys
import inspect


def log(string, level="info", give_caller=True):
    """
    Prints a colorful and fancy string and logs the same string.

    :param level: the log level
        'info' gives a green output and no level text
        'verbose' gives a YELLOW output and 'VERBOSE' text
        'debug' gives a magenta output and 'DEBUG' text
        Everything else gives a red color and capitalized level text.

    Copied from Lion Krischer's hypoDDpy
    """
    logging.info(string)
    caller_str = ''
    if give_caller:
        caller = inspect.stack()[1]
        # caller_str = caller[1] + ':' + caller[3] + ': '
        caller_str = caller[3] + '()'
    # Info is green
    if level == "info":
        print(bc.BRIGHTGREEN + f">>> {string} (in {caller_str})" + bc.RESET)
    else:
        level = level.upper()
        if level == "DEBUG":
            color = bc.BRIGHTMAGENTA
        elif level == "VERBOSE":
            color = bc.BRIGHTYELLOW
        else:
            color = bc.BRIGHTRED
        print(color + f">>> {level}: {string} (in {caller_str})" + bc.RESET)
    sys.stdout.flush()


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
