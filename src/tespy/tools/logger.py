# -*- coding: utf-8

"""Module for logging specification.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/logger.py

SPDX-License-Identifier: MIT
"""

import logging
import os
import sys
import warnings
from logging import handlers

import tespy

TESPY_LOGGER_ID = "TESPyLogger"
TESPY_PROGRESS_LOG_LEVEL = logging.INFO + 1  # 21
TESPY_RESULT_LOG_LEVEL = logging.INFO + 2  # 22

logging._levelToName[TESPY_PROGRESS_LOG_LEVEL] = 'PROGRESS'
logging._nameToLevel['PROGRESS'] = TESPY_PROGRESS_LOG_LEVEL
logging._levelToName[TESPY_RESULT_LOG_LEVEL] = 'RESULT'
logging._nameToLevel['RESULT'] = TESPY_RESULT_LOG_LEVEL

# Capture warnings globally instead of per file
logging.captureWarnings(True)
logger = logging.getLogger(TESPY_LOGGER_ID)
logger.setLevel(logging.DEBUG)


class FutureWarningHandler:
    def __init__(self, logger):
        self.logger = logger

    def __call__(self, message, category, filename, lineno, file=None, line=None):
        self.logger.warning(
            f"FutureWarning: {message}",
            stacklevel=2  # Adjust the stack level accordingly
        )

# Register the custom warning handler for FutureWarnings
warnings.showwarning = FutureWarningHandler(logger)


# Create a bunch of shorthand functions, this is mostly
# copied straight from the logging module.
def get_logger():
    return logger


def increment_stacklevel(kwargs):
    """"Method to force the logging framework to trace past this file"""
    if "stacklevel" not in kwargs:
        kwargs["stacklevel"] = 1
    return kwargs["stacklevel"] + 1


def log(level, msg, *args, **kwargs):
    """
    Log 'msg % args' with the integer severity 'level'.

    To pass exception information, use the keyword argument exc_info with
    a true value, e.g.

    log(level, "We have a %s", "mysterious problem", exc_info=1)
    """
    logger = get_logger()
    # We force the logging framework to trace past this file
    kwargs["stacklevel"] = increment_stacklevel(kwargs)
    # Last exit for Python < 3.8
    if (
        sys.version_info.major < 3
        or (sys.version_info.major == 3 and sys.version_info.minor < 8)
    ):
        kwargs.pop("stacklevel")
    logger.log(level, msg, *args, **kwargs)


def debug(msg, *args, **kwargs):
    """
    Log 'msg % args' with severity 'DEBUG'.

    To pass exception information, use the keyword argument exc_info with
    a true value, e.g.

    debug("Houston, we have a %s", "thorny problem", exc_info=1)
    """
    # We force the logging framework to trace past this file
    kwargs["stacklevel"] = increment_stacklevel(kwargs)
    return log(logging.DEBUG, msg, *args, **kwargs)


def info(msg, *args, **kwargs):
    """
    Log 'msg % args' with severity 'INFO'.

    To pass exception information, use the keyword argument exc_info with
    a true value, e.g.

    info("Houston, we have a %s", "interesting problem", exc_info=1)
    """
    # We force the logging framework to trace past this file
    kwargs["stacklevel"] = increment_stacklevel(kwargs)
    return log(logging.INFO, msg, *args, **kwargs)


def warning(msg, *args, **kwargs):
    """
    Log 'msg % args' with severity 'WARNING'.

    To pass exception information, use the keyword argument exc_info with
    a true value, e.g.

    warning("Houston, we have a %s", "bit of a problem", exc_info=1)
    """
    # We force the logging framework to trace past this file
    kwargs["stacklevel"] = increment_stacklevel(kwargs)
    return log(logging.WARNING, msg, *args, **kwargs)


def error(msg, *args, **kwargs):
    """
    Log 'msg % args' with severity 'ERROR'.

    To pass exception information, use the keyword argument exc_info with
    a true value, e.g.

    error("Houston, we have a %s", "major problem", exc_info=1)
    """
    # We force the logging framework to trace past this file
    kwargs["stacklevel"] = increment_stacklevel(kwargs)
    return log(logging.ERROR, msg, *args, **kwargs)


def exception(msg, *args, exc_info=True, **kwargs):
    """
    Convenience method for logging an ERROR with exception information.
    """
    # We force the logging framework to trace past this file
    kwargs["stacklevel"] = increment_stacklevel(kwargs)
    error(msg, *args, exc_info=exc_info, **kwargs)


def critical(msg, *args, **kwargs):
    """
    Log 'msg % args' with severity 'CRITICAL'.

    To pass exception information, use the keyword argument exc_info with
    a true value, e.g.

    critical("Houston, we have a %s", "major disaster", exc_info=1)
    """
    # We force the logging framework to trace past this file
    kwargs["stacklevel"] = increment_stacklevel(kwargs)
    return log(logging.CRITICAL, msg, *args, **kwargs)


# Custom logging function that abuses log level TESPY_PROGRESS_LOG_LEVEL
# to report progress information programmatically.
def progress(value, msg, *args, **kwargs):
    """
    Report progress values between 0 and 100, you can also use
    the `extra` dict to modify the progress limits. Additionally,
    log 'msg % args' with severity 'TESPY_PROGRESS_LOG_LEVEL'.

    progress(51, "Houston, we have completed step %d of 100.", 51)
    progress(0.51, "Houston, we have completed %f percent of the mission.", 0.51*100, extra=dict(progress_min=0.0, progress_max=1.0))
    """
    if "extra" not in kwargs:
        kwargs["extra"] = {}
    if "progress_min" not in kwargs["extra"] and "progress_max" not in kwargs["extra"]:
        kwargs["extra"]["progress_min"] = 0
        kwargs["extra"]["progress_max"] = 100
    kwargs["extra"]["progress_val"] = value
    # We force the logging framework to trace past this file
    kwargs["stacklevel"] = increment_stacklevel(kwargs)
    return log(TESPY_PROGRESS_LOG_LEVEL, msg, *args, **kwargs)


# Custom reporting function that abuses log level TESPY_RESULT_LOG_LEVEL
# to report result information programmatically.
def result(msg, *args, **kwargs):
    """
    Report result values by logging 'msg % args' with severity 'TESPY_RESULT_LOG_LEVEL'.

    result("The result is %f", 1.23456)
    """
    kwargs["stacklevel"] = increment_stacklevel(kwargs)
    return log(TESPY_RESULT_LOG_LEVEL, msg, *args, **kwargs)


def add_console_logging(
        logformat=None, logdatefmt="%H:%M:%S", loglevel=logging.INFO,
        log_the_version=True
        ):
    r"""Initialise customizable console logger.

    Parameters
    ----------
    logformat : str
        Format of the screen output.
        Default: "%(asctime)s-%(levelname)s-%(message)s"

    logdatefmt : str
        Format of the datetime in the screen output. Default: "%H:%M:%S"

    loglevel : int
        Level of logging to stdout. Default: 20 (logging.INFO)

    log_the_version : boolean
        If True, version information is logged while initialising the logger.

    """
    # Prepare the log settings
    logformat_setting = "%(asctime)s-%(levelname)s-%(message)s"
    if logformat is not None:
        logformat_setting = logformat

    # Create the console handler and apply the settings
    loghandler = logging.StreamHandler(sys.stdout)
    loghandler.setFormatter(logging.Formatter(logformat_setting, logdatefmt))
    loghandler.setLevel(loglevel)

    # Get the logger object and register the handler
    log = get_logger()
    log.addHandler(loghandler)

    # Submit the first messages to the logger
    if log_the_version:
        info("Used TESPy version: {0}".format(get_version()))

    return None


def add_file_logging(
        logpath=None, logfile=None, logrotation=None,
        logformat=None, logdatefmt=None, loglevel=logging.DEBUG,
        log_the_version=True, log_the_path=True
        ):
    r"""Initialise customisable file logger.

    Parameters
    ----------
    logpath : str
        The path for log files. By default a ".tespy' folder is created in your
        home directory with subfolder called 'log_files'.

    logfile : str
        Name of the log file, default: tespy.log

    logrotation : dict
        Option to pass parameters to the TimedRotatingFileHandler.

    logformat : str
        Format of the file output.
        Default: "%(asctime)s - %(levelname)s - %(module)s - %(message)s"

    logdatefmt : str
        Format of the datetime in the file output. Default: None

    loglevel : int
        Level of logging to file. Default: 10 (logging.DEBUG)

    log_the_version : boolean
        If True, version information is logged while initialising the logger.

    log_the_path : boolean
        If True, the used file path is logged while initialising the logger.

    Returns
    -------
    file : str
        Place where the log file is stored.
    """
    # Prepare the log file settings
    logpath_setting = tespy.tools.helpers.extend_basic_path('log_files')
    if logpath is not None:
        logpath_setting = logpath

    logfile_setting = os.path.join(logpath_setting, 'tespy.log')
    if logfile is not None:
        logfile_setting = os.path.join(logpath_setting, logfile)

    if not os.path.isdir(logpath_setting):
        os.makedirs(logpath_setting)

    logrotation_setting = {'when': 'midnight', 'backupCount': 10}
    if logrotation is not None:
        logrotation_setting.update(logrotation)

    logformat_setting = ("%(asctime)s - %(levelname)s - %(module)s - %(message)s")
    if logformat is not None:
        logformat_setting = logformat

    # Create the file handler and apply the settings
    loghandler = handlers.TimedRotatingFileHandler(logfile_setting, **logrotation_setting)
    loghandler.setFormatter(logging.Formatter(logformat_setting, logdatefmt))
    loghandler.setLevel(loglevel)

    # Get the logger object and register the handler
    log = get_logger()
    log.addHandler(loghandler)

    # Submit the first messages to the logger
    if log_the_path:
        info("Path for logging: {0}".format(logfile_setting))

    if log_the_version:
        info("Used TESPy version: {0}".format(get_version()))

    return logfile_setting


def define_logging(
        logpath=None, logfile='tespy.log', file_format=None,
        screen_format=None, file_datefmt=None, screen_datefmt=None,
        screen_level=logging.INFO, file_level=logging.DEBUG,
        log_the_version=True, log_the_path=True, timed_rotating=None
        ):
    r"""Initialise customisable logger.

    Parameters
    ----------
    logpath : str
        The path for log files. By default a ".tespy' folder is created in your
        home directory with subfolder called 'log_files'.

    logfile : str
        Name of the log file, default: tespy.log

    file_format : str
        Format of the file output.
        Default: "%(asctime)s - %(levelname)s - %(module)s - %(message)s"

    screen_format : str
        Format of the screen output.
        Default: "%(asctime)s-%(levelname)s-%(message)s"

    file_datefmt : str
        Format of the datetime in the file output. Default: None

    screen_datefmt : str
        Format of the datetime in the screen output. Default: "%H:%M:%S"

    screen_level : int
        Level of logging to stdout. Default: 20 (logging.INFO)

    file_level : int
        Level of logging to file. Default: 10 (logging.DEBUG)

    log_the_version : boolean
        If True the actual version or commit is logged while initialising the
        logger.

    log_the_path : boolean
        If True the used file path is logged while initialising the logger.

    timed_rotating : dict
        Option to pass parameters to the TimedRotatingFileHandler.

    Returns
    -------
    file : str
        Place where the log file is stored.

    Notes
    -----
    By default the INFO level is printed on the screen and the DEBUG level
    in a file, but you can easily configure the logger.
    Every module that wants to create logging messages has to import the
    logger.

    Examples
    --------
    To define the default logger you have to import the python logging
    library and this function. The first logging message should be the
    path where the log file is saved to.

    >>> import logging
    >>> from tespy.tools import logger
    >>> mypath = logger.define_logging(
    ...     log_the_path=True, log_the_version=True, timed_rotating={'backupCount': 4},
    ...     screen_level=logging.ERROR, screen_datefmt = "no_date"
    ... )
    >>> mypath[-9:]
    'tespy.log'
    >>> logger.debug('Hi')
    """
    add_console_logging(screen_format, screen_datefmt, screen_level, False)
    return add_file_logging(
        logpath, logfile, timed_rotating,
        file_format, file_datefmt, file_level,
        log_the_version, log_the_path
    )


def get_version():
    """
    Return a string part of the used version.

    If the commit and the branch is available the commit and the branch will b
    returned otherwise the version number.

    Example
    -------
    >>> from tespy.tools import logger
    >>> v = logger.get_version()
    >>> type(v)
    <class 'str'>
    """
    try:
        return check_git_branch()
    except FileNotFoundError:
        return '{0}'.format(check_version())


def check_version():
    """
    Return the actual version number of the used TESPy version.

    Example
    -------
    >>> from tespy.tools import logger
    >>> v = logger.check_version()
    >>> int(v.split('.')[0])
    0
    """
    try:
        version = tespy.__version__
    except AttributeError:
        version = 'No version found due to internal error.'
    return version


def check_git_branch():
    """
    Pass the used branch and commit to the logger.

    The following test reacts on a local system different than on Travis-CI.
    Therefore, a try/except test is created.

    Example
    -------
    >>> from tespy import logger
    >>> try:
    ...    v = logger.check_git_branch()
    ... except FileNotFoundError:
    ...    v = 'dsfafasdfsdf'
    >>> type(v)
    <class 'str'>
    """
    path = os.path.join(os.path.dirname(
        os.path.realpath(__file__)), os.pardir, os.pardir, os.pardir, '.git')

    # Reads the name of the branch
    f_branch = os.path.join(path, 'HEAD')
    f = open(f_branch, "r")
    first_line = f.readlines(1)
    name_full = first_line[0].replace("\n", "")
    name_branch = name_full.replace("ref: refs/heads/", "")
    f.close()

    # Reads the code of the last commit used
    f_commit = os.path.join(path, 'refs', 'heads', name_branch)
    f = open(f_commit, "r")
    last_commit = f.read(8)
    f.close()

    return "{0}@{1}".format(last_commit, name_branch)
