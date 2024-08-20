'''
HIV-64148  Copyright (C) 2024  Sara Wattanasombat
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it.

This moule provide logger
'''
__all__ = ['logger',]
__version__ = '0.1'
__author__ = 'Sara Wattanasombat'

import logging
from .settings import settings
class CustomFormatter(logging.Formatter):

    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "%(asctime)s | %(levelname)s | %(message)s (%(filename)s:%(lineno)d)"

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, "%Y-%m-%d %H:%M:%S")
        return formatter.format(record)

ch = logging.StreamHandler()
ch.setFormatter(CustomFormatter())

logger = logging.getLogger('logger')
logger.addHandler(ch)
logger.setLevel(settings['general']['logger']['level'])
