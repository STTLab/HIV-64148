'''
This moule provide logger
'''
__all__ = ['logger',]
__version__ = '0.1'
__author__ = 'Sara Wattanasombat'

import logging
from .settings import settings

FORMAT = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s","%Y-%m-%d %H:%M:%S")
ch = logging.StreamHandler()
ch.setFormatter(FORMAT)

logger = logging.getLogger('logger')
logger.addHandler(ch)
logger.setLevel(settings['general']['logger']['level'])
