
import logging

FORMAT = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s","%Y-%m-%d %H:%M:%S")
ch = logging.StreamHandler()
ch.setFormatter(FORMAT)

logger = logging.getLogger('logger')
logger.addHandler(ch)
