#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging

LEVEL = logging.INFO
LOGFILE, LOGFILEMODE = '', 'w'

HANDLERS = [logging.StreamHandler(sys.stdout)]
if LOGFILE:
    HANDLERS.append(logging.FileHandler(filename=LOGFILE, mode=LOGFILEMODE))

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(name)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', handlers=HANDLERS, level=LEVEL)

logger = logging.getLogger('[iTOL]')
warn, info, error = logger.warning, logger.info, logger.error

if __name__ == '__main__':
    pass
