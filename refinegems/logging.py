import logging
import sys

def setup_logging(logfile):
    """Sets up logging to write to both console and file"""
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    # Create console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Create file handler
    file_handler = logging.FileHandler(logfile)
    file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Redirect stdout to logger
    sys.stdout = LoggerWriter(logger.debug)

class LoggerWriter:
    """Redirects stdout to logger"""

    def __init__(self, level):
        self.level = level

    def write(self, message):
        if message.strip():
            self.level(message.strip())

    def flush(self):
        pass
