import logging
import sys
import os

COLOR_RESET = "\033[0m"
COLOR_DEBUG = "\033[96m"
COLOR_INFO = "\033[92m"
COLOR_WARNING = "\033[93m"
COLOR_ERROR = "\033[91m"
COLOR_CRITICAL = "\033[97;41m"

# Define string constants for log levels
DEBUG = "DEBUG"
INFO = "INFO"
WARNING = "WARNING"
ERROR = "ERROR"
CRITICAL = "CRITICAL"

# Mapping of string constants to logging levels
LOG_LEVELS = {
    DEBUG: logging.DEBUG,
    INFO: logging.INFO,
    WARNING: logging.WARNING,
    ERROR: logging.ERROR,
    CRITICAL: logging.CRITICAL
}

class PrettyFormatterNoLib(logging.Formatter):
    LEVEL_COLORS = {
        logging.DEBUG: COLOR_DEBUG,
        logging.INFO: COLOR_INFO,
        logging.WARNING: COLOR_WARNING,
        logging.ERROR: COLOR_ERROR,
        logging.CRITICAL: COLOR_CRITICAL
    }

    def format(self, record):
        levelname_color = self.LEVEL_COLORS.get(record.levelno, COLOR_RESET)
        levelname = record.levelname
        colored_levelname = f"{levelname_color}{levelname}{COLOR_RESET}"
        timestamp = self.formatTime(record, self.datefmt)
        log_message = f"{COLOR_RESET}[{timestamp}] [{colored_levelname}] \033[1m{record.name} - {record.funcName}:{record.lineno}{COLOR_RESET} - {record.msg}{COLOR_RESET}"
        if record.exc_info:
            log_message += '\n' + self.formatException(record.exc_info)
        return log_message

def setup_logger(name="svhet", level=DEBUG, log_file=""):
    logger = logging.getLogger(name)
    if logger.hasHandlers():
        return logger
    
    # Convert string level to logging level
    log_level = LOG_LEVELS.get(level.upper(), logging.INFO)
    logger.setLevel(log_level)
    
    for handler in logger.handlers:
        handler.setLevel(log_level)
    
    if log_file:
        if os.path.exists(log_file):
            os.remove(log_file)
        handler = logging.FileHandler(log_file)
    else:
        handler = logging.StreamHandler(sys.stdout)
    
    formatter = PrettyFormatterNoLib(
        fmt='%(asctime)s [%(levelname)s] %(name)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger
