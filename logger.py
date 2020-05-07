import os
import multiprocessing
import logging
import logging.config


class processLogFilter(logging.Filter):
    """
    This filter only show log entries for specified process name
    """

    def __init__(self, process_name, *args, **kwargs):
        logging.Filter.__init__(self, *args, **kwargs)
        self.process_name = process_name

    def filter(self, record):
        return record.processName == self.process_name


def start_process_logging(logs_dir):
    """
    Add a log handler to separate file for current process
    """
    process_name = multiprocessing.current_process().name
    file_name = f"{process_name}.log"
    log_file = os.path.join(logs_dir, file_name)
    log_handler = logging.FileHandler(log_file)

    log_handler.setLevel(logging.DEBUG)

    formatter = logging.Formatter(
        "%(asctime)s | %(processName)s | %(levelname)s |  %(message)s"
    )
    log_handler.setFormatter(formatter)

    log_filter = processLogFilter(process_name)
    log_handler.addFilter(log_filter)

    logger = logging.getLogger()
    logger.addHandler(log_handler)

    return log_handler


def stop_process_logging(log_handler):
    logging.getLogger().removeHandler(log_handler)
    log_handler.close()


def config_root_logger(logs_dir):
    log_file = os.path.join(logs_dir, "main.log")

    formatter = "%(asctime)s | %(processName)-10s | %(levelname)s | %(message)s"

    logging.config.dictConfig(
        {
            "version": 1,
            "formatters": {
                "root_formatter": {"format": formatter, "datefmt": "%m-%d %H:%M:%S"}
            },
            "handlers": {
                "console": {
                    "level": "INFO",
                    "class": "logging.StreamHandler",
                    "formatter": "root_formatter",
                },
                "log_file": {
                    "class": "logging.FileHandler",
                    "level": "DEBUG",
                    "filename": log_file,
                    "formatter": "root_formatter",
                },
            },
            "loggers": {
                "": {
                    "handlers": ["console", "log_file",],
                    "level": "DEBUG",
                    "propagate": True,
                }
            },
        }
    )
