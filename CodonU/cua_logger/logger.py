import logging
import logging.config
import sys

import yaml


def _config_logger(config_file: str):
    with open(config_file, "r") as config_file:
        config = yaml.safe_load(config_file.read())
        logging.config.dictConfig(config)


_config_logger("logger_config.yaml")

console_logger = logging.getLogger("console_logger")
file_logger = logging.getLogger("file_logger")


def div(x, y):
    try:
        z = x / y
        console_logger.info("Division Successful")
        file_logger.info("Division Successful")
        return z
    except ZeroDivisionError as ze:
        console_logger.critical("Aborting Process. Division by 0")
        file_logger.exception(ze)
        sys.exit(1)
    except Exception as e:
        console_logger.critical("Some error occurred")
        file_logger.exception(e)
        sys.exit(1)


if __name__ == '__main__':
    a = div(8, 'b')
    print(a)
