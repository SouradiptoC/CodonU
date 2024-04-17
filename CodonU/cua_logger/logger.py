import logging
import logging.config

import yaml

# path to config file
logger_config_file = 'logger_config.yaml'


def _create_config_dict(config_file: str) -> dict:
    """
    Creates config dict for logger

    :param config_file: path to yaml file for config
    :return: config dict
    """
    with open(config_file, "r") as config_file:
        return yaml.safe_load(config_file.read())


class CodonuLogger(object):
    def __init__(self):
        self.config_dict = _create_config_dict(logger_config_file)

    def config_logger(self):
        """
        Configures the logger
        """
        logging.config.dictConfig(self.config_dict)
