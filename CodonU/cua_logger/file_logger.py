import datetime
import logging

from CodonU.cua_logger.logger import CodonuLogger
from datetime import date
import os


def modify_log_files(path: str):
    """
    Deletes log files more than 30 days old

    :param path: path of the file
    """
    abs_file_path = os.path.abspath(os.path.join('logs', path))
    time_stamp = os.path.getmtime(abs_file_path)
    last_modify_date = int(datetime.datetime.fromtimestamp(time_stamp).strftime("%Y%m%d"))
    if last_modify_date < (int(date.today().strftime("%Y%m%d")) - 30):
        os.remove(abs_file_path)
    else:
        pass


def initialize_log_dir():
    """
    Creates `logs` dir if not existing. If existing, deletes all records more than 30 days old
    """
    if os.path.isdir('logs'):
        tuple(map(modify_log_files, (f for f in os.listdir('logs'))))
    else:
        os.mkdir('logs')


class FileLogger(CodonuLogger):
    def __init__(self):
        """
        Creates new log file the date
        """
        self.file_name = date.today().strftime("%Y%m%d") + 'log.log'
        initialize_log_dir()
        super().__init__()
        self.config_dict['handlers']['file']['filename'] = os.path.join('logs', self.file_name)
        self.config_logger()

    def get_logger(self):
        """
        Returns the file logger

        :return: file logger
        """
        return logging.getLogger("file_logger")


if __name__ == '__main__':
    l = FileLogger().get_logger()
    l.info("this is a log msg")
