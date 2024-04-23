import logging.config
import yaml
import os

# Creating logs dir
if not os.path.isdir('logs'):
    os.mkdir('logs')

# Creating config file not existing
# A new logfile will be generated at midnight
# A total 7 log files will be there in 'logs' dir
# 1 for the current run, 6 back up files
if not os.path.isfile('CodonU_logger_config.yaml'):
    config_yaml_dict = \
        {'version': 1,
         'Author': 'SouradiptoC',
         'formatters': {
             'console': {'format': '%(asctime)s - %(levelname)s - %(message)s',
                         'datefmt': '%Y-%m-%d %H:%M:%S'},
             'file': {'format': '%(asctime)s - %(module)s - %(funcName)s - %(levelname)s - %(message)s',
                      'datefmt': '%Y-%m-%d %H:%M:%S'}},
         'handlers': {
             'console':
                 {'class': 'logging.StreamHandler',
                  'formatter': 'console',
                  'stream': 'ext://sys.stdout'},
             'rotating_file_handler':
                 {'class': 'logging.handlers.TimedRotatingFileHandler',
                  'formatter': 'file',
                  'filename': 'logs/CodonU_log',
                  'when': 'midnight',
                  'interval': 1,
                  'backupCount': 6,
                  'encoding': 'utf-8',
                  'delay': True,
                  }},
         'loggers': {
             'console_logger': {'handlers': ['console'], 'level': 'INFO'},
             'file_logger': {'handlers': ['rotating_file_handler'], 'level': 'INFO',
                             'propagate': False}}}
    with open('CodonU_logger_config.yaml', 'w') as wf:
        yaml.dump(config_yaml_dict, wf, default_flow_style=False, sort_keys=False)

# Reading data from CodonU_logger_config.yaml file
with open('CodonU_logger_config.yaml', 'r') as f:
    dic = yaml.safe_load(f)
    logging.config.dictConfig(dic)
    print("I've been here")

# Creating the loggers
file_log = logging.getLogger('file_logger')
console_log = logging.getLogger('console_logger')
