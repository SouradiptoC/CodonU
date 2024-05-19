import speedtest
import sys
import time
import os
import yaml
from datetime import datetime, timezone

from CodonU.cua_logger import *

# File to store the last execution time
LAST_EXECUTION_FILE = 'last_execution.yaml'


def read_last_execution_time():
    if os.path.exists(LAST_EXECUTION_FILE):
        with open(LAST_EXECUTION_FILE, 'r') as file:
            data = yaml.safe_load(file)
            if data and 'last_execution_time' in data:
                return datetime.fromisoformat(data['last_execution_time'])
    return None


def write_last_execution_time(timestamp):
    with open(LAST_EXECUTION_FILE, 'w') as file:
        yaml.dump({'last_execution_time': timestamp.isoformat()}, file)


def _test_speed():
    try:
        console_log.info("Initiating network speed testing")
        file_log.info("Initiating network speed testing")
        speed_test = speedtest.Speedtest()
        download_speed = speed_test.download()
        upload_speed = speed_test.upload()
        console_log.info(f'download speed is {round(download_speed / (1024 * 1024), 2)} MBPS')
        console_log.info(f'upload speed is {round(upload_speed / (1024 * 1024), 2)} MBPS')
        file_log.info(f'download speed is {download_speed / (1024 * 1024)} MBPS')
        file_log.info(f'upload speed is {upload_speed / (1024 * 1024)} MBPS')
    except speedtest.ConfigRetrievalError as ce:
        console_log.error("Network not available")
        file_log.exception(ce)
        sys.exit(-1)
    except Exception as e:
        console_log.error("Some error occurred. Check log file for details.")
        file_log.exception(e)
        sys.exit(-1)


def test_speed():
    """
    Tests network speed on a span of 10 minutes creates a file named last_execution_time for capturing the last exec time
    """
    current_time = datetime.now(timezone.utc)
    last_execution_time = read_last_execution_time()

    if last_execution_time is None or (current_time - last_execution_time).total_seconds() >= 600:
        _test_speed()
        write_last_execution_time(current_time)
    else:
        pass
