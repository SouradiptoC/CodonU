import speedtest
import sys
from CodonU.cua_logger import *


def test_speed():
    """
    Tests network speed
    """
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
