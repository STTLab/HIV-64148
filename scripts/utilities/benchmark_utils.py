'''
Program: HIV-64148
Module: benchmark_utils

This module provides utility functions for benchmarking algorithms and measuring resource usage.

Author: Sara Wattanasombat (Faculty of Medicine, Chiang Mai University, Thailand)
Version: 1.0
Date Created: June 3, 2023
Last Updated: June 3, 2023

The module includes the following functions:

- log_memory_cpu: A context manager that logs memory and CPU usage metrics for a given process.
- stream_output: Streams the output of a subprocess line by line.
- run_command_with_logging: Runs a command with logging of memory and CPU usage.

Usage:
Import the module using `import benchmark_utils` and then access the functions as 
`benchmark_utils.function_name`.
'''

import sys
import subprocess
from contextlib import contextmanager
import time
import threading
import psutil
import pandas as pd

@contextmanager
def log_memory_cpu(process: subprocess.Popen, interval: float):
    '''
    Context manager that logs memory, CPU, and disk usage metrics for a given process.

    Args:
        process (subprocess.Popen): Process for which memory, CPU, and disk usage will be logged.
        interval (float): The interval between each metric measurement in seconds.

    Yields:
        tuple: A tuple containing five lists - 
        log_time, memory_usage, cpu_usage, disk_read_bytes, and disk_write_bytes.

    Raises:
        None.

    Example:
        with log_memory_cpu(my_process, 0.5) as 
        (log_time, memory_usage, cpu_usage, disk_read_bytes, disk_write_bytes):
            # Do something with the log_time, memory_usage, cpu_usage, 
            disk_read_bytes, and disk_write_bytes lists.
    '''
    log_memory_usage = []
    log_cpu_usage = []
    log_time = []
    log_disk_read_bytes = []
    log_disk_write_bytes = []
    start_time = time.time()
    def log_metrics():
        '''
        Internal function that continuously logs memory and CPU usage metrics.

        Args:
            None.

        Returns:
            None.
        '''
        while process.poll() is None:
            if psutil.pid_exists(process.pid):
                memory = 0.0
                cpu_percent = 0.0
                disk_read_bytes = 0.0
                disk_write_bytes = 0.0
                # Take sum of all child processes
                for child in psutil.Process(process.pid).children(recursive=True):
                    cpu_percent += child.cpu_percent(interval=0.1)
                    memory += child.memory_info().rss / 1024 / 1024
                    _io = child.io_counters()
                    disk_read_bytes += _io.read_bytes / 1024 / 1024
                    disk_write_bytes += _io.write_bytes / 1024 / 1024
                log_time.append(int(time.time()-start_time))
                log_memory_usage.append(memory)
                log_cpu_usage.append(cpu_percent)
                log_disk_read_bytes.append(disk_read_bytes)
                log_disk_write_bytes.append(disk_write_bytes)
                # print(f'Memory: {memory_info.rss / 1024 / 1024} MB, CPU: {cpu_percent}%')
                time.sleep(interval)
            else:
                break

    logging_thread = threading.Thread(target=log_metrics)
    logging_thread.start()
    try:
        yield log_time, log_memory_usage, log_cpu_usage, log_disk_read_bytes, log_disk_write_bytes
    finally:
        logging_thread.join()

@contextmanager
def log_background_resources(interval):
    '''
    Context manager that logs background resource usage metrics including memory, CPU, and disk I/O.

    Args:
        interval (float): The interval between each metric measurement in seconds.

    Yields:
        tuple: A tuple containing four lists - bg_memory_usage, bg_cpu_usage, 
        bg_disk_read_bytes, and bg_disk_write_bytes.

    Raises:
        None.

    Example:
        with log_background_resources(1.0) as (mem_usage, cpu_usage, disk_read, disk_write):
            # Do something with the mem_usage, cpu_usage, disk_read, and disk_write lists.
    '''
    bg_memory_usage = []
    bg_cpu_usage = []
    bg_read_bytes = []
    bg_write_bytes = []

    def monitor_resources():
        '''
        Internal function that continuously monitors background resource usage metrics.

        Args:
            None.

        Returns:
            None.
        '''
        while True:
            memory_info = psutil.virtual_memory()
            cpu_percent = psutil.cpu_percent(interval=interval)
            _io = psutil.disk_io_counters()
            bg_memory_usage.append(memory_info.percent)
            bg_cpu_usage.append(cpu_percent)
            bg_read_bytes.append(_io.read_bytes)   # pyright: ignore reportGeneralTypeIssues=false
            bg_write_bytes.append(_io.write_bytes) # pyright: ignore reportGeneralTypeIssues=false
            time.sleep(interval)

    bg_thread = threading.Thread(target=monitor_resources)
    bg_thread.start()

    try:
        yield bg_memory_usage, bg_cpu_usage, bg_read_bytes, bg_write_bytes
    finally:
        bg_thread.join()

def stream_output(process: subprocess.Popen) -> None:
    '''
    Stream the output of a subprocess line by line.

    Args:
        process (subprocess.Popen): The subprocess for which the output will be streamed.

    Returns:
        None.

    Raises:
        None.

    Example:
        stream_output(my_process)
    '''
    for line in iter(process.stdout.readline, b''): # pyright: ignore reportGeneralTypeIssues=false
        if process.poll() is not None:
            return
        if len(line) > 0:
            if isinstance(line, bytes):
                sys.stdout.write(line.decode().strip())
            elif isinstance(line, str):
                sys.stdout.write(line.strip())
            sys.stdout.flush()

def run_command_with_logging(command, save_to=None):
    '''
    Run a command with logging of memory and CPU usage.

    Args:
        command (str or list): The command to be executed.

    Returns:
        tuple: A tuple containing three lists - log_time, memory_usage, and cpu_usage.

    Raises:
        None.

    Example:
        log_time, memory_usage, cpu_usage = run_command_with_logging('ls -l')
    '''
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) as process:
        interval = 1
        with log_memory_cpu(process, interval) \
            as (log_time, memory_usage, cpu_usage, disk_read_bytes, disk_write_bytes): \
                # pyright: ignore reportGeneralTypeIssues=false
            output_thread = threading.Thread(target=stream_output, args=(process,))
            output_thread.start()

            process.wait()
            output_thread.join()

            with log_background_resources(interval) \
                as (background_memory_usage, background_cpu_usage, \
                    background_disk_read_bytes, background_disk_write_bytes):
                pass  # Do nothing, just let it run in the background

    stats = (log_time, memory_usage, cpu_usage, disk_read_bytes, disk_write_bytes, \
            background_memory_usage, background_cpu_usage, background_disk_read_bytes, \
            background_disk_write_bytes)

    if save_to:
        pd.DataFrame(
            zip(*stats),
            columns=[
                'log_time', 'memory_usage', 'cpu_usage',
                'disk_read_bytes', 'disk_write_bytes',
                'background_memory_usage', 'background_cpu_usage',
                'background_disk_read_bytes', 'background_disk_write_bytes'
            ]
        ).to_csv(save_to)

    return process, stats
