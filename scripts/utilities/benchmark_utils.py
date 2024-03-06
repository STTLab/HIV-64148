'''
HIV-64148  Copyright (C) 2024  Sara Wattanasombat
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it.

This module provides utility functions for benchmarking algorithms and measuring resource usage.

The module includes the following functions:

- log_memory_cpu: A context manager that logs memory and CPU usage metrics for a given process.
- stream_output: Streams the output of a subprocess line by line.
- run_command_with_logging: Runs a command with logging of memory and CPU usage.

Usage:
Import the module using `import utilities.benchmark_utils` 
and then access the functions as `utilities.benchmark_utils.function_name`.
'''
__all__ = ['stream_output', 'run_command_with_logging', 'log_resource_usage']
__version__ = '0.1'
__author__ = 'Sara Wattanasombat'

import sys
import subprocess
import time
import threading
import functools
import psutil
import pandas as pd
from .logger import logger

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
        `stream_output(my_process)`
    '''
    for line in iter(process.stdout.readline, b''): # type: ignore
        if process.poll() is not None:
            return
        if len(line) > 0:
            if isinstance(line, bytes):
                sys.stdout.write(line.decode())
            elif isinstance(line, str):
                sys.stdout.write(line)
            sys.stdout.flush()


def run_command_with_logging(command):
    """
    Run a command with logging of memory and CPU usage.

    Args:
        command (str or list): The command to be executed.

    Returns:
        tuple: A tuple containing three lists - log_time, memory_usage, and cpu_usage.

    Raises:
        None.

    Example:
        `log_time, memory_usage, cpu_usage = run_command_with_logging('ls -l')`
    """
    with subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    ) as process:
        output_thread = threading.Thread(target=stream_output, args=(process,))
        output_thread.start()
        process.wait()
        output_thread.join()
    return process


def log_resource_usage(interval=0.1, output_file='resource_usage_log.csv'):
    '''
    Decorator that logs resource usage during the execution of a function.

    Args:
        interval (float, optional): The interval between each resource usage
                                    measurement in seconds. Defaults to 0.1.
        output_file (str, optional): The path to the output CSV file where resource usage 
                                    data will be saved. Defaults to "resource_usage_log.csv".

    Returns:
        function: The decorated function.

    Example:
        @log_resource_usage(interval=0.5, output_file="usage.csv")
        def my_function():
            # Function body

        my_function()  # Resource usage data will be logged during the execution 
        of my_function and saved to "usage.csv".
    '''
    def decorator(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            t_start = time.time()

            # Create an empty list to store resource usage data
            resource_usage_data = []

            def log_resources():
                current_process = psutil.Process()
                while True:
                    try:
                        for child in current_process.children():
                            ps_name = child.cmdline()
                            ps_cpu_usage = 0
                            ps_threads_usage = 0
                            ps_memory_usage = 0
                            ps_disk_read = 0
                            ps_disk_write = 0
                            for grandchild in child.children(recursive=True):
                                if grandchild.is_running():
                                    ps_cpu_usage += grandchild.cpu_percent(interval=interval)
                                    ps_threads_usage += grandchild.num_threads()
                                    ps_memory_usage += grandchild.memory_info().rss
                                    _io = child.io_counters()
                                    ps_disk_read += _io.read_bytes
                                    ps_disk_write += _io.write_bytes

                            # Get the CPU usage as a percentage
                            bg_cpu_usage = psutil.cpu_percent(interval=interval)
                            # Get the memory usage in kilobytes
                            bg_memory_usage = psutil.virtual_memory().used
                            # Get the disk read and write parameters
                            _io = psutil.disk_io_counters()
                            bg_disk_read = _io.read_bytes # type: ignore
                            bg_disk_write = _io.write_bytes # type: ignore
                            # Append resource usage to the list
                            resource_usage_data.append((
                                round(time.time() - t_start, 3), ps_name, ps_threads_usage,
                                ps_cpu_usage, ps_memory_usage, ps_disk_read, ps_disk_write,
                                bg_cpu_usage, bg_memory_usage, bg_disk_read, bg_disk_write
                            ))
                    except psutil.Error as error:
                        logger.debug('%s', error)

                    if not func_running.is_set():
                        break

            # Create a flag to indicate if the function has completed
            func_running = threading.Event()
            func_running.set()

            # Start the resource logging in a separate thread
            log_thread = threading.Thread(target=log_resources)
            log_thread.start()

            try:
                # Execute the function and its subprocess calls
                result = func(self, *args, **kwargs)

            finally:
                # Signal the completion of the function
                func_running.clear()

                # Wait for the logging thread to finish
                log_thread.join()

                # Create a pandas DataFrame from the resource usage data
                pd.DataFrame(
                    resource_usage_data,
                    columns=[
                        'time (seconds)', 'Process name', 'THREADS', 'Proc CPU Usage (%)',
                        'Proc Memory Usage (kilobytes)', 'Proc Disk Read (bytes)',
                        'Proc Disk Write (bytes)', 'BG CPU Usage (%)',
                        'BG Memory Usage (kilobytes)', 'BG Disk Read (bytes)',
                        'BG Disk Write (bytes)'
                    ],
                ).to_csv(f'{self.output_dir}/{output_file}', index=False)

            return result

        return wrapper

    return decorator
