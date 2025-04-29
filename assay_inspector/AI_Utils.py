#! /usr/bin/env python

__author__ = "Luca Menestrina"
__date__ = "20250130"
__copyright__ = "Copyright 2025, Chemotargets"
__license__ = ""
__credits__ = ["Data Science & Translational Research Group"]
__maintainer__ = "Luca Menestrina"
__version__ = "20250224"
__deprecated__ = False

import os
from rdkit import Chem
import pandas as pd
import numpy as np
from tqdm import tqdm
import psutil
import numpy as np
import pandas as pd
from datetime import datetime
from tqdm import tqdm
import ray
import threading
import cloudpickle

from rdkit import RDLogger
from contextlib import contextmanager

from functools import wraps

DATE_TAG = datetime.now().strftime("%Y%m%d")
TIME_TAG = datetime.now().strftime("%H%M%S")
DATETIME_TAG = f"{DATE_TAG}_{TIME_TAG}"

import logging
PROCESS_LOG_FILE = "logging.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s # %(levelname)s # %(message)s",
    datefmt="%d/%m/%Y %H:%M:%S",
    handlers=[
        logging.FileHandler(filename=PROCESS_LOG_FILE, mode="a"),
        logging.StreamHandler()
    ]
)

from rdkit import RDLogger
from contextlib import contextmanager

@contextmanager
def no_rdkit_log():

    """
    Context manager to suppress all rdkit loggings.
    """

    RDLogger.DisableLog("rdApp.*")
    yield
    RDLogger.EnableLog("rdApp.*")

from rdkit.Chem.MolStandardize import rdMolStandardize

###
def standardize(originalMol):

    """
    Returns the standardized version of the input molecule.
    """

    with no_rdkit_log():
        repairedMol = Chem.MolFromSmiles(Chem.MolToSmiles(originalMol))
        rdmd = rdMolStandardize.MetalDisconnector()
        noMetalsMol = rdmd.Disconnect(repairedMol)
        # Fragment the structure and select the biggest non metalic fragment
        rdfr = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
        biggestNonMetalMol = rdfr.choose(noMetalsMol)
        rdu = rdMolStandardize.Uncharger()
        unchargedMol = rdu.uncharge(biggestNonMetalMol)
        normalizedMol = Chem.MolFromSmiles(Chem.MolToSmiles(rdMolStandardize.Normalize(unchargedMol)), sanitize=True)

    return normalizedMol

###
def molFromSmiles(smiles):
    if smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
        except:
            mol = None
    else:
        mol = None
    
    return mol

###
def tqdm4ray(ids, *args, **kwargs):

    """
    A tqdm equivalent for ray parallelized processes.

    This function is a wrapper around the `tqdm` function that allows it to work with ray parallelized processes.

    Parameters:
        ids (list): A list of task IDs to monitor.

    Returns:
        A tqdm iterator that yields the results of the tasks.

    Example:
        import ray
        ray.init()
        ids = [ray.remote(lambda x: x)(i) for i in range(10)]
        for result in tqdm4ray(ids):
            print(result)
    """

    def to_iterator(ids):
        while ids:
            done, ids = ray.wait(ids)
            yield ray.get(done[0])
    
    return tqdm(to_iterator(ids), *args, **kwargs)

###
def parallel_apply(dataframe, function, n_processes=None, n_batches=None, desc=None, safe:bool=True, **apply_kwargs):

    """
    Applies a function in parallel to each row of a pandas dataframe using the ray library.

    Parameters:
        - dataframe: the dataframe on which to apply the function
        - function: the function to apply on each row
        - n_processes: the number of parallel processes to run. If None, defaults to the number of available CPU cores.
        - n_batches: the number of batches to split the dataframe into. If None, the dataframe is split into 
                     batches based on 10 times the number of processes (never more than 1/10 of the number of rows 
                     of the dataframe).
        - safe: If True, wraps the function in a try-except block to catch errors during execution, returning None 
            for rows where the function raises an exception. If False, exceptions are not caught, and the 
            function behaves normally. Default is True.
        - apply_kwargs: the parameters/kwargs to pass to the pandas apply function

    Returns:
        The result of the parallel apply operation.

    Notes:
        This function uses the ray library to parallelize the apply operation.
        It works only on axis=1 (i.e., applies the function to each row).
        The number of processes is automatically set to the minimum of the number of CPUs and the length of the dataframe.
    """

    if isinstance(apply_kwargs, dict):
        if "axis" in apply_kwargs.keys():
            if apply_kwargs["axis"] == 1:
                apply_kwargs.pop("axis")
            else:
                raise Exception("Only axis = 1 is implemented")
    else:
        raise Exception("Provide `apply_kwargs` as a dict")

    n_processes = min(n_processes, os.cpu_count()) if n_processes else os.cpu_count()
    if safe:
        def func(row):
            try:
                return function(row)
            except Exception as e:  # TODO: log error
                return None
    else:
        func = function

    # Limit overhead performing a normal apply
    if len(dataframe) < n_processes*10:  # TODO: choose threshold
        if isinstance(dataframe, pd.Series):
            return dataframe.apply(func, **apply_kwargs)
        elif isinstance(dataframe, pd.DataFrame):
            return dataframe.apply(func, axis=1, **apply_kwargs)

    n_batches = min(len(dataframe)//10, max(n_batches, n_processes)) if n_batches else min(len(dataframe), n_processes)
    df_memory = cloudpickle.dumps(dataframe).__sizeof__() / 2**20  # Convert bytes to MiB
    available_memory = psutil.virtual_memory().available / 2**20
    max_memory_per_process = available_memory * 0.8 // os.cpu_count()  # Safety limit 80%
    max_batches_by_memory = int(np.ceil(df_memory / max_memory_per_process))  # Should avoid OOM
    # print("df_memory", df_memory)
    # print("available_memory", available_memory)
    # print("max_memory_per_process", max_memory_per_process)
    # print("max_batches_by_memory", max_batches_by_memory)
    logging.debug(f"Max batches allowed by available memory: {max_batches_by_memory}")
    n_batches = max(n_batches, max_batches_by_memory)

    # print(f"First batch size: {cloudpickle.dumps(dataframe.loc[np.array_split(dataframe.index, n_batches)[0]]).__sizeof__() / 2**20}")
    batches = (dataframe.loc[indeces] for indeces in np.array_split(dataframe.index, n_batches))
    # batches = (ray.put(batch) for batch in batches)

    if isinstance(dataframe, pd.Series):
        @ray.remote
        def remote_func(batch, func):
            return batch.apply(func, **apply_kwargs)
    elif isinstance(dataframe, pd.DataFrame):
        @ray.remote
        def remote_func(batch, func):
            return batch.apply(func, axis=1, **apply_kwargs)
    else:
        raise Exception("`dataframe` is neither a pandas DataFrame nor a pandas Series")

    if not ray.is_initialized():
        running_path = os.path.dirname(os.path.abspath(__file__))
        ray.init(
            num_cpus=n_processes,
            runtime_env={
                 "working_dir": running_path,
                 "excludes":[file for file in os.listdir(running_path) if not file.endswith(".py")]
                },  # slows down a bit the execution (overhead of a few seconds) but it makes it runnable on clusters TODO: improve with a general parallelization class
            log_to_driver=False,
            logging_level=logging.WARNING,
            )
        needs_shutdown = True
    else:
        needs_shutdown = False
    func = ray.put(func) if cloudpickle.dumps(func).__sizeof__()/2**20 >= 95 else func
    tasks_id = [remote_func.remote(batch, func) for batch in batches]
    del func
    results = (batch for batch in tqdm4ray(tasks_id, total=n_batches, desc=desc))
    result = pd.concat(results, axis=0).loc[dataframe.index]

    if needs_shutdown:
        ray.shutdown()

    return result

# Create a thread-local storage instance to track Ray initialization state.
_thread_local = threading.local()

###
def ray_manager(n_processes=None):

    """
    A decorator to manage the Ray lifecycle around the decorated function. 
    It ensures Ray is initialized with the specified number of CPUs, manages error propagation, 
    and guarantees that Ray is shut down if it was initialized by this decorator.

    Parameters:
    -----------
    n_processes : int, optional
        The number of CPU cores to allocate for Ray (default is None). 
        If None, the maximum available CPUs will be used. If a specific number is provided, 
        Ray will use up to that many cores. The value will not exceed the number of available CPU cores.

    Returns:
    --------
    function
        A wrapped function with Ray initialization and cleanup behavior. 
        The function will be executed with Ray initialized if it wasn't already.

    Raises:
    ------
    Exception
        Any exception raised by the decorated function will propagate after Ray has been shut down.
    
    Notes:
    -----
    - This decorator initializes Ray only once per thread, ensuring that a Ray server is 
      started and stopped as needed, preventing multiple redundant initializations within 
      the same thread.
    - The Ray server is shut down **only** if it was initialized by the decorator in the same thread.
    """

    def accessory(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Check if Ray has been initialized in this thread
            if not hasattr(_thread_local, 'ray_initialized') or not _thread_local.ray_initialized:
                # Initialize Ray only if it has not been initialized
                # print(f"START RAY WRAPPER {threading.get_ident()}")
                num_cpus = min(n_processes, os.cpu_count()) if n_processes else os.cpu_count()
                running_path = os.path.dirname(os.path.abspath(__file__))
                ray.init(
                    num_cpus=num_cpus,
                    runtime_env={
                        "working_dir": running_path,
                        "excludes":[file for file in os.listdir(running_path) if not file.endswith(".py")]
                        },  # slows down a bit the execution (overhead of a few seconds) but it makes it runnable on clusters TODO: improve with a general parallelization class
                    log_to_driver=False,
                    logging_level=logging.WARNING,
                )
                _thread_local.ray_initialized = True  # Mark Ray as initialized in this thread
                needs_shutdown = True
            else:
                needs_shutdown = False
            try:
                return func(*args, **kwargs)  # Execute the function
            except Exception as e:
                ray.shutdown()  # Stop Ray on error
                raise e  # Re-raise the error for handling by the caller
            finally:
                if needs_shutdown:
                    # print(f"STOP RAY WRAPPER {threading.get_ident()}")
                    ray.shutdown()  # Ensure Ray is always shut down in the same thread
                    _thread_local.ray_initialized = False  # Reset state

        return wrapper
    
    return accessory

###
def classify_skewness(skewness):
    if skewness < -2.0:
        return 'Severe Left Skewed'
    elif skewness < -1.0:
        return 'Moderate Left Skewed'
    elif skewness > 1.0:
        return 'Moderate Right Skewed'
    elif skewness > 2.0:
        return 'Severe Right Skewed'
    else:
        return 'Non-Skewed'

###
def formatWarningTitle(text):
    base_text = f" {text} Warning "
    num_equals = (75 - len(base_text)) // 2
    equals_str = '=' * num_equals

    # In case the total length is odd, add one more '=' to the end
    if (len(equals_str) * 2 + len(base_text)) < 75:
        return f"\n{equals_str}{base_text}{equals_str}=\n"
    else:
        return f"\n{equals_str}{base_text}{equals_str}\n"