import os, sys, inspect, resource
from loguru import logger

def setLog(name=None, debug:bool=None, trackMem=None):
    logger.remove()
    debug = debug if debug is not None else os.getenv("DEBUG", 'false').lower()=='true'
    env = os.getenv("ENV_PATH", '.env')
    if name is None:
        try:
            frame = inspect.stack()[-1]
            name = frame.filename
        except Exception: name = '<unknown>'

    base_path = os.getenv(
        "BASE_PATH",
        os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    )
    logsDir = os.path.join(base_path, 'logs'); os.makedirs(logsDir, exist_ok=True)
    debugPath = os.path.join(logsDir, 'debug'); os.makedirs(debugPath, exist_ok=True)
    errorPath = os.path.join(logsDir, 'error'); os.makedirs(errorPath, exist_ok=True)

    # Add stderr handler with custom format
    logger.add(
        sys.stderr,
        format="<level>{level}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
        colorize=True, level="DEBUG" if debug else "INFO"
    )

    # Add file handler for debug level
    logger.add(
        os.path.join(logsDir, f"{name}_debug.log"),
        rotation="500 MB", retention="10 days",
        level="DEBUG", compression="zip"
    )

    # Add file handler for errors only
    logger.add(
        os.path.join(logsDir, f"{name}_error.log"),
        rotation="100 MB", retention="1 month",
        level="ERROR", backtrace=True, diagnose=True
    )

    logger.configure(extra={"environment": env})

    # if trackMem: logger.add(lambda msg: print(f'Memory Usage: {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss} KB - {msg}'), level="TRACE")
    if trackMem:
        logger.add(lambda msg: logger.opt(depth=1).debug(
                f"Memory: {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.1f} MB"
            ), level="TRACE"
        )

    return logger

if __name__ == "__main__":
    log = setLog()
