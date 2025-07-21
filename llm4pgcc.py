import os, argparse, subprocess
from pathlib import Path
from loguru import logger
from dotenv import load_dotenv, set_key, get_key

def log(logFile=get_key("LOGS_PATH"), lvl="TRACE"):
    verbosity = ['TRACE', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']    
    if isinstance(lvl, str) and lvl.upper() in verbosity: lvl = lvl.upper()
    if isinstance(lvl, int) and lvl in range(0,6): lvl = verbosity[lvl]
    if isinstance(lvl, bool): lvl='TRACE' if lvl else 'ERROR'

    logger.remove()
    logger.add(
        logFile, 
        colorize=True,
        format="<red>[{level}]</red> Message : <green>{message}</green> @ {time}",
        level=lvl
    )

if __name__ == "__main__": 
    # create .env
    base_path = os.getcwd()
    env_path = os.path.join(base_path, '.env')
    if not os.path.exists(env_path): subprocess('touch $(pwd)/.env')
    if not load_dotenv(): logger.critical('.env not found..')

    # set paths in .env
    set_key(env_path, "BASE_PATH", base_path)
    set_key(env_path, "ENV_PATH", env_path)    
    subdirs = ["config", "data", "docs", "logs", "resources", "results", "src", "tools"]
    for dir in subdirs: set_key(env_path, f'{dir.upper()}_PATH', os.path.join(base_path, dir.lower()), export=True)
    
    # set empty API secret variables in .env
    accepted_keys = ['huggingface', 'openrouter', 'openai', 'gemini', 'anthropic']
    for key in accepted_keys: set_key(env_path, f'{key.upper()}_API_KEY', value_to_set='', export=True)

    # setup conda env
    
    # TODO: 
    # 1) logging setup
    # file paths, levels, formatting, versatility

    # 2) CondaEnv
    # * path handling, creation/activation
    # ! eval "$(conda shell.bash hook)"
    # ! conda env create -p "{env_path}" -f "{env_file}" python=3.10.13
    # ! conda activate "{env_path}"