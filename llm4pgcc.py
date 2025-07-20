import os, argparse, subprocess
from pathlib import Path
from loguru import logger
from dotenv import load_dotenv, set_key, get_key

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
    for dir in subdirs: set_key(env_path, f'{dir.upper()}_PATH', os.path.join(base_path, dir.lower()))

    # set empty API secret variables in .env
    set_key(env_path, "HUGGINGFACE_API_KEY", value_to_set='', export=True)
    set_key(env_path, "OPENROUTER_API_KEY", value_to_set='', export=True)
    set_key(env_path, "OPENAI_API_KEY", value_to_set='', export=True)
    set_key(env_path, "GEMINI_API_KEY", value_to_set='', export=True)
    set_key(env_path, "ANTHROPIC_API_KEY", value_to_set='', export=True)

