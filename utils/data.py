import os, shutil
import polars as pl
from utils import setLog
from dotenv import set_key
from loguru import logger

class Data:
    def __init__(self):
        self.logger = setLog(
            name='utils/Data',
            debug=os.getenv("DEBUG", 'false').lower()=='true',
            trackMem=os.getenv("TRACK_MEM", 'false').lower()=='true'
        )
        self.base_path = os.getenv("BASE_PATH", os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
        self.data_path = os.path.join(self.base_path, 'data'); os.makedirs(self.data_path, exist_ok=True)

    @logger.catch
    def ConfigureRawData(self, sourcePath:str=None) -> pl.DataFrame:
        log = self.logger

        source_path = sourcePath if sourcePath else os.getenv("SOURCE_PATH")
        if source_path is None: raise ValueError(f'"SOURCE_PATH" not included in arguments or found in .env')
        if not os.path.exists(source_path): raise FileNotFoundError(f'source_path is not a valid file path: {source_path}')

        _,ext = os.path.splitext(source_path); ext = ext.lower()
        if ext=='.csv': df = pl.read_csv(source_path, has_header=True)
        elif ext in ['.xlsx', '.xls']: df = pl.read_excel(source_path, has_header=True)
        elif ext=='.json': df = pl.read_json(source_path)
        elif ext=='.parquet': df = pl.read_parquet(source_path)
        else: raise NotImplementedError(f'Source file path is of an unsupported type: [{ext}]')

        if isinstance(df, pl.DataFrame):
            try:
                req = pl.DataFrame(df).select(['compound_name', 'SMILES', 'score'])
                opt = pl.DataFrame(df.select([optCol for optCol in ['pathway', 'target', 'info'] if optCol in df.columns]))
                rawData = req if not opt else pl.concat([req, opt], how='horizontal')
                rawData.write_csv(os.path.join(self.data_path, 'raw_lib.csv'), include_header=True)
                set_key(os.getenv("ENV_PATH"), "RAWDATA_PATH", os.path.join(self.data_path, 'raw_lib.csv'))
                return True
            except Exception as e: raise Exception(f'Error in Data.createRawData(): {e}')

    @logger.catch
    def getRawData(self) -> pl.DataFrame:
        path = os.getenv("RAWDATA_PATH", os.path.join(self.data_path, 'raw_lib.csv'))
        if not os.path.exists(path): raise FileNotFoundError(f'Raw data file path not found: {path}')
        return pl.read_csv(path, has_header=True)