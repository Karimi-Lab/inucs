#!/usr/bin/env python
# coding: utf-8

import argparse
import collections.abc as abc
import configparser
import contextlib
import copy
import datetime
import functools
import gzip
import io
import itertools
import logging
import multiprocessing
import os
import platform
import pprint
import re
import shutil
import subprocess
import sys
from multiprocessing.shared_memory import SharedMemory
from pathlib import Path
from timeit import default_timer as timer
from typing import Optional

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError


class S:  # S for Settings
    """ Holds all the Settings and Constants used globally """
    # todo the following constants should be read from a config file
    TESTING_MODE = False
    COMMENT_CHAR = '#'
    FIELD_SEPARATOR = '\t'
    NORM_DISTANCE = 200  # normalization parameter specifying genomic distance
    OUTPUT_NUCS_COLS__MATRIX = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
    OUTPUT_COLS__MATRIX = \
        ['nuc_id1', 'nuc_id2', 'counts', 'dist_sum', 'dist_avg',
         # 'num_same_dist_interas', 'num_same_dist_nuc_pairs', 'expected_counts',
         'norm_'] + OUTPUT_NUCS_COLS__MATRIX
    OUTPUT_COLS__NUC_INTERAS = ['chrom1', 'pos1', 'chrom2', 'pos2', 'nuc_id1', 'nuc_id2']
    HEADER_LINE_BEGINNING = '#columns: '
    USED_COLS = pd.Series(['chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2'])
    LOCATION_COLS = USED_COLS[:4]  # ['chrom1', 'pos1', 'chrom2', 'pos2']
    STRAND_COLS = USED_COLS[4:]  # ['strand1', 'strand2']
    CHROM_COLS = USED_COLS.iloc[[0, 2]]  # ['chrom1', 'chrom2']
    POS_COLS = USED_COLS.iloc[[1, 3]]  # ['pos1', 'pos2']

    @staticmethod
    def get_char_name(char: chr) -> str:
        char_names = {'\t': '<tab>', ' ': '<space>'}
        return char_names.get(char, char)  # return the name, or the char itself if no name is found

    @classmethod
    def get_norm_col_name(cls) -> str:
        for col in cls.OUTPUT_COLS__MATRIX:
            if col.startswith('norm_'):
                return col
        raise RuntimeError("Cannot find a column starting with 'norm_'")

    @classmethod
    def set_norm_col_name(cls, norm_distance: int):
        cls.NORM_DISTANCE = norm_distance
        # update the 'norm_' col and keep all the rest unchanged
        cls.OUTPUT_COLS__MATRIX = [
            f'norm_{norm_distance}' if col.startswith('norm_') else col for col in cls.OUTPUT_COLS__MATRIX]


class InputFileError(SyntaxError):
    def __init__(self, filename):
        message = '\n' \
                  f'Please make sure that the {S.get_char_name(S.FIELD_SEPARATOR)} character is used as ' \
                  f'the field-separator used in file\n{filename}\n' \
                  f'Another possibility is that there are the wrong number of columns in this input file.'
        super().__init__(message)
        self.filename = filename


# todo make Files a service so that multiple processes can interact with it
# todo re-write methods to better handle zip/unzip, fullpath/relative-path, exist/non-exist, ... cases
class Files:
    S_INTER = 'state_inter'
    S_NUC_INTER = 'state_nuc_inter'
    S_MATRIX = 'state_matrix'
    S_PLOT = 'state_plot'
    ZIP_SUFFIX = '.gz'

    __IMP_COLS = pd.Index(  # important cols
        ['chrom', 'strands', 'orientation', 'file', 'file_zip', 'subdir', 'needs_refresh'])

    __EMPTY = pd.DataFrame(columns=__IMP_COLS)

    __STRAND2ORIENT = {tuple('+-'): 'inner', tuple('-+'): 'outer', tuple('++'): 'right', tuple('--'): 'left'}
    __ORIENT2STRAND = {'inner': ['+-'], 'outer': ['-+'], 'right': ['++'], 'left': ['--'],
                       'tandem': ['++', '--'], 'all': ['+-', '-+', '++', '--']}

    __STATES = pd.DataFrame({
        'state': [S_INTER, S_NUC_INTER, S_MATRIX, S_PLOT],
        'subdir': ['cache/interactions', 'cache', 'matrices', '.'],  # plots are stored directly in the working dir
        'data_type': ['cached', 'cached', 'intermediary', 'results'],
        'refresh_col': ['refresh_inter', 'refresh_nuc_inter', 'refresh_matrix', None]  # plots aren't refreshed
    }).set_index('state')

    __CONFIG_FILE = 'last_run.ini'
    __WORKING_DIR_SUFFIX = '.inucs'

    def __init__(self, base_file_name=None, working_dir=None, norm_distance=None,
                 keep_cache=False, n_processes=0, zipped=False, refresh=False):
        self.__working_files: pd.DataFrame = pd.DataFrame()  # will be reset by reset_working_files
        self.__base_file_name = base_file_name
        self.__working_dir = working_dir
        self.__norm_distance = norm_distance
        self.__keep_cache = keep_cache
        self.__n_processes = n_processes
        self.__zipped = zipped
        self.__config = configparser.ConfigParser()
        self.reset_all(base_file_name=base_file_name, working_dir=working_dir, norm_distance=norm_distance,
                       keep_cache=keep_cache, n_processes=n_processes, zipped=zipped, refresh=refresh)

    def reset_all(self, base_file_name=None, working_dir=None, norm_distance=None,
                  keep_cache=False, n_processes=0, zipped=False, refresh=False):

        self.__zipped = zipped
        self.__keep_cache = keep_cache
        self.__norm_distance = norm_distance
        if norm_distance is not None:
            S.set_norm_col_name(norm_distance)

        if n_processes < 1:  # default value is 0
            self.__n_processes = min(20, os.cpu_count() - 1)  # on a large server, limit num of cores used by default
        elif n_processes > os.cpu_count():
            self.__n_processes = os.cpu_count()
        else:
            self.__n_processes = n_processes

        if self.__n_processes < 3:
            raise ValueError(f'Num of processes cannot be {self.__n_processes}. It has to be at least 3.')

        wd_suffix = self.__WORKING_DIR_SUFFIX
        if base_file_name and working_dir:
            self.__base_file_name = Path(base_file_name)
            self.__working_dir = Path(working_dir)
        elif base_file_name:  # then working_dir is None, so determine it
            self.__base_file_name = Path(base_file_name)
            self.__working_dir = Path('.') / (self.__base_file_name.name + wd_suffix)
        elif working_dir:  # then base_file_name is None, so determine it
            self.__working_dir = Path(working_dir)
            interas_file = self.input_interas_file_name
            if interas_file is None:
                wd = self.__working_dir
                self.__base_file_name = wd.stem if wd.suffix == wd_suffix else wd
            else:
                self.__base_file_name = interas_file
        else:
            self.__base_file_name = None
            self.__working_dir = None
            return  # no need to create working dir (so skip __init_working_dir)

        self.__init_working_dir(refresh)

    def reset_working_files(self, chrom_list, refresh: bool = False):
        self.__init_working_files(chrom_list)
        self.__decide_what_needs_refresh(refresh)

    def validate_working_dir(self):
        if self.working_dir is None or not self.working_dir.exists():
            raise RuntimeError(
                f'Working directory does not exist '
                f'{self.working_dir if self.working_dir else ""}')
        if not self.working_dir.is_dir():
            raise FileExistsError(f'Expecting a folder, but a file exists: {self.working_dir}')

        nucs_file = self.indexed_nucs_file()
        if nucs_file is None or not Path(self.working_dir / nucs_file).exists():
            raise RuntimeError(
                f'Cannot find the indexed nucleosome file in the working directory '
                f'{nucs_file if nucs_file else ""}\n'
                f'Please consider refreshing your working directory using "prepare --refresh"')

        matrix_subdir = self.__STATES.loc[self.S_MATRIX, 'subdir']
        matrix_subdir = self.working_dir / matrix_subdir
        if not matrix_subdir.exists() or not any(matrix_subdir.iterdir()):
            raise RuntimeWarning(
                f'WARNING: No matrices folder or files found {matrix_subdir}\n'
                f'Please consider refreshing your working directory using "prepare --refresh"'
            )

    @property
    def working_dir(self):
        return self.__working_dir

    @property
    def base_file_name(self):
        return self.__base_file_name

    def indexed_nucs_file(self, fullpath=False):
        file = None
        if self.__reload_config_file():
            file = self.__config['input_files']['nucs_file']
        if file is None:
            return None
        return self.fullpath(file) if fullpath else Path(file)

    @property
    def input_interas_file_name(self):
        file = None
        if self.__reload_config_file():
            file = self.__config['input_files']['interas_file']
        return Path(file) if file else None

    def __reload_config_file(self):
        config_file = self.working_dir / self.__CONFIG_FILE
        if config_file.exists():
            self.__config.read(config_file)
            return True
        else:
            return False

    def set_input_interas_file(self, interas_file=None):
        interas_file = Path(interas_file).name if interas_file else ''
        try:
            self.__config.add_section('input_files')
        except configparser.DuplicateSectionError:
            pass
        self.__config['input_files']['interas_file'] = interas_file

        with open(self.working_dir / self.__CONFIG_FILE, 'wt') as fh:
            print('# Please do not edit this file!', file=fh)
            self.__config.write(fh)

    def set_input_nucs_file(self, nucs_file=None):
        nucs_file = Path(nucs_file).name if nucs_file else ''
        try:
            self.__config.add_section('input_files')
        except configparser.DuplicateSectionError:
            pass
        self.__config['input_files']['nucs_file'] = nucs_file

        with open(self.working_dir / self.__CONFIG_FILE, 'wt') as fh:
            print('# Please do not edit this file!', file=fh)
            self.__config.write(fh)

    def reset_input_file_names(self):
        self.set_input_interas_file(None)
        self.set_input_nucs_file(None)

    @classmethod
    def is_file_zipped(cls, file) -> bool:
        return Path(file).suffix.endswith(cls.ZIP_SUFFIX)

    @property
    def zipped(self) -> bool:
        """Is output zipped?"""
        return self.__zipped

    @property
    def keep_cache(self) -> bool:
        """If False, the cache subdirectory will be removed after preparation of matrices."""
        return self.__keep_cache

    @property
    def norm_distance(self) -> int:
        """Genomic distance used for normalization."""
        return self.__norm_distance

    @property
    def n_processes(self) -> int:
        """Number of processes simultaneously analysing data."""
        return self.__n_processes

    @classmethod
    def orientations(cls) -> list:
        return list(cls.__STRAND2ORIENT.values())  # 'inner' 'outer' 'right' 'left'

    @classmethod
    def combined_orientations(cls) -> list:
        return [k for k in cls.__ORIENT2STRAND.keys() if k not in cls.orientations()]  # 'tandem' 'all'
        # return list(cls.__ORIENT2STRAND.keys())  # 'inner' 'outer' 'right' 'left' 'tandem' 'all'

    @classmethod
    def get_strands(cls, orientation) -> list:
        return cls.__ORIENT2STRAND.get(orientation)

    @classmethod
    def get_orientation(cls, strands) -> str:  # e.g. '++' -> 'right'
        return cls.__STRAND2ORIENT.get(tuple(strands))

    @classmethod
    def break_combined_orientation(cls, combined_orientation: str) -> list:  # e.g. 'tandem' -> ['right', 'left']
        strands_list = cls.get_strands(combined_orientation)  # 'tandem' -> ['++', '--']
        orients_list = [cls.get_orientation(strands) for strands in strands_list]  # e.g. ['right', 'left']
        return orients_list

    def fix_file_zip_suffix(self, file_name):

        file_name = Path(file_name)
        suffix = file_name.suffix
        if not self.zipped and suffix == self.ZIP_SUFFIX:
            suffix = ''
        elif self.zipped and suffix != self.ZIP_SUFFIX:
            suffix += self.ZIP_SUFFIX

        return file_name.stem + suffix

    def fullpath(self, file_name):
        return self.working_dir / file_name

    def set_needs_refresh_for_all_files(self, need_refresh: bool, state):
        refresh_col = self.__STATES.loc[state, 'refresh_col']
        self.__working_files.loc[:, refresh_col] = need_refresh

    def get_files_needing_refresh(self, state) -> pd.DataFrame:
        refresh_col = self.__STATES.loc[state, 'refresh_col']
        df = self.__working_files.copy().query(refresh_col)
        if df.empty:
            return self.__EMPTY
        df['needs_refresh'] = df[refresh_col]
        df['subdir'] = self.get_subdir(state)
        return df[self.__IMP_COLS]

    def get_files_needing_refresh__dict(self, state) -> dict:
        """returns working files in a dict, which is easier to handle for multiprocessing"""
        df = self.get_files_needing_refresh(state)
        if df.empty:
            return {}

        chroms_and_strands = \
            df[['chrom', 'strands']].apply(lambda row: (row['chrom'], *row['strands']), axis=1)  # e.g., (chrII, +, -)
        full_file_names = df['subdir'] / (df['file_zip'] if self.zipped else df['file'])

        working_files_dict = dict(zip(chroms_and_strands, full_file_names))  # e.g. {(chrII, +, -): 'path/to/outfile'}

        return working_files_dict

    def get_subdir(self, state) -> Path:
        subdir = self.working_dir / self.__STATES.loc[state, 'subdir']
        return subdir

    def get_working_files(self, chrom: str = None, orientation: str = 'all', state=None) -> pd.DataFrame:
        """
        :param chrom: e.g. 'II' or 'chrom14'. Default value is None which returns all chromosomes.
        :param orientation:  It can be any of ('+-', '-+', '++', '--') or ('inner', 'outer', 'right', 'left') or
        ('tandem', 'all'). Further, 'p' instead of '+', and 'n' or 'm' instead of '-' are acceptable.
        Default value is 'all' which returns all orientations.
        :param state: It can be any one of File.S_INTER, Files.S_NUC_INTER, Files.S_MATRIX, Files.S_PLOT
        :return:
        """
        working_files = self.__working_files
        if chrom:  # but if chrom is None, then return all chromosomes
            working_files = working_files.query(f'chrom == "{chrom}"')

        if working_files.empty:
            return self.__EMPTY

        if orientation != 'all':  # otherwise if it's equal to all, then don't do anything
            if orientation == 'tandem':
                orients = ['left', 'right']
            else:
                if isinstance(orientation, str):
                    orientation = orientation.lower()
                    if orientation not in self.__STRAND2ORIENT.values():
                        strands = ''.maketrans('mnpsa', '--++-')  # Minus, Negative, Positive (Plus), Sense, Anti-sense
                        orientation = tuple(orientation.translate(strands))

                if isinstance(orientation, tuple):
                    orientation = self.__STRAND2ORIENT.get(orientation, 'Unknown')

                orients = [orientation]

            working_files = working_files.query(f"orientation in {orients}")
            if working_files.empty:
                return self.__EMPTY

        working_files = working_files.copy()  # so that the following column manipulations is detached from original df

        if state:
            refresh_col = self.__STATES.loc[state, 'refresh_col']
            working_files['subdir'] = self.get_subdir(state)
            working_files['needs_refresh'] = working_files[refresh_col].squeeze()

        return working_files

    def __init_working_dir(self, refresh):
        try:
            if self.working_dir.exists():
                if not self.working_dir.is_dir():
                    raise FileExistsError(f'Expecting a folder, but a file exists: {self.working_dir}')

            self.working_dir.mkdir(parents=True, exist_ok=True)

            if refresh:  # then delete some of files and sub-folders
                try:
                    subdirs_to_remove = self.__STATES.query("data_type in ['cached', 'intermediary']")
                    for _, subdir_info in subdirs_to_remove.iterrows():
                        subdir = (self.working_dir / subdir_info.subdir)
                        if subdir.exists():
                            shutil.rmtree(subdir)  # delete subdir
                    if self.indexed_nucs_file():
                        nucs_file = self.working_dir / self.indexed_nucs_file()
                        nucs_file.unlink(missing_ok=True)  # delete file
                    self.reset_input_file_names()  # resets both arguments to None
                except Exception:
                    error_message = f'Not permitted to remove some of the content in the working directory.\n' \
                                    f'Please remove "{self.working_dir}" manually and try again.'
                    raise RuntimeError(error_message)

            if self.__norm_distance is not None:
                # delete all files containing a non-matching norm distance column
                # for matrix_file in FILES.get_working_files(state=Files.S_MATRIX):
                matrix_subdir = FILES.get_subdir(state=Files.S_MATRIX)
                error_message = ''
                for matrix_file in matrix_subdir.glob('*.txt'):
                    try:
                        header_only = pd.read_csv(matrix_file, sep=S.FIELD_SEPARATOR, nrows=0)
                    except:
                        extracted_norm_distance = None  # so that the matrix_file will be deleted next
                    else:
                        extracted_norm_distance = NucInteraMatrix.extract_norm_distance(header_only)
                    if extracted_norm_distance != self.norm_distance:
                        try:
                            Path(matrix_file).unlink()  # remove matrix_file as it's norm_distance is wrong
                        except Exception:
                            error_message += f'{self.working_dir}\n'

                if error_message:
                    raise RuntimeError('Not permitted to remove the following file(s). '
                                       'Please remove them manually and try again:\n' + error_message)

            # recreate sub-folders
            for _, subdir_info in self.__STATES.iterrows():
                subdir = (self.working_dir / subdir_info.subdir)
                subdir.mkdir(parents=True, exist_ok=True)

        except Exception as e:
            LOGGER.error(f'\nCannot create cache folders in: {self.working_dir}')
            LOGGER.exception(e)
            sys.exit(e)

    def __init_working_files(self, chrom_list):
        df = pd.DataFrame(  # add important cols except ['chrom', 'strands'] bcz they are added by reset_index()
            columns=self.__IMP_COLS.drop(['chrom', 'strands']),
            index=pd.MultiIndex.from_product([chrom_list, self.__STRAND2ORIENT.keys()], names=['chrom', 'strands']))
        df = df.reset_index()
        df['orientation'] = df.strands.map(self.__STRAND2ORIENT)
        suffix = '.txt'  # Regardless of what the original suffix was, replace it with .txt
        df['file'] = self.base_file_name.stem + '_' + df.chrom + '_' + df.orientation + suffix
        df['file_zip'] = df.file + self.ZIP_SUFFIX

        refresh_cols = self.__STATES.dropna()['refresh_col'].tolist()
        df[refresh_cols] = True  # makes three bool columns for file refreshing

        self.__working_files = df

    def __decide_what_needs_refresh(self, refresh):
        status = self.__STATES[::-1].dropna()  # __STATES with rows reversed and the 'plot' row dropped
        subdir_list = status['subdir'].to_list()  # ['matrices', 'cache', 'cache/interactions']
        refresh_col_list = status['refresh_col'].to_list()  # ['refresh_matrix', 'refresh_nuc_inter', 'refresh_inter']

        self.__working_files[refresh_col_list] = refresh  # init all columns to refresh

        if refresh:  # if all are set to True for refresh then nothing more to do
            return

        # Check all the following conditions:
        # Do all files (with or without '.gz') in 'matrices' subdir exist?
        # Do all files (with or without '.gz') in 'cache' subdir exist?
        # Do all files (with or without '.gz') in 'cache/interactions' subdir exist?

        def not_exists(file):
            return not file.exists()

        prev_needs_refresh = True
        for subdir, refresh_col in zip(subdir_list, refresh_col_list):
            subdir = self.working_dir / subdir
            file_names = subdir / self.__working_files.file
            file_names_zip = subdir / self.__working_files.file_zip
            needs_refresh = file_names.map(not_exists) & file_names_zip.map(not_exists)
            needs_refresh &= prev_needs_refresh  # means the two are "anded" together
            prev_needs_refresh = needs_refresh
            self.__working_files.loc[needs_refresh, refresh_col] = True
            no_refresh_needed = not self.__working_files[refresh_col].any()
            if no_refresh_needed:
                break

    @staticmethod
    def time_diff(start_time, end_time=None):
        seconds = round((timer() if end_time is None else end_time) - start_time)
        return f'{datetime.timedelta(seconds=seconds)} (h:m:s)'

    @staticmethod
    def tabulate_dict(dictionary: dict, keys_included: list = None, left_pad: str = '| ') -> str:
        if keys_included is None:
            keys_included = dictionary.keys()
        keys_width = max(map(len, keys_included))
        dictionary = {k: dictionary[k] for k in keys_included}  # reorder according to keys_included
        table = ''
        for k, v in dictionary.items():
            table += f'{left_pad}{k.rjust(keys_width)}  {v}\n'
        return table


class Chroms(abc.Sequence):
    """ Chromosomes """

    def __init__(self, chrom_list_file, comment: str = S.COMMENT_CHAR, name: str = None):
        chroms = pd.read_csv(chrom_list_file, names=['chrom'], comment=comment,
                             sep=S.FIELD_SEPARATOR, header=None, usecols=[0], squeeze=True)

        if len(chroms) > 0 and chroms[0] in ['chrom', 'chr']:  # remove the first row if it was a header row
            chroms = chroms.drop(index=0).reset_index(drop=True)

        chroms = chroms.astype(pd.CategoricalDtype(ordered=True))
        self.__chroms = chroms
        self.__name = name if name else str(Path(chrom_list_file).name)
        LOGGER.info(f'\nChroms file read: "{self.name}" containing {len(self):,} chromosomes:')
        LOGGER.info(', '.join(self.list) + '\n')

    @property
    def name(self) -> str: return self.__name

    @name.setter
    def name(self, name: str): self.__name = name

    @property
    def series(self) -> pd.Series: return self.__chroms.copy()

    @property
    def list(self) -> list: return self.__chroms.to_list()

    @property
    def categorical(self): return self.create_categorical(self.list)

    @staticmethod
    def create_categorical(chroms: list): return pd.CategoricalDtype(categories=chroms, ordered=True)

    def __len__(self): return self.__chroms.__len__()

    def __repr__(self): return self.__chroms.__repr__()

    def __iter__(self): return self.__chroms.__iter__()

    def __getitem__(self, index): return self.__chroms.__getitem__(index)

    def __contains__(self, item): return item in self.__chroms or self.__chroms.isin([item]).any()


class Nucs:
    """ Nucleosomes """

    def __init__(self, nucs_file=None, chroms_file=None, name: str = None, sep: chr = S.FIELD_SEPARATOR,
                 comment: str = S.COMMENT_CHAR, read_index: bool = None, load_nucs_file=True, verbose=True):

        self.__nucs: pd.DataFrame = self.schema()  # nucs is a table with cols: nuc_id,chrom,start,end
        self.__chrom_2_nucs: dict = {}

        if not load_nucs_file:
            return

        if nucs_file is None:
            if read_index is None:
                read_index = True  # nucs_file can only be None if previously been created in working dir with index
            FILES.validate_working_dir()
            nucs_file = FILES.working_dir / FILES.indexed_nucs_file()

        if read_index is None:
            read_index = False  # if read_index is not set so far, then set it to False as default

        try:
            df_nucs = self.read_csv(nucs_file, sep=sep, comment=comment, read_index=read_index)
        except pd.errors.ParserError as parser_error:
            raise InputFileError(nucs_file) from parser_error

        self.__name: str = name if name else str(Path(nucs_file).name)
        chrom_list = Chroms(chroms_file, comment=comment).list if chroms_file else None
        self.__update(df_nucs, chrom_list, self.name)

        if not read_index:
            nucs_file = self.to_csv(output_file=nucs_file)
            FILES.set_input_nucs_file(nucs_file=nucs_file)

        if verbose:
            LOGGER.info(f'\nNucs file read: "{self.name}" containing {len(self):,} nucleosomes:')
            LOGGER.info(self)

    @classmethod
    def from_dataframe(cls, id_chrom_start_end: pd.DataFrame, name: str = None):
        return cls(nucs_file=None, load_nucs_file=False).__update(id_chrom_start_end, name=name)

    def __update(self, id_chrom_start_end: pd.DataFrame, chrom_list: list = None, name: str = None):
        self.__name = name

        if chrom_list:  # keep only rows with chrom in chrom_list
            # id_chrom_start_end = id_chrom_start_end.query(f"chrom in {chrom_list}")  # too slow
            chrom_cat = Chroms.create_categorical(chrom_list)
            id_chrom_start_end['chrom'] = id_chrom_start_end['chrom'].astype(chrom_cat)
            id_chrom_start_end = id_chrom_start_end.dropna()  # drop any chroms not in chrom_list

        if 'nuc_id' in id_chrom_start_end.columns:
            id_chrom_start_end = id_chrom_start_end.set_index('nuc_id').sort_index()
        elif 'nuc_id' != id_chrom_start_end.index.name:  # index name is already 'nuc_id' when selected from regions
            chrom_cat = Chroms.create_categorical(id_chrom_start_end['chrom'].unique().tolist())
            id_chrom_start_end['chrom'] = id_chrom_start_end['chrom'].astype(chrom_cat)
            id_chrom_start_end = id_chrom_start_end.sort_values(['chrom', 'start', 'end'])
            id_chrom_start_end = id_chrom_start_end.reset_index(drop=True).rename_axis('nuc_id')

        # ensuring cols existence in correct order, and ignoring extra cols if any
        self.__nucs = id_chrom_start_end[['chrom', 'start', 'end']]

        # internal cache to reduce the lookup time for each chrom
        id_chrom_start_end = self.__nucs.reset_index()
        self.__chrom_2_nucs = \
            {k: v.set_index('nuc_id') for k, v in id_chrom_start_end.groupby('chrom')}

        return self

    @staticmethod
    def schema() -> pd.DataFrame:
        """:returns: an empty DataFrame table with appropriate columns and index name."""
        return pd.DataFrame(columns=['chrom', 'start', 'end']).rename_axis(['nuc_id'])

    @property
    def name(self) -> str:
        return self.__name if self.__name else 'nucleosomes'

    @name.setter
    def name(self, name: str):
        self.__name = name

    @property
    def chrom_list(self) -> list:
        return list(self.__chrom_2_nucs.keys())

    @property
    def chrom_categorical(self):
        return pd.CategoricalDtype(categories=self.chrom_list, ordered=True)

    @property
    def df_nucs(self):
        return self.__nucs  # .copy()  # not copying for efficiency

    def get_nucs(self, chrom) -> pd.DataFrame:
        return self.__chrom_2_nucs[chrom]

    def iter_nucs(self):
        return self.__nucs.iterrows()

    def iter_chroms(self):
        return iter(self.__chrom_2_nucs.keys())

    def iter_chrom_and_nucs(self):
        return iter(self.__chrom_2_nucs.items())

    def __len__(self):
        return len(self.__nucs)

    def __repr__(self):
        return f'{self.name}:\n' \
               f'{self.__nucs.head(3)}\n...\n' \
               f'({self.__nucs.shape[0]:,} nucleosomes)'

    @functools.lru_cache
    def find_nuc_id(self, chrom, pos):
        df = self.__chrom_2_nucs[chrom]
        i = df.start.searchsorted(pos)
        if i == len(df) or (i > 0 and pos != df.at[df.index[i], 'start']):
            i -= 1
        nuc = df.iloc[i, :]
        found = nuc.start <= pos <= nuc.end
        return nuc.name, found  # nuc.name is nuc_id

    @functools.lru_cache
    def find_nuc(self, chrom, pos):
        id_, found = self.find_nuc_id(chrom, pos)
        return self.__nucs.iloc[id_, :] if found else None

    @functools.lru_cache
    def find_nucs_in_region(self, chrom, start_region, end_region) -> Optional['Nucs']:
        if chrom not in self.__chrom_2_nucs:
            return None

        if start_region > end_region:  # swap if in the wrong order
            start_region, end_region = end_region, start_region

        start_id, _ = self.find_nuc_id(chrom, start_region)
        end_id, _ = self.find_nuc_id(chrom, end_region)
        nucs_in_region = self.__nucs.iloc[start_id:end_id + 1, :]

        region_name = f'{self.name} ({chrom}:{start_region}->{end_region})'
        return Nucs.from_dataframe(nucs_in_region, region_name)

    @staticmethod
    def read_csv(nucs_file, sep=S.FIELD_SEPARATOR, comment: str = S.COMMENT_CHAR, read_index=False) -> pd.DataFrame:
        nucs_cols_dict = {0: 'chrom', 1: 'start', 2: 'end'}
        nucs_cols = list(nucs_cols_dict.values())  # ['chrom', 'start', 'end']
        first_row = pd.read_csv(nucs_file, sep=sep, comment=comment, header=None, nrows=1)
        first_row = first_row.iloc[0, :].tolist()
        if nucs_cols == first_row[0:len(nucs_cols)] or \
                nucs_cols == first_row[1:len(nucs_cols) + 1]:  # needed because 1st column is optional
            header = 0  # header on the first row
        else:
            header = None  # no header line in the file

        try:
            df = pd.read_csv(nucs_file, sep=sep, comment=comment, header=header)  # read nucs info: chrom, start, end
        except pd.errors.ParserError as parser_error:
            raise InputFileError(nucs_file) from parser_error

        if header is None:
            df = df.rename(nucs_cols_dict, axis=1)

        df = df.rename({'chr': 'chrom'}, axis=1)
        if read_index:
            if 'nuc_id' in df:
                df = df.set_index('nuc_id')  # .sort_index()  # It is assumed here that nuc_id is already sorted
            else:
                raise RuntimeError(f"Cannot read index as no column is called 'nuc_id' in {nucs_file}")
        else:
            df = df.sort_values(nucs_cols)  # sorts and ensures cols exist
            df = df.reset_index(drop=True).rename_axis('nuc_id')
            df = df[nucs_cols]  # reordering cols and ignoring extra cols if any

        df['chrom'] = df['chrom'].astype(pd.CategoricalDtype(ordered=True))

        return df

    def to_csv(self, output_file=None, sep=S.FIELD_SEPARATOR, index: bool = True):
        output_file = FILES.fix_file_zip_suffix(output_file if output_file else self.name)
        self.__nucs.to_csv(
            FILES.fullpath(output_file), sep=sep, index=index, header=True)
        return output_file


class InterasFileSplitter:
    __DONE = 'DONE!'

    @classmethod
    def split_interas_file(cls, outfiles: dict, interas_file, chroms_list):
        manager = multiprocessing.Manager()
        q1_read = manager.Queue()
        q2_split = manager.Queue()
        q3_write = manager.Queue()

        header_line_num, used_cols = NucInteras.validate_interas_file(interas_file)

        n_splitter_processes = FILES.n_processes - 2  # two processes are used for reading from and writing to files
        args = (interas_file, header_line_num, n_splitter_processes, q1_read, q2_split)
        reader_process = multiprocessing.Process(target=cls._worker1_read_input_file, args=args)
        reader_process.start()

        args = (q3_write, n_splitter_processes)
        writer_process = multiprocessing.Process(target=cls._worker3_merge_split_files, args=args)
        writer_process.start()

        subdir_inter = next(iter(outfiles.values())).parent  # in process: no access to FILES.get_subdir(Files.S_INTER)
        # n_splitter_processes = min(n_splitter_processes, len(outfiles))  # not needed here; more processes are better
        args = [(outfiles, chroms_list, used_cols.to_dict(), q1_read, q2_split, q3_write)] * n_splitter_processes
        LOGGER.info(f'In total, creating at most {len(outfiles)} split files using {FILES.n_processes} '
                    f'parallel processes:\n\t{subdir_inter}\nThis may take a while. Please wait')
        with(multiprocessing.Pool(n_splitter_processes)) as splitter_pool:
            splitter_pool.starmap(cls._worker2_split_input, args)

        reader_process.join()
        writer_process.join()

        cls.Buffer.release_all_buffers()

    @classmethod
    def _worker1_read_input_file(cls, interas_file, header_line_num, n_splitter_processes, q1_read, q2_split):

        n_buffers = round(n_splitter_processes * 1.3)  # more buffers than workers to help prevent worker waiting time
        for _ in range(n_buffers):
            buffer_name = cls.Buffer.new_buffer_name()
            q1_read.put(buffer_name)  # put the new buffer into empty queue

        buffer = cls.Buffer()
        with io.FileIO(interas_file, 'rb') as in_f:  # FileIO supports "readinto" predefined buffer space
            cls.__skip_comments_and_header_line(in_f, header_line_num)
            for i in itertools.count():  # infinite loop with counter
                buffer_name = q1_read.get(block=True)
                with buffer.reset(buffer_name):
                    if 0 == buffer.read_into_buf_from(in_f, i):  # stop if the returned num of bytes read is zero
                        q1_read.put(buffer.name)  # put back buffer to free up its shared_memory later
                        q2_split.put(cls.__DONE)  # signal split workers to finish
                        break
                    q2_split.put(buffer.name)

    @classmethod
    def __skip_comments_and_header_line(cls, interas_file_handle, header_line_num, comment=S.COMMENT_CHAR):
        """Skips the leading comment lines, assuming no other comment lines occur later in the file"""
        comment = bytearray(comment, 'utf-8')
        last_pos = interas_file_handle.tell()
        line = interas_file_handle.readline()
        while line != '':
            if line.startswith(comment):
                last_pos = interas_file_handle.tell()
                line = interas_file_handle.readline()
            else:
                if header_line_num is None:  # so no_header_line_exists
                    interas_file_handle.seek(last_pos)  # put back the last line which is the first non-comment line
                break

    @classmethod
    def _worker2_split_input(cls, outfiles: dict, chroms_list, used_cols: dict, q1_read, q2_split, q3_write):
        # outfiles: e.g. {(chrII, +, -): 'path/to/file', ...}

        chroms_categorical = pd.CategoricalDtype(categories=chroms_list)

        chrom1, chrom2 = [col_num for col_num, col_name in used_cols.items() if col_name in S.CHROM_COLS.values]  # 1, 3
        strand1, strand2 = [col_num for col_num, col_name in used_cols.items() if col_name in S.STRAND_COLS.values]

        buffer = cls.Buffer()
        while True:
            buffer_name = q2_split.get(block=True)
            if buffer_name == cls.__DONE:
                q2_split.put(cls.__DONE)  # put it back in queue to signal other splitting workers to finish
                q3_write.put(cls.__DONE)  # also signal writer worker to finish
                break

            with buffer.reset(buffer_name):
                with io.BytesIO(buffer.buf) as buf:
                    df = pd.read_csv(  # low_memory=False is set to suppress a warning
                        buf, usecols=used_cols.keys(), sep=S.FIELD_SEPARATOR, header=None, low_memory=False)
                buffer_id = buffer.id
                q1_read.put(buffer.name)  # now buffer can be reused by the reader worker

            df[[chrom1, chrom2]] = df[[chrom1, chrom2]].astype(chroms_categorical)  # for faster filtering and comparing
            df = df.dropna()  # drops any rows not containing chroms stated in chroms_categorical
            df = df[df[chrom1] == df[chrom2]]  # fastest row filtering approach when used with categorical columns
            output_chunks = dict()
            for chrom_and_strands, chunk_df in df.groupby([chrom1, strand1, strand2]):
                outfile = outfiles.get(chrom_and_strands)  # e.g., {(chrII, +, -): 'path/to/outfile'}
                if outfile is None:
                    continue
                chunk_file = outfile.parent / f'{outfile.name}-chunk{buffer_id}-pid{os.getpid()}'
                chunk_df = chunk_df.iloc[:, :len(S.LOCATION_COLS)]  # keeping: ch1 pos1 ch2 pos2,  dropping the strands
                chunk_df.to_csv(chunk_file, sep=S.FIELD_SEPARATOR, index=False, header=False)
                output_chunks.update({outfile: chunk_file})

            q3_write.put(output_chunks)

    @classmethod
    def _worker3_merge_split_files(cls, q3_write, n_splitter_processes):

        n_finished_split_processes = 0
        while True:
            output_chunks = q3_write.get(block=True)
            if output_chunks == cls.__DONE:
                n_finished_split_processes += 1
                if n_finished_split_processes < n_splitter_processes:
                    continue
                else:
                    break

            for outfile, chunk_file in output_chunks.items():
                with open(outfile, 'ab', buffering=10_000_000) as fd_write:  # open for append
                    with open(chunk_file, 'rb') as fd_read:
                        shutil.copyfileobj(fd_read, fd_write, length=-1)  # read the entire file at once and copy it
                    try:
                        chunk_file.unlink()  # delete chunk file
                    except:
                        pass

    @classmethod
    def __worker3_write_split_strings(cls, q3_write, n_splitter_processes):

        n_finished_split_processes = 0
        while True:
            output_chunks = q3_write.get(block=True)
            if output_chunks == cls.__DONE:
                n_finished_split_processes += 1
                if n_finished_split_processes < n_splitter_processes:
                    continue
                else:
                    break

            for outfile, chunk_csv in output_chunks.items():
                with open(outfile, 'a') as out_f:  # open for append
                    out_f.write(chunk_csv)

    @classmethod
    def __concatenate(cls, outfile):
        """e.g., concatenates all 'outfile-*' files into one file 'outfile' and save the file"""

        import shutil
        with open(outfile, 'wb', buffering=10_000_000) as fd_write:
            for outfile_chunk in sorted(outfile.parent.glob(f'{outfile.name}-*')):
                with open(outfile_chunk, 'rb') as fd_read:
                    shutil.copyfileobj(fd_read, fd_write, length=-1)  # read the entire file at once and copy it
                try:
                    outfile_chunk.unlink()  # delete chunk file
                except:
                    pass

    class Buffer:
        """Must be used within a "with" statement"""

        BYTES_FOR_ID = 3  # Number of bytes reserved at the end of buffer to hold the id for each buffer content
        BYTES_FOR_LENGTH = 5  # Number of bytes reserved at the end of buffer to hold the length of the buffer content
        SIZE_IN_MEGA_BYTES = 10  # Estimated size of buffer in bytes. E.g. 100MB fits ~ 1 million lines, each 100 chars

        __buffer_names = []
        __should_unregister_shared_memory = 'Linux' == platform.system()

        def __init__(self):
            self.__shared_memory = None
            self.__buffer_name = None  # a name for shared_memory waiting to be assigned in next __enter__ call

        @classmethod
        def __get_shared_memory(cls, shared_memory_name: str = None):
            if shared_memory_name:
                shm = SharedMemory(create=False, name=shared_memory_name)
            else:
                shm = SharedMemory(create=True, size=cls.SIZE_IN_MEGA_BYTES * 1_000_000)
                cls.__buffer_names.append(shm.name)

            return shm

        @classmethod
        def new_buffer_name(cls):
            """creates shared_memory and returns its name"""
            return cls.__get_shared_memory().name

        def reset(self, buffer_name: str):
            """Must to be used with a "with" statement: e.g. "with buffer.reset('buf_name'):" """
            self.__buffer_name = buffer_name  # __buffer_name only needed to enforce use of "with" statement
            return self

        def __enter__(self):  # enter the context for with statement
            """find and reconnect to a previously created shared memory, using buffer_name previously set by reset()"""
            self.__shared_memory = self.__get_shared_memory(self.__buffer_name)

        def __exit__(self, exc_type, exc_val, exc_tb):  # exit the context for with statement
            if self.__should_unregister_shared_memory:  # depends on platform
                # unregistering shared memory is very slow, and is used here as a workaround for python bug:
                # https://bugs.python.org/issue38119
                multiprocessing.resource_tracker.unregister(self.__shared_memory._name, 'shared_memory')
            self.__shared_memory = None

        @classmethod
        def release_all_buffers(cls):
            for buffer_name in cls.__buffer_names:
                shm = cls.__get_shared_memory(buffer_name)
                shm.unlink()

        def read_into_buf_from(self, infile: io.RawIOBase, content_id: int = 0) -> int:
            shm_buf = self.__shared_memory.buf
            n_bytes_read = infile.readinto(shm_buf)  # overwrites id and length variables, so they must to be set later

            self.set_id(content_id)  # MUST be after reading in data, otherwise gets overwritten

            length = n_bytes_read
            if length <= self.max_length:  # including length == 0, when no data is read in
                self.__set_length(length)
                return length

            length = self.max_length
            lf, cr = ord('\n'), ord('\r')  # Unix: \n,  Windows \r\n,  Classic Mac \r
            for i in range(length - 1, 0, -1):  # reverse search for \n or \r
                if shm_buf[i] == lf or shm_buf[i] == cr:  # checks for \n first because searching in reverse direction
                    length = i + 1
                    break

            pos = infile.tell() - (n_bytes_read - length)
            infile.seek(pos)  # put back the last incomplete line to be read in again next time

            self.__set_length(length)
            return length

        @property
        def buf(self):
            return self.__shared_memory.buf[0:self.length]

        @property
        def name(self) -> str:
            return self.__shared_memory.name

        @property
        def id(self):
            start = - (self.BYTES_FOR_ID + self.BYTES_FOR_LENGTH)
            end = - self.BYTES_FOR_LENGTH
            int_bytes = self.__shared_memory.buf[start:end]  # read the last bytes before length bytes from buf
            return int.from_bytes(int_bytes, byteorder='big', signed=False)

        def set_id(self, content_id: int):
            """Writes id at the end of shared_memory.buf just before length"""
            max_content_id = 2 ** (8 * self.BYTES_FOR_ID) - 1
            if not 0 <= content_id <= max_content_id:
                raise ValueError(f'content_id must be between 0 and {max_content_id}')
            start = - (self.BYTES_FOR_ID + self.BYTES_FOR_LENGTH)
            end = - self.BYTES_FOR_LENGTH
            self.__shared_memory.buf[start:end] = \
                content_id.to_bytes(self.BYTES_FOR_ID, byteorder='big', signed=False)

        def __set_length(self, length: int):
            """Writes length at the end of shared_memory.buf just after id"""
            if not 0 <= length <= self.max_length:
                raise ValueError(f'length must be between 0 and {self.max_length}')
            start = - self.BYTES_FOR_LENGTH
            self.__shared_memory.buf[start:] = \
                length.to_bytes(self.BYTES_FOR_LENGTH, byteorder='big', signed=False)

        @property
        def length(self) -> int:
            int_bytes = self.__shared_memory.buf[-self.BYTES_FOR_LENGTH:]  # read the last bytes from buf
            return int.from_bytes(int_bytes, byteorder='big', signed=False)

        @property
        def max_length(self) -> int:
            return self.__shared_memory.buf.nbytes - self.BYTES_FOR_ID - self.BYTES_FOR_LENGTH


class NucInteras:
    """ Nucleosome Interactions """

    def __init__(self, interas_file=None, nucs=None, name: str = None,
                 create: bool = True, refresh: bool = False, verbose: bool = True):

        if interas_file is None:
            create = False  # if interas_file is not given, cannot create nuc interas files
            interas_file = FILES.input_interas_file_name

        self.__name = name if name else str(Path(interas_file).name)
        self.__len = 0  # length or number of nuc interactions

        if nucs is None:
            nucs = Nucs()  # loads nucs from the working directory

        # # The following code was moved to __create_nuc_interas to support multiprocessing
        # # Make a MultiIndex ['chrom', 'pos'] for Nucs to make them similar to Interactions, simplifying their merger
        # df_nucs = nucs.df_nucs.copy()
        # df_nucs['pos'] = df_nucs.start  # copy 'start' to 'pos' to make df_nucs and later df_inter alike
        # df_nucs = df_nucs.reset_index().set_index(['chrom', 'pos'])  # .sort_index()
        # self.__nucs__id_chrom_start_end = df_nucs

        self.__chrom_list = nucs.chrom_list

        if create:
            FILES.reset_working_files(chrom_list=self.__chrom_list, refresh=refresh)
            FILES.set_input_interas_file(interas_file=interas_file)

            self.__interas_file = Path(interas_file).resolve(strict=True)

            if verbose:
                LOGGER.info('\nSplitting interaction file based on chrom and orientation')

            start = timer()
            if 'Windows' == platform.system():
                if verbose:
                    LOGGER.info('NOTE: Due to a bug in Python for Windows (https://bugs.python.org/issue38119), '
                                'the splitting step here is much slower in Windows than in macOS or Linux.')
                self.__split_interas_file__single_process()
            else:
                self.__split_interas_file__multiprocessing()

            if verbose:
                LOGGER.info(f'Done! Splitting task finished in {Files.time_diff(start)}')

                LOGGER.info('\nFinding nucleosomes for interactions')

            start = timer()
            self.__create_nuc_interas_files__multiprocessing()

            if verbose:
                LOGGER.info(f'Done! Finding nucleosomes finished in {Files.time_diff(start)}')

    @property
    def name(self) -> str:
        return self.__name

    @name.setter
    def name(self, name: str):
        self.__name = name

    @staticmethod
    def schema() -> pd.DataFrame:
        """:returns: an empty DataFrame table with appropriate columns and index name."""
        return NucInteras.__EMPTY

    __EMPTY = pd.DataFrame(columns=pd.MultiIndex.from_product(
        [['side1', 'side2'], ['chrom', 'pos']], names=['side', 'loc'])).rename_axis(index='inter_id')

    @staticmethod
    def get_nuc_interas_file(chrom: str, orientation: str) -> Optional[Path]:
        """
        :param chrom: e.g. 'II' or 'chrom14'
        :param orientation:  It can be any of ('+-', '-+', '++', '--') or ('Inner', 'Outer', 'Right', 'Left')
        Further, 'p' instead of '+', and 'n' or 'm' instead of '-' are acceptable.
        :return:
        """
        working_file = FILES.get_working_files(chrom, orientation, Files.S_NUC_INTER).squeeze()
        if working_file.empty:
            return None

        nuc_interas_file = working_file.subdir / working_file.file
        if not nuc_interas_file.exists():
            nuc_interas_file = working_file.subdir / working_file.file_zip  # check if zip version exists
            if not nuc_interas_file.exists():
                return None

        return nuc_interas_file

    @classmethod
    def read_nuc_interas_file(cls, nuc_interas_file) -> pd.DataFrame:

        try:
            df_nuc_interas = pd.read_csv(nuc_interas_file, sep=S.FIELD_SEPARATOR, comment=S.COMMENT_CHAR,
                                         usecols=S.OUTPUT_COLS__NUC_INTERAS)
        except pd.errors.ParserError as parser_error:
            raise InputFileError(nuc_interas_file) from parser_error

        # Not checking for the following condition for seek of efficiency; assumed as a precondition
        # if not df_nuc_interas.eval("chrom1 == chrom2").all():
        #     raise RuntimeError(f"Not all chrom1 and chrom2 values are equal in {nuc_interas_file}")

        df_nuc_interas['counts'] = 1  # initially every row is counted as 1
        df_nuc_interas['dist_sum'] = (df_nuc_interas['pos1'] - df_nuc_interas['pos2']).abs()  # calc location distance

        df_nuc_interas = df_nuc_interas[['nuc_id1', 'nuc_id2', 'counts', 'dist_sum']]  # drop extra columns

        df_nuc_interas = df_nuc_interas.apply(pd.to_numeric, downcast='unsigned')  # to save space

        return df_nuc_interas

    @classmethod
    def _convert_interas_to_nuc_interas(
            cls, nuc_interas_outfile, interas_infile, indexed_nucs_infile, chrom, keep_cache):
        """This method is meant to be used by multiprocessing.Pool"""

        df_interas = cls.__read_interas_file(interas_infile)
        df_nuc_interas = cls.__create_nuc_interas(df_interas, chrom, indexed_nucs_infile)
        df_nuc_interas.to_csv(nuc_interas_outfile, sep=S.FIELD_SEPARATOR, index=False, header=True)
        if not keep_cache:
            try:
                Path(interas_infile).unlink()  # remove the cached input_file to reduce runtime space requirements
            except Exception as e:
                # LOGGER.error(e)  # The LOGGER object is not available processes, potentially running on other machines
                print(e, file=sys.stderr)  # used instead of LOGGER.error(e) due to multiprocessing

    @classmethod
    def _convert_interas_to_nuc_interas__in_chunks(
            cls, nuc_interas_outfile, interas_infile, indexed_nucs_infile, chrom, keep_cache):
        """This method is meant to be used by multiprocessing.Pool"""
        # this method is slower than _convert_interas_to_nuc_interas

        schema = cls.schema()

        try:
            interas_chunks = pd.read_csv(interas_infile, sep=S.FIELD_SEPARATOR, comment=S.COMMENT_CHAR,
                                         usecols=range(schema.columns.size), chunksize=10_000_000, header=None)
        except pd.errors.ParserError as parser_error:
            raise InputFileError(interas_infile) from parser_error

        df_nuc_interas_chunks = []
        for df_interas_chunks in interas_chunks:
            df_interas_chunks[0] = df_interas_chunks[0].astype(pd.CategoricalDtype(ordered=True))  # chrom1
            df_interas_chunks[2] = df_interas_chunks[2].astype(pd.CategoricalDtype(ordered=True))  # chrom2
            df_interas_chunks.columns = schema.columns
            df_interas_chunks = df_interas_chunks.rename_axis(schema.index.name)

            df_nuc_interas_chunk = cls.__create_nuc_interas(df_interas_chunks, chrom, indexed_nucs_infile)
            df_nuc_interas_chunks.append(df_nuc_interas_chunk)

        df_nuc_interas = pd.concat(df_nuc_interas_chunks, axis=0, ignore_index=True, sort=False)
        df_nuc_interas = df_nuc_interas.sort_values(['nuc_id1', 'nuc_id2'], ignore_index=True)

        df_nuc_interas.to_csv(nuc_interas_outfile, sep=S.FIELD_SEPARATOR, index=False, header=True)
        if not keep_cache:
            try:
                Path(interas_infile).unlink()  # remove the cached input_file to reduce runtime space requirements
            except Exception as e:
                # LOGGER.error(e)  # The LOGGER object is not available processes, potentially running on other machines
                print(e, file=sys.stderr)  # used instead of LOGGER.error(e) due to multiprocessing

    def __create_nuc_interas_files__multiprocessing(self):

        subdir_inter = FILES.get_subdir(Files.S_INTER)
        subdir_nuc_inter = FILES.get_subdir(Files.S_NUC_INTER)

        # decide if need to continue
        files_needing_refresh = FILES.get_files_needing_refresh(Files.S_NUC_INTER)
        if files_needing_refresh.empty:
            LOGGER.info('No nucleosome interaction files were updated')
            return

        indexed_nucs_file = FILES.indexed_nucs_file(fullpath=True)
        to_do_list = []

        for _, file_info in files_needing_refresh.iterrows():
            # if not file_info.needs_refresh:  # get_files_needing_refresh has already checked for this
            #     continue
            if file_info.chrom not in self.__chrom_list:
                continue
            input_file = subdir_inter / file_info.file
            if not input_file.exists():
                input_file = subdir_inter / file_info.file_zip  # check if the zip version exists
                if not input_file.exists():
                    # LOGGER.warning(f'Skipping cache file. File not found: {file_info.file}')
                    continue
            output_file = file_info.file_zip if FILES.zipped else file_info.file
            output_file = subdir_nuc_inter / output_file

            to_do_list.append((output_file, input_file, indexed_nucs_file, file_info.chrom, FILES.keep_cache))
            # LOGGER.info(f'Queued: creating nuc interaction file: {output_file}')

        if len(to_do_list) > 0:
            n_processes = min(FILES.n_processes, len(to_do_list))
            LOGGER.info(f'In total, creating {len(to_do_list)} nuc interaction files using {n_processes} '
                        f'parallel processes:\n\t{subdir_nuc_inter}\nThis may take a while. Please wait')
            with multiprocessing.Pool(n_processes) as pool:
                pool.starmap(self._convert_interas_to_nuc_interas, to_do_list)

        FILES.set_needs_refresh_for_all_files(False, Files.S_NUC_INTER)

    @classmethod
    def __read_interas_file(cls, input_file) -> pd.DataFrame:

        schema = cls.schema()

        try:
            df_inter = pd.read_csv(input_file, sep=S.FIELD_SEPARATOR, comment=S.COMMENT_CHAR, header=None, )
            # usecols=range(schema.columns.size))  # not needed as strands are not saved anymore
        except pd.errors.ParserError as parser_error:
            raise InputFileError(input_file) from parser_error

        df_inter[0] = df_inter[0].astype(pd.CategoricalDtype(ordered=True))  # chrom1
        df_inter[2] = df_inter[2].astype(pd.CategoricalDtype(ordered=True))  # chrom2
        df_inter.columns = schema.columns
        df_inter = df_inter.rename_axis(schema.index.name)

        return df_inter

    @staticmethod
    def __create_nuc_interas(df_inter: pd.DataFrame, chrom, indexed_nucs_file=None) -> pd.DataFrame:

        if indexed_nucs_file is None:
            df_nucs = Nucs().df_nucs  # loads indexed nucs from the working directory
        else:  # indexed_nucs_file is given when loading from a separate process
            df_nucs = Nucs.read_csv(indexed_nucs_file, read_index=True)

        # The following code for Step 1 was moved from __init__() to support multiprocessing.
        # To avoid data sharing between processes, each process needs to do some redundant job for Step 1

        # Step 1:
        # Both Nucs and Interactions will use the same MultiIndex ['chrom', 'pos']
        # Make a MultiIndex ['chrom', 'pos'] for Nucs to make them similar to Interactions, simplifying their merger
        df_nucs['pos'] = df_nucs.start  # copy 'start' to 'pos' to make df_nucs and later df_inter alike
        # df_nucs = df_nucs.query(f"chrom == '{chrom}'")  # too slow
        df_nucs = df_nucs[df_nucs['chrom'] == chrom]
        df_nucs = df_nucs.reset_index().set_index(['chrom', 'pos'])  # .sort_index()

        df_inter_stack = df_inter.stack(level='side').reset_index()  # 'side' becomes a col, with 'side1' or 'side2'
        df_inter_stack = df_inter_stack.set_index(['chrom', 'pos'])  # .sort_index()

        # Step 2:
        # Combine Nucs and Interactions (which use the same MultiIndex ['chrom', 'pos'] now)
        index_names = ['orig'] + df_inter_stack.index.names  # 'orig' column will be 'nuc' or 'inter'
        df_nuc_interas = pd.concat((df_nucs, df_inter_stack), keys=['nuc', 'inter'], names=index_names, sort=False)
        df_nuc_interas = df_nuc_interas.reset_index(level='orig').sort_index(level='pos')  # chrom is unique, so sorted
        df_nuc_interas = df_nuc_interas.fillna(method='ffill')  # Copy nuc_id,start,end for each inter from prev nuc
        df_nuc_interas = df_nuc_interas.query('pos <= end')  # nuc.start is already matched; nuc.end needs to match too
        df_nuc_interas = df_nuc_interas.drop(columns=['start', 'end'])  # they are not needed anymore

        # Step 3:
        # Recreate __nuc_interas, but now with nuc_id for side1 and side2
        df_nuc_interas = df_nuc_interas.query('orig == "inter"').drop(columns=['orig'])  # keep 'inter' rows only
        df_nuc_interas = df_nuc_interas.reset_index().set_index(['inter_id', 'side'])  # .sort_index()
        df_nuc_interas = df_nuc_interas.unstack(level='side').dropna().reset_index(drop=True)
        df_nuc_interas = df_nuc_interas.apply(pd.to_numeric, downcast='unsigned', errors='ignore')  # to save space

        # Replace MultiIndex columns with single level, and reorder the columns
        df_nuc_interas.columns = df_nuc_interas.columns.map(
            lambda col_level: ''.join((col_level[0], col_level[1].lstrip('side'))))
        df_nuc_interas = df_nuc_interas[S.OUTPUT_COLS__NUC_INTERAS]
        # ['chrom1', 'pos1', 'chrom2', 'pos2', 'nuc_id1', 'nuc_id2']

        df_nuc_interas = df_nuc_interas.sort_values(['nuc_id1', 'nuc_id2'], ignore_index=True)

        return df_nuc_interas

    def __split_interas_file__multiprocessing(self):
        # decide if need to continue
        outfiles = FILES.get_files_needing_refresh__dict(Files.S_INTER)
        if not outfiles:
            LOGGER.info('No interaction cache files were updated')
            return

        InterasFileSplitter.split_interas_file(outfiles, self.__interas_file, self.__chrom_list)

        self.__warn_for_missing_interas_files()

        FILES.set_needs_refresh_for_all_files(False, Files.S_INTER)

    def __split_interas_file__single_process(self):
        # This is slow! Recommended to use __split_interas_file__multiprocessing() when possible.

        # decide if need to continue
        outfile_names = FILES.get_files_needing_refresh__dict(Files.S_INTER)
        if not outfile_names:
            LOGGER.info('No interaction cache files were updated')
            return

        chroms_categorical = pd.CategoricalDtype(categories=self.__chrom_list)

        header_line_num, used_cols = self.validate_interas_file(self.__interas_file)

        chrom1, chrom2 = S.CHROM_COLS.values  # ['chrom1', 'chrom2']
        strand1, strand2 = S.STRAND_COLS.values  # ['strand1', 'strand2']

        try:
            chunks = pd.read_csv(self.__interas_file, sep=S.FIELD_SEPARATOR, comment=S.COMMENT_CHAR,
                                 usecols=used_cols.index, chunksize=1_000_000, header=header_line_num)
        except pd.errors.ParserError as parser_error:
            raise InputFileError(self.__interas_file) from parser_error

        open_files = {}
        with contextlib.ExitStack() as cm:  # allows to open an unknown number of files
            for df_chunk in chunks:
                df_chunk.columns = S.USED_COLS
                df_chunk[[chrom1, chrom2]] = df_chunk[[chrom1, chrom2]].astype(chroms_categorical)  # for fast filtering
                df_chunk = df_chunk.dropna()  # drops any rows not containing chroms stated in chroms_categorical
                df_chunk = df_chunk[df_chunk[chrom1] == df_chunk[chrom2]]  # fast row filtering when used with cat cols

                groups = df_chunk.groupby([chrom1, strand1, strand2])
                for chrom_and_strands, df in groups:
                    outfile_handler = open_files.get(chrom_and_strands)
                    if outfile_handler is None:
                        outfile_name = outfile_names.get(chrom_and_strands)  # e.g., {(chrII, +, -): 'path/to/outfile'}
                        if outfile_name is None:
                            continue
                        opener = gzip.open if Files.is_file_zipped(outfile_name) else open
                        outfile_handler = cm.enter_context(opener(outfile_name, 'wt'))
                        open_files[chrom_and_strands] = outfile_handler

                    df = df.iloc[:, :len(S.LOCATION_COLS)]  # keeping: ch1 pos1 ch2 pos2,  dropping the strands
                    df.to_csv(outfile_handler, sep=S.FIELD_SEPARATOR, index=False, header=False)

        self.__warn_for_missing_interas_files()
        FILES.set_needs_refresh_for_all_files(False, Files.S_INTER)

    @staticmethod
    def __warn_for_missing_interas_files():
        # if any files still need refresh, it means they have not been created, so no interactions existed for them
        outfiles = FILES.get_working_files(state=Files.S_INTER)
        if outfiles.empty:
            return  # no files expected!

        chroms_with_all_strands_missing = []
        chroms_with_some_strands_missing = []
        missing_files = outfiles[~ ((outfiles.subdir / outfiles.file).apply(Path.exists) |
                                    (outfiles.subdir / outfiles.file_zip).apply(Path.exists))].copy()
        if len(missing_files) == 0:
            return  # no more actions necessary as no files are missing
        missing_files['strand1'], missing_files['strand2'] = zip(*missing_files['strands'])  # split tuple into two cols
        missing_files = missing_files[['chrom', 'strand1', 'strand2']]  # ignore other columns
        for chrom, strands in missing_files.groupby(['chrom'], as_index=False):
            if len(strands) < 4:
                chroms_with_some_strands_missing.append(strands)
            else:
                chroms_with_all_strands_missing.append(chrom)

        if len(chroms_with_all_strands_missing) > 0 or len(chroms_with_some_strands_missing) > 0:
            LOGGER.warning('NOTE:')
            if len(chroms_with_all_strands_missing) > 0:
                LOGGER.warning('No interactions have been found for the following chromosome(s) (with any strands):')
                LOGGER.warning(', '.join(chroms_with_all_strands_missing))
            if len(chroms_with_some_strands_missing) > 0:
                LOGGER.warning('No interactions have been found for the following chromosome(s) and strands:')
                LOGGER.warning(pd.concat(chroms_with_some_strands_missing).sort_index().to_string())

    @classmethod
    def validate_interas_file(cls, interas_file, num_lines_to_check=1000):

        # skip the comment lines at the beginning of the file
        with open(interas_file, 'r') as in_f:
            for line_num, line in enumerate(in_f):
                if not line.startswith(S.COMMENT_CHAR):
                    break

        used_cols = S.USED_COLS.copy()
        line = pd.Series(line.strip().split(S.FIELD_SEPARATOR))  # convert line to series for easier comparison for cols
        if line.size < used_cols.size:
            raise RuntimeError(f'Not enough number of columns in: {interas_file}')

        found_header = line[used_cols.index].equals(used_cols)
        if not found_header and line.size > used_cols.size:
            used_cols.index += 1  # shift index to check if there is one extra column at the beginning
            found_header = line[used_cols.index].equals(used_cols)

        header_line_num = 0 if found_header else None  # Use 0 for header on 1st non-comment line and None for no header

        try:
            df = pd.read_csv(interas_file, sep=S.FIELD_SEPARATOR, comment=S.COMMENT_CHAR,
                             nrows=num_lines_to_check, header=header_line_num)
        except pd.errors.ParserError as parser_error:
            raise InputFileError(interas_file) from parser_error

        df = df.rename(used_cols, axis=1)  # this will only have an effect if headers are numbers
        df = df[used_cols]  # ensures used_cols exist with no extra cols

        checks = {'strand1': ['+', '-'], 'strand2': ['+', '-']}  # it could also check for chrom_list in the future
        dtypes = {S.CHROM_COLS.iloc[0]: np.dtype('O'), S.POS_COLS.iloc[0]: np.dtype('int64'),
                  S.CHROM_COLS.iloc[1]: np.dtype('O'), S.POS_COLS.iloc[1]: np.dtype('int64'),
                  S.STRAND_COLS.iloc[0]: np.dtype('O'), S.STRAND_COLS.iloc[1]: np.dtype('O')}

        try:
            df.astype(dtypes)
        except ValueError:
            type_test_failed = True
        else:
            type_test_failed = False

        if type_test_failed or not df[S.STRAND_COLS].isin(checks).all().all():
            raise RuntimeError(f'Invalid columns in: {interas_file}')

        return header_line_num, used_cols

    def copy(self):
        return copy.copy(self)

    def __len__(self):
        return self.__len  # raise NotImplementedError()

    def __repr__(self):
        return self.name + ' (nucleosome interactions)'


class NucInteraMatrix:
    """ Nucleosome-Nucleosome Interaction Matrix

        Using a sparse matrix where elements are recorded as ijv with
        i for rows, j for cols, and v for values
    """

    def __init__(self, nuc_interas: NucInteras, nucs: Nucs,
                 name: str = None, create=True, keep_nucs_cols=True, verbose=True):

        self.__nucs = nucs
        self.__nuc_interas = nuc_interas
        self.__name = name

        if not create:
            return

        if verbose:
            LOGGER.info('\nPreparing Nuc Interaction Matrices')

        start = timer()
        refreshed_files = self.__create_nuc_intera_matrix_files__multiprocessing(keep_nucs_cols)

        if verbose:
            # if refreshed_files == 0:
            #     LOGGER.info('No Nuc Interaction Matrix files were updated')
            LOGGER.info(f'Done! Preparing Nuc Interaction Matrices finished in {Files.time_diff(start)}')

            LOGGER.info(f'\nNucleosome interaction matrices ready with total size: '
                        f'{len(self):,} nucs * {len(self):,} nucs')
            # LOGGER.info(self)

    def read_submatrices_for_orient(
            self, chrom: str, start_region: int, end_region: int, orientation='all') -> Optional[pd.DataFrame]:
        """
        :param chrom: e.g. 'II' or 'chrom14'
        :param start_region: positive integer
        :param end_region: positive integer
        :param orientation:  It can be any of ('+-', '-+', '++', '--') or ('inner', 'outer', 'right', 'left') or
        ('tandem', 'all'). Further, 'p' instead of '+', and 'n' or 'm' instead of '-' are acceptable.
        Default value is 'all' which returns all orientations.
        :return: a submatrix for the region in a DataFrame
        """

        nuc_intera_matrix_files = self.get_nuc_intera_matrix_files(chrom, orientation)
        ijv = self.read_submatrices_from_files(nuc_intera_matrix_files, chrom, start_region, end_region)
        return ijv

    def read_submatrices_from_files(
            self, nuc_intera_matrix_files, chrom, start_region: int, end_region: int) -> Optional[pd.DataFrame]:

        nucs = self.__nucs.find_nucs_in_region(chrom, start_region, end_region)
        if nucs is None:
            return None

        nucs = nucs.df_nucs

        # nucs index is sorted
        min_nuc_id = nucs.index[0]
        max_nuc_id = nucs.index[-1]

        ijv = self.read_nuc_intera_matrix(nuc_intera_matrix_files)
        # ijv = ijv.query(f'{min_nuc_id} <= nuc_id1 <= {max_nuc_id} and '
        #                 f'{min_nuc_id} <= nuc_id2 <= {max_nuc_id}')   # this is very slow!
        idx = (min_nuc_id <= ijv['nuc_id1']) & (ijv['nuc_id1'] <= max_nuc_id) & \
              (min_nuc_id <= ijv['nuc_id2']) & (ijv['nuc_id2'] <= max_nuc_id)
        ijv = ijv[idx]

        ijv = self.__append_nucs_cols_if_missing(ijv, nucs)

        return ijv

    @classmethod
    def write_submatrix_for_orients(cls, output_file: Path, nuc_intera_matrix_files: list, interas_file, nucs_file,
                                    chrom: str, start_region: int, end_region: int) -> bool:
        """This method is meant to be used by multiprocessing.Pool"""

        nucs = Nucs(nucs_file=nucs_file, read_index=True, verbose=False)
        nuc_interas = NucInteras(interas_file=interas_file, nucs=nucs, create=False, refresh=False, verbose=False)
        nuc_intera_matrix = NucInteraMatrix(nuc_interas=nuc_interas, nucs=nucs, create=False, verbose=False)

        submatrix_ijv = nuc_intera_matrix.read_submatrices_from_files(
            nuc_intera_matrix_files, chrom, start_region, end_region)

        if submatrix_ijv is None or len(submatrix_ijv) == 0:
            return False

        submatrix_ijv.to_csv(output_file, S.FIELD_SEPARATOR, header=True, index=False)

        return True

    @classmethod
    def write_submatrices_for_all_orients__multiprocessing(
            cls, output_file: Path, chrom: str, start_region: int, end_region: int) -> dict:

        nucs_file = FILES.working_dir / FILES.indexed_nucs_file()
        interas_file = FILES.input_interas_file_name

        output_files = {orient: output_file.with_stem(f'{output_file.stem}_{orient}') for orient in
                        Files.orientations() + Files.combined_orientations()}  # inner outer right left tandem all

        to_do_list = []
        for orient in Files.orientations():  # 'inner', 'outer', 'right', 'left'
            nuc_intera_matrix_files = cls.get_nuc_intera_matrix_files(chrom, orient)
            to_do_list.append((output_files[orient], nuc_intera_matrix_files, interas_file, nucs_file,
                               chrom, start_region, end_region))

        n_processes = min(FILES.n_processes, len(to_do_list))
        LOGGER.info(f'\nSelecting {len(to_do_list)} nuc interaction submatrices for plotting using {n_processes} '
                    f'parallel processes:\n\t{output_file.parent}\nThis may take a while. Please wait')
        with multiprocessing.Pool(n_processes) as pool:
            success = pool.starmap(cls.write_submatrix_for_orients, to_do_list)

        orient_2_submatrix = dict()
        for i, orient in enumerate(Files.orientations()):  # 'inner', 'outer', 'right', 'left'
            output_with_orient, *_ = to_do_list[i]
            if success[i]:
                submatrix = pd.read_csv(output_with_orient, sep=S.FIELD_SEPARATOR)
                orient_2_submatrix[orient] = submatrix
            else:
                LOGGER.info(f'No records found for {orient}:  {output_with_orient}')

        for combined_orient in Files.combined_orientations():  # 'tandem', 'all'
            orient_list = Files.break_combined_orientation(combined_orient)  # e.g. 'tandem' -> ['right', 'left']
            submatrices = [orient_2_submatrix[orient] for orient in orient_list if orient in orient_2_submatrix]
            if len(submatrices) == 0:
                LOGGER.info(f'No records found for {combined_orient}:  {output_files[combined_orient]}')
            else:
                submatrix = cls.combine_nuc_intera_matrix(submatrices)
                submatrix.to_csv(output_files[combined_orient], S.FIELD_SEPARATOR, header=True, index=False)
                orient_2_submatrix[combined_orient] = submatrix

        if len(orient_2_submatrix) > 0:
            width = len('Interactions')
            bar = '-' * (width * 3)
            LOGGER.info('\n' + bar)
            LOGGER.info('Non-Zero Nuc')
            LOGGER.info('Interactions  Submatrix')
            LOGGER.info(bar)
            for orient, submatrix in orient_2_submatrix.items():
                LOGGER.info(f'{len(submatrix):>{width - 3}}     {output_files[orient]}')
            LOGGER.info(bar + '\n')

        return orient_2_submatrix

    @staticmethod
    def extract_norm_distance(matrix_ijv: pd.DataFrame) -> Optional[int]:
        """Extract norm_distance value from the column name starting with 'norm_'"""
        norm_col = matrix_ijv.filter(regex='^norm_', axis=1).columns.to_list()  # find 'norm_' column
        if len(norm_col) != 1:
            return None
        norm_col = norm_col[0]  # now there is exactly one col
        norm_distance = int(re.search('^norm_(.*)$', norm_col).group(1))
        return norm_distance

    @staticmethod
    def get_nuc_intera_matrix_files(chrom: str, orientation='all') -> list:
        """
        :param chrom: e.g. 'II' or 'chrom14'
        :param orientation:  It can be any of ('+-', '-+', '++', '--') or ('inner', 'outer', 'right', 'left') or
        ('tandem', 'all'). Further, 'p' instead of '+', and 'n' or 'm' instead of '-' are acceptable.
        Default value is 'all' which returns all orientations.
        :return: a list of matrix files
        """

        working_files = FILES.get_working_files(chrom, orientation, Files.S_MATRIX)
        if working_files.empty:
            return []

        matrix_files = []
        for _, file_info in working_files.iterrows():
            matrix_file = file_info.subdir / file_info.file
            if not matrix_file.exists():
                matrix_file = file_info.subdir / file_info.file_zip  # check if zip version exists
                if not matrix_file.exists():
                    continue

            matrix_files.append(matrix_file)

        return matrix_files

    @classmethod
    def read_nuc_intera_matrix(cls, nuc_intera_matrix_files: list) -> pd.DataFrame:
        """Reads and combines multiple matrix files and returns a single ijv sparse matrix"""

        if not nuc_intera_matrix_files:  # if None or empty list
            return cls.__EMPTY

        nuc_intera_matrices = []
        for input_file in nuc_intera_matrix_files:
            try:
                ijv = pd.read_csv(input_file, sep=S.FIELD_SEPARATOR, comment=S.COMMENT_CHAR)
            except pd.errors.ParserError as parser_error:
                raise InputFileError(input_file) from parser_error

            nuc_intera_matrices.append(ijv)

        return cls.combine_nuc_intera_matrix(nuc_intera_matrices)

    @classmethod
    def combine_nuc_intera_matrix(cls, nuc_intera_matrices: list) -> pd.DataFrame:

        ijv = pd.concat(nuc_intera_matrices, axis=0, ignore_index=True, sort=False)

        norm_distance = cls.extract_norm_distance(ijv)
        if norm_distance is None:
            raise RuntimeError(' '.join(f"""
                There most be exactly one norm column. 
                Please consider refreshing your working directory using "prepare --refresh" """.split()))
        S.set_norm_col_name(norm_distance)

        # for groupby aggregates: set 'first' for all cols, then overwrite that for some cols
        aggregates = dict(zip(ijv.columns, ['first'] * len(ijv.columns)))
        aggregates['counts'] = 'sum'
        aggregates['dist_sum'] = 'sum'
        aggregates['dist_avg'] = 'mean'
        aggregates[S.get_norm_col_name()] = 'sum'

        ijv = ijv.groupby(['nuc_id1', 'nuc_id2'], as_index=False).agg(aggregates)

        ijv = ijv.apply(pd.to_numeric, downcast='unsigned', errors='ignore')  # for space efficiency

        return ijv

    @classmethod
    def _convert_nuc_interas_to_matrices(
            cls, output_file, nuc_interas_file, indexed_nucs_file, norm_distance, keep_cache):
        """This method is meant to be used by multiprocessing.Pool"""

        df_nuc_intera_matrix = cls.__create_nuc_intera_matrix(nuc_interas_file, indexed_nucs_file, norm_distance)
        # print('.', end='')  # print a dot to show progress between to time consuming tasks
        df_nuc_intera_matrix.to_csv(output_file, sep=S.FIELD_SEPARATOR, index=False, header=True)
        if not keep_cache:
            try:
                Path(nuc_interas_file).unlink()  # remove the cached input_file to reduce runtime space requirements
            except Exception as e:
                # LOGGER.error(e)  # The LOGGER object is not available processes, potentially running on other machines
                print(e, file=sys.stderr)  # used instead of LOGGER.error(e) due to multiprocessing

    def __create_nuc_intera_matrix_files__multiprocessing(self, keep_nucs_cols) -> int:
        # decide if need to continue
        files_needing_refresh = FILES.get_files_needing_refresh(Files.S_MATRIX)
        if files_needing_refresh.empty:
            return 0
        if FILES.norm_distance is None:  # so we are not here via 'prepare' but some files are missing
            raise RuntimeError('Some matrix files are missing. '
                               'Please consider refreshing your working directory using "prepare --refresh"')

        indexed_nucs_file = FILES.indexed_nucs_file(fullpath=True) if keep_nucs_cols else None
        to_do_list = []

        for _, file_info in files_needing_refresh.iterrows():
            if file_info.chrom not in self.chrom_list:
                continue
            output_file = file_info.file_zip if FILES.zipped else file_info.file
            output_file = file_info.subdir / output_file

            nuc_interas_file = NucInteras.get_nuc_interas_file(file_info.chrom, file_info.orientation)
            if nuc_interas_file is None:
                continue

            to_do_list.append((output_file, nuc_interas_file, indexed_nucs_file, FILES.norm_distance, FILES.keep_cache))
            # LOGGER.info(f'Queued: creating nuc interaction Matrix file: {output_file}')

        if len(to_do_list) > 0:
            subdir_matrix = FILES.get_subdir(Files.S_MATRIX)
            n_processes = min(FILES.n_processes, len(to_do_list))
            LOGGER.info(f'In total, creating {len(to_do_list)} Matrix files using {n_processes} parallel processes:'
                        f'\n\t{subdir_matrix}\nThis may take a while. Please wait')
            with multiprocessing.Pool(n_processes) as pool:
                pool.starmap(self._convert_nuc_interas_to_matrices, to_do_list)

        FILES.set_needs_refresh_for_all_files(False, Files.S_MATRIX)

        n_files_refreshed = len(to_do_list)
        return n_files_refreshed

    @classmethod
    def __create_nuc_intera_matrix(
            cls, nuc_interas_file, indexed_nucs_file=None, norm_distance=S.NORM_DISTANCE) -> pd.DataFrame:

        df_nuc_intera = NucInteras.read_nuc_interas_file(nuc_interas_file)

        # Our goal is to create a ijv sparse matrix: i:nuc_id1, j:nuc_id2, v:counts
        ijv = df_nuc_intera

        # make ijv symmetric by swapping cols i,j to j,i:  'nuc_id2', 'nuc_id1'
        swapped = pd.DataFrame(ijv.values, columns=['nuc_id2', 'nuc_id1', 'counts', 'dist_sum'])
        ijv = pd.concat([ijv, swapped], sort=False)  # append ij rows with their swapped ji rows

        ijv = ijv.groupby(['nuc_id1', 'nuc_id2'], as_index=False).sum()  # calc counts and dist_sum

        diagonal_ids = (ijv.nuc_id1 == ijv.nuc_id2)
        ijv.loc[diagonal_ids, 'counts'] //= 2  # divide by 2 to fix duplicated diagonal values
        ijv.loc[diagonal_ids, 'dist_sum'] //= 2  # divide by 2 to fix duplicated diagonal values

        ijv = ijv.apply(pd.to_numeric, downcast='unsigned')  # for space efficiency

        ijv = cls.__norm_nuc_intera_matrix(ijv, norm_distance)  # add columns for normalized values

        if indexed_nucs_file:
            df_nucs = Nucs.read_csv(indexed_nucs_file, read_index=True)
            ijv = cls.__append_nucs_cols_if_missing(ijv, df_nucs)

        # keep existing cols but put them in correct order
        col_order = pd.Index(S.OUTPUT_COLS__MATRIX).intersection(ijv.columns)
        ijv = ijv[col_order]  # same columns, just in the order specified in S.OUTPUT_COLS__MATRIX

        return ijv

    @staticmethod
    def __append_nucs_cols_if_missing(ijv: pd.DataFrame, nucs: pd.DataFrame) -> pd.DataFrame:
        if set(S.OUTPUT_NUCS_COLS__MATRIX).issubset(ijv.columns):
            return ijv  # columns already exist

        ijv = pd.merge(ijv, nucs, how='left', left_on='nuc_id1', right_on='nuc_id')
        ijv = pd.merge(ijv, nucs, how='left', left_on='nuc_id2', right_on='nuc_id', suffixes=('1', '2'))

        return ijv

    @staticmethod
    def __norm_nuc_intera_matrix(ijv: pd.DataFrame, norm_distance) -> pd.DataFrame:
        # Adds columns for normalized values
        # Normalized_values = Observed / Expected_nonzero
        # Expected values are the sum per genomic distance divided by sum of non-zero contacts.
        # Only non-zero values are used to compute the expected values per genomic distance.
        # https://hicexplorer.readthedocs.io/en/latest/content/tools/hicTransform.html#observed-expected-non-zero

        # # Alternative expressions of the high-level goal:
        # Find all nucs with about same genomic distances, sum up the counts, divide by number of rows (m 6.57)
        # That is, for each genomic distance rolling window, calculate the number of interactions per nucleosome contact
        # That is, for interactions with about the same genomic distance, sum counts and divide by row (contacts)
        # That is, average of interaction counts between nuc-pairs whose distance is within the query nuc-pair average
        # distance +/- a given distance

        # Algorithm:
        # Datetime rolling for a ragged (meaning not-a-regular frequency), time-indexed DataFrame can be used
        # for the above high-level goal. First, convert dist_avg to datetime to make use of ragged rolling feature
        # But as of Pandas 1.3.1, this ragged rolling feature only works in one direction and not both ways.
        # To tackle this limitation, we use rolling twice, once in ascending order of dist_avg, and for a second
        # time we negate the dist_avg values and rolling them in reverse order.
        # Finally, combining the two directions, we get rolling results in both directions.

        ijv['dist_avg'] = ijv['dist_sum'] / ijv['counts']
        ijv = ijv.sort_values(['dist_avg', 'nuc_id1', 'nuc_id2'])  # sort to enable the ragged rolling
        ijv['for_rolling'] = pd.to_datetime(ijv['dist_avg'], unit='s')
        ijv = ijv.set_index('for_rolling')  # make the sorted ijv time-indexed to enable the ragged rolling
        same_dist_nucs = ijv['counts'].rolling(f'{norm_distance}s')  # size in 'seconds' for ragged rolling
        ijv['num_same_dist_interas'] = same_dist_nucs.sum()
        ijv['num_same_dist_nuc_pairs'] = same_dist_nucs.count()

        # Now, do the same in the reverse order and combine the sum and count values to the previous ones
        ijv = ijv[::-1]  # reverse the row order
        ijv['for_rolling'] = pd.to_datetime(-1 * ijv['dist_avg'], unit='s')  # negate dist_avg values
        ijv = ijv.set_index('for_rolling')  # re-set the negated dist_avg values
        same_dist_nucs = ijv['counts'].rolling(f'{norm_distance}s')  # same as above
        ijv['num_same_dist_interas'] += same_dist_nucs.sum() - ijv['counts']  # combine sum with previous direction
        ijv['num_same_dist_nuc_pairs'] += same_dist_nucs.count() - 1  # combine counts with previous direction

        ijv['expected_counts'] = ijv['num_same_dist_interas'] / ijv['num_same_dist_nuc_pairs']

        S.set_norm_col_name(norm_distance)
        ijv[S.get_norm_col_name()] = ijv['counts'] / ijv['expected_counts']

        # ijv = ijv[::-1]  # reverse the row order
        ijv = ijv.sort_values(['nuc_id1', 'nuc_id2'], ignore_index=True)  # back to original order

        return ijv

    @staticmethod
    def schema() -> pd.DataFrame:
        """:returns: an empty DataFrame table with appropriate columns."""
        return NucInteraMatrix.__EMPTY

    __EMPTY = pd.DataFrame(columns=S.OUTPUT_COLS__MATRIX)

    @property
    def chrom_list(self) -> list:
        return self.__nucs.chrom_list

    @property
    def name(self) -> str:
        if self.__name is None:
            name = {
                'nucleosomes': self.__nucs.name if self.__nucs.name else 'unknown',
                'interactions': self.__nuc_interas.name if self.__nuc_interas.name else 'unknown',
            }
            self.__name = str(name)

        return self.__name

    @name.setter
    def name(self, name: str):
        self.__name = name

    def __len__(self):
        return len(self.__nucs)

    def __repr__(self):
        return 'Nucleosome interaction matrices created using: \n' + self.name


class CLI:
    """Command Line Interface"""

    @classmethod
    def handler_pairtools_command(cls, command, pairtools_arguments):
        pairtools_command = command + ' ' + ' '.join(pairtools_arguments)
        print(pairtools_command)

        try:  # first make sure that bash is available
            process = subprocess.Popen(
                '/bin/bash', stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        except FileNotFoundError:
            print('Could not find "bash" on this system, so cannot proceed. '
                  'Please run pairtools directly on the commandline. Type "pairtools -h" for more information.')
        else:
            out, err = process.communicate(pairtools_command)
            if err:
                raise RuntimeError(err)

            print(out)

    @classmethod
    def handler_prepare_command(cls, command, chroms_file, nucs_file, interas_file, working_dir,
                                norm_distance, keep_nucs_cols, keep_cache, n_processes, refresh, zipped=False):

        start = timer()
        start_time = datetime.datetime.now().isoformat(timespec='seconds')

        FILES.reset_all(base_file_name=interas_file, working_dir=working_dir, norm_distance=norm_distance,
                        keep_cache=keep_cache, n_processes=n_processes, zipped=zipped, refresh=refresh)

        log_file_name = FILES.working_dir / (FILES.base_file_name.stem + '.log')
        handler_log_file = logging.FileHandler(filename=log_file_name, mode='a')
        LOGGER.addHandler(handler_log_file)
        LOGGER.info(f'\n\n _ _ _ _ _ Started command: {command} _ _ _ _ _')
        vars_to_log = ['start_time', 'n_processes', 'norm_distance', 'keep_nucs_cols', 'refresh', 'keep_cache',
                       'zipped', 'chroms_file', 'nucs_file', 'interas_file', 'working_dir', 'log_file_name']
        if n_processes == 0:
            vars_to_log.remove('n_processes')
        LOGGER.info(Files.tabulate_dict(locals(), vars_to_log))  # log important arguments and variables
        # LOGGER.info(f'Logging to file {log_file_name}\n')

        nucs = Nucs(nucs_file=nucs_file, chroms_file=chroms_file)
        nuc_interas = NucInteras(interas_file=interas_file, nucs=nucs, refresh=refresh)
        nuc_intera_matrix = NucInteraMatrix(nuc_interas=nuc_interas, nucs=nucs, keep_nucs_cols=keep_nucs_cols)

        # Check integrity of the results
        missing_files = []
        working_files = FILES.get_working_files(Files.S_MATRIX)
        for _, f in working_files.iterrows():
            if not (f.subdir / f.file).exists() and not (f.subdir / f.file_zip).exists():
                missing_files.append(f.file)

        if len(missing_files) > 0:
            raise RuntimeError(
                'Error: Some expected matrix files are missing:'
                '\n'.join(missing_files) + '\n' +
                'Could not produce expected matrix files! \n'
                'Please consider refreshing your working directory using "prepare --refresh"')

        if len(nuc_intera_matrix) == 0:
            LOGGER.warning('\nWARNING: The produced matrix is empty!')

        cache_dir = FILES.working_dir / 'cache'
        LOGGER.info(f'\nPrepared working directory: {FILES.working_dir}')
        if FILES.keep_cache:
            LOGGER.info('\nTo save space, use the --no-cache flag to remove cache files automatically.')
            LOGGER.info(f'You can safely remove cache folder manually: {cache_dir}')
        else:
            try:
                if cache_dir.exists():
                    shutil.rmtree(cache_dir)  # remove the cache dir
            except Exception as e:
                LOGGER.error(e)
                LOGGER.info(f'You can safely remove cache folder manually: {cache_dir}')

        LOGGER.info(f"\n _ _ _ _ _ Finished command {command} in {Files.time_diff(start)}\n")

    @classmethod
    def handler_plot_command_testing(cls, command, working_dir, chrom, start_region, end_region, prefix, save_only):
        # For testing only. Assumes the plot submatrix files have already been generated by handler_plot_command
        print(pprint.pformat(locals()))  # prints all the arguments

        FILES.reset_all(working_dir=working_dir, norm_distance=S.NORM_DISTANCE)  # norm_distance can be diff here
        files = list(FILES.working_dir.glob('my_plot*txt'))

        output_file = None
        submatrices = dict()
        for f in files:
            match = re.search('^(.*)_(.*).txt$', f.name)
            output_file = match.group(1)
            orient = match.group(2)
            submatrices[orient] = pd.read_csv(f, S.FIELD_SEPARATOR)

        # reorder submatrices in the predefined order
        submatrices = {orient: submatrices[orient] for orient in Files.orientations() + Files.combined_orientations()}

        cls.make_heat_map(submatrices, chrom, start_region, end_region, Path(output_file), save_only)
        print(files)

    @classmethod
    def handler_plot_command(cls, command, working_dir, chrom, start_region, end_region, prefix, save_only):
        start = timer()
        start_time = datetime.datetime.now().isoformat(timespec='seconds')

        FILES.reset_all(working_dir=working_dir)
        FILES.validate_working_dir()
        FILES.reset_working_files(chrom_list=[chrom], refresh=False)

        log_file_name = FILES.working_dir / (FILES.base_file_name.stem + '.log')
        handler_log_file = logging.FileHandler(filename=log_file_name, mode='a')
        LOGGER.addHandler(handler_log_file)
        LOGGER.info(f'\n\n _ _ _ _ _ Started command: {command} _ _ _ _ _')
        vars_to_log = ['start_time', 'chrom', 'start_region', 'end_region', 'prefix', 'working_dir', 'log_file_name']
        LOGGER.info(Files.tabulate_dict(locals(), vars_to_log))  # log important arguments and variables

        output_file = FILES.working_dir / f'{prefix}_{chrom}_{start_region}_{end_region}.txt'
        orient_2_submatrix = NucInteraMatrix.write_submatrices_for_all_orients__multiprocessing(
            output_file, chrom, start_region, end_region)

        # create the Heat Map:
        if len(orient_2_submatrix) > 0:
            cls.make_heat_map(orient_2_submatrix, chrom, start_region, end_region, output_file, save_only)
        else:
            LOGGER.warning('\n\nNOTE: No plots have been produced because no data is available. '
                           'You many want to consider expanding the chromosome region of study.')

        LOGGER.info(f"\n _ _ _ _ _ Finished command {command} in {Files.time_diff(start)}\n")

    @classmethod
    def make_heat_map(cls, matrices, chrom, start_region, end_region, output_file: Path, save_only):
        from bokeh.io import show, save
        from bokeh.layouts import layout
        from bokeh.models import ColorBar, Panel, Tabs, Range1d, LinearColorMapper  # , LogColorMapper
        from bokeh.models import HoverTool, PanTool, WheelZoomTool, BoxZoomTool, ResetTool, SaveTool
        from bokeh.plotting import figure, output_file as save_file
        from bokeh.palettes import Greys, Oranges, Greens, Blues, Purples, Reds

        output_file = output_file.with_suffix('.html')
        save_file(str(output_file))

        # seaborn: "the best sequential palettes will be perceptually uniform, meaning that the relative
        # discriminability of two colors is proportional to the difference between the corresponding data values"

        def set_fig_attr(fig: figure):
            fig.toolbar.active_scroll = fig.select_one(WheelZoomTool)
            fig.grid.grid_line_color = None
            fig.axis.axis_line_color = None
            fig.axis.major_tick_line_color = None
            fig.axis.major_label_standoff = 0
            fig.axis.major_label_text_font_size = '12pt'
            fig.xaxis.axis_label_text_font_size = '12pt'
            fig.yaxis.axis_label_text_font_size = '12pt'
            fig.yaxis.major_label_orientation = 'vertical'
            fig.title.text_font_size = '12pt'

        def r(colors, num_colors=256):  # Reverse the order going from light to dark, and select max num of colors
            return list(reversed(colors[num_colors]))

        palettes = dict(
            all=r(Oranges), tandem=r(Greys), right=r(Reds), left=r(Blues), inner=r(Purples), outer=r(Greens))
        # alphas = dict(all=.8, tandem=.4, right=.6, left=.6, inner=.5, outer=.5)

        full_ijv = pd.concat(matrices.values(), keys=matrices.keys(), names=['orient', None], sort=False)
        full_ijv['arcsinh_counts'] = np.arcsinh(full_ijv.counts)

        nuc_id_min = full_ijv[['nuc_id1', 'nuc_id2']].min().min()
        nuc_id_max = full_ijv[['nuc_id1', 'nuc_id2']].max().max()
        max_min = nuc_id_max - nuc_id_min
        nuc_id_min -= max_min * .05  # making min 5% smaller
        nuc_id_max += max_min * .05  # making max 5% larger

        x_range = Range1d(nuc_id_min, nuc_id_max, bounds='auto')
        y_range = Range1d(nuc_id_max, nuc_id_min, bounds='auto')  # flipped order: max, min

        tools = [PanTool(), WheelZoomTool(), BoxZoomTool(match_aspect=True), ResetTool()]

        norm_col = S.get_norm_col_name()
        norm_dist = S.NORM_DISTANCE

        tooltips = [('Nuc1', '@nuc_id1:  @chrom1  @start1-@end1'),
                    ('Nuc2', '@nuc_id2:  @chrom2  @start2-@end2'),
                    ('Counts', '@counts'),
                    (f'Norm ({norm_dist})', f'@{norm_col}'), ]

        # Bokeh Bug?
        # Setting match_aspect=True, aspect_scale=1 does not work for figure.rect, so setting width and height manually
        figure_args = dict(
            plot_width=600, plot_height=670, x_axis_location='above', toolbar_location='below',
            x_range=x_range, y_range=y_range, tools=tools, x_axis_label='Nucleosome 1', y_axis_label='Nucleosome 2', )

        p_overlaid = figure(
            title=f'Nuc-Nuc Interaction Counts for Chrom {chrom} from pos {start_region} to {end_region}',
            **figure_args)
        set_fig_attr(p_overlaid)
        p_overlaid.height -= 30

        p_overlaid_norm = figure(title="Normalized interaction counts", **figure_args)
        set_fig_attr(p_overlaid_norm)
        p_overlaid_norm.height -= 30

        lyt = layout([[p_overlaid, p_overlaid_norm]])

        tabs = list()
        tabs.append(Panel(child=lyt, title='Overlaid', closable=False))
        for orient, ijv in full_ijv.groupby('orient'):  # matrices.items():
            orient_label = f"{orient} ({','.join(Files.get_strands(orient))})".capitalize()
            title = f"Nuc-Nuc Interactions for range: {chrom},{start_region}-{end_region} and orient: {orient_label}"

            min_arcsinh_counts, max_arcsinh_counts = ijv['arcsinh_counts'].describe()[['min', 'max']]

            mapper = LinearColorMapper(  # LogColorMapper
                # palette=palettes[orient], low=full_ijv['counts'].min(), high=full_ijv['counts'].max())
                palette=palettes[orient], low=min_arcsinh_counts, high=max_arcsinh_counts)

            p = figure(title=title, **figure_args)

            r = p.rect(source=ijv, x='nuc_id1', y='nuc_id2', width=1, height=1, alpha=1,  # alphas[orient],
                       fill_color={'field': 'arcsinh_counts', 'transform': mapper}, line_color=None)

            p.add_tools(HoverTool(renderers=[r], tooltips=[(orient.title(), '')] + tooltips))
            p.add_tools(SaveTool())

            color_bar_mapper = LinearColorMapper(  # adjust low and high linearly
                palette=palettes[orient], low=np.sinh(min_arcsinh_counts), high=np.sinh(max_arcsinh_counts))
            color_bar = ColorBar(color_mapper=color_bar_mapper, label_standoff=10, height=20)
            p.height += color_bar.label_standoff + color_bar.height
            p.add_layout(color_bar, 'below')

            set_fig_attr(p)

            mapper_norm = LinearColorMapper(palette=palettes[orient], low=ijv[norm_col].min(), high=ijv[norm_col].max())

            p_norm = figure(title="Normalized interaction counts", **figure_args)

            r = p_norm.rect(source=ijv, x='nuc_id1', y='nuc_id2', width=1, height=1, alpha=1,  # alphas[orient],
                            fill_color={'field': norm_col, 'transform': mapper_norm}, line_color=None)

            p_norm.add_tools(HoverTool(renderers=[r], tooltips=[(orient.title(), '')] + tooltips))
            p_norm.add_tools(SaveTool())

            color_bar = ColorBar(color_mapper=mapper_norm, label_standoff=10, height=20)
            p_norm.height += color_bar.label_standoff + color_bar.height
            p_norm.add_layout(color_bar, 'below')

            set_fig_attr(p_norm)

            lyt = layout([[p, p_norm]])  # , sizing_mode='fixed')
            tabs.append(Panel(child=lyt, title=orient_label, closable=False))

            # Overlaid tab
            legend_label = '' if orient == 'all' else ' (' + ','.join(Files.get_strands(orient)) + ')'
            legend_label = orient.title() + legend_label
            r = p_overlaid.rect(source=ijv, x='nuc_id1', y='nuc_id2', width=1, height=1, alpha=.9,
                                fill_color={'field': 'arcsinh_counts', 'transform': mapper},
                                line_color=None, legend_label=legend_label)

            p_overlaid.add_tools(HoverTool(renderers=[r], tooltips=[(orient.title(), '')] + tooltips))
            if orient != 'all':
                r.visible = False

            r = p_overlaid_norm.rect(source=ijv, x='nuc_id1', y='nuc_id2', width=1, height=1, alpha=.9,
                                     fill_color={'field': 'arcsinh_counts', 'transform': mapper_norm},
                                     line_color=None, legend_label=legend_label)

            p_overlaid_norm.add_tools(HoverTool(renderers=[r], tooltips=[(orient.title(), '')] + tooltips))
            if orient != 'all':
                r.visible = False

        p_overlaid.add_tools(SaveTool())
        p_overlaid.legend.location = 'bottom_left'
        p_overlaid.legend.click_policy = 'hide'
        p_overlaid_norm.add_tools(SaveTool())
        p_overlaid_norm.legend.location = 'bottom_left'
        p_overlaid_norm.legend.click_policy = 'hide'

        last_tab = len(tabs) - 1
        tabs = Tabs(tabs=tabs)
        tabs.active = last_tab

        if save_only:
            save(tabs)  # only save the plot
            LOGGER.info(f'\nSince --save flag is used, the results plot is saved only and not displayed.')
        else:
            show(tabs)  # save and show the plot
            LOGGER.info(f'\nIf the results plot is not automatically displayed, '
                        f'please open the following plot file manually in a web browser.')

        LOGGER.info(f'\nOutput plot: \n'
                    f'   {output_file} \n'
                    f'or with full path:\n'
                    f'   {output_file.resolve()}\n')

    __help = {  # HELP texts used in different parts of the program
        'main': {
            'description': 'Creates and plots selected sections of nuc-nuc interaction matrices.',
            'title': 'possible commands',
            'help': 'Use one of the commands. Use -h option after any command for more information.',
            'args': {
                '--quiet': """
                    Suppress the progress output. Still, a log file will be saved in the working folder.""",
            }
        },
        'pairtools': {
            'help': """
                This command allows to create files in "pair" format from sam files or bam files. 
                All your commandline arguments will passed to pairtools for processing. 
                Please see instructions on pairtools website https://pairtools.readthedocs.io """,
            'epilog': 'NOTE: This command requires pairtools.',
            'args': {
                'pairtools_arguments': """
                All arguments presented here will be passed to pairtools for processing. To see the available option,
                type 'pairtools --help' on the commandline.""",
            },
        },
        'prepare': {
            'help': """
                Creates nuc-nuc interaction matrices, given a chromosomes file,
                a nucleosomes file, and an interactions file.""",
            'epilog': f"""
                All the following apply for all input files: 
                (1) Lines starting with "{S.COMMENT_CHAR}" will be ignored as comments.
                (2) The first non-comment line is interpreted as the header, containing 
                column names, and must be exactly as specified for each input type. If the first column is 
                indicated as [id], it means the column is optional and it can have any name 
                (e.g., "readID" or "nuc_id").
                (3) The optional [id] column may or may not be present as the first column. The values in the 
                optional [id] column will be ignored unless the are unique integers.
                (4) Values of columns on each line must be separated by a 
                {S.get_char_name(S.FIELD_SEPARATOR)} character.
                (5) Any columns after the mandatory columns will be ignored.
            """,
            'args': {
                'chroms_file': """
                    Chromosomes file with one mandatory column "chrom". Any other chromosome names that may appear
                    in the <nucs_file> or <interas_file> will be ignored if they not present in <chroms_file>.""",
                'nucs_file': """
                    Nucleosomes file with columns: "[id] chrom start end" with the [id] column being optional
                    (see below).""",
                'interas_file': """
                    Interactions file with columns: "[id] chrom1 pos1 chrom2 pos2 strand1 strand2" with the [id] column
                    being optional (see below). NOTE: The output of the pairtools program is accepted as <interas_file>
                    without requiring any changes (e.g., the extra columns will be ignored).""",
                '--dir': """
                    Optional output folder to store the nuc-nuc interaction matrices (and other files)
                    into. Without this option, <working_dir> is derived from the name of the
                    interaction input file, <interas_file>.""",
                '--norm': f"""
                    The normalized value for each nucleosome pair (NP) is calculated as the observed interaction
                    count for that NP divided by the expected count.
                    The expected count for an NP is the average of interaction count between nucleosome pairs whose
                    distance is within the NP's average distance +/- a given <distance>.
                    The default value for the parameter <distance> is {S.NORM_DISTANCE}.""",
                '--no-nucs': """
                    Avoid saving extra columns for nucleosomes (chrom1 start1 end1 chrom2 start2 end2) in the generated
                    nuc-nuc interaction matrices. Instead, it will count on nuc_id1 and nuc_id2 (which are generated
                    automatically) to identify nucleosomes. This option reduces the size of the produced matrix files,
                    which can improve efficiency in memory and storage usage and speed.""",
                '--no-cache': """
                    While computing the nuc-nuc matrices, using the prepare command, there are many intermediate files
                    that are stored in the <working_dir>/cache subdirectory and can optionally be removed at the end.
                    Keeping these intermediate files can speed up the subsequent executions of the prepare command
                    (for example with different --norm parameters). The cache files might also prove useful for some
                    downstream analysis by the users. However, keeping cache files increases the amount of disk space
                    requirements significantly.
                    The option --no-cache instructs the prepare command to remove the cache files as soon as the
                    matrix files are generated. The cache files are not needed for the subsequent executions of the
                    plot command and can be deleted manually if desired.""",
                '--multiprocessing': """
                    Maximum number of processes (one per CPU core) that simultaneously analyse the data. Typically, this
                    depends on the number of CPU cores and the amount of memory available on the system. For very large 
                    input files, the recommended number of processes is one per 5 GB of RAM available (e.g., use -m 12 
                    on a system with 60 GB of free RAM). The default value is 0 (zero), which means the program attempts
                    to guess a reasonable number of processes.""",
                '--zip': """
                    Compress intermediary and cache files using gzip. If these files already exist, you need to 
                    use the --refresh flag to regenerate them in zipped format. Zipping files saves space but requires 
                    more time to write and read files. """
                         'BUG: zipping does not currently work under macOS!' if platform.system() == 'Darwin' else '',
                '--refresh': """
                    This will start fresh, re-computing all the intermediate files such as the nuc-nuc matrices.
                    WARNING: Refreshing starts from scratch and may take a long time to complete,
                    depending on the size of the input data and your hardware specifications.""",
            }
        },
        'plot': {
            'help': 'Select and plot a subset (or window) of a nuc-nuc interaction matrix.',
            'epilog': None,
            'args': {
                'working_dir': 'This is the folder created by the "prepare" command.',
                'chrom': """
                    A chromosome, e.g. III or chr21, on which a region is being selected.
                    NOTE: chrom has to be from those mentioned in <chroms_file> provided
                    to the "prepare" command.""",
                'start_region': 'A positive integer, representing the start of the region being selected.',
                'end_region': 'A positive integer, representing the end of the region being selected.',
                '--prefix': """
                    A prefix used for output file names generated for the selected subset of the nuc-nuc
                    interaction matrix. The generated files using this prefix include submatrices for the
                    selected regions and orientations and an html file containing the plots.""",
                '--orientation': """
                    Optional parameter to specify the orientation of interaction. 
                    If this option is not present, all orientations and their combinations will be considered.""",
                '--save': """
                    Only save the output files and do not attempt to show them. Without this flag, the results are
                    both saved and shown. This flag can be useful in scripts, for example.""",
            }
        },
    }

    @classmethod
    def create_commandline_parser(cls):
        h = cls.__help

        parser_main = argparse.ArgumentParser(description=h['main']['description'])
        subparsers = parser_main.add_subparsers(title=h['main']['title'], help=h['main']['help'], dest='command')

        # # # Create the UI commands
        commands = dict()
        for c in ['pairtools', 'prepare', 'plot']:
            commands[c] = subparsers.add_parser(c, help=h[c]['help'], description=h[c]['help'], epilog=h[c]['epilog'])

        # # # arguments for main
        c = 'main'
        a = '--quiet'  # and '-q'
        parser_main.add_argument(a[1:3], a, help=h[c]['args'][a], default=False, action='store_true')

        # # # arguments for commands
        c = 'pairtools'  # c for Command
        a = 'pairtools_arguments'  # a for Argument
        commands[c].add_argument(a, help=h[c]['args'][a], nargs=argparse.REMAINDER)

        c = 'prepare'  # c for Command
        a = 'chroms_file'  # a for Argument
        commands[c].add_argument(a, help=h[c]['args'][a], metavar=f'<{a}>')
        a = 'nucs_file'
        commands[c].add_argument(a, help=h[c]['args'][a], metavar=f'<{a}>')
        a = 'interas_file'
        commands[c].add_argument(a, help=h[c]['args'][a], metavar=f'<{a}>')
        a = '--dir'  # and '-d'
        commands[c].add_argument(a[1:3], a, help=h[c]['args'][a], metavar='<working_dir>', dest='working_dir')
        a = '--norm'  # and '-n'
        commands[c].add_argument(a[1:3], a, help=h[c]['args'][a], metavar='<distance>', dest='norm_distance',
                                 type=int, default=S.NORM_DISTANCE)
        a = '--no-nucs'  # and '-N'
        commands[c].add_argument(
            '-N', a, help=h[c]['args'][a], default=True, action='store_false', dest='keep_nucs_cols')
        a = '--no-cache'  # and '-C'
        commands[c].add_argument('-C', a, help=h[c]['args'][a], default=True, action='store_false', dest='keep_cache')
        a = '--multiprocessing'  # and '-m'
        commands[c].add_argument(
            a[1:3], a, help=h[c]['args'][a], metavar='<processes>', dest='n_processes', type=int, default=0)
        a = '--refresh'  # don't offer -r option to avoid accidental refreshing
        commands[c].add_argument(a, help=h[c]['args'][a], default=False, action='store_true')
        if platform.system() != 'Darwin':  # todo BUG: zipping does not currently work under macOS!
            a = '--zip'  # and '-z'
            commands[c].add_argument(a[1:3], a, help=h[c]['args'][a], default=False, action='store_true', dest='zipped')

        c = 'plot'  # c for Command
        a = 'working_dir'  # a for Argument
        commands[c].add_argument(a, help=h[c]['args'][a], metavar=f'<{a}>')
        a = 'chrom'
        commands[c].add_argument(a, help=h[c]['args'][a], metavar=f'<{a}>')
        a = 'start_region'
        commands[c].add_argument(a, help=h[c]['args'][a], metavar=f'<{a}>', type=int)
        a = 'end_region'
        commands[c].add_argument(a, help=h[c]['args'][a], metavar=f'<{a}>', type=int)
        a = '--save'  # and '-s'
        commands[c].add_argument(a[1:3], a, help=h[c]['args'][a], default=False, action='store_true', dest='save_only')
        a = '--prefix'  # and '-p'
        commands[c].add_argument(
            a[1:3], a, help=h[c]['args'][a], metavar='<outfile_prefix>', default='plot', dest='prefix')

        return parser_main

    @staticmethod
    def test_parse_args(parser):
        assert S.TESTING_MODE
        LOGGER.debug('WARNING: Using TESTING mode!')
        # # # for development only:
        refresh = ""  # "--refresh"
        chroms_file = r"data/yeast_chromosomes.txt"
        nucs_file = r"data/yeast_nucleosomes.txt"
        interas_file = r"data/yeast_interactions.txt"
        working_dir = r"wd/yeast_debug"
        # chroms_file = r"data/human_chromosomes.txt"
        # nucs_file = r"data/human_nucleosomes.txt"
        # interas_file = r"data_local/4DNFI1O6IL1Q.pairs.gz"
        # working_dir = r"wd/human_debug"

        # testing_args = ""
        # testing_args = "--help"
        # testing_args = "prepare --help"
        # testing_args = "plot --help"
        testing_args = f"prepare {refresh} {chroms_file} {nucs_file} {interas_file} --dir {working_dir}"
        # testing_args = f"plot {working_dir} chrIII 237000 262000 --prefix my_plot"
        LOGGER.debug(f"inucs {testing_args}")
        arguments = parser.parse_args(testing_args.split())  # for debugging

        return arguments


def main():
    # logging.basicConfig(filename=Path(__file__).stem + '.log', level=logging.DEBUG, filemode='w')
    # LOGGER = logging.getLogger()
    # formatter = logging.Formatter('%(asctime)p %(levelname)p %(message)p')

    commandline_parser = CLI.create_commandline_parser()

    # S.TESTING_MODE = True  # for debugging
    if S.TESTING_MODE:
        args = CLI.test_parse_args(commandline_parser)  # for testing only
    else:
        args = commandline_parser.parse_args()  # the actual line

    if not args.quiet:
        handler_stdout = logging.StreamHandler(sys.stdout)
        # handler_stdout.setFormatter(formatter)
        LOGGER.addHandler(handler_stdout)
        LOGGER.info(f'\nUse flag --quiet to suppress printing to standard output.')

    del args.quiet  # argument not needed anymore

    # using args.command, select a command handler function
    command_handler = {
        None: lambda **_: commandline_parser.print_help(),  # prints help, ignoring any arguments passed to it
        'pairtools': CLI.handler_pairtools_command,
        'prepare': CLI.handler_prepare_command,
        'plot': CLI.handler_plot_command,
    }.get(args.command)

    command_handler(**vars(args))  # call the selected handler with all the arguments from commandline


if __name__ == "__main__":
    logging.getLogger('numexpr').setLevel(logging.WARNING)  # to suppress an irrelevant message from NumExpr
    LOGGER = logging.getLogger()
    LOGGER.setLevel(level=logging.DEBUG)  # use logging.INFO or logging.DEBUG
    FILES = Files()  # each command handler will reset FILES with appropriate parameters

    main()

# done Min Python version was reduced to 3.8
# done input checking, e.g.: validate_chroms_file, validate_nucs_file, validate_interas_file
# done check for zipped Nucs file
# done reverse y axis
# done input/output error checking
# done add a --no-zip flag!
# done doc: bulleted list. Efficiency; scalability (breaking done: input, matrix, submatrix)
# done doc: graphs:  1) usage pipe line  2) how efficiency is achieved
# done github
# done add Tandem
# done add to Legend '+-', '-+', '++', '--'
# done Colors, Font size
# done Bokeh, default tool? on scroll zoom
# todo Bokeh, Full GUI app?
# todo Bokeh, multiple HoverTools (tooltips)
# todo Bokeh, keep aspect ratio the same as zooming in and out
# todo Bokeh, resize plot?
# todo Efficiency: use parallelism to utilize multicore CPUs
# todo Efficiency, Mem: use categorical data type in DataFrames to reduce the size
# todo Efficiency, Mem: use memory profiling to see where is highest memory usage
# todo Efficiency, Mem: check lru_cache memory usage in find_nucs_in_region, for example.
# TODO doc: renew running time measurements (Yeast, Human, PC, HPC)
# TODO deployment methods? https://conda.io/projects/conda-build/en/latest/user-guide/tutorials/index.html
