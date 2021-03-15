#!/usr/bin/env python
# coding: utf-8

import argparse
import collections.abc as abc
import configparser
import contextlib
import copy
import functools
import gzip
import logging
import mimetypes
import pprint
import shutil
import subprocess
import sys
from pathlib import Path
from timeit import default_timer as timer
from typing import Optional

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError


class S:  # S for Settings
    """ Holds all the Settings and Constants used globally """
    # todo the following constants should be read from a config file
    OUTPUT_COLS__MATRIX = ['nuc_id1', 'nuc_id2', 'counts']
    OUTPUT_COLS__NUC_INTERS = ['chrom1', 'pos1', 'chrom2', 'pos2', 'nuc_id1', 'nuc_id2']
    COMMENT_CHAR = '#'
    FIELD_SEPARATOR = '\t'
    INPUT_CHUNK_LINES_READ = 1000 * 1000
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


class InputFileError(SyntaxError):
    def __init__(self, filename):
        message = '\n' \
                  f'Please make sure that the {S.get_char_name(S.FIELD_SEPARATOR)} character is used as ' \
                  f'the field-separator used in file\n{filename}\n' \
                  f'Another possibility is that there are the wrong number of columns in this input file.'
        super().__init__(message)
        self.filename = filename


class Files:
    S_INTER = 'state_inter'
    S_NUC_INTER = 'state_nuc_inter'
    S_MATRIX = 'state_matrix'
    S_PLOT = 'state_plot'

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
    __OUT_SUFFIX = '.inucs'

    def __init__(self, base_file_name=None, working_dir=None, refresh: bool = False, zipped: bool = True):
        self.__working_files: pd.DataFrame = pd.DataFrame()  # will be reset by reset_working_files
        self.__base_file_name = base_file_name
        self.__working_dir = working_dir
        self.__zipped = zipped
        self.__config = configparser.ConfigParser()
        self.reset_all(base_file_name, working_dir, refresh, zipped)

    def reset_all(self, base_file_name=None, working_dir=None, refresh: bool = False, zipped: bool = True):

        self.__zipped = zipped

        suffix = self.__OUT_SUFFIX
        if base_file_name and working_dir:
            self.__base_file_name = Path(base_file_name)
            self.__working_dir = Path(working_dir)
        elif base_file_name:  # then working_dir is None, so determine it
            self.__base_file_name = Path(base_file_name)
            self.__working_dir = Path('.') / (self.__base_file_name.name + suffix)
        elif working_dir:  # then base_file_name is None, so determine it
            self.__working_dir = Path(working_dir)
            interas_file = self.input_interas_file_name
            if interas_file is None:
                wd = self.__working_dir
                self.__base_file_name = wd.stem if wd.suffix == suffix else wd
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
        if FILES.working_dir is None or not FILES.working_dir.exists():
            raise RuntimeError(
                f'Working directory does not exist '
                f'{FILES.working_dir if FILES.working_dir else ""}')
        if not self.working_dir.is_dir():
            raise FileExistsError(f'Expecting a folder, but a file exists: {self.working_dir}')

        nucs_file = self.input_nucs_file_name
        if nucs_file is None or not Path(FILES.working_dir / nucs_file).exists():
            raise RuntimeError(
                f'Cannot find the indexed nucleosome file in the working directory '
                f'{nucs_file if nucs_file else ""}\n'
                'Please consider using the --refresh flag.')

        matrix_subdir = self.__STATES.loc[Files.S_MATRIX, 'subdir']
        matrix_subdir = FILES.working_dir / matrix_subdir
        if not matrix_subdir.exists() or not any(matrix_subdir.iterdir()):
            raise RuntimeWarning(
                f'WARNING: No matrices folder or files found {matrix_subdir}\n'
                'Please consider using the --refresh flag.'
            )

    def iter_working_files(self):
        df = self.__working_files.copy()
        df = df[self.__IMP_COLS.drop(['needs_refresh', 'subdir'])]  # dropping cols that depend on state and unspecified
        for _, file_info in df.iterrows():
            yield file_info

    @property
    def working_dir(self):
        return self.__working_dir

    @property
    def base_file_name(self):
        return self.__base_file_name

    @property
    def input_nucs_file_name(self):
        file = None
        if self.__reload_config_file():
            file = self.__config['input_files']['nucs_file']
        return Path(file) if file else None

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

    @property
    def zipped(self) -> bool:
        """Is output zipped?"""
        return self.__zipped

    @classmethod
    def orientations(cls) -> list:
        return list(cls.__STRAND2ORIENT.values())  # 'inner' 'outer' 'right' 'left'

    @classmethod
    def combined_orientations(cls) -> list:
        return list(cls.__ORIENT2STRAND.keys())  # 'inner' 'outer' 'right' 'left' 'tandem' 'all'

    @classmethod
    def get_strands(cls, orientation) -> list:
        return cls.__ORIENT2STRAND.get(orientation)

    def fix_file_zip_suffix(self, file_name):

        file_name = Path(file_name)
        suffix = file_name.suffix
        if not self.zipped and suffix == '.gz':
            suffix = ''
        elif self.zipped and suffix != '.gz':
            suffix += '.gz'

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
                    subdirs_to_remove = Files.__STATES.query("data_type in ['cached', 'intermediary']")
                    for _, subdir_info in subdirs_to_remove.iterrows():
                        subdir = (self.working_dir / subdir_info.subdir)
                        if subdir.exists():
                            shutil.rmtree(subdir)  # delete subdir
                    if self.input_nucs_file_name:
                        nucs_file = self.working_dir / self.input_nucs_file_name
                        nucs_file.unlink(missing_ok=True)  # delete file
                    self.reset_input_file_names()  # resets both arguments to None
                except Exception as e:
                    message = 'Cannot remove some of the content in the working directory.\n' \
                              f'Please remove "{self.working_dir}" manually and try again.'
                    # LOGGER.error(message)
                    raise RuntimeError(message) from e

            # recreate sub-folders
            for _, subdir_info in Files.__STATES.iterrows():
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
        # df['file'] = self.base_file_name.stem + '.' + df.chrom + '.' + df.orientation + self.base_file_name.suffix
        # df['file_zip'] = df.file + df.file.map(lambda f: '' if Path(f).suffix == '.gz' else '.gz')
        suffix = '' if self.base_file_name.suffix == '.gz' else self.base_file_name.suffix  # remove .gz
        df['file'] = self.base_file_name.stem + '.' + df.chrom + '.' + df.orientation + suffix
        df['file_zip'] = df.file + '.gz'  # todo remove: after adding --zip flag, file_zip is not needed anymore
        df[Files.__STATES.dropna()['refresh_col'].to_list()] = True  # makes three bool columns for file refreshing

        self.__working_files = df

    def __decide_what_needs_refresh(self, refresh):
        status = Files.__STATES[::-1].dropna()  # __STATES with rows reversed and the 'plot' row dropped
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


class Chroms(abc.Sequence):
    """ Chromosomes """

    def __init__(self, chrom_list_file, comment: str = S.COMMENT_CHAR, name: str = None):
        chroms = pd.read_csv(chrom_list_file, names=['chrom'], comment=comment,
                             sep=S.FIELD_SEPARATOR, header=None, usecols=[0], squeeze=True)

        if len(chroms) > 0 and chroms[0] in ['chrom', 'chr']:  # remove the first row if it was a header row
            chroms = chroms.iloc[1:, :].reset_index(drop=True)

        self.__chroms = chroms
        self.__name = name if name else str(Path(chrom_list_file).name)
        LOGGER.info(f'\nChroms file read: "{self.name}" containing {len(self)} chromosomes:')
        LOGGER.info(', '.join(self.list) + '\n')

    @property
    def name(self) -> str: return self.__name

    @name.setter
    def name(self, name: str): self.__name = name

    @property
    def series(self) -> pd.Series: return self.__chroms.copy()

    @property
    def list(self) -> list: return self.__chroms.to_list()

    def __len__(self): return self.__chroms.__len__()

    def __repr__(self): return self.__chroms.__repr__()

    def __iter__(self): return self.__chroms.__iter__()

    def __getitem__(self, index): return self.__chroms.__getitem__(index)

    def __contains__(self, item): return item in self.__chroms or self.__chroms.isin([item]).any()


class Nucs:
    """ Nucleosomes """

    def __init__(self, nucs_file=None, chroms_file=None, name: str = None,
                 sep: chr = S.FIELD_SEPARATOR, comment: str = S.COMMENT_CHAR, load_nucs_file=True):
        self.__chrom_2_id_chrom_start_end: dict = {}
        self.__id_chrom_start_end: pd.DataFrame = self.schema()

        if not load_nucs_file:
            return

        read_index = False
        if nucs_file is None:
            read_index = True  # nucs_file can be None only when it's previously been created in working dir with index
            FILES.validate_working_dir()
            nucs_file = FILES.working_dir / FILES.input_nucs_file_name

        try:
            df_input = self.read_csv(nucs_file, sep=sep, comment=comment, read_index=read_index)
        except pd.errors.ParserError as parser_error:
            raise InputFileError(nucs_file) from parser_error

        self.__name: str = name if name else str(Path(nucs_file).name)
        chrom_list = Chroms(chroms_file, comment=comment).list if chroms_file else None
        self.__update(df_input, chrom_list, self.name)

        if not read_index:
            nucs_file = self.to_csv(output_file=nucs_file)
            FILES.set_input_nucs_file(nucs_file=nucs_file)

        LOGGER.info(f'\nNucs file read: "{self.name}" containing {len(self)} nucleosomes.')
        LOGGER.info(self)

    @classmethod
    def from_dataframe(cls, id_chrom_start_end: pd.DataFrame, name: str = None):
        return cls(nucs_file=None, load_nucs_file=False).__update(id_chrom_start_end, name=name)

    def __update(self, id_chrom_start_end: pd.DataFrame, chrom_list: list = None, name: str = None):
        self.__name = name
        if 'nuc_id' in id_chrom_start_end.columns:
            id_chrom_start_end = id_chrom_start_end.set_index('nuc_id')
        elif 'nuc_id' != id_chrom_start_end.index.name:
            id_chrom_start_end = id_chrom_start_end.reset_index(drop=True)
            id_chrom_start_end = id_chrom_start_end.rename_axis('nuc_id')

        # ensuring cols existence and order, and ignoring extra cols if any
        self.__id_chrom_start_end = id_chrom_start_end[['chrom', 'start', 'end']]

        if chrom_list:  # keep only rows with chrom in chrom_list
            self.__id_chrom_start_end = self.__id_chrom_start_end.query(f"chrom in {chrom_list}")

        # internal cache to reduce the lookup time for each chrom
        id_chrom_start_end = self.__id_chrom_start_end.reset_index()
        self.__chrom_2_id_chrom_start_end = \
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
        return list(self.__chrom_2_id_chrom_start_end.keys())

    @property
    def id_chrom_start_end(self):
        return self.__id_chrom_start_end.copy()

    def get_nucs(self, chrom) -> pd.DataFrame:
        return self.__chrom_2_id_chrom_start_end[chrom]

    def iter_nucs(self):
        return self.__id_chrom_start_end.iterrows()

    def iter_chroms(self):
        return iter(self.__chrom_2_id_chrom_start_end.keys())

    def iter_chrom_and_nucs(self):
        return iter(self.__chrom_2_id_chrom_start_end.items())

    def __len__(self):
        return len(self.__id_chrom_start_end)

    def __repr__(self):
        return self.name + '\n' + \
               str(self.__id_chrom_start_end.head(2)) + '\n...\n(' + \
               str(self.__id_chrom_start_end.shape[0]) + ' nucleosomes)'

    @functools.lru_cache
    def find_nuc_id(self, chrom, pos):
        df = self.__chrom_2_id_chrom_start_end[chrom]
        i = df.start.searchsorted(pos)
        if i == len(df) or (i > 0 and pos != df.at[df.index[i], 'start']):
            i -= 1
        nuc = df.iloc[i, :]
        found = nuc.start <= pos <= nuc.end
        return nuc.name, found  # nuc.name is nuc_id

    @functools.lru_cache
    def find_nuc(self, chrom, pos):
        id_, found = self.find_nuc_id(chrom, pos)
        return self.__id_chrom_start_end.iloc[id_, :] if found else None

    @functools.lru_cache
    def find_nucs_in_region(self, chrom, start_region, end_region) -> Optional['Nucs']:
        if chrom not in self.__chrom_2_id_chrom_start_end:
            return None

        if start_region > end_region:  # swap if in the wrong order
            start_region, end_region = end_region, start_region

        start_id, _ = self.find_nuc_id(chrom, start_region)
        end_id, _ = self.find_nuc_id(chrom, end_region)
        nucs_in_region = self.__id_chrom_start_end.iloc[start_id:end_id + 1, :]

        region_name = f'{self.name} ({chrom}:{start_region}->{end_region})'
        return Nucs.from_dataframe(nucs_in_region, region_name)

    @staticmethod
    def read_csv(nucs_file, sep=S.FIELD_SEPARATOR, comment: str = S.COMMENT_CHAR, read_index=False):
        try:
            df = pd.read_csv(nucs_file, sep=sep, comment=comment)  # read nucs info: chrom, start, end
        except pd.errors.ParserError as parser_error:
            raise InputFileError(nucs_file) from parser_error

        df = df.rename({'chr': 'chrom'}, axis=1)
        if read_index:
            if 'nuc_id' in df:
                df = df.set_index('nuc_id')  # this is just resetting the same nuc_id read from file
        else:
            df = df.sort_values(['chrom', 'start', 'end'])  # sorts and ensures cols exist
            df = df.reset_index(drop=True).rename_axis('nuc_id')
            df = df[['chrom', 'start', 'end']]  # reordering cols and ignoring extra cols if any

        return df

    def to_csv(self, output_file=None, sep=S.FIELD_SEPARATOR, index: bool = True):
        output_file = FILES.fix_file_zip_suffix(output_file if output_file else self.name)
        self.__id_chrom_start_end.to_csv(
            FILES.fullpath(output_file), sep=sep, index=index, header=True)
        return output_file


class NucInteras:
    """ Nucleosome Interactions """

    def __init__(self, interas_file=None, nucs=None, name: str = None, refresh: bool = False):

        should_create_nuc_interas_files = True
        if interas_file is None:
            should_create_nuc_interas_files = False
            interas_file = FILES.input_interas_file_name

        self.__name = name if name else str(Path(interas_file).name)
        self.__len = 0  # length or number of nuc interactions

        if nucs is None:
            nucs = Nucs()

        # Make a MultiIndex ['chrom', 'pos'] for Nucs to make them similar to Interactions, simplifying their merger
        df_nucs = nucs.id_chrom_start_end
        df_nucs['pos'] = df_nucs.start  # copy 'start' to 'pos' to make df_nucs and later df_inter alike
        df_nucs = df_nucs.reset_index().set_index(['chrom', 'pos'])  # .sort_index()
        self.__nucs__id_chrom_start_end = df_nucs

        self.__chrom_list = nucs.chrom_list

        FILES.reset_working_files(chrom_list=self.__chrom_list, refresh=refresh)

        FILES.set_input_interas_file(interas_file=interas_file)

        if should_create_nuc_interas_files:
            self.__interas_file = Path(interas_file).resolve(strict=True)

            LOGGER.info('\nSplitting interaction file based on chrom and orientation:')
            start = timer()
            self.__split_interas_file()
            LOGGER.info(f'Done! Splitting task finished in {round((timer() - start) / 60.0, 1)} minutes')

            LOGGER.info('\nFinding nucleosomes for interactions:')
            start = timer()
            self.__create_nuc_interas_files()
            LOGGER.info(f'Done! Finding nucleosomes finished in {round((timer() - start) / 60.0, 1)} minutes')

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

    def iter_nuc_interas(self):
        for file_info in FILES.iter_working_files():
            nuc_interas = self.read_nuc_inter(file_info.chrom, file_info.orientation)
            yield nuc_interas  # todo test to see it's ok where nuc_interas.empty is True

    def read_nuc_inter(self, chrom: str, orientation: str) -> pd.DataFrame:
        """
        :param chrom: e.g. 'II' or 'chrom14'
        :param orientation:  It can be any of ('+-', '-+', '++', '--') or ('Inner', 'Outer', 'Right', 'Left')
        Further, 'p' instead of '+', and 'n' or 'm' instead of '-' are acceptable.
        :return:
        """  # todo Default orientation should be None for all.

        working_file = FILES.get_working_files(chrom, orientation, Files.S_NUC_INTER).squeeze()
        if working_file.empty:
            return self.__EMPTY

        input_file = working_file.subdir / working_file.file
        if not input_file.exists():
            input_file = working_file.subdir / working_file.file_zip  # check if zip version exists
            if not input_file.exists():
                return self.__EMPTY

        try:
            df_nuc_interas = pd.read_csv(input_file, sep=S.FIELD_SEPARATOR, comment=S.COMMENT_CHAR,
                                         usecols=S.OUTPUT_COLS__NUC_INTERS[-2:])  # ['nuc_id1', 'nuc_id2']
        except pd.errors.ParserError as parser_error:
            raise InputFileError(input_file) from parser_error

        df_nuc_interas['counts'] = 1  # initially every row is counted as 1

        return df_nuc_interas  # so now: ['nuc_id1', 'nuc_id2', 'counts']

    def __create_nuc_interas_files(self):

        subdir_inter = FILES.get_subdir(Files.S_INTER)
        subdir_nuc_inter = FILES.get_subdir(Files.S_NUC_INTER)

        # decide if need to continue
        files_needing_refresh = FILES.get_files_needing_refresh(Files.S_NUC_INTER)
        if files_needing_refresh.empty:
            LOGGER.info('\nNo nucleosome interaction files were updated')
            return

        for _, file_info in files_needing_refresh.iterrows():
            # if not file_info.needs_refresh:  # get_files_needing_refresh has already checked for this
            #     continue
            if file_info.chrom not in self.__chrom_list:
                continue
            input_file = subdir_inter / file_info.file
            if not input_file.exists():
                input_file = subdir_inter / file_info.file_zip  # check if the zip version exists
                if not input_file.exists():
                    LOGGER.warning(f'\nSkipping cache file. File not found: {file_info.file}.')
                    continue
            output_file = file_info.file_zip if FILES.zipped else file_info.file
            output_file = subdir_nuc_inter / output_file
            LOGGER.info(f'Creating nuc interaction file: {output_file}')
            df_interas = self.__read_interas_file(input_file)
            df_nuc_interas = self.__create_nuc_interas(df_interas)
            df_nuc_interas.to_csv(output_file, sep=S.FIELD_SEPARATOR, index=False, header=True)
            try:
                Path(input_file).unlink()  # remove the cached input_file
            except Exception as e:
                LOGGER.error(e)

        FILES.set_needs_refresh_for_all_files(False, Files.S_NUC_INTER)
        # files.__working_files[refresh_col] = False  # flag all as done

    @classmethod
    def __read_interas_file(cls, input_file) -> pd.DataFrame:

        schema = cls.schema()

        try:
            df_inter = pd.read_csv(input_file, sep=S.FIELD_SEPARATOR, comment=S.COMMENT_CHAR,
                                   usecols=range(schema.columns.size), header=None)
        except pd.errors.ParserError as parser_error:
            raise InputFileError(input_file) from parser_error

        df_inter.columns = schema.columns
        df_inter = df_inter.rename_axis(schema.index.name)

        return df_inter

    def __create_nuc_interas(self, df_inter: pd.DataFrame) -> pd.DataFrame:

        # Step 1:
        # Both Nucs and Interactions will use the same MultiIndex ['chrom', 'pos']
        # See __init__() for preprocessing of self.__nucs__id_chrom_start_end
        df_nucs = self.__nucs__id_chrom_start_end  # already has .set_index(['chrom', 'pos'])

        df_inter_stack = df_inter.stack(level='side').reset_index()  # 'side' becomes a col, with 'side1' or 'side2'
        df_inter_stack = df_inter_stack.set_index(['chrom', 'pos'])  # .sort_index()

        # Step 2:
        # Combine Nucs and Interactions (which use the same MultiIndex ['chrom', 'pos'] now)
        index_names = ['orig'] + df_inter_stack.index.names  # 'orig' column will be 'nuc' or 'inter'
        df_nuc_interas = pd.concat((df_nucs, df_inter_stack), keys=['nuc', 'inter'], names=index_names)
        df_nuc_interas = df_nuc_interas.reset_index(level='orig').sort_index()  # sort on ['chrom', 'pos']
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
        df_nuc_interas = df_nuc_interas[S.OUTPUT_COLS__NUC_INTERS]
        # ['chrom1', 'pos1', 'chrom2', 'pos2', 'nuc_id1', 'nuc_id2']

        # todo efficiency: sorting is useful for output files but not needed. Should keep it?
        df_nuc_interas = df_nuc_interas.sort_values(['nuc_id1', 'nuc_id2']).reset_index(drop=True)

        return df_nuc_interas

    def __split_interas_file(self):
        # decide if need to continue
        if FILES.get_files_needing_refresh(Files.S_INTER).empty:
            LOGGER.info('\nNo interaction cache files were updated')
            return

        bash_available = self.__split_interas_file__bash()
        if not bash_available:
            LOGGER.warning(
                '\nNOTE: Could not find "bash" on this system, so falling back to an alternative '
                'approach, which is slower. Using a system with bash (e.g., Linux or macOS) can '
                'speed up the process of preparing internal cache files.')
            self.__split_interas_file__python()

    def __split_interas_file__bash(self) -> bool:
        """
        :returns: false if "bash" is not available and true if the bash script has successfully executed
        """

        try:  # first make sure that bash is available
            process = subprocess.Popen(  # ..., text=True)  # python 3.7
                '/bin/bash', stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        except FileNotFoundError:
            return False

        header_line_num, used_cols = self.__validate_interas_file(self.__interas_file)

        sep = S.FIELD_SEPARATOR

        # start with files that need to be (re)created
        df = FILES.get_files_needing_refresh(Files.S_INTER)
        if df.empty:
            return True

        # create subcommands for filtering input records:
        df['strands_match'] = df.strands.map(lambda strand_tuple: '\\' + (sep + '\\').join(strand_tuple))
        df['awk'] = "awk '/^" + df.chrom + sep + ".*" + sep + df.strands_match + "$/'"
        if FILES.zipped:
            df['filtering'] = '>(' + df.awk + ' | gzip > ' + df.file_zip + ') \\'
        else:
            # files_with_no_zip_suffix = df.file.map(FILES.fix_file_zip_suffix)
            # df['filtering'] = '>(' + df.awk + ' > ' + files_with_no_zip_suffix + ') \\'
            df['filtering'] = '>(' + df.awk + ' > ' + df.file + ') \\'
        filtering_subcommands = '\n'.join(df['filtering'].to_list()) + '\n'

        # subdir = FILES.get_subdir(Files.S_INTER)
        subdir = df.head(1).subdir.squeeze()  # read subdir from first row, which is the same here for all rows

        # create the complete bash command lines to split the input interaction file
        command_lines = f'cd "{subdir}" \n'
        command_lines += f'cat "{self.__interas_file}" \\\n'
        command_lines += ' | gunzip ' if mimetypes.guess_type(self.__interas_file)[1] == 'gzip' else ''
        command_lines += f" | grep -v '^{S.COMMENT_CHAR}'"
        command_lines += '' if header_line_num is None else f" | tail -n +{header_line_num + 2}"
        command_lines += f" | cut -f{','.join((1 + used_cols.index).map(str))}"
        command_lines += f" | awk '${'==$'.join((1 + S.CHROM_COLS.index).map(str))}'"
        command_lines += ' | tee \\\n' + filtering_subcommands
        command_lines += '> /dev/null'

        # COMMENT: In filtering_subcommands above, using "awk" instead of "grep" for filtering
        # made almost no difference in terms of performance in our benchmarking, but the regex
        # code for awk seems cleaner than grep (e.g., how '+' is treated)

        LOGGER.info(f'\nSplitting input file: {self.__interas_file} \n')
        LOGGER.debug(command_lines)

        out, err = process.communicate(command_lines)
        if err:
            raise RuntimeError(err)

        LOGGER.info(out)

        FILES.set_needs_refresh_for_all_files(False, Files.S_INTER)

        # delete any empty csv files which may have been created by the bash script above
        for file in subdir.iterdir():
            if file.is_file():
                try:  # an easy trick to detect zipped empty files (which have a size > 0 so are hard to detect)
                    pd.read_csv(file, nrows=1)  # try to read 1 line
                except EmptyDataError:
                    file.unlink()  # remove file if it has no data

        return True

    def __split_interas_file__python(self):

        header_line_num, used_cols = self.__validate_interas_file(self.__interas_file)

        try:
            chunks = pd.read_csv(self.__interas_file, sep=S.FIELD_SEPARATOR, comment=S.COMMENT_CHAR,
                                 usecols=used_cols.index, chunksize=S.INPUT_CHUNK_LINES_READ, header=header_line_num)
        except pd.errors.ParserError as parser_error:
            raise InputFileError(self.__interas_file) from parser_error

        groupby_cols = [S.CHROM_COLS.iloc[0]] + S.STRAND_COLS.to_list()  # ['chrom1', 'strand1', 'strand2']
        output_files = {}
        with contextlib.ExitStack() as cm:  # allows to open an unknown number of files
            for df_chunk in chunks:
                df_chunk.columns = S.USED_COLS
                df_chunk = df_chunk.query('=='.join(S.CHROM_COLS))  # rows where chrom1==chrom2
                df_chunk = df_chunk.query(  # todo efficiency: this approach is slow if all chrom_list are acceptable
                    f'{S.CHROM_COLS.iloc[0]} in {self.__chrom_list}')  # rows with listed chrom
                groups = df_chunk.groupby(groupby_cols)
                for chrom_and_orientation, df in groups:
                    chrom = chrom_and_orientation[0]  # e.g., 'chrom21'
                    orientation = chrom_and_orientation[1:]  # e.g., ('+', '-')
                    working_file = FILES.get_working_files(chrom, orientation, Files.S_INTER).squeeze()
                    if working_file.empty or not working_file.needs_refresh:
                        continue
                    file = output_files.get(chrom_and_orientation)
                    if file is None:
                        open_file = gzip.open if FILES.zipped else open
                        file = working_file.file_zip if FILES.zipped else working_file.file
                        file = working_file.subdir / file
                        file = cm.enter_context(open_file(file, 'wt'))
                        output_files[chrom_and_orientation] = file

                    # Here header=False because the code for bash is not writing out header (yet)
                    df.to_csv(file, sep=S.FIELD_SEPARATOR, index=False, header=False)

        FILES.set_needs_refresh_for_all_files(False, Files.S_INTER)

    @classmethod
    def __validate_interas_file(cls, interas_file, num_lines_to_check=1000):
        # check the header names if there are any in the interas_file
        df = pd.read_csv(interas_file, sep='\n', nrows=num_lines_to_check, names=['line'])
        comment_lines = df[df['line'].str.startswith(S.COMMENT_CHAR)]
        data_is_in_pairs_format = False
        if not comment_lines.empty:
            header_line = comment_lines.iloc[-1, 0]  # last line starting with '#'
            if header_line.startswith(S.HEADER_LINE_BEGINNING):
                data_is_in_pairs_format = True
                file_header_names = header_line.lstrip(S.HEADER_LINE_BEGINNING).split()
                if file_header_names[1:(1 + len(S.USED_COLS))] != S.USED_COLS.to_list():  # ignoring any extra columns
                    raise RuntimeError(f'Invalid input file format: {interas_file}')

        header_line_num = None if data_is_in_pairs_format else 0  # Use 0 for header on 1st line and None for no header
        try:
            df = pd.read_csv(interas_file, sep=S.FIELD_SEPARATOR, comment=S.COMMENT_CHAR,
                             nrows=num_lines_to_check, header=header_line_num)
        except pd.errors.ParserError as parser_error:
            raise InputFileError(interas_file) from parser_error

        used_cols = S.USED_COLS.copy()
        if df.columns.size < used_cols.size:
            raise RuntimeError(f'Not enough number of columns in: {interas_file}')

        # TODO this test is not enough to shift the index (also can't delete validate method as it decides on used_cols)
        if df.columns.size > used_cols.size:
            used_cols.index += 1  # there is one extra column at the beginning

        df = df.rename(used_cols, axis=1)  # this will have no effect if not data_is_in_pairs_format
        df = df[used_cols]  # ensures used_cols exist with no extra cols

        # if used_cols.to_list() == df.iloc[0, :].to_list():  # remove the first row if it was a header row
        #     df = df.iloc[1:, :].reset_index(drop=True)

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
    #     return (self.__name if self.__name else '') + '\n' + \
    #            'At most,' + str(files.__working_files.shape[0]) + ' nucleosome interaction files)'
    #     # str(self.__working_files.head(2).file) + '\n...\n(' + \


class NucInteraMatrix:
    """ Nucleosome-Nucleosome Interaction Matrix

        Using a sparse matrix where elements are recorded as ijv with
        i for rows, j for cols, and v or series
    """

    def __init__(self, nuc_interas: NucInteras, nucs: Nucs, name: str = None):
        self.__nucs = nucs
        self.__nuc_interas = nuc_interas
        self.__name = name

        LOGGER.info('\nPreparing Nuc Interaction Matrices:')
        start = timer()
        refreshed_files = self.__create_nuc_intera_matrix_files()
        LOGGER.info(
            f'Done! Preparing Nuc Interaction Matrices finished in {round((timer() - start) / 60.0, 1)} minutes')

        LOGGER.info(f'\nNumber of updated Nuc Interaction Matrices: {refreshed_files}')
        LOGGER.info(f'\nNuc Interaction Matrices ready with total size: {len(self)}nucs * {len(self)}nucs.')
        LOGGER.info(self)

    def read_nuc_intera_matrix_region(self, chrom: str, start_region: int, end_region: int,
                                      orientation='all', extra_cols=False) -> Optional[pd.DataFrame]:
        """
        :param chrom: e.g. 'II' or 'chrom14'
        :param start_region: positve integer
        :param end_region: positve integer
        :param orientation:  It can be any of ('+-', '-+', '++', '--') or ('inner', 'outer', 'right', 'left') or
        ('tandem', 'all'). Further, 'p' instead of '+', and 'n' or 'm' instead of '-' are acceptable.
        Default value is 'all' which returns all orientations.
        :param extra_cols: additionally return 'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'
        :return: a dense matrix in a DataFrame
        """

        nucs = self.__nucs.find_nucs_in_region(chrom, start_region, end_region)
        if nucs is None:
            return None
        nucs = nucs.id_chrom_start_end.sort_index()

        min_nucs_id = nucs.index[0]
        max_nucs_id = nucs.index[-1]

        ijv = self.read_nuc_intera_matrix(chrom, orientation)
        ijv = ijv.query(f'{min_nucs_id} <= nuc_id1 <= {max_nucs_id}')
        ijv = ijv.query(f'{min_nucs_id} <= nuc_id2 <= {max_nucs_id}')

        # # now converting sparse matrix ijv to a dense matrix:
        # submatrix = pd.DataFrame(0, index=nuc_ids, columns=nuc_ids)
        # submatrix = submatrix.rename_axis(index='nuc_id', columns='nuc_id')
        #
        # for _, r in ijv.iterrows():
        #     submatrix.loc[r.nuc_id1, r.nuc_id2] += r.counts

        if extra_cols:
            ijv = pd.merge(ijv, nucs, how='left', left_on='nuc_id1', right_on='nuc_id')  # suffixes won't work here
            ijv = pd.merge(ijv, nucs, how='left', left_on='nuc_id2', right_on='nuc_id', suffixes=(None, '2'))
            ijv = ijv.rename(columns={'chrom': 'chrom1', 'start': 'start1', 'end': 'end1'})  # fix col names

        return ijv

    def read_nuc_intera_matrix(self, chrom: str, orientation='all') -> pd.DataFrame:
        """
        :param chrom: e.g. 'II' or 'chrom14'
        :param orientation:  It can be any of ('+-', '-+', '++', '--') or ('inner', 'outer', 'right', 'left') or
        ('tandem', 'all'). Further, 'p' instead of '+', and 'n' or 'm' instead of '-' are acceptable.
        Default value is 'all' which returns all orientations.
        :return: a sparse ijv matrix
        """

        working_files = FILES.get_working_files(chrom, orientation, Files.S_MATRIX)
        if working_files.empty:
            return self.__EMPTY

        schema = self.schema()
        matrices = []
        for _, file_info in working_files.iterrows():
            input_file = file_info.subdir / file_info.file
            if not input_file.exists():
                input_file = file_info.subdir / file_info.file_zip  # check if zip version exists
                if not input_file.exists():
                    continue

            try:
                matrix_ijv = pd.read_csv(
                    input_file, sep=S.FIELD_SEPARATOR, comment=S.COMMENT_CHAR, usecols=range(schema.columns.size))
            except pd.errors.ParserError as parser_error:
                raise InputFileError(input_file) from parser_error

            matrices.append(matrix_ijv)

        if len(matrices) == 0:
            return self.__EMPTY

        matrix_ijv = pd.concat(matrices, axis=0, ignore_index=True)

        matrix_ijv = matrix_ijv.groupby(['nuc_id1', 'nuc_id2']).sum().reset_index()  # calculate counts

        matrix_ijv = matrix_ijv.apply(pd.to_numeric, downcast='unsigned')  # for space efficiency

        return matrix_ijv

    def __create_nuc_intera_matrix_files(self) -> int:
        # decide if need to continue
        files_needing_refresh = FILES.get_files_needing_refresh(Files.S_MATRIX)
        if files_needing_refresh.empty:
            return 0

        # subdir_matrix = files_needing_refresh.head(1).subdir.squeeze()  # read subdir from first row

        refreshed_files = 0
        for _, file_info in files_needing_refresh.iterrows():
            if file_info.chrom not in self.chrom_list:
                continue
            output_file = file_info.file_zip if FILES.zipped else file_info.file
            output_file = file_info.subdir / output_file
            LOGGER.info(f'Creating nuc interaction Matrix file: {output_file}')
            df_nuc_inter = self.__nuc_interas.read_nuc_inter(file_info.chrom, file_info.orientation)
            df_nuc_intera_matrix = self.__create_nuc_intera_matrix(df_nuc_inter)
            df_nuc_intera_matrix.to_csv(output_file, sep=S.FIELD_SEPARATOR, index=False, header=True)
            refreshed_files += 1
            # todo decide about when to keep or remove the cache files

        FILES.set_needs_refresh_for_all_files(False, Files.S_MATRIX)

        return refreshed_files

    @staticmethod
    def __create_nuc_intera_matrix(df_nuc_inter: pd.DataFrame) -> pd.DataFrame:

        # Goal is to create a ijv sparse matrix: i:nuc_id1, j:nuc_id2, v:counts
        matrix_ijv = df_nuc_inter

        # make matrix_ijv symmetric
        swapped = pd.DataFrame(matrix_ijv.values, columns=['nuc_id2', 'nuc_id1', 'counts'])  # swapped cols 2,1 or j,i
        matrix_ijv = pd.concat([matrix_ijv, swapped])  # append ij rows with their swapped ji rows

        matrix_ijv['counts'] = np.uint16(1)  # set all rows to 1, but using unsigned type for space efficiency
        matrix_ijv = matrix_ijv.groupby(['nuc_id1', 'nuc_id2']).sum().reset_index()  # calculate counts

        matrix_ijv = matrix_ijv.apply(pd.to_numeric, downcast='unsigned')  # for space efficiency

        diagonal_ids = (matrix_ijv.nuc_id1 == matrix_ijv.nuc_id2)
        matrix_ijv.loc[diagonal_ids, 'counts'] //= 2  # divide by 2 to fix duplicated diagonal counts

        return matrix_ijv

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
        return 'Nucleosome interaction matrix created from: \n' + self.name


class CLI:
    """Command Line Interface"""

    @classmethod
    def handler_prepare_command(cls, command, chroms_file, nucs_file, interas_file, working_dir, refresh, zipped):
        start = timer()

        FILES.reset_all(base_file_name=interas_file, working_dir=working_dir, refresh=refresh, zipped=zipped)

        log_file_name = FILES.working_dir / (FILES.base_file_name.stem + '.log')
        handler_log_file = logging.FileHandler(filename=log_file_name, mode='a')
        LOGGER.addHandler(handler_log_file)
        LOGGER.info(f"\n\n\n = = = = = Arguments for command '{command}':")
        LOGGER.info(pprint.pformat(locals()))  # log all the arguments
        LOGGER.info(f'Logging to file {log_file_name}.\n')

        nucs = Nucs(nucs_file=nucs_file, chroms_file=chroms_file)
        nuc_interas = NucInteras(interas_file=interas_file, nucs=nucs, refresh=refresh)
        nuc_intera_matrix = NucInteraMatrix(nuc_interas=nuc_interas, nucs=nucs)

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
                'Could not produce expected matrix files!')

        if len(nuc_intera_matrix) == 0:
            LOGGER.warning('\nWARNING: The produced matrix is empty!')

        LOGGER.info(f'\nPrepared working directory: {FILES.working_dir}')
        LOGGER.info(f"\n = = = = = Finished command '{command}' in {round((timer() - start) / 60.0, 1)} minutes")

    @classmethod
    def handler_plot_command(cls, command, working_dir, chrom, start_region, end_region, prefix, save_only):
        start = timer()

        FILES.reset_all(working_dir=working_dir)

        log_file_name = FILES.working_dir / (FILES.base_file_name.stem + '.log')
        handler_log_file = logging.FileHandler(filename=log_file_name, mode='a')
        LOGGER.addHandler(handler_log_file)
        LOGGER.info(f"\n\n\n = = = = = Arguments for command '{command}':")
        LOGGER.info(pprint.pformat(locals()))  # log all the arguments
        LOGGER.info(f'Logging to file {log_file_name}.\n')

        nucs = Nucs()  # loads nucs from the working directory

        nuc_interas = NucInteras(nucs=nucs, refresh=False)

        matrix = NucInteraMatrix(nuc_interas=nuc_interas, nucs=nucs)

        # output_name = Path(Path(output_name).name)  # dismiss the path to just save in working dir
        # output_name = Path(f'{output_name.stem}_{chrom}_{start_region}_{end_region}{output_name.suffix}')
        output_name = f'{prefix}_{chrom}_{start_region}_{end_region}'
        submatrices = dict()
        for orient in Files.combined_orientations():  # includes 'all' and 'tandem'
            # output_with_orient = Path(f'{output_name.stem}_{orient}{output_name.suffix}')
            output_with_orient = Path(f'{output_name}_{orient}')
            if output_with_orient.suffix != '.tsv':
                output_with_orient = Path(output_with_orient.name + '.tsv')
            output_with_orient = FILES.working_dir / output_with_orient.name

            submatrix_ijv = matrix.read_nuc_intera_matrix_region(
                chrom, start_region, end_region, orient, extra_cols=True)

            if submatrix_ijv is None or len(submatrix_ijv) == 0:
                LOGGER.info(f'\nNo records found for , {output_with_orient}')
            else:
                submatrix_ijv.to_csv(output_with_orient, S.FIELD_SEPARATOR, header=True, index=True)
                submatrices[orient] = submatrix_ijv
                LOGGER.info(f'{len(submatrix_ijv)} lines saved to file:  {output_with_orient}')

        # create the Heat Map:
        if len(submatrices) > 0:
            cls.make_heat_map(submatrices, chrom, start_region, end_region, output_name, save_only)
        else:
            LOGGER.warning('\n\nNOTE: No plots have been produced because no data is available. '
                           'You many want to consider expanding the chromosome region of study.')

        LOGGER.info(f"\n = = = = = Finished command '{command}' in {round((timer() - start), 0)} seconds")

    @classmethod
    def make_heat_map(cls, matrices, chrom, start_region, end_region, output_name, save_only):
        from bokeh.io import show, save
        from bokeh.models import LinearColorMapper, HoverTool, WheelZoomTool  # BasicTicker, ColorBar
        from bokeh.plotting import figure, output_file
        from bokeh.palettes import Greys, Oranges, Greens, Blues, Purples, Reds

        output_name = Path(output_name)
        suffix = '' if output_name.suffix == '.html' else '.html'
        output_name = FILES.working_dir / (output_name.name + suffix)
        output_file(output_name)

        # seaborn: "the best sequential palettes will be perceptually uniform, meaning that the relative
        # discriminability of two colors is proportional to the difference between the corresponding data values"

        def r(colors, num_colors=256):  # Reverse the order going from light to dark, and select max num of colors
            return list(reversed(colors[num_colors]))

        palettes = dict(
            all=r(Oranges), tandem=r(Greys), right=r(Reds), left=r(Blues), inner=r(Purples), outer=r(Greens))
        # alphas = dict(all=.8, tandem=.4, right=.6, left=.6, inner=.5, outer=.5)

        tools = 'box_zoom, wheel_zoom, pan, reset, save'
        tooltips = [('Nuc1', '@nuc_id1:  @chrom1  @start1-@end1'),
                    ('Nuc2', '@nuc_id2:  @chrom2  @start2-@end2'),
                    ('Counts', '@counts'), ]  # ('log(counts+1)', '@log_counts')])

        p = figure(title=f'Nuc-Nuc Interaction Counts for Chrom {chrom} from pos {start_region} to {end_region}',
                   x_axis_location='above', plot_width=900, plot_height=900,
                   tools=tools, toolbar_location='below')  # tooltips=tooltips,

        full_ijv = pd.concat(matrices.values(), keys=matrices.keys(), names=['orient', None])
        full_ijv['log_counts'] = np.log2(full_ijv.counts + 1)  # Plus 1 to avoid log(0)
        # log_counts_stats = full_ijv.query(f'orient in {Files.orientations()}')['log_counts'].describe()
        log_counts_stats = full_ijv['log_counts'].describe()  # e.g., min and max
        for orient, ijv in full_ijv.groupby('orient'):  # matrices.items():
            legend_label = '' if orient == 'all' else ' (' + ','.join(FILES.get_strands(orient)) + ')'
            legend_label = orient.title() + legend_label
            mapper = LinearColorMapper(
                palette=palettes[orient], low=log_counts_stats['min'], high=log_counts_stats['max'])
            r = p.rect(source=ijv, x='nuc_id1', y='nuc_id2', width=1, height=1,  # alpha=alphas[orient],
                       fill_color={'field': 'log_counts', 'transform': mapper},
                       line_color=None, legend_label=legend_label)

            p.add_tools(HoverTool(renderers=[r], tooltips=[(orient.title(), '')] + tooltips))

            if orient != 'all':
                r.visible = False

        p.toolbar.active_scroll = p.select_one(WheelZoomTool)
        p.legend.location = "bottom_left"
        p.legend.title = 'Orientation'
        p.legend.click_policy = 'hide'
        p.xaxis.axis_label = 'Nucleosome 1'
        p.yaxis.axis_label = 'Nucleosome 2'
        p.y_range.flipped = True
        p.grid.grid_line_color = None
        p.axis.axis_line_color = None
        p.axis.major_tick_line_color = None
        p.axis.major_label_standoff = 0
        p.axis.major_label_text_font_size = "14px"
        p.title.text_font_size = '14pt'
        p.xaxis.axis_label_text_font_size = "14pt"
        p.yaxis.axis_label_text_font_size = "14pt"
        # p.xaxis.major_label_orientation = pi / 2
        # p.match_aspect = True

        # color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="7px",
        #                      ticker=BasicTicker(desired_num_ticks=len(colors)),
        #                      formatter=PrintfTickFormatter(format="%d%%"),
        #                      label_standoff=6, border_line_color=None, location=(0, 0))
        # p.add_layout(color_bar, 'right')

        # # Set up widgets
        # from bokeh.models import TextInput, Slider
        # text_title = TextInput(title="Plot Title", value='')
        # text_xaxis = TextInput(title="X Axis Title", value='')
        # text_yaxis = TextInput(title="Y Axis Title", value='')
        # # offset = Slider(title="offset", value=0.0, start=-5.0, end=5.0, step=0.1)
        #
        # def update_title(attrname, old, new):
        #     p.title.text = text_title.value
        #
        # text_title.on_change('value', update_title)
        #
        # # Set up layouts and add to document
        # from bokeh.layouts import column, row
        # from bokeh.io import curdoc
        # inputs = column(text_title, text_xaxis, text_yaxis)
        #
        # curdoc().add_root(row(p, inputs, width=1400))
        # curdoc().title = "Sliders"

        # noinspection PyBroadException
        if save_only:
            save(p)  # only save the plot
            LOGGER.info(f'\n\nSince --save flag is used, the results plot is saved only and not displayed.')
        else:
            show(p)  # save and show the plot
            LOGGER.info(f'\n\nIf the results plot is not automatically displayed, '
                        f'please open the following plot file manually in a web browser.')

        LOGGER.info(f'\nOutput plot: \n'
                    f'   {output_name} \n'
                    f'or with full path:\n'
                    f'   {Path(output_name).resolve()}\n')

    __help = {  # HELP texts used in different parts of the program
        'main': {
            'description': 'Creates and plots selected sections of nuc-nuc interaction matrices.',
            'title': 'possible commands',
            'help': 'Use one of the commands. Use -h option after any command for more information.',
            'args': {  # todo --quiet not implemented; see help text below
                '--quiet': """
                    Suppress the progress output. Still, a log file will be saved in the working folder.""",
            }
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
            # If the exact expected names are not found in the first non-comment line,
            # then it will be assumed that there is no header row, and the first row
            # is read in as data just like the next rows.
            'args': {
                'chroms_file': """
                    Chromosomes file with one mandatory column "chrom". Any other chromosome names that may appear 
                    in the <nucs_file> or <interas_file> will be ignored if they not present in <chroms_file>.""",
                'nucs_file': """
                    Nucleosomes file with columns: "[id] chrom start end" with the [id] column being optional 
                    (see below).""",
                'interas_file': """
                    Interactions file with columns: "[id] chrom1 pos1 chrom2 pos2 strand1 strand2" with the [id] column 
                    being optional (see below). NOTE: The output of the "Pairs" program is accepted as <interas_file> 
                    without requiring any changes (e.g., the extra columns will be ignored).""",
                # todo what is the appropriate name to use instead of the "Pairs" program?
                '--dir': """
                    Optional output folder to store the nuc-nuc interaction matrices (and other files) 
                    into. Without this option, <working_dir> is derived from the name of the 
                    interaction input file, <interas_file>.""",
                '--refresh': """
                    This will start fresh, re-computing all the intermediate files such as the nuc-nuc matrices. 
                    WARNING: Refreshing starts from scratch and may take a long time to complete, 
                    depending on the size of the input data and your hardware specifications.""",
                '--zip': """
                    Compress intermediary and cache files using gzip. Zipping files saves space but requires more time 
                    to write and read files.""",
            }
        },
        'plot': {
            'help': 'Select and plot a subset (or window) of a nuc-nuc interaction matrix.',
            'epilog': None,
            'args': {
                'working_dir': 'This is the folder created by the "prepare" command.',
                'chrom': """
                    A chromosome, e.g. III or chr21, on which a region is being selected.  
                    Note that the chrom has to be from those mentioned in <chroms_file> provided  
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
        for c in ['prepare', 'plot']:  # disabling some commands for now ['prepare', 'select', 'plot', 'selplot']:
            commands[c] = subparsers.add_parser(c, help=h[c]['help'], description=h[c]['help'], epilog=h[c]['epilog'])

        # # # arguments for main
        c = 'main'
        a = '--quiet'  # and '-q'   # todo Implement --quiet
        parser_main.add_argument(a[1:3], a, help=h[c]['args'][a], default=False, action='store_true')

        # # # arguments for commands
        c = 'prepare'  # c for Command
        a = 'chroms_file'  # a for Argument
        commands[c].add_argument(a, help=h[c]['args'][a], metavar=f'<{a}>')
        a = 'nucs_file'
        commands[c].add_argument(a, help=h[c]['args'][a], metavar=f'<{a}>')
        a = 'interas_file'
        commands[c].add_argument(a, help=h[c]['args'][a], metavar=f'<{a}>')
        a = '--dir'  # and '-d'
        commands[c].add_argument(a[1:3], a, help=h[c]['args'][a], metavar='<working_dir>', dest='working_dir')
        a = '--refresh'  # don't offer -r option to avoid accidentally refreshing
        commands[c].add_argument(a, help=h[c]['args'][a], default=False, action='store_true')
        a = '--zip'  # not zipped
        commands[c].add_argument(a[1:3], a, help=h[c]['args'][a], dest='zipped', default=False, action='store_true')

        c = 'plot'  # c for Command
        a = 'working_dir'  # a for Argument
        commands[c].add_argument(a, help=h[c]['args'][a], metavar=f'<{a}>')
        a = 'chrom'
        commands[c].add_argument(a, help=h[c]['args'][a], metavar=f'<{a}>')
        a = 'start_region'
        commands[c].add_argument(a, help=h[c]['args'][a], metavar=f'<{a}>', type=int)
        a = 'end_region'
        commands[c].add_argument(a, help=h[c]['args'][a], metavar=f'<{a}>', type=int)
        a = '--save'
        commands[c].add_argument(a[1:3], a, help=h[c]['args'][a], dest='save_only', default=False, action='store_true')
        a = '--prefix'
        commands[c].add_argument(
            a[1:3], a, help=h[c]['args'][a], metavar='<outfile_prefix>', dest='prefix', default='plot')
        # a = 'orientation'  # '-o' and '--orientation'
        # commands[c].add_argument(f'-{a[0]}', f'--{a}', metavar=f'<{a}>', dest=a, help=h[c]['args'][f'--{a}'])

        return parser_main

    @staticmethod
    def test_parse_args(parser):
        LOGGER.debug('WARNING: Using TESTING mode!!')
        # # # for development only:
        # chroms_file = 'Chromosomes_Human.txt'
        chroms_file = r'last_version\Chromosomes_Yeast.txt'
        nucs_file = r'last_version\Nucleosomes_Yeast.txt'
        # interas_file = r'C:\Temp\Orientation\cat\A1-3NC-972-form_S7_R1_2_001_ALL.txt'
        interas_file = r'last_version\Yeast.tsv.gz'
        working_dir = 'Yeast-win'

        # testing_args = ''
        # testing_args = '--help'
        # testing_args = 'prepare --help'
        # testing_args = 'plot --help'
        testing_args = f'prepare {chroms_file} {nucs_file} {interas_file} --dir {working_dir}'
        # testing_args = f'plot {working_dir} II 1 50000 --prefix my_plot'
        arguments = parser.parse_args(testing_args.split())  # for debugging

        return arguments


def main():
    # logging.basicConfig(filename=Path(__file__).stem + '.log', level=logging.DEBUG, filemode='w')
    # LOGGER = logging.getLogger()
    # formatter = logging.Formatter('%(asctime)p %(levelname)p %(message)p')

    commandline_parser = CLI.create_commandline_parser()

    # for testing only
    import platform
    if platform.system() == 'Windows':
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
        'prepare': CLI.handler_prepare_command,
        'plot': CLI.handler_plot_command,
    }.get(args.command)

    command_handler(**vars(args))  # call the selected handler with all the arguments from commandline

    # # before exiting the program, remove the log file if it's empty
    # @atexit.register
    # def cleanup_tasks():
    #     for log_handler in iter(logger.handlers):
    #         log_handler.flush()
    #         log_handler.close()
    #
    #     if Path(log_file_name).stat().st_size == 0:
    #         Path(log_file_name).unlink()
    #         # todo Should remove working dir IF empty?

    pass


if __name__ == "__main__":
    LOGGER = logging.getLogger()
    LOGGER.setLevel(level=logging.DEBUG)
    FILES = Files()  # each command handler will reset FILES with appropriate parameters

    main()

# done Min Python version was reduced to 3.6
# done input checking, e.g.: validate_chroms_file, validate_nucs_file, validate_interas_file
# done check for zipped Nucs file
# done reverse y axis
# done input/output error checking
# done add a --no-zip flag!
# done doc: bulleted list. Efficiency; scalability (breaking done: input, matrix, submatix)
# TODO doc: graphs:  1) usage pipe line  2) how efficiency is achieved
# TODO doc: running time measurements (Yeast, Human, PC, HPC)
# done github
# TODO deployment methods? https://conda.io/projects/conda-build/en/latest/user-guide/tutorials/index.html
# done add Tandem
# done add to Legend '+-', '-+', '++', '--'
# done Colors, Font size
# todo Bokeh, Full GUI app?
# done Bokeh, default tool? on scroll zoom
# todo Bokeh, multiple HoverTools (tooltips)
# todo Bokeh, keep aspect ratio the same as zooming in and out
# todo Bokeh, resize plot?
# todo Efficiency: use parallelism to utilize multicore CPUs
