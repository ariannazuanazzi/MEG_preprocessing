#!/usr/bin/env python 3.7.6
# -*- coding: utf-8 -*-
"""
part of MEG_processing_general
"""
'''
script_version == 1
Numpy == 1.18.1
MNE == 0.20.8
python == 3.7.6
matplotlib == 3.2.2
'''

__author__ = "Arianna Zuanazzi, originally from Hao Zhu and Laura Gwilliams"
__copyright__ = "Copyright 2019"
__date__ = '2020/09/17'
__license__ = "MIT"
__version__ = "1.2.0"
__email__ = "az1864@nyu.edu"


#################################################################################################################
# IMPORT NEEDED MODULES
#################################################################################################################
import numpy as np
from matplotlib import pyplot as plt
import os
import mne
from PyQt5 import QtWidgets
import sys
import pandas as pd
from scipy.io import loadmat
mne.set_log_level(False)
from general_info_general import *

#################################################################################################################
# SET BASIC FUNCTIONS
#################################################################################################################
# setup dialog boxes
class Dialog (QtWidgets.QDialog) :
    from PyQt5 import QtWidgets

    def __init__ (self , dinput) :
        super (Dialog , self).__init__ ()
        self.createFormGroupBox (dinput)

        buttonBox = QtWidgets.QDialogButtonBox (QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttonBox.accepted.connect (self.accept)
        buttonBox.rejected.connect (self.reject)

        mainLayout = QtWidgets.QVBoxLayout (self)
        mainLayout.addWidget (self.formGroupBox)
        mainLayout.addWidget (buttonBox)

        self.setWindowTitle ("Input Box")

    def createFormGroupBox (self , dinput) :
        layout = QtWidgets.QFormLayout ()
        self.linedit1 = QtWidgets.QLineEdit ('')
        self.combox1 = QtWidgets.QComboBox ()
        self.combox1.setToolTip ('Choose Your Status')
        self.combox1.addItems ([ '0' ])
        self.spinbox1 = QtWidgets.QSpinBox ()

        for text , w in zip (dinput , (self.linedit1 , self.combox1 , self.spinbox1)) :
            layout.addRow (text , w)

        self.formGroupBox = QtWidgets.QGroupBox ("")
        self.formGroupBox.setLayout (layout)

    def accept (self) :

        self._output = self.linedit1.text ()  # , self.combox1.currentText(), self.spinbox1.value()
        super (Dialog , self).accept ()

    def get_output (self) :
        return self._output

# setup dialog box to select ICA components
# def call_msg_box2proceed (title = "Proceed?") :
#
#     app = QtWidgets.QApplication (sys.argv)
#     # Staring Functions for Execution
#     dinput = [ title ]
#     # Call the UI and get the inputs
#     dialog = Dialog (dinput)
#     if dialog.exec_ () == Dialog.Accepted :
#         proceed = dialog.get_output ()
#         if proceed == 1 :
#             return


# setup dialog box to select ICA components
def call_msg_box (title = 'ICA Component') :

    app = QtWidgets.QApplication (sys.argv)
    # Staring Functions for Execution
    dinput = [ title ]
    # Call the UI and get the inputs
    dialog = Dialog (dinput)
    if dialog.exec_ () == Dialog.Accepted :
        ICAs = dialog.get_output ()
        return ICAs


# setup function to use eyelink data
def convert_sacevents (onlyevent_file) :

    file = open (onlyevent_file , 'r')
    data = file.readlines ()

    # chunk sac&blk out of rawdata
    pick = [ ]
    for i in data :
        if 'START' in i :
            pick = [ ]
            pick.append (i.split ())
        elif 'ESACC' in i :
            pick.append (i.split () [ :5 ])
        elif 'EBLINK' in i :
            pick.append (i.split ())
        elif 'END' in i :
            pick.append (i.split ())
            break
        else :
            continue

    # filter sac&blk if duration is too short
    f_pick = [ ]
    min_duration = 30

    for p in pick :
        if p [ 0 ] == 'ESACC' :
            if int (p [ -1 ]) < min_duration :
                pass
            else :
                f_pick.append (p)
        elif p [ 0 ] == 'EBLINK' :
            if int (p [ -1 ]) < min_duration :
                pass
            else :
                f_pick.append (p)
        else :
            f_pick.append (p)

    # create eog array
    start = int (f_pick [ 0 ] [ 1 ])
    end = int (f_pick [ -1 ] [ 1 ])

    SAC = np.zeros (end - start ,)
    BLK = np.zeros (end - start ,)

    blinks = [ ]
    # create SAC&BLK array
    for p in f_pick [ 1 :-1 ] :
        if p [ 0 ] == 'ESACC' :
            ssacc = int (p [ 2 ]) - start
            esacc = int (p [ 3 ]) - start
            SAC [ ssacc :esacc ] = 1
        else :
            sblink = int (p [ 2 ]) - start
            eblink = int (p [ 3 ]) - start
            blinks.append ([ sblink , eblink ])
            BLK [ sblink :eblink ] = 4

    # filter sac&blk overlapping: remove saccade that is overlapped with blink
    for p in f_pick [ 1 :-1 ] :
        if p [ 0 ] == 'ESACC' :
            ssacc = int (p [ 2 ]) - start
            esacc = int (p [ 3 ]) - start
            for pair in blinks :
                if (ssacc < pair [ 0 ] < esacc) or (pair [ 0 ] < ssacc < pair [ 1 ]) :
                    SAC [ ssacc :esacc ] = 0
                else :
                    continue
    return [ SAC , BLK ]

# find eye events
def get_eyelink (raw , events , eyelink_file , eog_channel , high_pass , low_pass , save_plot, figure_path) :

    eyelink = convert_sacevents (eyelink_file)
    start = events [ 0 ] [ 0 ]

    raw.set_channel_types ({ eog_channel : 'eog' })
    # set eye tracking channel as eog channel

    eog_data = raw.copy ().pick_types (meg = False , eog = True)
    eog = eog_data.get_data () [ 0 ] [ start : ]
    eog = eog - eog.mean ()

    blinks = np.squeeze (np.where (eog < -0.5))
    blinks_time = [ blinks [ 0 ] ]
    for b in blinks :
        if b - blinks_time [ -1 ] > 1000 :
            blinks_time.append (b)

    b_eyelink = np.squeeze (np.where (eyelink [ 1 ] > 1))
    b_eyelink_time = [ b_eyelink [ 0 ] ]
    for b in b_eyelink :
        if b - b_eyelink_time [ -1 ] > 400 :
            b_eyelink_time.append (b)

    time_difference = np.array (b_eyelink_time [ 0 ]) - np.array (blinks_time [ 0 ])

    delay = int (time_difference.mean ())

    print (f'Delay of eyelink data and recording MISC 019 is {delay} ms')
    print (f'Fixing delay ......')

    eyelink [ 1 ] = eyelink [ 1 ] [ delay : ]
    eyelink [ 0 ] = eyelink [ 0 ] [ delay : ]

    print (f'Delay fixed')


    eyep = plt.figure (figsize = (20 , 6))
    plt.subplot (211)
    plt.plot (eog , label = 'MISC 019')
    plt.plot (eyelink [ 0 ] , label = 'SACCADE')
    plt.legend ()
    plt.title ('SACCADE')

    plt.subplot (212)
    plt.plot (eog , label = 'MISC 019')
    plt.plot (eyelink [ 1 ] , label = 'BLINK')
    plt.legend ()
    plt.title ('BLINK')

    plt.tight_layout ()
    plt.show ()

    if save_plot :
        eyep.savefig (f'{figure_path}eye_events.png')
    else :
        pass

    # raw_filtered.add_channels([eog_data])
    sac_info = eog_data.info.copy ()
    sac_info [ 'ch_names' ] = [ 'HEOG' ]
    sac_info [ 'chs' ] [ 0 ] [ 'ch_name' ] = 'HEOG'
    sac_info [ 'highpass' ] = high_pass
    sac_info [ 'lowpass' ] = low_pass

    blk_info = eog_data.info.copy ()
    blk_info [ 'ch_names' ] = [ 'VEOG' ]
    blk_info [ 'chs' ] [ 0 ] [ 'ch_name' ] = 'VEOG'
    blk_info [ 'highpass' ] = high_pass
    blk_info [ 'lowpass' ] = low_pass

    saccde_data = np.concatenate (
        [ np.zeros (start ,) , eyelink [ 0 ] , np.zeros (raw.n_times - len (eyelink [ 0 ]) - start ,) ])
    blink_data = np.concatenate (
        [ np.zeros (start ,) , eyelink [ 1 ] , np.zeros (raw.n_times - len (eyelink [ 1 ]) - start ,) ])
    sac_channel = mne.io.RawArray (saccde_data.reshape (1 , raw.n_times) , sac_info)
    blk_channel = mne.io.RawArray (blink_data.reshape (1 , raw.n_times) , blk_info)

    return sac_channel , blk_channel

# ensure that the triggers match
def match_list(A, B, on_replace='delete'):
    """Match two lists of different sizes and return corresponding indice

    created by JR king

    Parameters
    ----------
    A: list | array, shape (n,)
        The values of the first list
    B: list | array: shape (m,)
        The values of the second list

    Returns
    -------
    A_idx : array
        The indices of the A list that match those of the B
    B_idx : array
        The indices of the B list that match those of the A
    """
    from Levenshtein import editops #pip install python-Levenshtein

    A = np.nan_to_num(np.squeeze(A))
    B = np.nan_to_num(np.squeeze(B))
    assert A.ndim == B.ndim == 1

    unique = np.unique(np.r_[A, B])
    label_encoder = dict((k, v) for v, k in enumerate(unique))

    def int_to_unicode(array):
        return ''.join([str(chr(label_encoder[ii])) for ii in array])

    changes = editops(int_to_unicode(A), int_to_unicode(B))
    B_sel = np.arange(len(B)).astype(float)
    A_sel = np.arange(len(A)).astype(float)
    for type, val_a, val_b in changes:
        if type == 'insert':
            B_sel[val_b] = np.nan
        elif type == 'delete':
            A_sel[val_a] = np.nan
        elif on_replace == 'delete':
            # print('delete replace')
            A_sel[val_a] = np.nan
            B_sel[val_b] = np.nan
        elif on_replace == 'keep':
            # print('keep replace')
            pass
        else:
            raise NotImplementedError
    B_sel = B_sel[np.where(~np.isnan(B_sel))]
    A_sel = A_sel[np.where(~np.isnan(A_sel))]
    assert len(B_sel) == len(A_sel)
    return A_sel.astype(int), B_sel.astype(int)


import numpy as np


# least square
def least_square_reference (inst , empty_room = None , max_times_samples = 2000 ,
                             bad_channels = None , scaler = None , mrk = None ,
                             elp = None , hsp = None) :
    """
    # downloaded function least_square_reference from https://github.com/kingjr/jr-tools/blob/master/jr/meg/kit.py and added to base_funcs
    Fits and applies Least Square projection of the reference channels
    (potentially from an empty room) and removes the corresponding component
    from the recordings of a subject.
    Parameters
    ----------
        inst : Raw | str
            Raw instance or path to raw data.
        empty_room : str | None
            Path to raw data acquired in empty room.
        max_times_samples : int
            Number of time sample to use for pinv. Defautls to 2000
        bad_channels : list | array, shape (n_chans) of strings
            Lists bad channels
        scaler : function | None
            Scaler functions to normalize data. Defaults to
            sklearn.preprocessing.RobustScaler.
    Returns
    -------
        inst : Raw
    adapted from Adeen Flinker 6/2013 (<adeen.f@gmail.com>) LSdenoise.m
    Main EHN
        - Automatically detects channel types.
        - Allows flexible scaler; Robust by default.
        - The data is projected back in Tesla.
        - Allows memory control.
    TODO:
        - Allow other kind of MNE-Python inst
        - Allow baseline selection (pre-stim instead of empty room)
        - Clean up memory
        - Allow fancy solver (l1, etc)
    """
    from scipy.linalg import pinv
    from mne.io import read_raw_fif
    from mne.io import BaseRaw

    # Least square can be fitted on empty room or on subject's data
    if empty_room is None :
        if not isinstance (inst , BaseRaw) :
            raw = read_raw_fif (inst , preload = True)
        else :
            raw = inst
    else :
        if not isinstance (empty_room , BaseRaw) :
            raw = read_raw_fif (empty_room , preload = True)
        else :
            raw = empty_room

    # Parameters
    n_chans , n_times = raw._data.shape
    chan_info = raw.info [ 'chs' ]

    # KIT: axial gradiometers (equiv to mag)
    ch_mag = np.where ([ ch [ 'coil_type' ] == 6001 for ch in chan_info ]) [ 0 ]
    # KIT: ref magnetometer
    ch_ref = np.where ([ ch [ 'coil_type' ] == 6002 for ch in chan_info ]) [ 0 ]
    # Other channels
    ch_misc = np.where ([ ch [ 'coil_type' ] not in [ 6001 , 6002 ]
                           for ch in chan_info ]) [ 0 ]

    # check if refs is included
    assert len(ch_ref) != 0, "MEG refs are not among the channels! They are needed for denoise!"

    # Bad channel
    ch_bad = np.empty (0)
    if (bad_channels is not None) and len (bad_channels) :
        if np.all ([ isinstance (ch , int) for ch in bad_channels ]) :
            bad_channels = np.array (bad_channels)
        elif np.all ([ isinstance (ch , str) for ch in bad_channels ]) :
            bad_channels = [ ii for ii , ch in enumerate (raw.ch_names)
                             if ch in bad_channels ]
        else :
            raise ValueError ('bad_channels needs array of int or array of str')
    else :
        bad_channels = [ ]
    default_bad_channels = [ ii for ii , ch in enumerate (raw.ch_names)
                             if ch in raw.info [ 'bads' ] ]
    bad_channels = np.array (default_bad_channels + bad_channels , int)

    print ('bad channels:' , [ raw.ch_names [ bad ] for bad in bad_channels ])
    # To avoid memory error, let's subsample across time
    sel_times = slice (0 , n_times , int (np.ceil (n_times // max_times_samples)))

    # Whiten data
    if scaler is None :
        from sklearn.preprocessing import RobustScaler
        scaler = RobustScaler ()
    data_bsl = scaler.fit_transform (raw._data.T)

    # Fit Least Square coefficients on baseline data
    empty_sensors = data_bsl [ : , ch_mag ]
    if len (ch_bad) :
        empty_sensors [ : , ch_bad ] = 0  # remove bad channels
    coefs = np.dot (pinv (data_bsl [ sel_times , ch_ref ]) ,
                     empty_sensors [ sel_times , : ])
    empty_sensors , data_bsl = None , None  # clear memory

    # Apply correction on subject data
    if empty_room is not None :
        del raw
        raw = read_raw_fif (inst , preload = True)

    data_subject = scaler.transform (raw._data.T)
    subject_sensors = (data_subject [ : , ch_mag ] -
                       np.dot (data_subject [ : , ch_ref ] , coefs))

    # Remove bad channels
    if len (ch_bad) :
        subject_sensors [ : , ch_bad ] = 0

    # Reproject baseline
    new_ref = np.dot (subject_sensors , pinv (coefs))

    # Un-whiten data to get physical units back
    data = np.concatenate ((subject_sensors , new_ref ,
                             raw._data [ ch_misc , : ].T) , axis = 1)
    data = scaler.inverse_transform (data)

    # Output
    raw._data = data.T
    return raw


# bad channels for non-epoched signal
def find_bad_channels(raw, events, window_size=200, epoch_tmax=60.,
                      low_percentile=4., high_percentile=97.,
                      freq_selected=.4):
    '''
    created by Laura Gwilliams (lg5@nyu.edu)  and Arianna Zuanazzi (az1864@nyu.edu) - March2020

    INPUT:

    raw: instance of mne raw object
    events: trigger events
    window_size: how many time samples to include in each computation
    epoch_tmax: (s) how much data to cut around each event
    low_percentile: mark channels that have a range less than this percentile
    high_percentile: mark channels that have a range larger than this percentile
    freq_selected: pick channels that are marked at least this proportion of the time


    OUTPUT:

    list of identified bad channels
    '''

    # make some epochs
    epochs = mne.Epochs(raw, events=events, tmin=0., tmax=epoch_tmax,
                        preload=True, baseline=None)

    # concatenate epoch data and get stds
    epoch_concat = np.concatenate(np.transpose(epochs.pick_types(meg=True)._data,
                                               [0, 2, 1]),
                                                axis=0)

    # now we take 1 second chunks, compute the range, and see how this fits with the distribution
    n_times, n_chs = epoch_concat.shape
    n_windows = int(n_times/window_size)
    picked_chs = np.zeros([n_windows, n_chs])
    wind_ranges = np.zeros([n_windows, n_chs])

    for wi in range(n_windows):

        # time dims
        tstart = window_size*wi
        tstop = tstart+window_size

        # get data
        data_chunk = epoch_concat[tstart:tstop, :]

        # get standard deviation
        ch_ranges = np.abs(np.min(data_chunk, axis=0)) + np.abs(np.max(data_chunk, axis=0))
        wind_ranges[wi, :] = ch_ranges

        # remove flats
        thresh = np.percentile(ch_ranges, low_percentile)
        picked_flats = np.where(ch_ranges <= thresh)
        picked_chs[wi, picked_flats] = 1

        # remove crazies
        thresh = np.percentile(ch_ranges, high_percentile)
        picked_flats = np.where(ch_ranges >= thresh)
        picked_chs[wi, picked_flats] = 1

    # define selected channels
    picked_chs_bool = picked_chs.mean(0) > freq_selected
    raw.info['bads'] = (np.array(raw.info['ch_names'][0:n_chs])[picked_chs_bool]).tolist()

    print("Identified %s bad channels: %s" % (len(raw.info['bads']), raw.info['bads']))

    return raw.info['bads']

# bad channels for epoched signal
def find_bad_channels_epochs(epochs, megchs=157, low_percentile=3., high_percentile=97.,
                             freq_selected=.2):
    '''
    created by Laura Gwilliams (lg5@nyu.edu) and Arianna Zuanazzi (az1864@nyu.edu) - March2020

    INPUT:

    epochs: instance of mne epochs object, sliced to have only meg and not ref_meg
    low_percentile: mark channels that have a range less than this percentile
    high_percentile: mark channels that have a range larger than this percentile
    freq_selected: pick channels that are marked at least this proportion of the time

    OUTPUT:

    list of identified bad channels
    '''

    # params
    n_trials, n_chs, n_times = epochs._data[:, 0:megchs, :].shape # only use meg channels, not refs
    epoch_data = epochs._data[:, 0:megchs, :]

    # now we take 1 second chunks, compute the range, and see how this fits with the distribution
    picked_chs = np.zeros([n_trials, n_chs])
    wind_ranges = np.zeros([n_trials, n_chs])

    for trial_n in range(n_trials):

        # get data
        data_chunk = epoch_data[trial_n, ...].T

        # get standard deviation
        ch_ranges = np.abs(np.min(data_chunk, axis=0)) + np.abs(np.max(data_chunk, axis=0))
        wind_ranges[trial_n, :] = ch_ranges

        # remove flats
        thresh = np.percentile(ch_ranges, low_percentile)
        picked_flats = np.where(ch_ranges <= thresh)
        picked_chs[trial_n, picked_flats] = 1

        # remove crazies
        thresh = np.percentile(ch_ranges, high_percentile)
        picked_flats = np.where(ch_ranges >= thresh)
        picked_chs[trial_n, picked_flats] = 1

    # define selected channels
    picked_chs_bool = picked_chs.mean(0) > freq_selected
    bads = (np.array(epochs.info['ch_names'][0:n_chs])[picked_chs_bool]).tolist()

    print("Identified %s bad channels: %s" % (len(bads), bads))

    return bads


def save_bad_channel_list(bad_channels, fname):
    # save list of bad channels
    df = pd.DataFrame()
    df['bads'] = bad_channels
    df.to_csv(fname)
    pass
