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
import matplotlib
matplotlib.use ('TkAgg')
import glob as glob
import re as re
import mne
import os

#################################################################################################################
# BASIC INFO TO SET: SUBJECT, MARKERS, RUN NUMBER, BLOCK NUMBER, EYE FILE, ER FILE, SAVE, VIEW PLOT
#################################################################################################################
subject = 'R1551'
# set subject name

marker_1 = 0
marker_2 = 1
# select markers that were saved correctly

blockonset_recorded = 0
# whether block onset trigger was recorded

concat_blocks = False
# if analysis done on single blocks and epochs concatenated later

#blockn = [1, 2, 3, 4, 5, 6, 7, 8];
blockn = [1] #if only one block
# specify run numbers

if concat_blocks==True:
    blockn_folder = 'all/'
else:
    blockn_folder = str(blockn).strip('[]') + '/'
# select the name of the save folder based on whether data are concat or not

blocks = [f'_Twotones_block{block}' for block in blockn]
blocks_behav = ''
# blocks_behav = blocks #if only one block
# define block number

eye_link_events = False  #if eyedata exist you will have to concat the blocks in main preprocessing script unless
# a different eye file was saved for each block
# eyedata have to be separately converted from edf to asc and saved as _onlyevents
# Eyedata from eyelink in asc format

eye = '_onlyevents'
# eye file name

emptyroom = '_EmptyRoom'
# emptyroom name

save_plot , dpi = True , 360
# save plot into figure directory, default False
view_plot = True
# view plot

#################################################################################################################
# MAIN FOLDERS AND MAIN FILE NAMES
#################################################################################################################
mainFolder = 'D:\\PostDoc_NYU\\Experiments\\MEG_Checks\\1stlevel'
# main folder

sharedFolder = 'D:\\PostDoc_NYU\\Experiments\\MEG_Checks\\shared'
# folder with shared files

behavioralfolder = 'behav'
#behavioralfolder = 'behav\\blocks' #if only one block
# folder where behav results are saved

megfolder = 'meg'
# folder where MEG data are saved

anatomyFolder = 'anat'
# folder where T1 is saved

megformat_old = '.sqd'
megformat_new = '.fif'
# format of meg data

eyeformat = '.asc'
# format of eye data

elpformat = '*.elp'
hspformat = '*.hsp'
mrkformat = '*Marker*'
# formats for conversion of maindata to fif

txtformat = '.p'
# text format

textformat = '.txt'
# text format

behavformat = '.mat'
# behavioral folder

### DEFINE EACH BLOCK'S NAME ###
maindata_sqd = []
maindata_fif = []

for block_idx, block in enumerate(blockn):
    maindata_sqd.append(subject + "*" + blocks[block_idx] + "*" + megformat_old)
    maindata_fif.append(subject + "*" + blocks[block_idx] + "*" + megformat_new)
    # set main data file name

    #behavdata_mat.append(subject + "*" + blocks[block_idx] + behavformat) #if only one block
    # set behavioral file names
##################################

behavdata_mat = "*" + subject + "*" + blocks_behav + behavformat
# set behavioral file names

emptyroomdata_sqd = subject + "*" + emptyroom + "*" + megformat_old
emptyroomdata_fif = subject + "*" + emptyroom + "*" + megformat_new
# set emptyroom file name

eyedata = subject + eye + "*" + eyeformat
# set eyedata file name

alldata_fif = subject + '_allblocks' + megformat_new
# set all MEG blocks file name

kitpos = 'KIT-157_new2.lout'
# layout file that contains the MEG electrodes' position, from the Japanese engineer group (new system)

kitsel = 'KIT.sel'
# selection file that marks the corresponding electrode to different anatomical locations, from Hao (new system)

#################################################################################################################
# NAMES OF ADDITIONAL FOLDERS
#################################################################################################################
subjectfolder = 'subjects'
# where all subjects are saved

saveFolder = '_saveddata'
# folder where to save data

figureFolder = '_figures/'
# folder where to save data

evokedfigureFolder = 'evoked/'
# folder where to save evoked results

#################################################################################################################
# SET UP DIRECTORIES
#################################################################################################################

dataFolder = os.path.join (mainFolder , subjectfolder, subject, megfolder)
# MEG data folder

behavFolder =  os.path.join(mainFolder, subjectfolder, subject, behavioralfolder, subject, 'MEG', 'main', str(blockn)[1])
# behavioral data folder

figure_path = os.path.join (mainFolder , subjectfolder, subject , figureFolder, blockn_folder)
if not os.path.exists (figure_path) :
    os.makedirs (figure_path)
# set figure directory

evokedfigure_path = os.path.join (mainFolder , subjectfolder, subject , figureFolder, blockn_folder, evokedfigureFolder)
if not os.path.exists (evokedfigure_path) :
    os.makedirs (evokedfigure_path)
# set evoked figure directory

save_path = os.path.join (mainFolder , subjectfolder, subject , saveFolder, blockn_folder)
if not os.path.exists (save_path) :
    os.makedirs (save_path)
# set saved file directory

anat_path = os.path.join (mainFolder , subjectfolder, subject , anatomyFolder)
if not os.path.exists (anat_path) :
    os.makedirs (anat_path)
# set up anat (MRI) path

log_path = os.path.join (mainFolder , subjectfolder, subject , saveFolder)
if not os.path.exists (log_path) :
    os.makedirs (log_path)
# set saved file directory

#################################################################################################################
# SET PARAMETERS
#################################################################################################################
ica_method = 'fastica'
# ica method to use

ica_components = 20  # amount of variance explained
'''
Controls the number of PCA components from the pre-ICA PCA entering the ICA
decomposition in the ICA.fit() method. If None (default), all PCA components
will be used (== max_pca_components). If int, must be <= max_pca_components.
If float, value between 0 and 1 to select n_components based on the percentage
of variance explained by the PCA components.
pre-ICA PCA can be potentially opted out in fastica. For reasons why to opt out see https://doi.org/10.1016/j.neuroimage.2018.03.016.
'''

ica_decim = 3
# save time to set decim bigger than 1

n_max_ecg , n_max_eog = 3 , 1
# maximum number of components to reject

stim_channel = 'STI 014'
# set trigger channel

'''
define triggers
triggers = [1]; %all the triggers

trig.singletone = triggers(1); %block starts

'''

eog_channels = 'MISC 019' #'MISC 017'
# set eye tracking channel

badchannelsUI = ['none']#['MEG 122']#[ 'MEG 023', 'MEG 087', 'MEG 115', 'MEG 134', 'MEG 152' ]
# some sensors that are always noisy - user input

## NEW SYSTEM
#bad_channels = [ f'MEG 0{m}' for m in [ 13 , 19 , 21 , 26 , 32 , 39 , 44 , 52 , 61 , 70 , 86 ] ] +
#                [ f'MEG {n}' for in [ 103 , 110 , 115 , 141 ] ]
# some sensors that are always noisy

tmin = -0.2
tmax = 0.8
# set epoch duration


low_pass = 40
high_pass = 1.0
filter_type = 'bandpass'
# choose filter type: either 'bandpass' for direct bandpass filter or 'hl' for highpass and lowpass filter in order
# filer also depends on the line of the country
'''
ICA is sensitive to low-frequency drifts and therefore requires the data to be high-pass filtered prior to fitting.
Typically, a cutoff frequency of 1 Hz is recommended.
'''

aud = 0
# if using StimTracker sending pulses for actual auditory onsets

baseline = (tmin , 0)
# set baseline

# in some cases the block onset trigger was not recorded so this takes care of the case when it was not
# if blockonset_recorded:
event_id = { 'tone1' : 1, 'tone2' : 2}
#event_id = { 'tone1' : 1}

scaling = 0.000000000005
# scalings for plot
#################################################################################################################
# SET NAMES FOR FILES AND CONVERT SQD
#################################################################################################################

### DEFINE EACH BLOCK'S FILE ###
raw_file_sqd = []
raw_file_fif = []
for block_idx, block in enumerate(blockn):
    raw_file_sqd.append(os.path.join (dataFolder, maindata_sqd[block_idx]))
    raw_file_fif.append(os.path.join (dataFolder, maindata_fif[block_idx]))
    # raw MEG data file

ER_file_sqd = os.path.join (dataFolder, emptyroomdata_sqd)
ER_file_fif = os.path.join (dataFolder, emptyroomdata_fif)
# empty room recording file / interpolated file

if eye_link_events: #if eyedata exist you will have to concat the blocks in main preprocessing script
    eyelink_file = glob.glob(os.path.join (dataFolder,  eyedata))[0]
    # eyelink data file

#behav_file_mat = glob.glob(os.path.join (behavFolder, behavdata_mat))[0]
# behavioral data file

######## CONVERT SQD to FIF#################
# if it does not already exist
filesinfolder = glob.glob(os.path.join (dataFolder, '*' + megformat_old))#list of sqd files
if any(blocks[block_idx] in string for string in filesinfolder):#only if there is an sqd
    raw_file_2convert = []
    raw_file_converted = []
    for block_idx, block in enumerate(blockn):
        raw_file_2convert.append(glob.glob(raw_file_sqd[block_idx])[0])
        # names of files to convert (str)
        raw_file_converted.append(re.sub(str(megformat_old), str(megformat_new), str(raw_file_2convert[block_idx])))
        # names of converted (str)
        # raw_file_converted = ['D:\\PostDoc_NYU\\Experiments\\Trajectories\\1stlevel\\subjects\\R1435\\meg\\R1435_block1_9.26.18.fif']

ER_file_2convert = glob.glob(ER_file_sqd)[0]
# names of files to convert (str)
ER_file_converted = re.sub(str(megformat_old), str(megformat_new), str(ER_file_2convert))
# names of converted (str)

markers = glob.glob(os.path.join (dataFolder, mrkformat))
# markers

elpfile = glob.glob(os.path.join (dataFolder, elpformat))[0]
# location of fiducials

hspfile = glob.glob(os.path.join (dataFolder, hspformat))[0]
# headshape info

## CONVERTS
for block_idx, block in enumerate(blockn):
    filesinfolder = glob.glob(os.path.join (dataFolder, '*' + megformat_new))#list of all fif files
    if not any(blocks[block_idx] in string for string in filesinfolder):#only if there isn't a fif
        raw2save = mne.io.read_raw_kit(raw_file_2convert[block_idx], [markers[marker_1], markers[marker_2]], elpfile, hspfile , allow_unknown_format = True)
        raw2save.save(raw_file_converted[block_idx])
    # converts raw from sqd to fif

if not (glob.glob(os.path.join (dataFolder, emptyroomdata_fif))):
    ER2save = mne.io.read_raw_kit(ER_file_2convert, [markers[marker_1], markers[marker_2]], elpfile, hspfile , allow_unknown_format = True)
    ER2save.save(ER_file_converted)
# converts ER from sqd to fif
else:
    ER_file_converted = glob.glob(os.path.join (dataFolder, emptyroomdata_fif))[0]

### CONCATENATE ALL BLOCKS ###########
if not os.path.exists (glob.glob(os.path.join (dataFolder, maindata_fif[block_idx]))[0]):
    if concat_blocks:
        raw_alldata_fileName = os.path.join(dataFolder, alldata_fif)
        if not os.path.exists (raw_alldata_fileName) :
            raw_file_converted_all = []
        for block_idx, block in enumerate(blockn):
            raw_file_converted_all.append(mne.io.read_raw_fif (raw_file_converted[block_idx], preload=True))

        raw_alldata = mne.io.concatenate_raws(raw_file_converted_all)
        # concatenate all blocks
        raw_alldata.save(raw_alldata_fileName, overwrite=True) #saves data split in 2 files if one is too big
    else:
        raw_alldata_fileName = raw_file_converted[block_idx]
else:
    raw_alldata_fileName = glob.glob(os.path.join (dataFolder, maindata_fif[block_idx]))[0]
###########################################

KIT_layout = mne.channels.read_layout (os.path.join (sharedFolder , kitpos))
# load KIT system layout

KIT_selection = os.path.join (sharedFolder , kitsel)
# selection file for KIT system, partially modified (temporal), need update

#################################################################################################################
# SET NAMES FOR DERIVED FILES
#################################################################################################################

trial_info_file = os.path.join (save_path , subject + '_trial_info' + txtformat)
# saved trial info file

event_file = os.path.join (save_path , subject + '_eve' + megformat_new)
# event file

raw_lsd_file = os.path.join (save_path , subject + '_lsd' + megformat_new)
# saved LSD raw file

raw_filt_file = os.path.join (save_path , subject + '_filt_' + str (high_pass) + '_' + str (low_pass) + megformat_new)
# saved filtered raw dataset

raw_ica_file = os.path.join (save_path , subject + '_ica' + megformat_new)
# saved ica raw dataset

raw_interpolate_file = os.path.join (save_path , subject + '_interpolate' + megformat_new)
# saved interpolated raw dataset

ER_interpolate_file = os.path.join (save_path , subject + '_ERinterpolate' + megformat_new)
# saved ER interpolated file

epoch_file = os.path.join (save_path , subject + '_epo' + megformat_new)
# saved epoch dataset

epoch_rej =  os.path.join (save_path , subject + '_epochs_rej' + textformat)
# saved text file with rejected epochs
