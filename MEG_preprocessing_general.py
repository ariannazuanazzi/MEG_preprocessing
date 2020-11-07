#!/usr/bin/env python 3.7.6
# -*- coding: utf-8 -*-
"""
general MEG_processing Pipeline: adapted for singletone or two tones.
calls the functions: general_info_general and base_funcs_general
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
__date__ = '2020/11/07'
__license__ = "MIT"
__version__ = "1.2.0"
__email__ = "az1864@nyu.edu"


#########################################################
# IMPORT NEEDED MODULES----------------------------------------------------


import matplotlib
matplotlib.use('TkAgg')
import autoreject
import pickle
import eelbrain
from eelbrain.gui import select_epochs
eelbrain.configure(prompt_toolkit=False)
# disable running the GUI main loop in parallel to the Terminal

from general_info_general import *
from base_funcs_general import *


#########################################################
# READ RAW -------------------------------------------------------
raw = mne.io.read_raw_fif (raw_alldata_fileName , preload = True)
# read fif image
print(raw.info)
# print raw info
raw.plot(block=True, n_channels=40, title='RawData', scalings = scaling)

# check for STIM and eyedata channels
raw.set_channel_types ({ eog_channels : 'eog' })
# set eye tracking channel as eog channel

if view_plot:
    plot_movement (raw, save_plot, figure_path)
    # plot head movement in coordinates x,y,z

if save_plot:
    matplotlib.pyplot.savefig(f'{figure_path}headmovement.png',dpi=dpi)

if view_plot:
    raw.plot_psd (tmax = np.inf , fmax = 250 , average = False)
    # psd plots

# SAVE TRIAL INFORMATION -----------------------------------------------------------------------------------------

if os.path.isfile (trial_info_file) :
    trial_info = pickle.load (open (trial_info_file , 'rb'))
# load trial info

else :
    trial_info = { }
# save every process into a dictionary

#########################################################
# FIND EVENTS -----------------------------------------------------------------------------------------------------

if os.path.isfile (event_file) :
    events = mne.read_events (event_file)
    # read events from event file

else :

    events = mne.find_events(raw, stim_channel = stim_channel, shortest_event = 1, min_duration = 0.003) #min_duration allows to define an event as one and not two (if trigger is of a different duration)
    # find events based on stim channel or on trigger
    #events = events[events[:,1]==1]  #if necessary to filter events cause trigger box sent multiple triggers

    trial_info [ 'events' ] = events
    trial_info [ 'event_id' ] = event_id
    pickle.dump (trial_info , open (trial_info_file , 'wb'))
    # save events into trial info

    mne.write_events (event_file , events)
    # save events to event fif file

    # plot events
    infoEvents = mne.viz.plot_events(events, event_id=event_id, sfreq=raw.info['sfreq'],
                          first_samp=raw.first_samp, show=False)

    if save_plot:
        infoEvents.savefig(f'{figure_path}infoEvents.png', dpi=dpi)
    if view_plot:
        plt.show()


    # ensure that number of events from the triggers is the same as the number of trials in the block*6triggers per block
    # i.e., event1, event2, event3, probe, response, feedback + block number
    # throws an assertion error if they do not
    if blockonset_recorded:
        # if the blockn does not exist, it is = to 0
        trig_blocksn =  sum(events[:,2]==event_id['blockOnset'])
    else:
        trig_blocksn = 0
    # number of block's triggers in the block (should be 1, but sometimes it is 2 - if 2 blocks together - or 3 - because the 3rd block started and stopped)
    #assert events.shape[0] == log['trigger'].values[log_sel].shape[0]*6+trig_blocksn, "trigger numbers and conditions in log file do not correspond!"


#########################################################
# REMOVE BAD CHANNELS (Autoreject: Reject a trial only if most sensors “agree” that the trial is bad, otherwise interpolate as many sensors as possible)-------

if os.path.isfile (raw_interpolate_file) :
    raw_inter = mne.io.read_raw_fif (raw_interpolate_file , preload = True)
    # read interpolated data

else :
    raw_inter = raw.copy()
    # load raw data for rejection and interpolation

    raw_ER = mne.io.read_raw_fif (ER_file_converted , preload = True)
    # load empty room data to update interpolated channels there as well

    raw_inter.pick_types (meg = True)
    raw_ER.pick_types (meg = True)
    # pick channel type: only MEG otherwise signal of other channel (eg sound) is used for interpol

    epoch_for_autoreject = mne.Epochs (raw_inter , events , event_id = event_id , tmin = tmin , tmax = tmax ,
                                        proj = True , baseline = baseline , preload = True ,
                                        reject = None, reject_by_annotation=False, flat=None)
    # create epoch for automatic bad channel detection

    epoch_for_autoreject.average()
    #epoch_for_autoreject.average().plot_joint();

    # #Ransac
    # ransac = autoreject.Ransac (verbose = 'progressbar' , min_corr=0.75, n_jobs = 1 , unbroken_time = 0.4)
    # # Ransac algorithm to autodetect bad channels before deleting user defined bad channels (arXiv:1612.08194)
    # ransac.fit (epoch_for_autoreject)
    # # fit epoch data
    # print (f'Bad channels detected by Ransac: {ransac.bad_chs_}, number = {len (ransac.bad_chs_)}')
    # #get bad channel detected
    # bad_channels = list (set (badchannelsUI + ransac.bad_chs_))
    # # bad channels = user predetermined set + autoreject detected set

    #Alternative BadChannels identification
    badchannelsRG = find_bad_channels_epochs(epoch_for_autoreject)
    # algorithm to autodetect bad channels based on range (LGwilliams)
    #bad_channels = list (set (badchannelsUI + badchannelsRG))
    bad_channels = list (set (badchannelsRG))
    # bad channels = user predetermined set + autoreject detected set

    print (f'Bad channels in total: {bad_channels}, number = {len (bad_channels)}')
    # check bad channels

    raw_inter.info [ 'bads' ] = bad_channels
    raw_ER.info [ 'bads' ] = bad_channels
    # set bad channels

    raw_inter.plot (block = True, duration = 100, n_channels = 40, title = 'Raw_selectBads',  scalings = scaling)
    # visualise bad channels, click the bad channel to make new / cancel selection (if ERP, the more bad channels are deleted the better)

    trial_info [ 'bads' ] = raw_inter.info [ 'bads' ]
    # save the bad channel info
    raw_ER.info [ 'bads' ] = raw_inter.info [ 'bads' ]
    # set bad channels consistent in ER as well

    print (f'Bad channels final: {trial_info [ "bads" ]}, number = {len (trial_info [ "bads" ])}')
    # show the final decision of bad channels

    badsensorsp = raw_inter.plot_sensors (ch_type = 'mag' , show_names = 1, show=False)
    # show the bad channel in 2D topomap

    if save_plot :
        badsensorsp.savefig (f'{figure_path}badsensors.png' , dpi = dpi)
    if view_plot:
        plt.show()

    pickle.dump (trial_info , open (trial_info_file , 'wb'))
    # save channels to interpolate

    raw_inter.interpolate_bads (reset_bads = True)
    raw_ER.interpolate_bads (reset_bads = True)
    # interpolate bads (MNE function interpolate_bads) and reset so that we have same number of channels for all blocks/subjects

    raw_inter.save (raw_interpolate_file , overwrite = True)
    raw_ER.save (ER_interpolate_file , overwrite = True)
    # overwrite w/  bad channel info/interpolated bads

if view_plot:
    raw_inter.plot (block = True, duration = 100, n_channels = 40, title = 'Raw_interpolated',  scalings = scaling)
    # plot new list of channels after bad channels have been interpolated


#########################################################
# GET AND SAVE HEARTBEAT PLOT --------------------------------------------------------------------------------------

ecg_epochs = mne.preprocessing.create_ecg_epochs (raw_inter)
avg_ecg_epochs = ecg_epochs.average()
# average across epochs
heartbeatp = avg_ecg_epochs.plot_topomap (times = np.linspace (-0.05 , 0.05 , 11), title='RawData_heartbeat', layout=KIT_layout, show = False)
#avg_ecg_epochs.plot_joint (times=[-0.25, -0.025, 0, 0.025, 0.25], title='RawData_heartbeat',  show = False)
# topos
if save_plot:
    heartbeatp.savefig (f'{figure_path}heartbeat.png',dpi=dpi)
if view_plot:
    plt.show()

# GET AND SAVE SACCADE AND BLINK FROM EYELINK DATA ---------------------------------------------------------------------
if eye_link_events :
    sac_channel , blk_channel = get_eyelink (raw , events , eyelink_file , eog_channels , high_pass, low_pass, save_plot, figure_path)
else :
    raw_eog = raw.copy()
    raw_eog._data[0:157, :] = raw_inter._data # put excluded channels (included EOG) back
    eog_channel = raw_eog.copy().pick_types (meg = False , eog = True)
    eog_channel.info [ 'highpass' ] = high_pass
    eog_channel.info [ 'lowpass' ] = low_pass


#########################################################
# EMPTY ROOM DENOISE with LSD ------------------------------------------------------------------------------------

if os.path.isfile (raw_lsd_file) :
    raw_lsd = mne.io.read_raw_fif (raw_lsd_file , preload = True)
    # read denoised data

else :
    raw_lsd = least_square_reference (raw_interpolate_file , ER_interpolate_file)
    # apply LSD

    raw_lsd.pick_types (meg = True , misc = False)
    # get meg channel only since information in MISC channel is lost after LSD

    raw_lsd.save (raw_lsd_file , overwrite = True)
    # save denoised raw data

if view_plot:
    raw_lsd.plot (block = True, duration = 100, n_channels = 40, title = 'Raw_denoised',  scalings = scaling)
    # plot new list of channels after LS denoise


#########################################################
# FILTER -----------------------------------------------------------------------------------------

if os.path.isfile (raw_filt_file) :
    raw_filtered = mne.io.read_raw_fif (raw_filt_file , preload = True)
    # read pre-filtered data
else :

    if filter_type == 'bandpass' :
        raw_filtered = raw_lsd.filter (l_freq = high_pass , h_freq = low_pass , n_jobs = 2)
    # apply band pass filter

    else :
        # apply high & low pass filter in order
        raw_filtered = raw_lsd.filter (l_freq = high_pass , h_freq = None , n_jobs = 2)
        # highpass filter

        raw_filtered = raw_lsd.filter (l_freq = None , h_freq = low_pass , n_jobs = 2)
        # lowpass filter

    if eye_link_events :
        raw_filtered.add_channels ([ sac_channel , blk_channel ])
    else :
        raw_filtered.add_channels ([ eog_channel ])
    # add back eye tracking channel for ICA stage

# plot psd before&after filter to check
psd_before_filter = raw_lsd.plot_psd (area_mode = 'range' , average = False , show = False)
psd_after_filter = raw_filtered.plot_psd (area_mode = 'range' , average = False , fmax = 60. , show = False)

if save_plot :
    psd_before_filter.savefig(f'{figure_path}psd_before_filter.png', dpi=dpi)
    psd_after_filter.savefig (f'{figure_path}psd_after_filter.png' , dpi = dpi)
if view_plot:
    plt.show()


trial_info [ 'filter' ] = [ high_pass , low_pass , filter_type ]
pickle.dump (trial_info , open (trial_info_file , 'wb'))
# save filter

raw_filtered.save (raw_filt_file , overwrite = True)
# save filtered raw dataset

if view_plot:
    raw_filtered.plot (block = True, duration = 100, n_channels = 40, title = 'Raw_filtered', scalings = scaling)
    # plot new list of channels after filtering


#########################################################
# ICA -----------------------------------------------------------------------------------------
if os.path.isfile (raw_ica_file) :
    raw_ica = mne.io.read_raw_fif (raw_ica_file , preload = True)
# read pre-ica data
else :
    ica = mne.preprocessing.ICA (n_components = ica_components , method = ica_method , random_state = 0)
    # generate ica

    ica.fit (raw_filtered , picks = 'meg', decim = ica_decim)
    # fit ica to dataset

    eog_inds , eog_scores = ica.find_bads_eog (raw_filtered)
    eogp = ica.plot_scores (eog_scores , exclude = eog_inds , title = 'eog' , labels = 'eog', show=False)
    #detect EOG by correlation: automatically find eog component & plot the scores

    if save_plot :
        eogp.savefig (f'{figure_path}ica_correlation_eog.png' , dpi = dpi)
    if view_plot:
        plt.show()

    print (f'The found EOG component is:{eog_inds}')
    # show the found eog component

    ecg_inds, ecg_scores = ica.find_bads_ecg(raw_filtered)
    ecgp = ica.plot_scores (ecg_scores , exclude = ecg_inds , title = 'ecg' , labels = 'ecg', show=False)
    # detect EOG by correlation: find which ICs match the ECG pattern

    if save_plot :
        ecgp.savefig (f'{figure_path}ica_correlation_ecg.png' , dpi = dpi)
    if view_plot:
        plt.show()

    print (f'The found ECG component is:{ecg_inds}')
    # show the found eog component

    raw_filtered.pick_types (meg = True , eog = False)
    # keep only meg channels

    while True :
        # loop for ICA component selection.
        ica_components_plot = ica.plot_components (layout = KIT_layout, show=False)
        # plot all components
        if view_plot:
            plt.show()

        ica.exclude = [ int (i) for i in call_msg_box().split (',') ]
        # exclude components

        print(ica.exclude)

        icap = ica.plot_sources (raw_filtered, block=True, show = False)
        # plot continuous ICA component data


        if save_plot :
            icap.savefig (f'{figure_path}ica_timeseries.png' , dpi = dpi)
        if view_plot:
            plt.show()

        ica_property = ica.plot_properties (raw_filtered , picks = ica.exclude, topomap_args = {'layout' : KIT_layout }, show = False) #
        # plot all components with property, take couple minutes, use with patience
        if save_plot :
            for ind , fig in enumerate (ica_property): fig.savefig (f'{figure_path}ica{ind}_property.png', dpi = dpi)
        if view_plot:
            plt.show()


        #ica.exclude = [ int (i) for i in call_msg_box().split (',') ]
        # exclude components

        trial_info [ 'ica' ] = ica.exclude
        pickle.dump (trial_info , open (trial_info_file , 'wb'))
        # save ica component into trial info

        ica_overlay = ica.plot_overlay (raw_filtered, show=False)
        # check ica result before&after artifacts rejection
        if save_plot :
            ica_overlay.savefig (f'{figure_path}ica_overlay.png' , dpi = dpi)
            # save ica plots
        if view_plot:
            plt.show()

        check = call_msg_box ('CLEAR?')
        # check final call

        if check == '1' :
            break

        else :
            continue
        # iterative ICA selection put 1 to end

    raw_ica = raw_filtered.copy()

    raw_ica.pick_types (meg = True , eog = False)
    # get rid of eog channel

    ica.apply (raw_ica)
    # apply ica to raw data

    raw_ica.save (raw_ica_file , overwrite = True)
    # save raw dataset after ICA

if view_plot:
    raw_ica.plot (block = True, duration = 100, n_channels = 40, title = 'Raw_ICAd',  scalings = scaling)
    # plot new list of channels after ICA


#########################################################
# EPOCH DATA -------------------------------------------------------------------------
if os.path.isfile (epoch_file) :
    epoch_def = mne.read_epochs (epoch_file)
# read pre-existed epoch data

else :
    epoch_def = mne.Epochs (raw_ica , events , event_id = event_id , tmin = tmin , tmax = tmax ,
                         proj = True , baseline = baseline , preload = True ,
                         reject = None, reject_by_annotation=False, flat=None)


    picks = mne.pick_types (epoch_def.info , meg = True , eog = False)
    # pick channel types

    # ar = autoreject.AutoReject (picks = picks , n_jobs = 4)
    # # generate autoreject object
    #
    # reject = autoreject.get_rejection_threshold (epoch_def)
    # # get autoreject threshold
    #
    # print (reject)
    #
    # epoch_def.drop_bad (reject = reject)
    # # reject according to the threshold
    # epochdrop = epoch_def.plot_drop_log(show=False)
    # if save_plot :
    #     epochdrop.savefig(f'{figure_path}epochs_reasondrop.png', dpi=dpi)
    # if view_plot:
    #     plt.show()

    epochsp = epoch_def.plot (block = True, n_channels = 40 , title = 'Epochs', events=events, show=False)
    # manual drop bad epoch by clicking
    if view_plot:
        plt.show()

    trial_info [ 'drop_log' ] = [ i for i in epoch_def.drop_log if i != [ 'IGNORED' ] and i != [ ] ]
    pickle.dump (trial_info , open (trial_info_file , 'wb'))
    # save reject info into trial info

    epoch_def.save (epoch_file , overwrite = True)
    # save epoch dataset

    if save_plot :
        epochsp.savefig(f'{figure_path}epochs.png', dpi=dpi)

    epochdrop_overlap = select_epochs(epoch_def, data='meg', accept='accept',
                                      vlim=3e-12, path=epoch_rej)


#########################################################
# PLOT ERP DATA ------------------------------------------------------------------------------

# if concat_blocks==True:
    # plots each condition butterfly and M1 topography
    for condition in event_id.keys() : #['probe_trig5', 'probe_trig6']:# ['event1', 'event2', 'event3']:
        evoked = epoch_def [ condition ].average()
        # average each condition epoch

        # selection = mne.read_selection ('temporal' , fname = KIT_selection , info = epoch.info)
        # select temporal lobe to pick up auditory cortex response

        butterfly = evoked.plot ( gfp = True , time_unit = 'ms' ,
                                  window_title = f'Condition:{condition}', show=False) #picks = selection
        # plot butterfly&GFP

        topomap = evoked.plot_topomap ( time_unit = 'ms' , ch_type = 'mag' ,
                                        average = 0.025 , colorbar = True, layout = KIT_layout, show=False) #times = [ 0 , 0.05 , 0.1 , 0.17 ], layout = KIT_layout
        # plot topography

        joint = evoked.plot_joint (topomap_args = { 'layout' : KIT_layout }, show=False) #times = [ 0 , 0.08 , 0.15 ], topomap_args = { 'layout' : KIT_layout })
        # plot joint graph

        if save_plot:
            butterfly.savefig(f'{evokedfigure_path}{condition}_butterfly.png',dpi=dpi)
            topomap.savefig(f'{evokedfigure_path}{condition}_topomap.png',dpi=dpi)
            joint.savefig(f'{evokedfigure_path}{condition}_joint.png',dpi=dpi)
            #save evoked plots
        # if view_plot:
        #     plt.show()
