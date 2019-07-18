# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 12:56:37 2019

updated on  Jul 15, anntation changed

@author: gansheng.tan
"""

# # Below are personalised functions -> to py script and import


import warnings

def events_time_forward(events,forward_time):
    ### forward_time in second
    for event in events[0]:
        event[0]=event[0]-forward_time*512
    return events

def removeItem_from_dict(d,key):
    r = dict(d)
    del r[key]
    return r

def dict_getValue(dictionary, search_value):
    for key,value in dictionary.items():
        if value==search_value:
            return key

def mne_annotation_empty(annotations):
    index_dlt = 0
    for i in range(annotations.__len__()):
        annotations.delete(i-index_dlt)
        index_dlt+=1
    return annotations

def mne_annotation_postpone (pptime, annotations):
    onset = []
    duration = []
    description = []
    index_dlt = 0
    for i in range(annotations.__len__()):
        onset.append(annotations.__getitem__(i)['onset']+pptime)
        duration.append(0.0)
        description.append(annotations.__getitem__(i)['description'])
    for i in range(annotations.__len__()):
        annotations.delete(i-index_dlt)
        index_dlt+=1
    annotations.append(onset=onset,duration=duration,description=description)
    print ('annotation time shift succeed')
    
        

def mne_annotation_recode_by_adding(session,state,annotations):
    onset = []
    duration = []
    description = []
    for i in range(annotations.__len__()):
        if annotations.__getitem__(i)['description'] in ['131.0','132.0']:
            onset,duration,description = mne_annotation_recode_info_extract(session=session,state=state,
                                                                        original_annotation = 
                                                                        annotations.__getitem__(i),
                                                                       onset=onset,duration=duration,
                                                                        description=description)
        else:
            continue
    index_dlt = 0
    for i in range(annotations.__len__()):
        if annotations.__getitem__(i-index_dlt)['description'] in ['131.0','132.0']:
            annotations.delete(i-index_dlt)
            index_dlt+=1
        else:
            continue
    onset.append(0.0)
    duration.append(0.0)
    description.append(mne_annotation_add_baseline(session=session,state=state))
    annotations.append(onset=onset,duration=duration,description=description)
    print ('annotation engineering succeed')
    return True

def mne_annotation_add_baseline(session,state):
    if session == '1':
        if state == 'VD':
            return '111.0'
        elif state == 'FA':
            return '211.0'
        elif state == 'OP':
            return '311.0'
        else:
            warnings.warn("unknown state detected", DeprecationWarning)
    elif session == '2':
        if state == 'VD':
            return '112.0'
        elif state == 'FA':
            return '212.0'
        elif state == 'OP':
            return '312.0'
        else:
            warnings.warn("unknown state detected", DeprecationWarning)
    else:
        warnings.warn("add baseline function only apply on rawfile having 2 sessions", DeprecationWarning)
    return '999.0'
        

def mne_annotation_recode_info_extract(session,state,original_annotation,onset,duration,description):
    if session =='1':
        if state == 'VD':
            if original_annotation['description']=='131.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('121.0')
                
            elif original_annotation['description']=='132.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('131.0')
            else:
                print('this function only detect safe and threat period, please check original annotations')
        elif state == 'FA':
            if original_annotation['description']=='131.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('221.0')

            elif original_annotation['description']=='132.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('231.0')
            else:
                print('this function only detect safe and threat period, please check original annotations')
        elif state == 'OP':
            if original_annotation['description']=='131.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('321.0')
            elif original_annotation['description']=='132.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('331.0')
            else:
                print('this function only detect VD, FA, OP states, please check original annotations')
    elif session =='2':
        if state == 'VD':
            if original_annotation['description']=='131.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('122.0')
            elif original_annotation['description']=='132.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('132.0')
            else:
                print('this function only detect safe and threat period, please check original annotations')
        elif state == 'FA':
            if original_annotation['description']=='131.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('222.0')
            elif original_annotation['description']=='132.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('232.0')
            else:
                print('this function only detect safe and threat period, please check original annotations')
        elif state == 'OP':
            if original_annotation['description']=='131.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('322.0')
            elif original_annotation['description']=='132.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('332.0')
            else:
                print('this function only detect VD, FA, OP states, please check original annotations')
    else:
        print('3rd session dected, please check annotations')
    return(onset,duration,description)
        
