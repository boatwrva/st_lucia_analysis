# -*- coding: utf-8 -*-
"""
Spyder Editor
Lin Hou
Extracting Metadata from Cruise info matching ifcb files
This is a temporary script file.
"""
import numpy as np
import pandas as pd
import os
from os import listdir
from os.path import isfile, join

def compile_meta(ifcbfiles,metadata,updates):
    if ifcbfiles==[]:
        return metadata,updates
    else:
        file=ifcbfiles[0]
        filename=file[:-4]
        if metadata[metadata['filename'].isin([filename])].empty:
            filenum=int(input(f'how many samples are with {filename}:' ))
            La=float(input('input the Latitude:'))
            Lon=float(input('input the Longitude:'))
            depth=input('Listing depths in order (i.e 1,2,3):')
            depth=[int(s) for s in depth.split(',')]
            sampletype=input('Listing Sample types in order (i.e Surface,DCM):')
            sampletype=[s for s in sampletype.split(',')]
            niskin=input('Listing the niskin in order:')
            niskin=[int(s) for s in niskin.split(',')]
            subset_ifcbfiles=ifcbfiles[0:filenum]
            for i in range(len(subset_ifcbfiles)):
                file=subset_ifcbfiles[i]
                df2 = pd.DataFrame({'filename': file[:-4], 'Latitude': [La], 'Longitude':[Lon],'Depth':[depth[i]],
                                    'Cruise':'CC2404','Sampletype':sampletype[i],'Niskin':[niskin[i]]})
                metadata = pd.concat([metadata,df2], ignore_index = True) 
                updates+=1
            return compile_meta(ifcbfiles[filenum:],metadata,updates)
        else:
            return compile_meta(ifcbfiles[1:],metadata,updates)
def getting_meta(start):
    if start=='y':
        shipinfo=('')# enter the ship metadata path here
        # shipinfo='/Users/allenlab/Documents/ifcbdb/ifcbdb/ifcb_data/P2402_underway/P2402_Met/' # address of the ship info, connected to the server and shared folder
        dataname=input('the experiment name calcofi2404_underway/discrete:')
        ifcbfolder=f'/Users/allenlab/Documents/ifcbdb/ifcbdb/ifcb_data/calcofi2404_{dataname}/'
        if os.path.isdir(ifcbfolder):
            print('folder path:' + ifcbfolder)
        else:
            print('Error: folder no exist')            
        metadatapath=f'/Users/allenlab/Documents/ifcbdb/ifcbdb/ifcb_data/calcofi2404_{dataname}/metadata.csv'
        if not os.path.exists(metadatapath):
            if dataname=='underway':
                metadata=pd.DataFrame({'filename':[],'Latitude':[],'Longitude':[]})
            else:
                metadata=pd.DataFrame({'filename':[],'Latitude':[],'Longitude':[],'Depth':[],'Cruise':[],'Sampletype':[],'Niskin':[]})
        else:
            metadata=pd.read_csv(metadatapath)

        ifcbfiles = [f for f in listdir(ifcbfolder) if isfile(join(ifcbfolder, f)) and '.adc' in f]
        update,lasttime=0,0
        if dataname=='underway':
            for file in ifcbfiles:
                filename=file[:-4]
                if metadata[metadata['filename'].isin([filename])].empty:
                    time=int(file[10:16]) # convert the string to interger
                    day=file[3:9]
                    shipdata=np.loadtxt(shipinfo+day+'.MET',usecols=(0,61,62))
                    if shipdata[1,2]>0 or shipdata[1,2]<-100:
                        shipdata=np.loadtxt(shipinfo+day+'.MET',usecols=(0,60,61))
                    metatime=min(shipdata[:,0], key=lambda x:abs(x-time))
                    if metatime!=lasttime:
                        lasttime=metatime
                        subset=shipdata[shipdata[:,0]==metatime,:]
                        La,Lon=subset[0,1:3]
                        df2 = pd.DataFrame({'filename': filename, 'Latitude': [La], 'Longitude':[Lon]})
                        metadata = pd.concat([metadata,df2], ignore_index = True) 
                        update+=1
                    else:
                        print('no latest metadata, please update the CalCOFI_MET file')
        else:
           metadata,updates=compile_meta(ifcbfiles,metadata,0)
           update+=updates
        metadata.to_csv(metadatapath, index=False)
        print('updates '+str(update)+' files, total '+str(len(ifcbfiles))+' files')
        print(f'Save to {metadatapath}')
        start=input('Do you want to do another dataset?[y/n]:')
        return getting_meta(start)
    else:
        print('End getting_metadata.py')

getting_meta('y')