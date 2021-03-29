import numpy as np
import pandas as pd
import os

# build a dataframe from .npz files
indir = 'detections/all_timepoints_threshold_150'
infiles = [x for x in os.listdir(indir) if x.endswith('.npz')]

# prepare buckets for the data
file_name = []
ch = []
strain = []
time = []
treatment = []
total_vRNAs = []
cytosolic_RNAs = []
repSites = []
mol_per_repSite = []


# parse .npz files
times = ['2h', '6h', '8h', '24h']
for file in infiles:
    for t in times:
        data = np.load(os.path.join(indir, file))
        file_name.append(file)
        total_vRNAs.append(data['rna_coord'].shape[0]+np.sum(data['foci'][:,4]))
        cytosolic_RNAs.append(data['rna_coord'].shape[0])
        repSites.append(data['foci'].shape[0])
        mol_per_repSite.append(np.mean(data['foci'][:,4]))

        if 'MOCK' in file:
            strain.append('None')
            time.append(t.split('h')[0])
            treatment.append('MOCK')
        if 'B117_INF' in file:
            strain.append('B117')
            time.append(t.split('h')[0])
            treatment.append('INF')
        if 'B117_RDV' in file:
            strain.append('B117')
            time.append(t.split('h')[0])
            treatment.append('RDV')
        if 'Vic_INF' in file:
            strain.append('Vic')
            time.append(t.split('h')[0])
            treatment.append('INF')  
        if 'Vic_RDV' in file:
            strain.append('Vic')
            time.append(t.split('h')[0])
            treatment.append('RDV')
        if 'ch3' in file:
            ch.append(3)
        if 'ch4' in file:
            ch.append(4)

data_tall = pd.DataFrame(list(zip(file_name, ch, strain, time, treatment, total_vRNAs,
                            cytosolic_RNAs, repSites, mol_per_repSite)),
            columns =['file_name', 'ch', 'strain', 'time', 'treatment', 'total_vRNAs',
                            'cytosolic_RNAs', 'repSites', 'mol_per_repSite'])
data_tall.to_csv('Data_summary.csv')
print('finished')

