#!/usr/bin/python3
#
# Download data from an EcoTaxa project
#
# (c) 2022 Jean-Olivier Irisson, GNU General Public License v3

'''
Download data from Ecotaxa projects 7544 (test) and 7545 (learn) and store in appropriate directories.
'''

import os
import yaml

from tqdm import tqdm
import pandas as pd
import pyarrow as pa
from pyarrow import parquet as pq

import ecotaxa_py_client

from ecotaxa_py_client.api import authentification_api
from ecotaxa_py_client.model.login_req import LoginReq

from ecotaxa_py_client.api import objects_api

from ecotaxa_py_client.api import projects_api
from ecotaxa_py_client.model.project_filters import ProjectFilters

from ecotaxa_py_client.api import samples_api
from ecotaxa_py_client.api import acquisitions_api

from ecotaxa_py_client.api import taxonomy_tree_api
from ecotaxa_py_client.model.taxon_model import TaxonModel


# read list of samples to extract
#wanted_sam = pd.read_csv('data/morphocluster_samples_science.csv')
#wanted_sam = wanted_sam.loc[0, 'samples']


# read config
with open(r'config.yaml') as config_file:
    cfg = yaml.safe_load(config_file)
print('### Download data for {}'.format(cfg['dataset']))

# prepare storage
data_dir = os.path.expanduser(cfg['base_dir'])
os.makedirs(data_dir, exist_ok=True)

# authenticate in EcoTaxa
with ecotaxa_py_client.ApiClient() as client:
    api = authentification_api.AuthentificationApi(client)
    token = api.login(LoginReq(
      username=cfg['ecotaxa_user'],
      password=cfg['ecotaxa_pass']
    ))

config = ecotaxa_py_client.Configuration(
    access_token=token, discard_unknown_keys=True)

# get validated objects and their metadata from project
with ecotaxa_py_client.ApiClient(config) as client:
    # list free fields for project
    projects_instance = projects_api.ProjectsApi(client)
    proj = projects_instance.project_query(cfg['proj_id'])
    obj_fields = ['fre.'+k for k in proj['obj_free_cols'].keys()]
    
    # get objects
    objects_instance = objects_api.ObjectsApi(client)
    # only validated
    filters = ProjectFilters()
    #filters = ProjectFilters(statusfilter="V")
    # get taxonomic name and image file name and all free fields
    fields = 'txo.id,txo.display_name,img.file_name,obj.depth_min,obj.latitude,obj.longitude,obj.classif_qual,obj.orig_id' #fre.cnn_score
    fields = fields + ',' + ','.join(obj_fields)
    
    # fetch one object to get the total number of objects to fetch
    objs = objects_instance.get_object_set(cfg['proj_id'], filters,
      fields=fields, window_start=0, window_size=1)

    ## fetch per sample
    samples_instance = samples_api.SamplesApi(client)
    samples = samples_instance.samples_search(
      project_ids=str(cfg['proj_id']),
      id_pattern='%'
    )
    # keep only samples from back transects
    #samples = [sam for sam in samples if sam['orig_id'] in wanted_sam]
    
    # get acquisitions
    acqs_instance = acquisitions_api.AcquisitionsApi(client)
    acqs = acqs_instance.acquisitions_search(
      project_id=str(cfg['proj_id']),
    )
    acquisitions = {
        'acq_id': [],
        'img_name': [],
    }
    for acq in acqs:
        #break
        acquisitions['acq_id'].append(acq['acquisid'])
        acquisitions['img_name'].append(acq['orig_id'])
    
    
    # prepare storage
    objs_dfs = []

    with tqdm(total=objs['total_ids']) as pbar:
        for sam in samples:
        
            # update filters to add sampleid
            filters.samples = str(sam['sampleid'])
            #filters.aquisitions = str(acq['acquisid'])
            
            # fetch a batch of objects
            objs = objects_instance.get_object_set(cfg['proj_id'], filters,
              fields=fields)
            n_fetched_objs = len(objs.details)
            
            # format retrieved data as a DataFrame
            objs_df = pd.DataFrame(objs['details'], columns=fields.split(','))
            # add object id as an identifier
            objs_df['objid'] = objs['object_ids']
            # and sample ID
            objs_df['sample_id'] = sam['orig_id']
            
            # store with the previous batches
            objs_dfs.append(objs_df)
            # and update progress bar
            ok = pbar.update(n_fetched_objs)

# combine all batches in a single DataFrame
df = pd.concat(objs_dfs, ignore_index=True)
# fix id column types
df['txo.id'] = df['txo.id'].astype('int32')
df['objid'] = df['objid'].astype('int32')

# get all unique taxa ids
taxo_ids = list(set(df['txo.id']))
# get lineage for each
with ecotaxa_py_client.ApiClient(config) as client:
    taxo_instance = taxonomy_tree_api.TaxonomyTreeApi(client)
    taxa = [taxo_instance.query_taxa(t) for t in taxo_ids]

lineages = ['/' + '/'.join(t['lineage'][::-1]) for t in taxa]

# add lineages to the DataFrame
taxo = pd.DataFrame({'lineage': lineages}, index=taxo_ids)
df = df.join(taxo, on='txo.id')

## Reformat table
# Keep only relevant columns
#df = df[['obj.orig_id', 'objid', 'img.file_name', 'txo.display_name', 'obj.classif_qual', 'obj.depth_min', 'obj.latitude', 'obj.longitude', 'sample_id', 'fre.area']]

df['datetime'] = df['obj.orig_id']

# Rename columns
df = df.rename(columns={
    'obj.depth_min': 'depth', 
    'obj.orig_id': 'orig_id',
    'obj.latitude': 'latitude', 
    'obj.longitude': 'longitude',
    'txo.display_name': 'taxon', 
    'obj.classif_qual': 'classif_qual', 
    'img.file_name': 'path_to_img',
})

colnames = df.columns
colnames = [x.replace('fre.', '') for x in colnames]
df.columns = colnames

# Reorder a few columns
cols_to_order = ['objid', 'orig_id', 'sample_id', 'depth', 'datetime', 'latitude', 'longitude', 'taxon', 'lineage', 'classif_qual', 'path_to_img']
new_columns = cols_to_order + (df.drop(cols_to_order, axis = 1).columns.tolist())
df = df[new_columns]  

# update image path
#df['path_to_img'] = [os.path.join('data/images', x) for x in df['path_to_img']]


print('  downloaded {} objects. Saving them'.format(df.shape[0]))

# write to disk
pq.write_table(
  pa.Table.from_pandas(df, preserve_index=False),
  where=os.path.join(data_dir, '01.uvp6_objects.parquet'),
  compression='NONE' # for compatibility with R
)

