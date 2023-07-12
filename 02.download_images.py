#!/usr/bin/python3
#
# Download images from EcoTaxa
#
# (c) 2022 Jean-Olivier Irisson, GNU General Public License v3

'''
Download images from Ecotaxa projects 7544 (test) and 7545 (learn) and store in appropriate directories.
'''

import os
import yaml
import pandas as pd
import pyarrow as pa
from pyarrow import parquet as pq

import shutil
import urllib.request

from tqdm import tqdm
from pyarrow.parquet import read_table

# read config
with open(r'config.yaml') as config_file:
    cfg = yaml.safe_load(config_file)
print('### Download images for {}'.format(cfg['dataset']))

# prepare storage
data_dir = os.path.expanduser(cfg['base_dir'])
img_dir = os.path.join(data_dir, 'images')
os.makedirs(img_dir, exist_ok=True)

# read data from EcoTaxa
df = read_table(os.path.join(data_dir, 'uvp6_objects.parquet')).to_pandas()

# detect image extension
file, ext = os.path.splitext(df['path_to_img'][0])
# TODO extract file extension for each file, like in step 4.
# name image according to ecotaxa image paths to avoid having all images in one folder
# do not use taxonomy subfolders as these are going to change with every classification iteration
#df['dest_path'] = [os.path.join(img_dir, str(this_id)+ext) for this_id in df['objid']]
df['dest_path'] = [os.path.join(img_dir, x) for x in df['path_to_img']]


# download images (with a nice progress bar)
vault_path = '/remote/ecotaxa/vault'
for i in tqdm(range(df.shape[0])):
    # if the file has not been copied already
    if not os.path.isfile(df['dest_path'][i]):
        # create directory
        os.makedirs(os.path.dirname(df['dest_path'][i]), exist_ok=True)
        # copy from vault
        if os.path.isdir(vault_path):
            res = shutil.copyfile(
                src=os.path.join(vault_path, df['path_to_img'][i]),
                dst=df['dest_path'][i]
            )
        # or copy through the internet
        else:
            res = urllib.request.urlretrieve(
                url='https://ecotaxa.obs-vlfr.fr/vault/'+df['path_to_img'][i],
                filename=df['dest_path'][i]
            )
                
n_imgs = len(os.listdir(img_dir))
print('  {} images in {}'.format(n_imgs, img_dir))


## write to disk
#df = df.drop('dest_path', axis = 1)
#pq.write_table(
#  pa.Table.from_pandas(df),
#  where=os.path.join(data_dir, 'all_copepods.parquet'),
#  compression='NONE' # for compatibility with R
#)
