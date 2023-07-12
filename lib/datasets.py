import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler, MultiLabelBinarizer
from tensorflow.keras import utils
import lycon
import random
import matplotlib.pyplot as plt
import imgaug as ia
from imgaug import augmenters as iaa
import tarfile
import cv2 


def read_data_cnn(path, frac=1, random_state=None):
    """
    Read a csv file containing data to train the cnn
    
    Args:
        path (str): path to the file
        frac (float, int): fraction of dataset to use
        random_state (int or RandomState): controls the randomness of shuffling and augmentation; default=None
        valid_plankton_prop (float, int, None): desired proportion of plankton objects in validation set. 'None' does not subsample non plankton objects.
    
    Returns:
        df_train (DataFrame): training data containing path to image and classif_id
        df_valid (DataFrame): validation data containing path to image and classif_id
        df_calib (DataFrame): calibration data (all validation data) containing path to image and classif_id
        df_test (DataFrame): testing data containing path to image and classif_id
        df_classes (DataFrame): classes with their plankton attribute
        df_comp (DataFrame): dataset composition
    """
    
    # Read CSV file
    df = pd.read_csv(path).rename(columns = {'taxon':'classif_id'})
    
    # Extract classes ('classif_id_1') and plankton attribute
    df_classes = df[['classif_id','plankton']].drop_duplicates().sort_values('classif_id').reset_index(drop=True)
    
    # The classifier is a CNN, keep 'classif_id_1', 'path_to_img' and 'set' split
    df = df[['path_to_img', 'classif_id', 'set', 'plankton']]
    
    # Fraction subsample 
    if frac < 1:
        df = df.groupby(['classif_id', 'set', 'plankton'], group_keys=False).apply(lambda x: x.sample(frac=frac, random_state=random_state)).reset_index(drop=True)
        
    
    # Extract training, validation and test splits
    df_train = df[df['set'] == 'train'].drop(['set', 'plankton'], axis = 1).reset_index(drop=True)
    df_valid = df[df['set'] == 'val'].drop(['set', 'plankton'], axis = 1).reset_index(drop=True)
    
    # Compute dataset composition
    df_comp = df.groupby(['classif_id','set']).size().unstack(fill_value=0)

    return df_train, df_valid, df_classes, df_comp


## Define a data generator 
class DataGenerator(utils.Sequence):
    """
    Generate batches of data for CNN.
    
    Args:
        df (DataFrame): dataframe with path to images and classif_id
        classes (list, array): name of classes
        data_dir (str): directory containing data
        batch_size (int): number of images per batch
        image_dimensions (tuple): images dimensions for CNN input
        shuffle (bool): whether to shuffle data
        augment (bool): whether to augment (zoom in/out, flip, shear) data
        px_del (int): number of pixels to delete at bottom of images (e.g. to remove a scale bar)
        preserve_size (bool): whether to preserve size of small images.
            If False, all image are rescaled. If True, large images are rescaled and small images are padded to CNN input dimension.
        random_state (int or RandomState): controls the randomness
    
    Returns:
        A batch of `batch_size` images (4D ndarray) and one-hot encoded labels (2D ndarray)
    """
    
    def __init__(self, df, classes, data_dir, batch_size=32, image_dimensions = (224, 224, 3), shuffle=True, augment=False, px_del = 0, preserve_size=False, random_state=None):
        self.df               = df  
        self.classes          = classes
        self.data_dir         = data_dir            
        self.batch_size       = batch_size          
        self.image_dimensions = image_dimensions    
        self.shuffle          = shuffle             
        self.augment          = augment             
        self.px_del           = px_del               
        self.preserve_size    = preserve_size   
        self.random_state     = random_state
        
        # initialise the one-hot encoder
        mlb = MultiLabelBinarizer(classes=classes)
        self.class_encoder = mlb
        
        self.on_epoch_end()

    def __len__(self):
        'Denotes the number of batches per epoch'
        return int(np.ceil(len(self.df) / self.batch_size))

    def on_epoch_end(self):
        'Updates indexes after each epoch'
        self.indexes = np.arange(len(self.df))
        # shuffle data if chosen
        if self.shuffle:
            self.df = self.df.sample(frac=1, random_state=self.random_state).reset_index(drop=True)
             
    def get_padding_value(self, img):
        'Compute value to use to pad an image, as the median value of border pixels'
        # get height and width of image
        h = img.shape[0]
        w = img.shape[1]
        
        # concatenate border pixels in an array
        borders = np.concatenate((
            img[:, 0],         # left column
            img[:, w-1],       # right column
            img[0, 1:w-2],     # top line without corners
            img[h-1, 1:w-2],   # bottom line without corners        
        ), axis=0)
        
        # compute the median
        pad_value = np.median(borders)
        
        return pad_value
    
    def augmenter(self, images):
        'Define a data augmenter which doses horizontalf flip (50% chance), vertical flip (50% chance), zoom and shear'
        
        if self.random_state is not None:
            ia.seed(self.random_state) # set seed for randomness
        
        seq = iaa.Sequential(
            [
                iaa.Fliplr(0.5),  # horizontally flip 50% of all images
                iaa.Flipud(0.5),  # vertically flip 20% of all images
                iaa.Affine(
                    scale={"x": (0.8, 1.2), "y": (0.8, 1.2)}, # scale images to 80-120% of their size, individually per axis
                    shear=(-15, 15),  # shear by -15 to +15 degrees
                    mode='edge', # pad images with border picels
                ),
            ],
            random_order=True # apply these transformations in random order
        )
        return seq(images=images)

    def __getitem__(self, index):
        'Generate one batch of data'
        # selects indices of data for next batch
        indexes = self.indexes[index * self.batch_size : (index + 1) * self.batch_size]

        # select data and load images
        paths = [os.path.join(self.data_dir, self.df.path_to_img[k]) for k in indexes]
        images = [lycon.load(p)/255 for p in paths]
        
        batch_prepared_images = []
        output_size = self.image_dimensions[0]
        # resize images to proper dimension
        for img in images:
            h,w = img.shape[0:2]
            
            # delete scale bar of px_del px at bottom of image
            img = img[0:h-self.px_del,:]
            h = img.shape[0]
            
            # compute largest dimension (hor or ver)
            dim_max = int(max(h, w))
            
            # if size is not preserved or image is larger than output_size, resize image to output_size
            if not(self.preserve_size) or (dim_max > output_size):
                # Resize image so that largest dim is now equal to output_size
                img = lycon.resize(
                    img, 
                    height = max(h*output_size//max(h,w),1), 
                    width = max(w*output_size//max(h,w),1), 
                    interpolation=lycon.Interpolation.AREA
                )
                h = img.shape[0]
                w = img.shape[1]  
            
            # create a square, empty output, of desired dimension, filled with padding value
            pad_value = self.get_padding_value(img)
            img_square = np.full(self.image_dimensions, pad_value)
            
            # compute number of pixels to leave blank 
            offset_ver = int((output_size-h)/2) # on top and bottom of image
            offset_hor = int((output_size-w)/2) # on left and right of image
            
            # replace pixels in output by input image
            img_square[offset_ver:offset_ver+img.shape[0], offset_hor:offset_hor+img.shape[1]] = img
            batch_prepared_images.append(img_square)
        
        # convert to array of images        
        batch_prepared_images = np.array([img for img in batch_prepared_images], dtype='float32')
        
        # data augmentation
        if self.augment == True:
            batch_prepared_images = self.augmenter(batch_prepared_images)
        
        ## Labels
        # extract the labels corresponding to the selected indexes, when they are provided
        if self.classes is not None:
            batch_labels = [self.df.classif_id[i] for i in indexes]
            batch_encoded_labels = self.class_encoder.fit_transform([[l] for l in batch_labels])
        else :
            batch_encoded_labels = None
            


        
        # Return reshaped images with labels
        return batch_prepared_images, batch_encoded_labels

    
## Create a data generator
class PredictionDataGenerator(utils.Sequence):
    """
    Generate batches of data from tarfile for CNN prediction.
    
    Args:
        tarname (str): name of particles tarfile 
        df (DataFrame): dataframe with path to images and classif_id
        classes (list, array): name of classes
        data_dir (str): directory containing data
        batch_size (int): number of images per batch
        image_dimensions (tuple): images dimensions for CNN input
        shuffle (bool): whether to shuffle data
        augment (bool): whether to augment (zoom in/out, flip, shear) data
        px_del (int): number of pixels to delete at bottom of images (e.g. to remove a scale bar)
        preserve_size (bool): whether to preserve size of small images.
            If False, all image are rescaled. If True, large images are rescaled and small images are padded to CNN input dimension.
        random_state (int or RandomState): controls the randomness
    
    Returns:
        A batch of `batch_size` images (4D ndarray) and one-hot encoded labels (2D ndarray)
    """
    def __init__(self, tarname, batch_size=256, image_dimensions = (224, 224, 3), shuffle=False, augment=False, px_del = 0, preserve_size=False):
        self.tar              = tarfile.open(tarname, mode='r')# tar file
        self.batch_size       = batch_size                     # batch size
        self.image_dimensions = image_dimensions               # image_dimensions
        self.shuffle          = shuffle                        # shuffle bool
        self.augment          = augment                        # augment bool
        self.px_del           = px_del                         # pixels to delete at bottom of images  
        self.preserve_size    = preserve_size                  # preserve_size bool    
        self.images_paths     = [x for x in self.tar.getnames() if 'png' in x] # paths for png files if tar file
        self.on_epoch_end()
        
    def __len__(self):
        'Denotes the number of batches per epoch'
        return int(np.ceil(len(self.images_paths) / self.batch_size))

    def on_epoch_end(self):
        'Updates indexes after each epoch'
        self.indexes = np.arange(len(self.images_paths))
        if self.shuffle:
            np.random.shuffle(self.indexes)
    
    def get_image_from_tar(self, name):
        'Gets a image as numpy array from file tar name'
        
        image = self.tar.extractfile(name)
        image = image.read()
        image = cv2.imdecode(
            np.asarray(
                bytearray(
                    image
                ), dtype=np.uint8
            ), 0
        )/255

        return image
    
    def __getitem__(self, index):
        'Generate one batch of data'
        # selects indices of data for next batch
        indexes = self.indexes[index * self.batch_size : (index + 1) * self.batch_size]
        
        #list_files            = self.tar.getnames()            # paths in tar file
        #self.images_paths     = [x for x in list_files if 'png' in x] # paths for png files
        
        #self.images_paths     = [x for x in self.tar.getnames() if 'png' in x] # paths for png files
        
        # select data and load images
        #images = [cv2.imread(self.images_paths[k])/255 for k in indexes]
        images = [self.get_image_from_tar(self.images_paths[k]) for k in indexes]
        
        square_images = []
        output_size = self.image_dimensions[0]
        # resize images to proper dimension
        for img in images:
            # Create a thrid axis
            #img = np.stack([img, img, img], axis=2)
            
            h = img.shape[0]
            w = img.shape[1] 
            
            # delete scale bar of 31px at bottom of image
            img = img[0:h-self.px_del,:]
            h = img.shape[0]
            
            # compute largest dimension (hor or ver)
            dim_max = int(max(h, w))
            
            # if size is not preserved or image is larger than output_size, resize image to output_size
            if not(self.preserve_size) or (dim_max > output_size):
                # Resize image so that largest dim is now equal to output_size
                img = lycon.resize(
                    img, 
                    height = max(h*output_size//max(h,w),1), 
                    width = max(w*output_size//max(h,w),1), 
                    interpolation=lycon.Interpolation.AREA
                )
                h = img.shape[0]
                w = img.shape[1]  
            
            # create a square, blank output, of desired dimension
            #img_square = np.ones(output_shape)
            img_square = np.ones(self.image_dimensions[0:2])
            
            # compute number of pixels to leave blank 
            offset_ver = int((output_size-h)/2) # on top and bottom of image
            offset_hor = int((output_size-w)/2) # on left and right of image
            
            # replace pixels in output by input image
            img_square[offset_ver:offset_ver+img.shape[0], offset_hor:offset_hor+img.shape[1]] = img
            
            # Make it 3 dimensional
            img_square = np.stack([img_square, img_square, img_square], axis=2)
            square_images.append(img_square)
        
        # convert to array of images        
        square_images = np.array([img for img in square_images])
        
        #object_ids = [os.path.splitext(os.path.split(self.images_paths[k])[1])[0] for k in indexes]
        #object_ids = np.array([obj for obj in object_ids])
        
        return square_images#, object_ids
    
    
def batch_glimpse(batches, classes, n=1):
    """
    Randomly select an image from a batch and display it with its label
    
    Args:
        batches (DataGenerator): data generator to glimpse at
        classes (array): array of taxonomic classes
        n(int): numer of images to look at
    
    Returns:
        nothing
        
    """
    
    for _ in range(n):
        b = random.randint(0, len(batches)-1)
        image_batch, label_batch = batches[b]
        i = random.randint(0, len(label_batch)-1)
        plt.imshow(image_batch[i][:,:,0], cmap='gray')
        plt.title(classes[np.argmax(label_batch[i])])
        plt.show()
        
    pass

