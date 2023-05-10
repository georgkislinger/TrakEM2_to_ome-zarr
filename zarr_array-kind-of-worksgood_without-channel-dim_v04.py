# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 13:40:59 2023

@author: kislingerg



"""


import os
import re
import glob
import math
import shutil
import tifffile as tif
from tifffile import imread
import numpy as np
#import h5py
from skimage import io
import zarr
from skimage.data import binary_blobs
from ome_zarr.io import parse_url
#from ome_zarr.writer import write_image
from datetime import date
import tkinter as tk
from tkinter import filedialog
import ast
from tqdm import tqdm

#SOMETHING WRONG SINCE UPDATE CHECK!!!!
#%% Get Input output dirs and Voxelsize


def choose_input_dir():
    input_dir = filedialog.askdirectory(title="Choose Input Directory")
    if os.path.isdir(input_dir):
        input_dir_var.set(input_dir)
    else:
        input_dir_var.set("")
        print("Invalid Input Directory")

def choose_output_dir():
    output_dir = filedialog.askdirectory(title="Choose Output Directory")
    if os.path.isdir(output_dir):
        output_dir_var.set(output_dir)
    else:
        output_dir_var.set("")
        print("Invalid Output Directory")

def create_output_dir():
    output_dir = output_dir_var.get()
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    else:
        i = 1
        while os.path.isdir(output_dir + "_" + str(i)):
            i += 1
        output_dir += "_" + str(i)
        os.makedirs(output_dir)
    output_dir_var.set(output_dir)

def validate_entry():
    pixelsize = pixelsize_entry.get()
    try:
        pixelsize_tuple = tuple(float(p) for p in pixelsize.split(","))
        pixelsize_var.set(pixelsize_tuple)
    except:
        pixelsize_var.set("")
        print("Invalid Pixelsize")

def close_dialog():
    root.destroy()

root = tk.Tk()

# Input directory selection
input_dir_var = tk.StringVar()
input_dir_label = tk.Label(root, text="Input Directory:")
input_dir_entry = tk.Entry(root, textvariable=input_dir_var)
input_dir_button = tk.Button(root, text="Choose Input Directory", command=choose_input_dir)
input_dir_label.grid(row=0, column=0, padx=5, pady=5)
input_dir_entry.grid(row=0, column=1, padx=5, pady=5)
input_dir_button.grid(row=0, column=2, padx=5, pady=5)

# Output directory selection
output_dir_var = tk.StringVar()
output_dir_label = tk.Label(root, text="Output Directory:")
output_dir_entry = tk.Entry(root, textvariable=output_dir_var)
output_dir_button = tk.Button(root, text="Choose Output Directory", command=choose_output_dir)
output_dir_create_button = tk.Button(root, text="Create Output Directory", command=create_output_dir)
output_dir_label.grid(row=1, column=0, padx=5, pady=5)
output_dir_entry.grid(row=1, column=1, padx=5, pady=5)
output_dir_button.grid(row=1, column=2, padx=5, pady=5)
output_dir_create_button.grid(row=1, column=3, padx=5, pady=5)

# Pixelsize entry
pixelsize_var = tk.StringVar()
pixelsize_label = tk.Label(root, text="Pixelsize (X,Y,Z):")
pixelsize_entry = tk.Entry(root, textvariable=pixelsize_var)
pixelsize_validate_button = tk.Button(root, text="Validate", command=validate_entry)
pixelsize_label.grid(row=2, column=0, padx=5, pady=5)
pixelsize_entry.grid(row=2, column=1, padx=5, pady=5)
pixelsize_validate_button.grid(row=2, column=2, padx=5, pady=5)

# Create the widgets for dataset name input
dataset_name_var = tk.StringVar()
dataset_name_label = tk.Label(root, text="Dataset name:")
dataset_name_label.grid(row=3, column=0, padx=5, pady=5)
dataset_name_entry = tk.Entry(root, textvariable=dataset_name_var, width=50)
dataset_name_entry.grid(row=3, column=1, padx=5, pady=5)

# Done button
done_button = tk.Button(root, text="Done", fg="blue", cursor="hand2")
done_button.grid(row=4, column=1, padx=5, pady=5)
done_button.bind("<Button-1>", lambda event: close_dialog())

root.mainloop()

input_dir = input_dir_var.get()
input_dir_base = os.path.dirname(input_dir)
reshape_dir = output_dir_var.get()
name_var = dataset_name_var.get()
dataset_name = os.path.normpath(os.path.join(input_dir_base,name_var+'.ome.zarr'))
# Create pixelsize tuple based on pixelsize_var after it has been validated
voxelsize = ast.literal_eval(pixelsize_var.get())



# %% Reshape the input (from TrakEM2) to better usable format

#input_dir = "D:\\test\\source"
#reshape_dir = "D:\\test\\reshaped"

#voxelsize = (4,4,50)

for z_dir in tqdm(os.listdir(input_dir)):
    z_path = os.path.join(input_dir, z_dir)
    if not os.path.isdir(z_path):
        continue

    for mip_dir in os.listdir(z_path):
        mip_path = os.path.join(z_path, mip_dir)
        if not os.path.isdir(mip_path):
            continue

        for filename in os.listdir(mip_path):
            if not filename.endswith('.tif'):
                continue

            y_x = os.path.splitext(filename)[0]
            y, x = map(int, y_x.split('_'))

            # Compute the output path based on the original directory structure
            output_mip_dir = f'{mip_dir.zfill(4)}'
            output_z_dir = f'{z_dir.zfill(4)}'
            output_yx_dir = f'{y:04d}_{x:04d}'
            output_filename = f'{y:04d}_{x:04d}.tif'

            output_path = os.path.join(reshape_dir, output_mip_dir, output_z_dir, output_filename)
            os.makedirs(os.path.dirname(output_path), exist_ok=True)

            shutil.copyfile(os.path.join(mip_path, filename), output_path)
            #print('Currently copying:',os.path.join(mip_path, filename))
            
#%%Multithreaded copy STILL NOT WORKING TODO
'''
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import os
import shutil

def copy_file(filename, input_dir, reshape_dir):
    try:
        if not filename.endswith('.tif'):
            return

        y_x = os.path.splitext(filename)[0]
        y, x = map(int, y_x.split('_'))

        z_path, mip_dir = os.path.split(os.path.split(os.path.join(input_dir, filename))[0])

        # Compute the output path based on the original directory structure
        output_mip_dir = f'{mip_dir.zfill(4)}'
        output_z_dir = f'{os.path.basename(z_path).zfill(4)}'
        output_yx_dir = f'{y:04d}_{x:04d}'
        output_filename = f'{y:04d}_{x:04d}.tif'

        output_path = os.path.join(reshape_dir, output_mip_dir, output_z_dir, output_filename)
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        shutil.copy(os.path.join(input_dir, filename), output_path)
    except Exception as e:
        print(f"Error copying {filename}: {e}")
        
def copy_files(input_dir, reshape_dir):
    files = [(filename, input_dir, reshape_dir) for filename in filter(lambda f: f.endswith('.tif'), [os.path.join(dp, f) for dp, dn, filenames in os.walk(input_dir) for f in filenames])]
    with Pool() as p:
        list(tqdm(p.imap_unordered(copy_file, files, chunksize=25), total=len(files)))

copy_files(input_dir, reshape_dir)
'''
# %% Getting Information on number of mipmaps, dataset shape etc


data_dir = reshape_dir

mip_levels = []
z_range = []
y_range = []
x_range = []

for root, dirs, files in os.walk(data_dir):
    for file in files:
        if file.endswith(".tif"):
            path = os.path.join(root, file)
            # Extract the mip, z, y, and x values from the path using a regular expression
            m = re.search(r"\\(\d+)\\(\d+)\\(\d+)_(\d+).tif$", path)
            if m:
                mip = int(m.group(1))
                z = int(m.group(2))
                y = int(m.group(3))
                x = int(m.group(4))
                mip_levels.append(mip)
                z_range.append(z)
                y_range.append(y)
                x_range.append(x)

# Create arrays of all found numbers for mip, z, y, and x
mip_array = np.array(mip_levels)
z_array = np.array(z_range)
y_array = np.array(y_range)
x_array = np.array(x_range)

# Get unique values for each array and sort them
unique_mip = np.sort(np.unique(mip_array))
unique_z = np.sort(np.unique(z_array))
unique_y = np.sort(np.unique(y_array))
unique_x = np.sort(np.unique(x_array))

print("Unique MIP values:", unique_mip)
print("Unique Z values:", unique_z)
print("Unique Y values:", unique_y)
print("Unique X values:", unique_x)

#%%

max_x_values = {}
max_y_values = {}

#Loop through every unique mip level
for mip in unique_mip:
    # Get the indices of all files in the mip level
    mip_indices = np.where(mip_array == mip)[0]
    # Get the unique z, y, and x values for this mip level
    mip_z = np.unique(z_array[mip_indices])
    mip_y = np.unique(y_array[mip_indices])
    mip_x = np.unique(x_array[mip_indices])
    # Store the maximum unique x and y values for this mip level
    max_x_values[mip] = np.max(mip_x)
    max_y_values[mip] = np.max(mip_y)

# Print the maximum unique values for x and y for every mip level
print("Maximum unique X values for each MIP level:", max_x_values)
print("Maximum unique Y values for each MIP level:", max_y_values)
min_y_lowres = list(max_y_values.items())[-1]
min_x_lowres = list(max_x_values.items())[-1]
min_y_tiles = min_y_lowres[1]+1
min_y_mip = min_y_lowres[0]
min_x_tiles = min_x_lowres[1]+1
min_x_mip = min_x_lowres[0]

#%%

print(len(unique_mip))
#load an image to get dtype and dims
for root, dirs, files in os.walk(input_dir):
    for file in files:
        if file.endswith('.tif'):
            datatile_path = os.path.join(root, file)
            datatile = tif.imread(datatile_path)
            break


# %% Setting the parameters for the zarr array

chunk_height = datatile.shape[0]
chunk_width = datatile.shape[1]
dtype = datatile.dtype
chunks = (1, chunk_height, chunk_width)
#dataset_shape = (1,len(unique_z), len(unique_y)*chunk_height, len(unique_x)*chunk_width)
fill_value = 0
def round_mp(x,num=chunk_height):
    return num * math.ceil(x/num)
#%% Creating the metadata/.zattrs file for the root lvl

root = zarr.open(dataset_name, mode='w')
root.attrs['description'] = 'My Dataset'
root.attrs['author'] = 'Me'
root.attrs['date'] = str(date.today())
      
metadata = [
     
        {   "name": "0",
            "version": "0.4",
            "axes": [
             #   {"name": "c", "type": "channel"},
                {"name": "z", "type": "space", "unit": "nanometer"},
                {"name": "y", "type": "space", "unit": "nanometer"},
                {"name": "x", "type": "space", "unit": "nanometer"},
            ],
            "datasets": [
                {
                    #"axes": [
                    #    {"name": "c", "type": "channel"},
                    #   {"name": "z", "type": "space", "unit": "nanometer"},
                    #   {"name": "y", "type": "space", "unit": "nanometer"},
                    #    {"name": "x", "type": "space", "unit": "nanometer"}],
                    "path": f"{i}".format(i=i),
                    "coordinateTransformations": [{
                        "type": "scale",
                        "scale": [1*voxelsize[2], voxelsize[0]*2**i, voxelsize[1]*2**i],
                        #"translate": [0,(voxelsize[0] * 2**i)*sum(y_shift[:i]),(voxelsize[1] * 2**i)*sum(x_shift[:i])]
                    }]
                } for i in unique_mip
            ],
        }
    ]


print(metadata)

root.attrs["multiscales"] = metadata
root.create_group('OME')
root.OME.attrs['metadata'] = 'currently empty'


#%%
# Open the Zarr array
'''
for mip in unique_mip:
    m=sum(y_shift[:mip])
    print(m)
sum(y_shift)
'''
#store = zarr.DirectoryStore("D:\\test\\my.zarr")
#z_array = zarr.open(store, mode='w')

# Loop through each mip level

for mip in tqdm(unique_mip):
    #create groups
    #mip_group = root.create_group('s'+str(mip))
    #populate with metadata TODO
    #define shape of each level
    shape = (max(unique_z),((min_y_tiles*2**min_y_mip)/2**mip)*chunk_height, ((min_x_tiles*2**min_x_mip)/2**mip)*chunk_width)
    print(shape)
    #shape = (1,max(unique_z),round_mp(int(max((unique_y+1)*chunk_height)/2**mip)), round_mp(int(max((unique_x+1)*chunk_width)/2**mip)))
    z_array= root.create_dataset(str(mip), shape=shape, dtype=dtype, chunks=chunks, fill_value=fill_value,dimension_separator= '/', write_empty_chunks=False)

    #z_array = mip_group.zeros('mip'+str(mip), shape=shape, chunks=(1,1,chunk_height, chunk_width), dtype=dtype)
    for z in tqdm(unique_z):
        try:
            tile_path = os.path.join(data_dir, f'{mip:04d}', f'{z:04d}')
            for tile_file in os.listdir(tile_path):
                if tile_file.endswith(".tif"):
                        
                    # Get the x and y indices from the filename
                    x_y = os.path.splitext(tile_file)[0].split("_")
                    x = int(x_y[1])
                    y = int(x_y[0])
                    tif_path = os.path.join(tile_path, tile_file)
                    print(tif_path)
                    tif_data = tif.imread(tif_path)
                #tif_data = tif.imread(os.path.join(tile_path, tile_file))
                # Write the TIFF data to the Zarr array
                z_array[z-1,y*chunk_height:(y+1)*chunk_height, x*chunk_width:(x+1)*chunk_width] = tif_data
        except FileNotFoundError:
            print(f"Skipping File: {tif_path}")
            continue
          
                
#%%working
#final metadata STILL TODO
'''
mip=4
print(z_array)
unique_z
print(range(unique_z))
range(0,int(max(unique_y)/2**mip))
for mip in unique_mip:
    for z in unique_z:
        for y in range(int(max(unique_y)/2**mip)):
            for x in range(int(max(unique_x)/2**mip)):
                #print(int(max((unique_y)*chunk_height)/2**mip))
                #print(int(max((unique_x)*chunk_height)/2**mip))
                #print(int(max(unique_y)/2**mip))
                print(int(max(unique_x)/2**mip))

            
for mip in range(2,5):
    for z in unique_z:
        for y in unique_y:
            for x in unique_x:
                # Load the TIFF data for this y-x position
                tif_path = os.path.join(data_dir, f'{mip:04d}', f'{z:04d}', f"{y:04d}_{x:04d}.tif")
                try:
                    tif_data = tif.imread(tif_path)
                except FileNotFoundError:
                    print(f"File not found: {tif_path}")
                    continue
                print(tif_path)



for mip in unique_mip:


    for z in unique_z:
        for y in range((math.ceil(max(unique_y)/2**mip))):
            for x in range((math.ceil(max(unique_x)/2**mip))):




for mip in unique_mip:
    mip_group = root.create_group('mip'+str(mip))
    groups.append(mip_group)
    shape = (1,max(unique_z),(math.ceil(max(unique_y)/2**mip))*chunk_height, (math.ceil(max(unique_x)/2**mip))*chunk_width)
    z_array= root['mip'+str(mip)].create_dataset('data', shape=shape, dtype=dtype, chunks=chunks, fill_value=fill_value)
    #z_array = mip_group.zeros('mip'+str(mip), shape=shape, chunks=(1,1,chunk_height, chunk_width), dtype=dtype)
    for z in unique_z:
        for y in range((math.ceil(max(unique_y)/2**mip))):
            for x in range((math.ceil(max(unique_x)/2**mip))):
                # Load the TIFF data for this y-x position
                tif_path = os.path.join(data_dir,f'{mip:04d}', f'{z:04d}', f"{y:04d}_{x:04d}.tif")
                print(tif_path)
                try:
                    tif_data = tif.imread(tif_path)
                except FileNotFoundError:
                    print(f"File not found: {tif_path}")
                    continue
                # Write the TIFF data to the Zarr array
                z_array[0,z-1,y*chunk_height:(y+1)*chunk_height, x*chunk_width:(x+1)*chunk_width] = tif_data

#%%



                
for mip in unique_mip:
        # Calculate the shape of the Zarr array for this mip level and z level

store = zarr.DirectoryStore("D:\\test\\array_out")
z_array = zarr.open(store, mode='w')
groups = []
# Loop through each mip level
for mip in unique_mip:
    mip_group = root.create_group('mip'+str(mip))
    groups.append(mip_group)

root = zarr.group()
root.attrs['description'] = 'My Zarr dataset'
root.attrs['author'] = 'Me'
root.attrs['date'] = '2023-04-21'




for level in groups:
    zarray = level.zeros('data', shape=(1,1,1,1), dtype=dtype)
    
    # Loop through each z level within each mip level
    for z in unique_z:
        # Calculate the shape of the Zarr array for this mip level and z level
        shape = (len(unique_z),(math.ceil(max(unique_y)/2**mip))*chunk_height, (math.ceil(max(unique_x)/2**mip))*chunk_width)
        z_array = mip_group.create_array(
            'section'+str(z),
            shape=shape, 
            chunks=(1,chunk_height, chunk_width), 
            dtype=dtype,
            )

        
        # Loop through each y-x position in the grid
        for y in range((math.ceil(max(unique_y)/2**mip))):
                       
            for x in range((math.ceil(max(unique_x)/2**mip))):
                
                # Load the TIFF data for this y-x position
                tif_path = os.path.join(data_dir,f'{mip:04d}', f'{z:04d}', f"{y:04d}_{x:04d}.tif")
                print(tif_path)
                try:
                    tif_data = tif.imread(tif_path)
                except FileNotFoundError:
                    print(f"File not found: {tif_path}")
                    continue
                # Write the TIFF data to the Zarr array
                z_array[z,y*chunk_height:(y+1)*chunk_height, x*chunk_width:(x+1)*chunk_width] = tif_data
#%%

math.ceil(max(unique_y)/2**6)
range((math.ceil(max(unique_y)/2**mip)))
y_ceil=[]
x_ceil=[]
y_actual=[]
x_actual=[]
y_shift=[]
x_shift=[]

for mip in unique_mip:
    #fill the arrays to calculate x and y shift
    y_c = round_mp(int(max((unique_y)*chunk_height)/2**mip))
    x_c = round_mp(int(max((unique_x)*chunk_width)/2**mip))
    y_a = max((unique_y)*chunk_height)/2**mip
    x_a = max((unique_x)*chunk_width)/2**mip
    y_s = y_c-y_a
    x_s = x_c-x_a
    y_ceil.append(y_c)
    x_ceil.append(x_c)
    y_actual.append(y_a)
    x_actual.append(x_a)
    y_shift.append(y_s)
    x_shift.append(x_s)
    
    
  
40158-39441 
9005-8772
y_shift
x_shift
for mip in unique_mip:
    print(mip)
max(unique_y)/2**6
root = zarr.group()
root.attrs['description'] = 'My Zarr dataset'
root.attrs['author'] = 'Me'
root.attrs['date'] = '2023-04-21'
compressor = None
fill_value = 0





root = zarr.group(store, overwrite=True)
root = zarr.open('D:\\test\\my.zarr')
for mip in unique_mip:
    mip_group = root.create_group('mip'+str(mip))
    groups.append(mip_group)
root.attr[]
for level in groups:
    
mip1 = root.create_group('mip1')
mip2 = root.create_group('mip2')

for group in [mip1, mip2]:
    group.attrs['multiscales'] = [{'version': '0.1', 'datasets': [{'path': '1', 'key': 'image', 'shape': (1, 3, 4, 4), 'dtype': 'uint8', 'compressor': None, 'fill_value': 0}]}]

    zarr_array = group.create_dataset(
        '1/image',
        shape=(1, 3, 4, 4),
        chunks=(1, 3, 4, 4),
        dtype='uint8',
        compressor=None,
        fill_value=0
    )
'''