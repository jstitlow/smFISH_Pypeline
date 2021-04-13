# smFISH_Pypeline
Processing and analysis of smFISH data using Python repos (cellpose and bigfish)

## TODO
* Edit the detect_spots() functions
* Edit the img_params() function
* Edit the calculate_psf function

## Setting up the Python environment
### How I set it up:    
1. Clone this environment into your source code directory   
```git clone https://github.com/jstitlow/smFISH_Pypeline.git```
2. Create a Conda environment from the cellpose .yml file (environment.yaml)   
```conda env create -f environment.yml```   
```conda activate cellpose```
3. Install bigfish   
```pip install big-fish```
4. Install OMERO and ICE   
```conda install -c ome zeroc-ice36-python omero-py```
5. Install other dependencies with conda as needed
* e.g., tifffile and pandas

### How you should set it up:
1. Create a .yml file from this repo   
```conda env export > smFISH_environment.yml```
2. Create an environment from the smFISH environment file    
```conda enc create -f smFISH_environment.yml```   
```conda activate cellpose```

## Overview of the code
### Processing on a server
#### Advantages of processing on a server
* can spread the analysis over multiple CPUs
* run the analysis in a screen session
** our data are on an OMERO server in Oxford, so running the code on a server in the same room is better than moving data back and forth to the US. I/O is a major bottleneck.
** screen also maintains the network connection

#### Disadvantages of processing on a server
* graphical interfaces are slow, at least on mac
#### How to run the code
1. Edit the CoV-FISH-10_B117_TIF_RAW_multiprocess.py file
* set OMERO login parameters   
```
PASS = '****'
conn = BlitzGateway('bioc1301', PASS,
        port=4064, group='COVID-19 smFISH', host='omero1.bioch.ox.ac.uk') #davisgroup
conn.connect()
conn.SERVICE_OPTS.setOmeroGroup(-1)
```   
* setup dictionary with datasetId and string to search for filenames   
```
dataset_id = {
    '25097':'2h',
    '25098':'6h'}
```
* setup separate functions to process the data in parallel
```
def func1():
    print ('started func1')
    for image in omero_dir:
        if list(dataset_id.values())[0] in image and 'B117_INF' in image:
            print ('processing ', image)
            detect_spots(image, omero_dir.get(image), 2, 'ch3')
            detect_spots(image, omero_dir.get(image), 3, 'ch4')
    print ('finished func1')
```
** In the function above, I am iterating through the list of files in the omero dataset and selecting files that have the search string from the dataset_id dictionary, and B117_INF.
** The function then runs the detect_spots() function for both smFISH channels
* Edit the detect_spots() functions
* Edit the img_params() function
* Edit the calculate_psf function
* Edit the cellpose code
```
# segment with cellpose
if 'R1' in image:
    seg_img = np.max(get_z_stack(img,2),0)
    seg_img = np.clip(seg_img,0,200)
else:
    seg_img = np.max(get_z_stack(img,1),0)
model = models.Cellpose(gpu=False, model_type='cyto')
channels = [0,0]
masks, flows, styles, diams = model.eval(seg_img, channels=channels, diameter=375, do_3D=False)
```
2. Call the multiprocess python script:
```
python CoV-FISH-10_B117_TIF_RAW_multiprocess.py file
```
This will generate .npz files for every cell that is detected by cellpose.
*
3. Call the numpy
