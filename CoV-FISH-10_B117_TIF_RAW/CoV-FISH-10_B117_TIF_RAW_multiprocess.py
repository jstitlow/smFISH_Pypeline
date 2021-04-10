import getpass
import omero
from omero.gateway import BlitzGateway
import os
import multiprocessing as mp
from multiprocessing import Process
import os
import threading
import time
from skimage.morphology import white_tophat, black_tophat, disk
import numpy as np
import tifffile
import bigfish
import bigfish.stack as stack
import bigfish.segmentation as segmentation
import bigfish.plot as plot
import bigfish.detection as detection
from cellpose import models, io

# setup OMERO connection
PASS = '****'
conn = BlitzGateway('bioc1301', PASS,
        port=4064, group='COVID-19 smFISH', host='omero1.bioch.ox.ac.uk') #davisgroup
conn.connect()
conn.SERVICE_OPTS.setOmeroGroup(-1)

# helper function to display OMERO objects
def print_obj(obj, indent=0):
    print ("""%s%s:%s  Name:"%s" (owner=%s)""" % (
        " " * indent,
        obj.OMERO_CLASS,
        obj.getId(),
        obj.getName(),
        obj.getOwnerOmeName()))

# keep connection to OMERO alive
def keep_connection_alive():
    while True:
        conn.keepAlive()
        time.sleep(60)

th_ka = threading.Thread(target = keep_connection_alive)
th_ka.daemon = True
th_ka.start()

# hard-code the output directory and datasetId
if not os.path.exists('detections'):
    os.makedirs('detections')
path_output = 'detections'

# create a dictionary with key = dataset_id and value = file search string
dataset_id = {
    '25097':'2h',
    '25098':'6h'}

#    '25099':'8h',
#    '25100':'24h'

# make a list of all files in directory
omero_dir = {}
for id in dataset_id:
    for image in conn.getObject('Dataset', id).listChildren():
        for orig_file in image.getImportedImageFiles():
            #omero_dir.append(orig_file.getName())
            omero_dir[orig_file.getName()] = conn.getObject('Image', image.getId())

# helper function to display OMERO objects
def print_obj(obj, indent=0):
    print ("""%s%s:%s  Name:"%s" (owner=%s)""" % (
        " " * indent,
        obj.OMERO_CLASS,
        obj.getId(),
        obj.getName(),
        obj.getOwnerOmeName()))

# write separate functions to be processed in parallel with multiprocessing
def func1():
    print ('started func1')
    for image in omero_dir:
        if list(dataset_id.values())[0] in image and 'B117_INF' in image:
            print ('processing ', image)
            detect_spots(image, omero_dir.get(image), 2, 'ch3')
            detect_spots(image, omero_dir.get(image), 3, 'ch4')
    print ('finished func1')

def func2():
    print ('started func2')
    for image in omero_dir:
        if list(dataset_id.values())[0] in image and 'B117_RDV' in image:
            print ('processing ', image)
            detect_spots(image, omero_dir.get(image), 2, 'ch3')
            detect_spots(image, omero_dir.get(image), 3, 'ch4')
    print ('finished func2')

def func3():
    print ('started func3')
    for image in omero_dir:
        if list(dataset_id.values())[0] in image and 'MOCK' in image:
            print ('processing ', image)
            detect_spots(image, omero_dir.get(image), 2, 'ch3')
            detect_spots(image, omero_dir.get(image), 3, 'ch4')
    print ('finished func3')

def func4():
    print ('started func4')
    for image in omero_dir:
        if list(dataset_id.values())[0] in image and 'Vic_INF' in image:
            print ('processing ', image)
            detect_spots(image, omero_dir.get(image), 2, 'ch3')
            detect_spots(image, omero_dir.get(image), 3, 'ch4')
    print ('finished func4')

def func5():
    print ('started func5')
    for image in omero_dir:
        if list(dataset_id.values())[0] in image and 'Vic_RDV' in image:
            print ('processing ', image)
            detect_spots(image, omero_dir.get(image), 2, 'ch3')
            detect_spots(image, omero_dir.get(image), 3, 'ch4')
    print ('finished func5')

def func6():
    print ('started func6')
    for image in omero_dir:
        if list(dataset_id.values())[1] in image and 'B117_INF' in image:
            print ('processing ', image)
            detect_spots(image, omero_dir.get(image), 2, 'ch3')
            detect_spots(image, omero_dir.get(image), 3, 'ch4')
    print ('finished func6')

def func7():
    print ('started func7')
    for image in omero_dir:
        if list(dataset_id.values())[1] in image and 'B117_RDV' in image:
            print ('processing ', image)
            detect_spots(image, omero_dir.get(image), 2, 'ch3')
            detect_spots(image, omero_dir.get(image), 3, 'ch4')
    print ('finished func7')

def func8():
    print ('started func8')
    for image in omero_dir:
        if list(dataset_id.values())[1] in image and 'MOCK' in image:
            print ('processing ', image)
            detect_spots(image, omero_dir.get(image), 2, 'ch3')
            detect_spots(image, omero_dir.get(image), 3, 'ch4')
    print ('finished func8')

def func9():
    print ('started func9')
    for image in omero_dir:
        if list(dataset_id.values())[1] in image and 'Vic_INF' in image:
            print ('processing ', image)
            detect_spots(image, omero_dir.get(image), 2, 'ch3')
            detect_spots(image, omero_dir.get(image), 3, 'ch4')
    print ('finished func9')

def func10():
    print ('started func10')
    for image in omero_dir:
        if list(dataset_id.values())[1] in image and 'Vic_RDV' in image:
            print ('processing ', image)
            detect_spots(image, omero_dir.get(image), 2, 'ch3')
            detect_spots(image, omero_dir.get(image), 3, 'ch4')
    print ('finished func10')

def detect_spots(image, imageId, chan, ch):
    start_time = time.time()
    #for file in conn.getObject('Dataset', dataset_id).listChildren():
    #    img = conn.getObject('Image', file.getId())
    img = imageId
    print ("data inport took ", time.time() - start_time, "sec")

    # segment with cellpose
    if 'R1' in image:
        seg_img = np.max(get_z_stack(img,2),0)
        seg_img = np.clip(seg_img,0,200)
    else:
        seg_img = np.max(get_z_stack(img,1),0)
    model = models.Cellpose(gpu=False, model_type='cyto')
    channels = [0,0]
    masks, flows, styles, diams = model.eval(seg_img, channels=channels, diameter=375, do_3D=False)

    # detect spots
    rna = get_z_stack(img, chan)

    # subtract background
    rna_no_bg = []
    for z in rna:
        z_no_bg = subtract_background(z)
        rna_no_bg.append(z_no_bg)
    rna = np.array(rna_no_bg)

    # calculate_psf(voxel_size_z, voxel_size_yx, 570, 610, 1.4, 1.364, 'confocal')[0]
    psf_z = calculate_psf(img_params(image)[0], img_params(image)[1],
                          570, 610, 1.4, 1.364, 'confocal')[0]
    psf_yx = calculate_psf(img_params(image)[0], img_params(image)[1],
                          570, 610, 1.4, 1.364, 'confocal')[1]
    sigma_z, sigma_yx, sigma_yx = detection.get_sigma(img_params(image)[0],
                                    img_params(image)[1], psf_z, psf_yx)
    sigma = (sigma_z, sigma_yx, sigma_yx)

    # LoG filter
    rna_log = stack.log_filter(rna, sigma)

    # local maximum detection
    mask = detection.local_maximum_detection(rna_log, min_distance=sigma)

    # thresholding
#    if 'R1' in image:
#        threshold = 100
#    else:
#        threshold = 30
    threshold = detection.automated_threshold_setting(rna_log, mask)
#    threshold = 150
    spots, _ = detection.spots_thresholding(rna_log, mask, threshold)

    # detect and decompose clusters
    spots_post_decomposition, clusters, reference_spot = detection.decompose_cluster(
    rna, spots,
    img_params(image)[0], img_params(image)[1], psf_z, psf_yx,
    alpha=0.7,  # alpha impacts the number of spots per cluster
    beta=1)   # beta impacts the number of detected clusters

    # separate spots from clusters
    radius = 350
    nb_min_spots = 4
    spots_post_clustering, foci = detection.detect_foci(spots_post_decomposition,
        img_params(image)[0], img_params(image)[1], radius, nb_min_spots)

    # extract cell level results
    image_contrasted = stack.rescale(rna, channel_to_stretch=0)
    image_contrasted = stack.maximum_projection(image_contrasted)
    rna_mip = stack.maximum_projection(rna)

    fov_results = stack.extract_cell(
        cell_label=masks.astype(np.int64),
        ndim=3,
        rna_coord=spots_post_clustering,
        others_coord={"foci": foci},
        image=image_contrasted,
        others_image={"smfish": rna_mip})

    # save results
    for i, cell_results in enumerate(fov_results):
        path = os.path.join(path_output, image.split('.')[0]+'_'+ch+'_results_cell_{0}.npz'.format(i))
        stack.save_cell_extracted(cell_results, path)

    print ("bigfish analysis took ", time.time() - start_time, "sec")

# get acquisition parameters
def img_params(file_name):
    if 'R1' in file_name:
        voxel_size_z = 200
        voxel_size_yx = 65
        return voxel_size_z, voxel_size_yx
    else:
        voxel_size_z = 200
        voxel_size_yx = 110
        return voxel_size_z, voxel_size_yx
#psf_z = 350
#psf_yx = 150

# calculate psf (thank you MK), with edit for consistent nomenclature
def calculate_psf(voxel_size_z, voxel_size_yx, Ex, Em, NA, RI, microscope):
    '''
    Use the formula implemented in Matlab version (sigma_PSF_BoZhang_v1)
    to calculate the theoretical PSF size.
    '''
    if microscope == 'widefield':
        psf_yx = 0.225*Em/NA
        psf_z = 0.78*RI*Em/(NA**2)
    elif microscope in {'confocal', 'nipkow'}:
        psf_yx = 0.225/NA*Ex*Em/np.sqrt(Ex**2 + Em**2)
        psf_z = 0.78*RI/NA**2*Ex*Em/np.sqrt(Ex**2 + Em**2)
        #original matlab code commented below:
        #widefield:
        #sigma_xy = 0.225 * lambda_em / NA ;
        #sigma_z  = 0.78 * RI * lambda_em / (NA*NA) ;
        #confocal/nipkow
        #sigma_xy = 0.225 / NA * lambda_ex * lambda_em / sqrt( lambda_ex^2 + lambda_em^2 ) ;
        #sigma_z =  0.78 * RI / (NA^2) *  lambda_ex * lambda_em / sqrt( lambda_ex^2 + lambda_em^2 ) ;
    else:
        print(f'microscope={microscope} is not a valid option')
        sys.exit()
    return psf_z, psf_yx

# convert OMERO Image object into np array
def get_z_stack(img, c, t=0):
    z_range = range(0, img.getSizeZ(), 1)
    zct_list = [(z, c, t) for z in z_range]
    pixels = img.getPrimaryPixels()
    return np.array(list(pixels.getPlanes(zct_list)))

# subtract background
def subtract_background(image, radius=5, light_bg=False):
    str_el = disk(radius) #you can also use 'ball' here to get a slightly smoother result at the cost of increased computing time
    if light_bg:
        return black_tophat(image, str_el)
    else:
        return white_tophat(image, str_el)

if __name__=='__main__':
    mp.set_start_method('spawn')
    p1 = Process(target = func1) # '2h_B117_INF'
    p1.start()
    p2 = Process(target = func2) # '2h_B117_RDV'
    p2.start()
    p3 = Process(target = func3) # '2h_MOCK'
    p3.start()
    p4 = Process(target = func4) # '2h_Vic_INF'
    p4.start()
    p5 = Process(target = func5) # '2h_Vic_RDV'
    p5.start()
    p6 = Process(target = func6) # '2h_B117_INF'
    p6.start()
    p7 = Process(target = func7) # '2h_B117_RDV'
    p7.start()
    p8 = Process(target = func8) # '2h_MOCK'
    p8.start()
    p9 = Process(target = func9) # '2h_Vic_INF'
    p9.start()
    p10 = Process(target = func10) # '2h_Vic_RDV'
    p10.start()
# based on this code: https://stackoverflow.com/questions/18864859/python-executing-multiple-functions-simultaneously
