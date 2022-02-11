import os
import sys
import subprocess

import os.path as osp
import numpy as np
import SimpleITK as sitk

def cip_segmentation(ct_image_file, output_file, metadata):
    input_image = sitk.ReadImage(ct_image_file)

    seeds = list()
    for i, x in metadata.iterrows():
        center = (
            int((x.end_x + x.start_x)/2),
            int((x.end_y + x.start_y)/2),
            int(input_image.GetSize()[2] - x.sliceIndex - 1)
        )
        center_phy = input_image.TransformIndexToPhysicalPoint(center)
        
        ## chest image platform (CIP) segmenation
        partSolid = '--echo'
        # if nodule.loc['Texture'] != 5:
        #    partSolid = '--partSolid'

        # seed generation
        seeds.append('--seeds')
        seeds.append('%f,%f,%f' % center_phy)


    levelset_file = output_file.replace("-label","-levelset")
    nodule_segmentation_cip = 'GenerateLesionSegmentation'
    p = subprocess.Popen(
        [nodule_segmentation_cip, '-i', ct_image_file, '-o', levelset_file, '--maximumRadius', str(metadata.iloc[0].length/2+3), partSolid, '--echo'] + seeds)
    p.wait()
    
    try:
        levelset_image = sitk.ReadImage(levelset_file)
    except:
        raise('Failed CIP')

    ## post processing
    cip_mask_image = sitk.BinaryThreshold(levelset_image, 0)
    sitk.WriteImage(cip_mask_image, output_file, True)
