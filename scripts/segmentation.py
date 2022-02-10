import os
import sys
import subprocess

import os.path as osp
import numpy as np
import SimpleITK as sitk

def cip_segmentation(ct_image_file, output_file, metadata):
    input_image = sitk.ReadImage(ct_image_file)

    center = (
        int((metadata.iloc[0].end_x + metadata.iloc[0].start_x)/2),
        int((metadata.iloc[0].end_y + metadata.iloc[0].start_y)/2),
        int(input_image.GetSize()[2] - metadata.iloc[0].sliceIndex - 1)
    )
    center_phy = input_image.TransformIndexToPhysicalPoint(center)
    
    ## chest image platform (CIP) segmenation
    partSolid = '--echo'
    # if nodule.loc['Texture'] != 5:
    #    partSolid = '--partSolid'

    # seed generation
    seeds = list()
    seeds.append('--seeds')
    seeds.append('%f,%f,%f' % center_phy)
    # try:
    #     np_mask_image = sitk.GetArrayFromImage(mask_image) > 0
    #     np_input_image = sitk.GetArrayFromImage(input_image)

    #     #points = np.matrix(np.where(np.bitwise_and(np_staple_image > 0.9, np_input_image > np.percentile(np_input_image[np_mask_image], 95))))
    #     points = np.matrix(np.where(np.bitwise_and(np_mask_image, np_input_image > np.percentile(np_input_image[np_mask_image], 90))))
    #     seed_points = points[:, np.random.choice(range(0, points.shape[1]), size=min(points.shape[1], 100), replace=False)]

    #     for seed in seed_points.transpose().tolist():
    #         seed_phy = mask_image.TransformIndexToPhysicalPoint(
    #             (seed[2], seed[1], seed[0]))
    #         seeds.append('--seeds')
    #         seeds.append('%f,%f,%f' % seed_phy)
    #         # print(seed_phy)
    # except:
    #     print(seeds)

    levelset_file = output_file.replace("-label","-levelset")
    nodule_segmentation_cip = osp.relpath('/home/wxc151/gitRepos/radiomics-tools/Tools/NoduleSegmentation_CIP/GenerateLesionSegmentation')
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
