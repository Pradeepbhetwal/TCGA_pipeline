import subprocess

import os.path as osp


feature_extractor = 'FeatureExtraction'
image_resample = '/radiomics-tools/PythonTools/image_resample.py'
iso_size = "1"

def feature_extraction(image_file, mask_file, output_file):
    mask_file_resize = mask_file.replace("-label.nrrd", "-" + iso_size + "mm-label.nrrd")
    print(mask_file[0], mask_file_resize)
    image_file_resize = image_file.replace(".nrrd", "-" + iso_size + "mm.nrrd")
    print(image_file, image_file_resize)
    p = subprocess.Popen(["python", image_resample, image_file, image_file_resize, iso_size, iso_size, iso_size])
    p.wait()
    p = subprocess.Popen(["python", image_resample, mask_file, mask_file_resize, image_file_resize, "1"])
    p.wait()
    p = subprocess.Popen([feature_extractor, image_file_resize, mask_file_resize, output_file, '1', '2igscr', '64'])
    p.wait()
