dataset = "TCGA-BRCA"
raw_data_loc = "/data/apps/users/wxc151/manifest-AmUKkZHx1564984923148296294/TCGA-BRCA"
les_loc = "/data/apps/users/wxc151/manifest-AmUKkZHx1564984923148296294/TCGA_Segmented_Lesions_UofC"
data_loc = "${projectDir}/data/${dataset}/"

Channel.fromPath( "${raw_data_loc}/TCGA*", type: 'dir').set{ cases_ch }


process dicom2nrrd {
    echo true
    publishDir "${data_loc}/${x.baseName}"
    container 'wookjinchoi/radiomics-tools:latest'
    containerOptions "--volume ${raw_data_loc}:${raw_data_loc}"

    input:
    path x from cases_ch
    path les_loc

    output:
    file "${x.baseName}_MR.nrrd" into mr_images
    file "${les_loc}/${x.baseName}*.les" into lesions
    
    script:
    """
    DICOM2NRRDConverter $x output
    patterns=(output/*/*MR*SAG*3D*POST*.nrrd output/*/*MR*Sag*post*.nrrd output/*/*MR*sar*.nrrd output/*/*MR*SAR*.nrrd output/*/*MR*Sag*+C*.nrrd output/*/*MR*106*.nrrd output/*/*MR*.nrrd)
    for p in "\${patterns[@]}"; do
        main_image_file=(\$p)
        if [ -f \${main_image_file} ]; then
            break
        fi
    done
    mv "\$main_image_file" "${x.baseName}_MR.nrrd" && rm -rf output
    """
}


"""
How to use the Segmentations

With regards to the naming structure, *S2-1.les: S2 means DCE-MRI sequence 2, lesion #1. Sometimes, there are multiple DCE-MRI sequences on TCIA data, and so the team used the sequence that corresponded to the one on which the radiologists annotated the truth.  Each of our tumor segmentation files is a binary file, consisting of the following format:

1. six uint16 values for the inclusive coordinates of the lesion’s cuboid , relative to the image:
y_start y_end
x_start x_end
z_start z_end

2. the N int8 on/off voxels (0 or 1) for the above specified cube, where N = (y_end y_start +1) * (x_end - x_start + 1) * (z_end - z_start + 1).

A voxel value of 1 denotes that it is part of the lesion, while a value of zero denotes it is not.

Please reference these data  extracted using version  V2010  of the UChicago MRI Quantitative Radiomics workstation.
"""
process segmentation {
    //echo true
    publishDir "${data_loc}/${mr.baseName.split("_")[0]}"

    //errorStrategy 'ignore'

    input:
    file mr from mr_images
    file les from lesions

    output:
    file "${mr.baseName}-label.nrrd" into labels
    file mr into mr_images1
    
    """
    #!/usr/bin/env python
    import SimpleITK as sitk
    import numpy as np

    def read_les(les_filename, MR_filename = None):
        with open(les_filename, "rb") as f:
            bbox = np.fromfile(f, np.uint16, 6)
            
            y_start, x_start, z_start = bbox[0:3] - 1
            y_end, x_end, z_end = bbox[3:6]

            mask_np = np.fromfile(f, np.uint8)
            mask_np = mask_np.reshape(z_end-z_start, x_end-x_start, y_end-y_start)
            mask_np = mask_np.swapaxes(1,2)        

            mask = sitk.GetImageFromArray(mask_np)
            mask.SetDirection(np.asarray([0,0,-1,1,0,0,0,-1,0], dtype=np.double)) # Sagittal LR

            if MR_filename is not None:
                image = sitk.ReadImage(MR_filename)
                origin = image.GetOrigin()
                spacing = image.GetSpacing()
                start_idx = np.asarray([z_start, -int(x_start), y_start])
                start = start_idx * spacing[::-1]

                mask.SetOrigin(origin-start)
                mask.SetSpacing(spacing)
            
            return mask

    
    mask = read_les(r"${les}","${mr}")
    sitk.WriteImage(mask, "${mr.baseName}-label.nrrd")
    """
}

process feature_extraction {
    //echo true
    publishDir "${data_loc}/${mr.baseName.split("_")[0]}"
    container 'wookjinchoi/radiomics-tools:latest'

    input:
    file mr from mr_images1
    file label from labels

    output:
    file "${label.baseName}.txt" into features

    when:
    mr.baseName.split("_")[0] == label.baseName.split("_")[0]

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.path.append("$projectDir")
    from scripts.feature_extraction import *

    feature_extraction(r"$mr", "$label", "${label.baseName}.txt")
    """
}

process feature_organization {
    //echo true
    publishDir "${data_loc}"
    container 'wookjinchoi/radiomics-tools:latest'

    input:
    file x from features.collect()

    output:
    file "features.csv"

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.path.append("$projectDir")
    from scripts.organization import *

    feature_organization("$x".split(" "), "features.csv")
    """
}
