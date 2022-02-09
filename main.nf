dataset = "TCGA-LUAD"
raw_data_loc = "/data/apps/users/wxc151/manifest-KTt2tScD7164745271364348431/TCGA-LUAD"
data_loc = "${projectDir}/data/${dataset}/"
annotations = "${projectDir}/data/${dataset}/annotations.csv" 
clinical_data = "${projectDir}/data/${dataset}/clinical.csv" 

new File(annotations).write(file("https://wiki.cancerimagingarchive.net/download/attachments/33948774/CrowdsCureCancer2017Annotations.csv?version=1&modificationDate=1526921134720&api=v2").text)
new File(clinical_data).write(file("https://wiki.cancerimagingarchive.net/download/attachments/33948774/ccc2017clinical.csv?version=1&modificationDate=1544564445857&api=v2").text)

Channel.fromPath("${raw_data_loc}/TCGA*", type: 'dir').set{ cases_ch }

process dicom2nrrd {
    publishDir "${data_loc}/${x.baseName}"

    input:
    path x from cases_ch

    output:
    file "${x.baseName}_*_CT.nrrd" into ct_images
    
    script:
    """
    /home/wxc151/gitRepos/radiomics-tools/bin/DICOM2NRRDConverter $x output
    mv output/*/*.nrrd . && rm -rf output
    """
}

process segmentation {
    echo true
    publishDir "${dataset}/${x.baseName.split("_")[0]}"

    input:
    file x from ct_images

    output:
    file "${x.baseName}-label.nrrd" into labels
    file x into ct_images1
    
    """
    #!/usr/bin/env python
    import os
    import pandas as pd

    anno = pd.read_csv("$annotations").query("patientID=='${x.baseName.split("_")[0]}'")
    print(anno[['length','start_x', 'start_y', 'end_x', 'end_y', 'sliceIndex', 'annotator','radiologist_status', '_id',]])
    os.system("touch ${x.baseName}-label.nrrd")
    """
}

process feature_extraction {
    echo true

    input:
    file ct from ct_images1
    file label from labels

    output:
    file "${label.baseName}.txt"

    script:
    """
    echo $ct $label
    touch "${label.baseName}.txt"
    """
}
