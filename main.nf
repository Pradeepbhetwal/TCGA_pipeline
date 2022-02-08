dataset = "TCGA-LUAD"
raw_data_loc = "/data/apps/users/wxc151/manifest-KTt2tScD7164745271364348431/TCGA-LUAD"
data_loc = "./data/${dataset}/"

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
    publishDir "${dataset}/${x.baseName.split("_")[0]}"

    input:
    file x from ct_images

    output:
    file "${x.baseName}-label.nrrd"
    
    script:
    """
    touch "${x.baseName}-label.nrrd"
    """
}
