Channel.fromPath("/data/apps/users/wxc151/manifest-KTt2tScD7164745271364348431/TCGA-LUAD/TCGA*", type: 'dir').set{ cases_ch }
println cases_ch

process dicom2nrrd {
    input:
    path x from cases_ch
    
    script:
    """
    /home/wxc151/gitRepos/radiomics-tools/bin/DICOM2NRRDConverter $x test
    """
}