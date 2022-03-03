dataset = "TCGA-LUAD"
raw_data_loc = "/data/apps/users/wxc151/manifest-KTt2tScD7164745271364348431/TCGA-LUAD"
data_loc = "${projectDir}/data/${dataset}/"
annotations = "${data_loc}/annotations.csv" 
clinical_data = "${data_loc}/clinical.csv" 

Channel.fromPath("${raw_data_loc}/TCGA*", type: 'dir').set{ cases_ch }

process metadata {
    publishDir "${data_loc}"

    output:
    file "annotations.csv" into ann
    file "clinical.csv" into cl

    exec:
    new File("${task.workDir}/annotations.csv").write(file("https://wiki.cancerimagingarchive.net/download/attachments/33948774/CrowdsCureCancer2017Annotations.csv?version=1&modificationDate=1526921134720&api=v2").text)
    new File("${task.workDir}/clinical.csv").write(file("https://wiki.cancerimagingarchive.net/download/attachments/33948774/ccc2017clinical.csv?version=1&modificationDate=1544564445857&api=v2").text)
}

process dicom2nrrd {
    publishDir "${data_loc}/${x.baseName}"
    container 'wookjinchoi/radiomics-tools:latest'
    containerOptions "--volume ${raw_data_loc}:${raw_data_loc}"

    input:
    path x from cases_ch

    output:
    file "${x.baseName}_CT.nrrd" into ct_images

    when:
    x.baseName != "TCGA-17-Z016" && x.baseName != "TCGA-17-Z042" && x.baseName != "TCGA-50-8459"
    
    script:
    """
    DICOM2NRRDConverter $x output
    mv output/*/*CT*.nrrd ${x.baseName}_CT.nrrd && rm -rf output
    """
}

process segmentation {
    //echo true
    publishDir "${data_loc}/${ct.baseName.split("_")[0]}"
    container 'acilbwh/chestimagingplatform:latest'

    input:
    file ct from ct_images
    file a from ann
    file c from cl

    output:
    file "${ct.baseName}-label.nrrd" into labels
    file ct into ct_images1
    
    """
    #!/usr/bin/env python
    import os
    import sys
    sys.path.append("$projectDir/scripts") # python2
    import pandas as pd
    from segmentation import *

    anno = pd.read_csv("$annotations") \
                .query("patientID=='${ct.baseName.split("_")[0]}'") # and radiologist_status=='radiologist'") \
                .sort_values(by=['length'], ascending=False)
    
    cip_segmentation("$ct", "${ct.baseName}-label.nrrd", anno)
    """
}

process feature_extraction {
    //echo true
    publishDir "${data_loc}/${ct.baseName.split("_")[0]}"
    container 'wookjinchoi/radiomics-tools:latest'

    input:
    file ct from ct_images1
    file label from labels

    output:
    file "${label.baseName}.txt" into features
    file "${ct.baseName}-1mm.nrrd"
    file "${ct.baseName}-1mm-label.nrrd"

    when:
    ct.baseName.split("_")[0] == label.baseName.split("_")[0]

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.path.append("$projectDir")
    from scripts.feature_extraction import *

    feature_extraction("$ct", "$label", "${label.baseName}.txt")
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
