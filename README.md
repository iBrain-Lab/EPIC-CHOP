# EPIC-CHOP
EPIlepsy-surgery Cavity CHaracterisation Pipeline
Usage notes:

====== ITEM 1

This algorithm was developed and documented in the following papers,
please cite if you use:

Cahill, V., Sinclair, B., Malpas, C. B., McIntosh, A. M., Chen, Z., Vivash, L. E., ... & O'Brien, T. J. (2019). 
Metabolic patterns and seizure outcomes following anterior temporal lobectomy. 
Annals of neurology, 85(2), 241-250.

for support email: ben.sinclair@monash.edu

====== ITEM 2

The steps you want to run can be specified:
0 - make a copy of images to avoid overwriting, and move them to a new directory. optional
1 - reorients via autoreoirent or center of mass
MRI_pipeline 
2 - steps_segment   
3 - skull strip
4 - steps_create_mask (to get image of brain tissue)
10 - linear registration
5.1 - longitudinal registration
5.2 - apply longitudinal registration to get postop images in preop space
5.3 - apply reverse longitudinal registration to get preop images in postop space
6 - subtract postop from preop to get resected tissue
7 - clean resection region
PET_pipeline
8 - normalise PET to MNI
Normalise to MNI
9.1 -normalise all images in preop space to MNI
9.2 - normalise all relevant postop images to MNI space (postop T1 and manual resection cavity)

====== ITEM 3
Calling the function

EPIC_CHOP(preop_anat_fname, postop_anat_fname,'preop_pet_fname',preop_pet_fname,'postop_rsctman_fname',postop_rsctman_fname,'out_path', out_path, 'subject_name', subject_name, 'options',options)

- Only preop_anat_fname, postop_anat_fname are mandatory
- Other arguments must be specified with preceding string to indicate which argument they are.
- You can specify an output folder (out_path), otherwise the output folder will default to the input folder.
- You can specify a subject name, creating a directory of that name to which the subjects files will be written, otherwise, the output folder will default to the name of preop_anat

====== ITEM 4
if you want to just do the preop pet to preop anat registration (step 8), and warping
to MNI (steps 2 + 9), you still need to specify mandatory argument postop_anat_pet_fname.
Suggest reusing preop_anat_pet_fname

====== ITEM 5
The difference map between preop_anat and postop_anat will contain
registration error and segmentation error. This is detached from the
resection cavity by erosion. The default errosion is set to 1, but you
can change to 2 or 3 if the resection cavity is still connected to
registration/segmentatio  error.

====== ITEM 6
If the resection region is small and some region hasnt been properly segmented the postop image, chop1 could be somethong other than the reseciton cavity.
In this case, check chop2 and chop3

====== ITEM 7
This algorithm depends on the resection cavity having the intensity of
CSF. Thus, it only works well for postop anatomical images acquired >
3 months after surgery, by which time surgical debris has mostly cleared. 

====== ITEM 8
Output Files

Note 1: processed images are kept in the directory of their unprocessed image, not in the directory corresponding to the image coordinates of the processed file 
Note 2: From the outset, all images' origin is set to their center of mass, for facilitating initial alignment needed for registrations. If you wish to move the processed images back to the space of the raw image, just need to copy the affine matrix from the raw image to the processed images.

Important output files:
SPACE: preop (centred to center of mass)
preop_anat		- /out_path/subject_name/preop/anat/com_[preop-anat-origname].nii 
postop_anat		- /out_path/subject_name/postop/anat/LRpre_com_[postop-anat-origname].nii
resection_cavity	- /out_path/subject_name/postop/anat/chop1_d2_clus1_e1_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii
preop_pet		- /out_path/subject_name/preop/anat/pet/cpre_com_[preop-pet-origname].nii
postop_rsctman		- /out_path/subject_name/postop/anat/LRpre_com_manual_[postop-anat-origname].nii
SPACE: postop (centred to center of mass)
preop_anat		- /out_path/subject_name/preop/anat/
postop_anat		- /out_path/subject_name/postop/anat/com_[postop-anat-origname].nii
resection_cavity	- /out_path/subject_name/postop/anat/thr_LRpost_chop1_d2_clus1_e1_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii
preop_pet		- NA
postop_rsctman		- /out_path/subject_name/postop/anat/com_manual_[postop-anat-origname].nii
SPACE: MNI
preop_anat		- /out_path/subject_name/preop/anat/wcom_[preop-anat-origname].nii
postop_anat		- /out_path/subject_name/postop/anat/wLRpre_com_[postop-anat-origname].nii
resection_cavity	- /out_path/subject_name/postop/anat/wchop1_d2_clus1_e1_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii
preop_pet		- /out_path/subject_name/preop/anat/pet/wcpre_com_[preop-pet-origname].nii
postop_rsctman		- /out_path/subject_name/postop/anat/wLRpre_com_manual_[postop-anat-origname].nii

Full list of output files:
______________________________
/out_path/subject_name/preop/anat:

STEP 1: 
com_[preop-anat-origname].mat          	- the affine matrix of preop_anat centered to centre of mass
com_[preop-anat-origname].nii          	- preop_anat centered to centre of mass
STEP 2: 
c1com_[preop-anat-origname].nii        	- GM
c2com_[preop-anat-origname].nii        	- WM
c3com_[preop-anat-origname].nii        	- CSF
c4com_[preop-anat-origname].nii        	- Tissue 4
c5com_[preop-anat-origname].nii        	- Tissue 5
com_[preop-anat-origname]_seg8.mat     	-
mcom_[preop-anat-origname].nii         	-
ySeg_com_[preop-anat-origname].nii     	- deformation field from preop_anat to MNI
STEP 3:
brain_com_[preop-anat-origname].nii    	- T1 masked with GM+WM+CSF, for better registration to PET
STEP 4:
c1c2_com_[preop-anat-origname].nii     	- GM+WM
STEP 5:
avg_com_[preop-anat-origname].nii      	- average of preop and postop anat for unbiased registration
dv_com_[preop-anat-origname]_com_[postop-anat-origname].nii 	- velocity field of longitudinal registration
y_LR_com_[preop-anat-origname].nii     	- longitudinal registration deformation field, from preop_anat to postop_anat
STEP 6:
thr_c1c2_com_[preop-anat-origname].nii 	- thresholded GM+WM
STEP 9:
wcom_[preop-anat-origname].nii         	- preop_anat warped to MNI
_________________________________
/out_path/subject_name/postop/anat:

STEP 1:
com_manual_[postop-anat-origname].nii		- manual resection centered to centre of mass
com_[postop-anat-origname].mat			- the affine matrix of postop_anat centered to centre of mass
com_[postop-anat-origname].nii			- postop_anat centered to centre of mass
STEP 2:
c1com_[postop-anat-origname].nii		- GM probability
c2com_[postop-anat-origname].nii		- WM probability
c3com_[postop-anat-origname].nii		- CSF probability
c4com_[postop-anat-origname].nii		- Tissue 4 probability
c5com_[postop-anat-origname].nii		- Tissue 5 probability
com_[postop-anat-origname]_seg8.mat		- 
mcom_[postop-anat-origname].nii		-
ySeg_com_[postop-anat-origname].nii		- deformation field from postop_anat to MNI
STEP 3:
brain_com_[postop-anat-origname].nii		- brain image masked with GM+WM+CSF
STEP 4:
c1c2_com_[postop-anat-origname].nii		- GM+WM probability
STEP 5:
LRpre_com_[postop-anat-origname].nii		- postop T1 lonigtudinaly registered to preop space
y_LR_com_[postop-anat-origname].nii		- longitudinal registration deformation field,  from postop to preop
LRpre_c1c2_com_[postop-anat-origname].nii	- GM+WM probability moved to preop space
LRpre_c1com_[postop-anat-origname].nii		- GM probability moved to preop space
LRpre_c2com_[postop-anat-origname].nii		- WM probability moved to preop space
LRpre_com_manual_[postop-anat-origname].nii	- manual segmentation moved to preop space
LRpost_d2_largest_e2_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii	- resection region moved from preop space to postop
STEP 6:
thr_LRpre_c1c2_com_[postop-anat-origname].nii			- threshold GM+WM probability at 0.1 to give brain mask
c1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii		- preop brain mask minus postop brain mask		
thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii	- retain the positive voxels of subtraction mask to give resection mask
STEP 7:
e1_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii			- erode resection mask to isolate resection cavity
clus1_e1_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii		- select cluster n
d2_clus1_e1_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii		- dilate cluster n to restor orginal size (plus a bit)
chop1_d2_clus1_e1_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii	- multiply restored image by subtraction mask to remove overdilations
chop2_d2_clus2_e1_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii	- second and third largest clusters also supplied
STEP 9:
wcom_manual_[postop-anat-origname].nii						- manual segmentation warped to MNI
wcom_[postop-anat-origname].nii							- postop_anat warped to MNI
wchop1_d2_clus1_e1_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii	- resection mask warped to MNI
wLRpre_com_[postop-anat-origname].nii						- postop_anat moved to preop_anat, then moved to MNI
wLRpre_com_manual_[postop-anat-origname].nii					- manual segmentation moved to preop_anat, then to MNI
STEP 5.3:
LRpost_chop1_d2_clus1_e1_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii	- resection mask moved to postop_anat space
thr_LRpost_chop1_d2_clus1_e1_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii- resection mask moved to postop_anat space and thresholded


_________________________________
/out_path/subject_name/preop/pet:

STEP 1:
com_[preop-pet-origname].mat		- the affine matrix of preop_pet centered to centre of mass		
com_[preop-pet-origname].nii		- postop_anat centered to centre of mass
STEP 8:
tmp_com_[preop-pet-origname].nii	- copy of preop_pet
cpre_com_[preop-pet-origname].nii	- preop_pet coregistered to preop_anat
brain_cpre_com_[preop-pet-origname].nii	- mask coregistered pet with GM+WM+CSF
STEP 9:
wcpre_com_[preop-pet-origname].nii	- coregistered pet moved to MNI space
wbrain_cpre_com_[preop-pet-origname].nii- masked coregistered pet moved to MNI space
