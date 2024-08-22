
options.nclus=3;
%options.RRRprefix='subthr';
%options.reorient='none';
options.reorient = 'none';

list_directory='';
list2run=fullfile(list_directory,'list2run.txt');
subjects=textread(list2run,'%s');

for i_sub=1:length(subjects)
    
    subject_name=subjects{i_sub};
    
    preop_anat_fname = fullfile(in_path,subject_name,'preop','anat',[subject_name,'_ses-preop_acq-T1w.nii']);
    postop_anat_fname = fullfile(in_path,subject_name,'postop','anat',[subject_name,'_ses-postop_acq-T1w.nii']);
    preop_pet_fname = fullfile(in_path,subject_name,'preop','pet',[subject_name,'_ses-preop_acq-pet.nii']);
    postop_rsctman_fname = fullfile(in_path,subject_name,'postop','anat',['manual_',subject_name,'_ses-postop_acq-T1w.nii']); % manual segmentation in postop space (Optional)


    try

        % MRI preop + MRI postop + PET 
        EPIC_CHOP(preop_anat_fname,postop_anat_fname,'preop_pet_fname',preop_pet_fname,'out_path',out_path,'subject_name',subject_name,'options',options) ; 
    
        % MRI preop + MRI_postop
        %EPIC_CHOP(preop_anat_fname,postop_anat_fname,'out_path',out_path,'subject_name',subject_name,'options',options); 
        % MRI with manual resection region in postop space
        %EPIC_CHOP(preop_anat_fname,postop_anat_fname,'postop_rsctman_fname',postop_rsctman_fname,'out_path',out_path,'subject_name',subject_name,'options',options) ; 
        % MRI preop + PET 
        %EPIC_CHOP(preop_anat_fname,preop_anat_fname,'preop_pet_fname',preop_pet_fname,'out_path',out_path,'subject_name',subject_name,'options',options) ;
        % MRI + PET + manual resection 
        %EPIC_CHOP(preop_anat_fname,postop_anat_fname,'preop_pet_fname',preop_pet_fname,'postop_rsctman_fname',postop_rsctman_fname,'out_path',out_path,'subject_name',subject_name,'options',options) ; 
        %EPIC_CHOP(preop_anat_fname,postop_anat_fname,'preop_pet_fname',preop_pet_fname,'postop_rsctman_fname',postop_rsctman_fname,'out_path',out_path,'options',options) ; 

    catch ME
        disp(ME.message);
    end
    
end
