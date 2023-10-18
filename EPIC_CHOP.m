function EPIC_CHOP(preop_anat_fname, postop_anat_fname, varargin)
% EPIC-CHOP the greatest chopper of them all.
% EPIlepsy surgery Cavity CHaracterisatiOn Pipeline 
%
% This program does 4 things:
% 1) Calculates a resection region in preoperative space
% 2) Coregisters a PET image to a preop MRI image
% 3) Moves all preoperative and postoperative and resection images to MNI
% space
% 4) if a manual resection mask in postoperative space is specified,
% registers this to preop space and MNI space
%
%
% The special features of this algorithm are:
% 1) corrects for postoperative tissue deformation by using Jhn Ashburner
% and Ged Ridgway's longitudinal registration in SPM
% 2) removes registration and segmentation error in resection region by
% eroding and dilating.
%__________________________________________________________________________
% 
% Inputs:
% 
% FORMAT EPIC_CHOP(out_path, subject_name, preop_anat_fname, postop_anat_fname,'preop_pet_fname',preop_pet_fname,'postop_rsctman_fname',postop_rsctman_fname,'options',options)
%
% Mandatory arguments:
% preop_anat_fname  - the full filename (including path and extension) of the preoperative T1w image. Must be a char between single quotation marks ''
% postop_anat_fname - the full filename (including path and extension) of the postoperative T1w image. Must be a char between single quotation marks ''
%
%
% Optional arguments:
% NB, these must be specified with ('argument_name',argument_value) format
% preop_pet_fname       - the full filename (including path and extension) of the preoperative pet image. Must be a char between single quotation marks ''
% postop_rsctman_fname  - the full filename (including path and extension) of a manually segmented resection mask in postop space. 
% out_path              - directory to which output files will be written in BIDS format: out_path/sub/sess/acq
% subject_name          - a string with the subjects name, used for outputing in BIDS format, deafult: preop_anat_fname without path or extension
% options               - a struct with the following fields:
% options.steps         - a vector with the steps you want to execute
% options.LRprefix      - a string with the prefix for images longitudinally registered to the preop space, default: 'LRpre_'
% options.LRparams      - a vector with the parameters for longitudinal registration, default: [0 0 100 25 100]
% options.RRRprefix     - a string specifying which order you want the threshhold and subtract postop from preop segmentations. default: 'thrsub';
% options.reorient      - a string specifying whether and how you wush to do initial reoritation of image, default: 'centreofmass';
%
% Outputs:
%
% Important output files:
% SPACE: preop (centred to center of mass)
% preop_anat		- /out_path/subject_name/preop/anat/com_[preop-anat-origname].nii 
% postop_anat		- /out_path/subject_name/postop/anat/LRpre_com_[postop-anat-origname].nii
% resection_cavity	- /out_path/subject_name/postop/anat/chop1_d2_clus1_e1_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii
% preop_pet         - /out_path/subject_name/preop/anat/pet/cpre_com_[preop-pet-origname].nii
% postop_rsctman	- /out_path/subject_name/postop/anat/LRpre_com_manual_[postop-anat-origname].nii
% SPACE: postop (centred to center of mass)
% preop_anat		- /out_path/subject_name/preop/anat/
% postop_anat		- /out_path/subject_name/postop/anat/com_[postop-anat-origname].nii
% resection_cavity	- /out_path/subject_name/postop/anat/thr_LRpost_chop1_d2_clus1_e1_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii
% preop_pet         - NA
% postop_rsctman	- /out_path/subject_name/postop/anat/com_manual_[postop-anat-origname].nii
% SPACE: MNI
% preop_anat		- /out_path/subject_name/preop/anat/wcom_[preop-anat-origname].nii
% postop_anat		- /out_path/subject_name/postop/anat/wLRpre_com_[postop-anat-origname].nii
% resection_cavity	- /out_path/subject_name/postop/anat/wchop1_d2_clus1_e1_thrsubc1c2loss_thr_LRpre_c1c2_com_[postop-anat-origname].nii
% preop_pet         - /out_path/subject_name/preop/anat/pet/wcpre_com_[preop-pet-origname].nii
% postop_rsctman	- /out_path/subject_name/postop/anat/wLRpre_com_manual_[postop-anat-origname].nii
% 
%__________________________________________________________________________
% Usage notes:
%
% ====== ITEM 1
%
% This algorithm was developed and documented in the following papers,
% please cite if you use:
%
% Cahill, V., Sinclair, B., Malpas, C. B., McIntosh, A. M., Chen, Z., Vivash, L. E., ... & O'Brien, T. J. (2019). 
% Metabolic patterns and seizure outcomes following anterior temporal lobectomy. 
% Annals of neurology, 85(2), 241-250.
%
% for support email: ben.sinclair@monash.edu
%
% ====== ITEM 2
%
% The steps you want to run can be specified:
% 0 - make a copy of images to avoid overwriting, and move them to a new directory. optional
% 1 - reorients via autoreoirent or center of mass
% MRI_pipeline 
% 2 - steps_segment   
% 3 - skull strip
% 4 - steps_create_mask
% 10 - linear registration
% 5.1 - steps_longitudinal registration
% 5.2 - apply longitudinal registration to get postop images in preop space
% 5.3 - apply reverse longitudinal registration to get preop images in postop space
% 6 - subtract postop from preop to get resected tissue
% 7 - clean resection region
% PET_pipeline
% 8 - normalise PET to MNI
% Normalise to MNI
% 9.1 -normalise all images in preop space to MNI
% 9.2 - normalise all relevant postop images to MNI space (postop T1 and manual resection cavity)
%
% ====== ITEM 3
% you can specify an output folder, otherwise the output folder will default to the input folder

% ====== ITEM 4
% if you want to just do the preop pet to preop anat registration (step 8), and warping
% to MNI (steps 2 + 9), you still need to specify mandatory argument postop_anat_pet_fname.
% Suggest reusing preop_anat_pet_fname
%
% ====== ITEM 5
% The difference map between preop_anat and postop_anat will contain
% registration error and segmentation error. This is detached from the
% resection cavity by erosion. The default errosion is set to 1, but you
% can change to 2 or 3 if the resection cavity is still connected to
% registration/segmentation  error.
%
% ====== ITEM 6
% If the resection region is small and some region hasnt been properly segmented the postop image, 
% chop1 could be somethong other than the reseciton cavity.
% In this case, check chop2 and chop3
%
% ====== ITEM 7
% This algorithm depends on the resection cavity having the intensity of
% CSF. Thus, it only works well for postop anatomical images acquired >
% 3 months after surgery, by which time surgical debris has mostly cleared. 
%
%__________________________________________________________________________
% version history
%
% v10: Cahill et al. 2019
% v11: 1) removes the multiple different pipeline streams in v10 and sticks
% with a single pipeline stream, in order to simplify code.
% 2) adds in transformation to MNI space, and adds back in different ordering of
% thresholding and subtracting mask (option RRR prefix) from v10.
% v12: removes the c_ from filenames and doesnt do reorient to centre of mass
% v13: does 3 upgrades: 1) adds back in reoient to centre of mass from v11
% 2) allows autoreorient and updates the names with ar_ without needing to
% specify this on every line.  3) outputs processed images to user chosen directory.
% v14: does not require BIDS format input for directories. This is helpful for when the
% manual resection mask is in a separate directory. Note, output is still
% in BIDS format.
% v15: uses structs rather than cells for input variable names, for easier
% referencing
% v16: moves job_functions to this script and allows generalised (non-BIDS)
% input format for filenames.
% v17: erodes by 2, dilates by 3, then multplies by difference mask, to
% salvage unrecovered eroded voxels. Also keeps top 3 clusters in case
% largest cluster is not the resection cavity
%__________________________________________________________________________


    

%==================================
% get input parameters
%==================================

% create input parser
p = inputParser;
n_req = 2;

% define parameters and their default values and validators
paramNames = {'preop_anat_fname','postop_anat_fname','preop_pet_fname','postop_rsctman_fname','out_path','subject_name','options'};
default_values = {"","","","","","", struct()};
validators = {@ischar,@ischar,@ischar,@ischar,@ischar,@ischar,@isstruct};

% initialize parameters struct
for i = 1:numel(paramNames)
    parameters.(paramNames{i}).default_value = default_values{i};
    parameters.(paramNames{i}).validator = validators{i};
end

% define required parameters
for i = 1:n_req  % first n_req parameters are required
    addRequired(p, paramNames{i}, parameters.(paramNames{i}).validator);
end

% define optional parameters
for i = n_req+1:numel(paramNames)  % remaining parameters are optional
    addParameter(p, paramNames{i}, parameters.(paramNames{i}).default_value, parameters.(paramNames{i}).validator);
end

parse(p, preop_anat_fname, postop_anat_fname, varargin{:}); 

% Assign parsed inputs to variables without using eval
for i = 1:numel(paramNames)
    assignin('caller', paramNames{i}, p.Results.(paramNames{i}));
end

% check which arguments were supplied
for i = 1:numel(paramNames)
    if ~isequal(p.Results.(paramNames{i}), parameters.(paramNames{i}).default_value) % if the 
        fprintf('Supplied: %s \n', paramNames{i});
        parameters.(paramNames{i}).wasSupplied = true;
    else
        fprintf('Not Supplied: %s \n', paramNames{i});
        parameters.(paramNames{i}).wasSupplied = false;
    end
end


% assign inread variables to variable names
for i = 1:numel(paramNames)
    if parameters.(paramNames{i}).wasSupplied       
        parameters.(paramNames{i}).origname = p.Results.(paramNames{i});
        if ischar(p.Results.(paramNames{i}))
            eval(paramNames{i}+"="+''''+p.Results.(paramNames{i})+'''')
        else
            eval(paramNames{i}+"="+"p.Results."+paramNames{i})
        end
    end
end

% Define the images
imgnames = {'preop_anat_fname','postop_anat_fname','preop_pet_fname','postop_rsctman_fname'};

% create a substruct of parameters containing the imgs supplied names
% it will have fields:
% ImgStruct.origname (inherited from parameters struct)
% ImgStruct.currentname: input to current step
% ImgStruct.outname: output of current step
ImgStruct = struct;
for i = 1:numel(imgnames)
    if parameters.(imgnames{i}).wasSupplied
        ImgStruct.(imgnames{i}) = parameters.(imgnames{i});
        ImgStruct.(imgnames{i}).currentfname = ImgStruct.(imgnames{i}).origname;
        ImgStruct.(imgnames{i}).outfname = ImgStruct.(imgnames{i}).origname;
    end
end
n_img = length(ImgStruct);
in_img = fieldnames(ImgStruct);

% check that images exist
for i_img=1:numel(in_img)  
    if ~exist(ImgStruct.(in_img{i_img}).origname,'file')
        error(['ERROR: ',ImgStruct.(in_img{i_img}).origname,' doesnt exist']);
    end
end

% if out_path wasnt specified, set to inpath of the preop_anat
if ~parameters.out_path.wasSupplied
    [out_path,dummy1,dummy2] = fileparts(ImgStruct.preop_anat_fname.origname);
end
% if subject_name wasnt supplied, set it to the image name of preop_anat
if ~parameters.subject_name.wasSupplied
    [dummy1,subject_name,dummy2] = fileparts(ImgStruct.preop_anat_fname.origname);
end    
out_path_sub = fullfile(out_path,subject_name); 

%=====================================
%=====================================
% Options 
%=====================================
%=====================================

%-------------------------
% Steps
%-------------------------

% MRI_pipeline %options.steps=[0:4,8,9.1,9.2]; % MRI preop + pet preop

% 0 - make a copy of images to avoid overwriting, and move them to a new directory. optional
% 1 - reorients via autoreoirent or center of mass
% 2 - steps_segment   
% 3 - skull strip
% 4 - steps_create_mask
% 10 - linear registration of postop anat to preop anat 
% 5.1 - steps_longitudinal registration
% 5.2 - apply longitudinal registration to get postop images in preop space
% 5.3 - apply reverse longitudinal registration to get preop images in postop space
% 6 - subtract postop from preop to get resected tissue
% 7 - clean resection region
% PET_pipeline
% 8 - normalise PET to MNI
% Normalise to MNI
% 9.1 -normalise all images in preop space to MNI
% 9.2 - normalise all relevant postop images to MNI space (postop T1 and manual resection cavity)


if isempty(fieldnames(p.Results.options))
    options = struct;
end
if ~isfield(options,'steps')
    options.steps=[0:8,5.1,5.2,5.3,9.1,9.2];
    %options.steps=[9.1,9.2];
end
if ~isfield(options,'reorient')
    %options.reorient = 'none';
    %options.reorient = 'autoreorient';
    options.reorient = 'centreofmass';
end
if ~isfield(options,'LRprefix')
    options.LRprefix='LRpre_';
    %options.LRprefix='LRpre501250_';
    %options.LRprefix='LRpre10025100_';
end
if ~isfield(options,'LRparams')
    options.LRparams=[0 0 100 25 100];
    %options.LRparams=[0 0 50 12 50];
end
if ~isfield(options,'RRRprefix')
    options.RRRprefix='thrsub'; % threshold brain masks, then subtract them
    %options.RRRprefix='subthr'; % subtract brain masks, then threshold them
end
if ~isfield(options,'nerode')
    options.nerode=1;
    options.ndilate=options.nerode+1;
    options.nclus=3; % number of clusters to extract
end


% END OF OPTIONS
%=====================================



%=====================================
%=====================================
% Run image processing
%=====================================
%=====================================
tstart_sub = tic; 
spm('defaults','PET');
disp(['Working on subject: ',subject_name]);
disp('Images: ');
disp(in_img);


%---------------------------------
% make a copy of images
%---------------------------------
% we create a copy of the input images because spm likes to overwrite the
% matrices of the images during reorientation and coregistration

% specify new names
for i_img=1:numel(in_img)  
    if strcmp(in_img{i_img},'preop_anat_fname')       ses='preop'; acq='anat'; end;
    if strcmp(in_img{i_img},'postop_anat_fname')      ses='postop'; acq='anat';end;
    if strcmp(in_img{i_img},'preop_pet_fname')        ses='preop'; acq='pet';  end;
    if strcmp(in_img{i_img},'postop_rsctman_fname')   ses='postop'; acq='anat';end;
    ImgStruct.(in_img{i_img}).outpath = fullfile(out_path_sub,ses,acq);
    ImgStruct.(in_img{i_img}).outfname = spm_file(ImgStruct.(in_img{i_img}).origname,'path', ImgStruct.(in_img{i_img}).outpath); % changes path to scratch_nobackup
end
% copy
for i_img=1:numel(in_img)
    if any(options.steps==0)
        mkdir(ImgStruct.(in_img{i_img}).outpath);
        disp(['copying to: ',ImgStruct.(in_img{i_img}).outpath]);
        status = copyfile(ImgStruct.(in_img{i_img}).origname,ImgStruct.(in_img{i_img}).outfname); % makes a copy of the file
    end
end
% update file names
for i_img=1:numel(in_img)
    ImgStruct.(in_img{i_img}).currentfname = ImgStruct.(in_img{i_img}).outfname;
end


%---------------------------------
% reorient images
%---------------------------------

% specify new names
% this is the silly multiple stages for file naming you have to do if you want
% to have multiple different pipelines of processing, but not always run each step.
% nipype fixes this.
for i_img=1:numel(in_img)
    if options.reorient=="none"
        ImgStruct.(in_img{i_img}).outfname = ImgStruct.(in_img{i_img}).currentfname;
    elseif options.reorient=="autoreorient"
        ImgStruct.(in_img{i_img}).outfname = spm_file(ImgStruct.(in_img{i_img}).currentfname,'prefix','aro_','ext','nii'); % declare new name
    elseif options.reorient=="centreofmass"
        ImgStruct.(in_img{i_img}).outfname = spm_file(ImgStruct.(in_img{i_img}).currentfname,'prefix','com_','ext','nii');
    end
end

% set origin to mni or centre of mass.
if any(options.steps==1) 

    for i_img=1:numel(in_img)

        if ~strcmp(in_img{i_img},'postop_rsctman_fname') 
            if options.reorient=="autoreorient"
                job_autoreorient(ImgStruct.(in_img{i_img}).currentfname, ImgStruct.(in_img{i_img}).outfname);
            elseif options.reorient=="centreofmass"
                job_cat_vol_set_com(ImgStruct.(in_img{i_img}).currentfname, ImgStruct.(in_img{i_img}).outfname);
            end
        % apply postop affine matrix to manual segmentation in postop space
        elseif strcmp(in_img{i_img},'postop_rsctman_fname')
            job_cat_vol_set_com(ImgStruct.(in_img{i_img}).currentfname, ImgStruct.(in_img{i_img}).outfname,...
                                spm_file(ImgStruct.postop_anat_fname.outfname,'ext','mat')); 
        end
    end
end

% update filenames
for i_img=1:numel(in_img)
    ImgStruct.(in_img{i_img}).currentfname = ImgStruct.(in_img{i_img}).outfname;
end


% the optional steps are done now, so we can abandon ImgStruct and use more
% interpretable variable names
for i_img=1:numel(in_img)
    eval( in_img{i_img} + "="  + '''' + ImgStruct.(in_img{i_img}).currentfname + '''' + ';'); 
end


%=====================================
% Process MRI
%=====================================  

%---------------------------------
% Segment
%---------------------------------
if any(options.steps==2)
    matlabbatch{1}.spm.spatial.preproc.channel.vols  = cellstr({spm_file(preop_anat_fname,'prefix','','ext','nii');spm_file(postop_anat_fname,'prefix','','ext','nii')});
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
    matlabbatch{1}.spm.spatial.preproc.warp.write    = [0 1];
    spm_jobman('run',matlabbatch); clear matlabbatch;
    % rename deformation field to avoid name clash with longitudinal registration
    system(['mv ',spm_file(preop_anat_fname,'prefix','y_','ext','nii'),' ',spm_file(preop_anat_fname,'prefix','ySeg_','ext','nii')]);
    system(['mv ',spm_file(postop_anat_fname,'prefix','y_','ext','nii'),' ',spm_file(postop_anat_fname,'prefix','ySeg_','ext','nii')]);    
end

%---------------------------------
% skullstrip T1, for better coregistration with PET
%---------------------------------
if any(options.steps==3) 
    job_imcalc( spm_file(preop_anat_fname,'prefix','','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','c1','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','c2','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','c3','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','brain_','ext','nii'),...
                'i1.*((i2+i3+i4)>0.5)');
    job_imcalc( spm_file(postop_anat_fname,'prefix','','ext','nii'),...
                spm_file(postop_anat_fname,'prefix','c1','ext','nii'),...
                spm_file(postop_anat_fname,'prefix','c2','ext','nii'),...
                spm_file(postop_anat_fname,'prefix','c3','ext','nii'),...
                spm_file(postop_anat_fname,'prefix','brain_','ext','nii'),...
                'i1.*((i2+i3+i4)>0.5)');        
end


% make mask 1
if any(options.steps==4) 
    job_imcalc( spm_file(preop_anat_fname,'prefix','c1','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','c2','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','c1c2_','ext','nii'),...
                'i1+i2');  
    job_imcalc( spm_file(postop_anat_fname,'prefix','c1','ext','nii'),...
                spm_file(postop_anat_fname,'prefix','c2','ext','nii'),...
                spm_file(postop_anat_fname,'prefix','c1c2_','ext','nii'),...
                'i1+i2');   
end

%---------------------------------
% Linear Registration
%---------------------------------
% register postop MRI to preop for visualisation purposes only 
% comparing with LRpre_postop_anat helps to see how much the longitudinal registration rescued any
% deformation of the tissue around the resection cavity.
% UNDER CONSTRUCTION
if any(options.steps==10) 
    disp(['Running step: Linear Registration']);

end

%---------------------------------
% Longitudinal Registration
%---------------------------------
if any(options.steps==5.1) % register postop MRI to preop  
    disp(['Running step: Longitudinal Registration']);

    matlabbatch{1}.spm.tools.longit{1}.pairwise.vols1 = {spm_file(preop_anat_fname,'prefix','','ext','nii')};
    matlabbatch{1}.spm.tools.longit{1}.pairwise.vols2 = {spm_file(postop_anat_fname,'prefix','','ext','nii')};
    matlabbatch{1}.spm.tools.longit{1}.pairwise.wparam = options.LRparams;
    matlabbatch{1}.spm.tools.longit{1}.pairwise.write_def = 1;
    spm_jobman('run',matlabbatch); clear matlabbatch;

    % rename longitudinal registration defotrmation fields
    system(['mv ',spm_file(preop_anat_fname,'prefix','y_','ext','nii'),' ',spm_file(preop_anat_fname,'prefix','y_LR_','ext','nii')]);   % has to be y_, cant have yLR_ 
    system(['mv ',spm_file(postop_anat_fname,'prefix','y_','ext','nii'),' ',spm_file(postop_anat_fname,'prefix','y_LR_','ext','nii')]);

end

% Longitudinal Registration: deform
if any(options.steps==5.2) % register postop MRI to preop  

    disp(['Running step: Longitudinal Registration Deform']);

    images2LRpre = {    spm_file(postop_anat_fname,'prefix','','ext','nii');...
                        spm_file(postop_anat_fname,'prefix','c1','ext','nii');...
                        spm_file(postop_anat_fname,'prefix','c2','ext','nii');...
                        spm_file(postop_anat_fname,'prefix','c1c2_','ext','nii')}                                             
    if any(contains(in_img,'postop_rsctman_fname')); images2LRpre{end+1} =  postop_rsctman_fname ; end % include the manual segmentation mask if provided

    % resample postop MRI to preop   
    matlabbatch{1}.spm.util.defs.comp{1}.def = {spm_file(postop_anat_fname,'prefix','y_LR_','ext','nii')}; 
    matlabbatch{1}.spm.util.defs.comp{2}.inv.comp{1}.def = {spm_file(preop_anat_fname,'prefix','y_LR_','ext','nii')};
    matlabbatch{1}.spm.util.defs.comp{2}.inv.space = {spm_file(preop_anat_fname,'prefix','','ext','nii')};
    %matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {spm_file(postop_anat_fname,'prefix','c_','ext','nii');};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = cellstr(images2LRpre);
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = options.LRprefix;
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    spm_jobman('run',matlabbatch); clear matlabbatch;

end



%=====================================
% Make resection mask
%=====================================

%---------------------------------
% threshold GMWM masks and subtract
%---------------------------------
if options.RRRprefix=='thrsub'
    resected_tissue_fname = spm_file(postop_anat_fname,'prefix',[options.RRRprefix,'c1c2loss_','thr_',options.LRprefix,'c1c2_'],'ext','nii');
elseif options.RRRprefix=="subthr"
    resected_tissue_fname = spm_file(postop_anat_fname,'prefix',[options.RRRprefix,'thr_','c1c2loss_',options.LRprefix,'c1c2_'],'ext','nii');
end

if any(options.steps==6)

    if options.RRRprefix=='thrsub'
        %threshold preop
        job_imcalc( spm_file(preop_anat_fname,'prefix','c1c2_','ext','nii'),...
                    spm_file(preop_anat_fname,'prefix','thr_c1c2_','ext','nii'),...
                    'i1>0.1');
        %threshold postop
        job_imcalc( spm_file(postop_anat_fname,'prefix',[options.LRprefix,'c1c2_'],'ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',['thr_',options.LRprefix,'c1c2_'],'ext','nii'),...
                    'i1>0.1'); 
        %subtract
        job_imcalc( spm_file(preop_anat_fname,'prefix','thr_c1c2_','ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',['thr_',options.LRprefix,'c1c2_'],'ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',['c1c2loss_','thr_',options.LRprefix,'c1c2_'],'ext','nii'),...
                    'i1-i2'); 
        %retain positive         
        job_imcalc( spm_file(postop_anat_fname,'prefix',['c1c2loss_','thr_',options.LRprefix,'c1c2_'],'ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',[options.RRRprefix,'c1c2loss_','thr_',options.LRprefix,'c1c2_'],'ext','nii'),...
                    'i1>0');    

    elseif options.RRRprefix=='subthr'
        %subtract


        job_imcalc( spm_file(preop_anat_fname,'prefix','c1c2_','ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',[options.LRprefix,'c1c2_'],'ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',['c1c2loss_',options.LRprefix,'c1c2_'],'ext','nii'),...
                    'i1-i2'); 

        %threshold 
        job_imcalc( spm_file(postop_anat_fname,'prefix',['c1c2loss_',options.LRprefix,'c1c2_'],'ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',['thr_','c1c2loss_',options.LRprefix,'c1c2_'],'ext','nii'),...
                    'i1>0.1'); 

        % retain positive 
        job_imcalc( spm_file(postop_anat_fname,'prefix',['thr_','c1c2loss_',options.LRprefix,'c1c2_'],'ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',[options.RRRprefix,'thr_','c1c2loss_',options.LRprefix,'c1c2_'],'ext','nii'),...
                    'i1>0');     
    end
end

%---------------------------------
% clean resected tissue
%---------------------------------
for i_clus=1:options.nclus
    resected_tissue_clean_fname{i_clus} = spm_file(resected_tissue_fname,'prefix',['chop',num2str(i_clus),'_d',num2str(options.ndilate),'_','clus',num2str(i_clus),'_e',num2str(options.nerode),'_'],'ext','nii');
end
if any(options.steps==7)
    disp(['Running step 7: Erode Resected Tissue']);
    for i_clus=1:options.nclus
        job_dilate_select_erode(resected_tissue_fname,resected_tissue_clean_fname{i_clus},options.nerode,options.ndilate,i_clus);  
    end
end




%=====================================
% Match up PET with MRI
%=====================================
%job_processing_matchPET2MRI(options, preop_anat_fname, preop_pet_fname);
if any(options.steps==8)
    if isfile(preop_pet_fname)
        disp(['Running step 8: Coreg']);

        system(['cp ',spm_file(preop_pet_fname,'prefix','','ext','nii'),' ',spm_file(preop_pet_fname,'prefix','tmp_','ext','nii')]); % make a copy because spm coreg likes to overwrite the original image 
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {spm_file(preop_anat_fname,'prefix','brain_','ext','nii')};
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {spm_file(preop_pet_fname,'prefix','tmp_','ext','nii')};
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'cpre_';
        spm_jobman('run',matlabbatch); clear matlabbatch;

        system(['mv ',spm_file(preop_pet_fname,'prefix','cpre_tmp_','ext','nii'),' ',spm_file(preop_pet_fname,'prefix','cpre_','ext','nii')]);

        % get the ICV extracted coregistered PET
        job_imcalc(spm_file(preop_pet_fname,'prefix','cpre_','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','c1','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','c2','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','c3','ext','nii'),...
                spm_file(preop_pet_fname,'prefix','brain_cpre_','ext','nii'),...
                'i1.*((i2+i3+i4)>0.5)');
    end
end

%=====================================
% Move all images in preop space to MNI space
%=====================================

files_to_mni = {};
    
if any(contains(in_img,'preop_anat_fname'));    files_to_mni{end+1} = spm_file(preop_anat_fname,'prefix','','ext','nii'); end
if any(contains(in_img,'postop_anat_fname'));   files_to_mni{end+1} = spm_file(postop_anat_fname,'prefix',options.LRprefix,'ext','nii'); end
if any(contains(in_img,'preop_pet_fname'));     files_to_mni{end+1} = spm_file(preop_pet_fname,'prefix','cpre_','ext','nii'); files_to_mni{end+1} = spm_file(preop_pet_fname,'prefix','brain_cpre_','ext','nii'); end
if any(contains(in_img,'postop_rsctman_fname'));files_to_mni{end+1} = spm_file(postop_rsctman_fname,'prefix',options.LRprefix,'ext','nii'); end

if any(options.steps==9.1)
    disp(['Running step 9.1: Normalise Write']);
    matlabbatch{1}.spm.spatial.normalise.write.subj.def      = cellstr(spm_file(preop_anat_fname,'prefix','ySeg_','ext','nii'));
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = files_to_mni';
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    spm_jobman('run',matlabbatch); clear matlabbatch;

    % do nearest neighbour interpolation for the resection region image
    matlabbatch{1}.spm.spatial.normalise.write.subj.def      = cellstr(spm_file(preop_anat_fname,'prefix','ySeg_','ext','nii'));
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = resected_tissue_clean_fname';
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
    spm_jobman('run',matlabbatch); clear matlabbatch;
end

%=====================================
% Move all relevant postop images to MNI space (postop T1 and manual resection cavity) 
% Note: this uses SPMs normalise function, and does not involve the
% longitudinal registration to account for postoperative tissue deformation
%=====================================

if any(options.steps==9.2)

    disp(['Running step 9.2: Normalise Write']);
    matlabbatch{1}.spm.spatial.normalise.write.subj.def      = cellstr(spm_file(postop_anat_fname,'prefix','ySeg_','ext','nii'));
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {spm_file(postop_anat_fname,'prefix','','ext','nii')};     
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    spm_jobman('run',matlabbatch); clear matlabbatch;

    if any(contains(in_img,'postop_rsctman_fname'))    
        matlabbatch{1}.spm.spatial.normalise.write.subj.def      = cellstr(spm_file(postop_anat_fname,'prefix','ySeg_','ext','nii'));
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(postop_rsctman_fname);
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
        spm_jobman('run',matlabbatch); clear matlabbatch;
    end
end

%=====================================
% Move resection cavity image to postop space, in case this is needed
% for comparison with other algorithms etc.
%=====================================

% Longitudinal Registration: deform
if any(options.steps==5.3) % register preop MRI to postop space  
    disp(['Running step: Longitudinal Registration Deform']);

    % resample postop MRI to preop   
    matlabbatch{1}.spm.util.defs.comp{1}.def = {spm_file(preop_anat_fname,'prefix','y_LR_','ext','nii')}; 
    matlabbatch{1}.spm.util.defs.comp{2}.inv.comp{1}.def = {spm_file(postop_anat_fname,'prefix','y_LR_','ext','nii')};
    matlabbatch{1}.spm.util.defs.comp{2}.inv.space = {spm_file(postop_anat_fname,'prefix','','ext','nii')};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = resected_tissue_clean_fname';
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'LRpost_';
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    spm_jobman('run',matlabbatch); clear matlabbatch;
    
    %threshold. The longitudinal registration algorithm doesnt allow
    %nearest neighbout interpolation, so threshold mask
    for i_clus=1:options.nclus
        job_imcalc( spm_file(resected_tissue_clean_fname{i_clus},'prefix','LRpost_','ext','nii'),...
                    spm_file(resected_tissue_clean_fname{i_clus},'prefix','thr_LRpost_','ext','nii'),...
                    'i1>0.5'); 
    end

end

t_sub = toc(tstart_sub);
disp(['time = ',num2str(t_sub/60)]);

end % end of function


%=====================================
% Job functions
%=====================================

function job_auto_reorient(img_fname,acq) 

    addpath(fullfile(spm('Dir'),'toolbox/OldNorm'));

    if strcmp(acq,'anat')
        vg = spm_vol(fullfile(spm('Dir'),'canonical','avg152T1.nii'));
    elseif strcmp(acq,'pet')
        vg = spm_vol(fullfile(spm('Dir'),'toolbox/OldNorm','PET.nii'));
    else
        disp('acq not recognised, must be anat or pet');
        error;
        return;
    end

    flags.regtype='rigid';
    out_fname = spm_file(img_fname,'prefix','r_','ext','nii');
    system(['cp ',img_fname,' ',out_fname]);
    %create temp smoothed image
    spm_smooth(img_fname,'temp.nii',[12 12 12]);
    vf=spm_vol('temp.nii');
    [M,scal] = spm_affreg(vg,vf,flags); %vg=template, vf=smoothed input file
    M3=M(1:3,1:3);
    [u s v]=svd(M3);
    M3=u*v';
    M(1:3,1:3)=M3;
    N=nifti(out_fname);
    N.mat=M*N.mat;
    create(N);
    system(['rm temp.nii']);

end
    

function job_cat_vol_set_com(img_fname,out_fname, varargin)
% use center-of-mass (COM) to roughly correct for differences in the
% position between image and template.
% This script changes the affine matrix stored in the image header.
% If a mat file is specified as the third argument, the code will simply
% replace the affine matrix with the one supplied.

    copyfile(img_fname,out_fname);
    V = spm_vol(out_fname);

    if ~isempty(varargin)  

        Affine = load(varargin{1});
        Affine = Affine.Affine;

    else

        % pre-estimated COM of MNI template
        com_reference = [0 -20 -15];

        fprintf('Correct center-of-mass for %s\n',V.fname);
        Affine = eye(4);
        vol = spm_read_vols(V);
        avg = mean(vol(:));
        avg = mean(vol(vol>avg));

        % don't use background values
        [x,y,z] = ind2sub(size(vol),find(vol>avg));
        com = V.mat(1:3,:)*[mean(x) mean(y) mean(z) 1]';
        com = com';

        Affine(1:3,4) = (com - com_reference)';
        % save affine to apply to other images (particularly the manual segmentation in the postop space)
        save(spm_file(out_fname,'ext','mat'),'Affine');

    end

    M = spm_get_space(V.fname);
    spm_get_space(V.fname,Affine\M); % this line overwrites the affine matrix

end



function job_imcalc(varargin)
    
    for i_img=1:nargin-2
        img_input{i_img,:} = varargin{i_img};    
    end

    
    [in_dir,in_name,in_ext] = fileparts(varargin{1});
    [out_dir,out_name,out_ext] = fileparts(varargin{end-1});
    expr = varargin{end};

    matlabbatch{1}.spm.util.imcalc.input = img_input;
    matlabbatch{1}.spm.util.imcalc.output = out_name;
    matlabbatch{1}.spm.util.imcalc.outdir = {out_dir};
    matlabbatch{1}.spm.util.imcalc.expression = expr;
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run', matlabbatch);



end


function job_dilate_select_erode(in_fname,out_fname,ne,nd,i_clus)

    kernel = cat(3,[0 0 0; 0 1 0; 0 0 0],[0 1 0; 1 1 1; 0 1 0],[0 0 0; 0 1 0; 0 0 0]);
    
    %------------------------------
    % read image
    %------------------------------
    fname = in_fname;
    in_V = spm_vol(fname);
    in_data_orig = spm_read_vols(in_V);
    in_data = in_data_orig;
    
    %------------------------------
    % Dilate
    %------------------------------
    fname = spm_file(fname,'prefix',['e',num2str(ne),'_']);
    
    % erode
    for ie = 1:ne
        in_data = imerode(in_data,kernel);
    end

    % save image
    out_V = in_V;
    out_V.fname = fname;
    out_V.private.dat.fname = fname;
    out_V = spm_write_vol(out_V,in_data);

    %------------------------------
    % Select
    %------------------------------
    fname = spm_file(fname,'prefix',['clus',num2str(i_clus),'_']);
    in_data = job_getclusters(in_data,i_clus);    
    % save image
    out_V = in_V;
    out_V.fname = fname;
    out_V.private.dat.fname = fname;
    out_V = spm_write_vol(out_V,in_data);

    %------------------------------
    % Dilate
    %------------------------------
    fname = spm_file(fname,'prefix',['d',num2str(nd),'_']);
    % erode
    for id = 1:nd
        in_data = imdilate(in_data,kernel);
    end
    % save image
    out_V = in_V;
    out_V.fname = fname;
    out_V.private.dat.fname = fname;
    out_V = spm_write_vol(out_V,in_data);

    %------------------------------
    % Mask with original subtraction image
    %------------------------------   
    %fname = spm_file(fname,'prefix',['cut','_'],'ext','nii');
    fname = out_fname;
    in_data = in_data .* in_data_orig;
    % save image
    out_V = in_V;
    out_V.fname = fname;
    out_V.private.dat.fname = fname;
    out_V = spm_write_vol(out_V,in_data);
   

end

function result = job_getclusters(image,i_clus)
    % function result = extentThreshold(image, indices, k)
    % Applies the extent threshold function of SPM5 on 'image' with given
    % 'indices' and cluster size 'k' 

    indices = find(image);
    [x, y, z] = ind2sub(size(image), indices);
    XYZ = [x y z];

    % from spm_getSPM line 676:

    %-Calculate extent threshold filtering (from spm_getSPM, line 676)
    %-------------------------------------------------------------------
    A     = spm_clusters(XYZ');

    tbl = tabulate(A);
    clusterindex = tbl(:,1);
    try
        clustersize = tbl(:,2);
    catch ME
        if length(tbl)==0
            msg=['no clusters found'];
            causeException = MException('MATLAB:myCode:clusters',msg);
            ME = addCause(ME,causeException);
        end
        rethrow(ME);
        return;
    end
    clusterpc = tbl(:,3);
    [clustersize_sort,inds_sort] = sort(clustersize,'descend');
    largestclusterindex = inds_sort(1);

    
    clusters_keep = inds_sort(i_clus);
    %{
    %if cluster is at least 0.3 of next largest cluster, keep.
    for i_ind = 1:length(clustersize_sort)-1
        if clustersize_sort(i_ind+1) > 0.3*clustersize_sort(i_ind)
            clusters_keep = [clusters_keep,inds_sort(i_ind+1)];
        else
            break;
        end
    end
    %}
    inds_keep = find(ismember(A,clusters_keep));

    % ...eliminate voxels
    %-------------------------------------------------------------------
    XYZ   = XYZ(inds_keep,:);

    result = zeros(size(image));
    inds = sub2ind(size(image), XYZ(:,1), XYZ(:,2), XYZ(:,3));
    result(inds) = image(inds);
end

