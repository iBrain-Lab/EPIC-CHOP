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
% postop_anat		- /out_path/subject_name/postop/anat/LRpre_[postop-anat-origname].nii
% resection_cavity	- /out_path/subject_name/postop/anat/chop1_d2_clus1_e1_thrsubp1p2loss_thr_LRpre_p1p2_[postop-anat-origname].nii
% preop_pet         - /out_path/subject_name/preop/anat/pet/cpre_[preop-pet-origname].nii
% postop_rsctman	- /out_path/subject_name/postop/anat/LRpre_manual_[postop-anat-origname].nii
% SPACE: postop (centred to center of mass)
% preop_anat		- /out_path/subject_name/preop/anat/
% postop_anat		- /out_path/subject_name/postop/anat/com_[postop-anat-origname].nii
% resection_cavity	- /out_path/subject_name/postop/anat/thr_LRpost_chop1_d2_clus1_e1_thrsubp1p2loss_thr_LRpre_p1p2_[postop-anat-origname].nii
% preop_pet         - NA
% postop_rsctman	- /out_path/subject_name/postop/anat/com_manual_[postop-anat-origname].nii
% SPACE: MNI
% preop_anat		- /out_path/subject_name/preop/anat/wcom_[preop-anat-origname].nii
% postop_anat		- /out_path/subject_name/postop/anat/wLRpre_[postop-anat-origname].nii
% resection_cavity	- /out_path/subject_name/postop/anat/wchop1_d2_clus1_e1_thrsubp1p2loss_thr_LRpre_p1p2_[postop-anat-origname].nii
% preop_pet         - /out_path/subject_name/preop/anat/pet/wcpre_[preop-pet-origname].nii
% postop_rsctman	- /out_path/subject_name/postop/anat/wLRpre_manual_[postop-anat-origname].nii
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
% 0 - make a copy of images to avoid overwriting (spm realign/spm coregister), and move them to a new directory. optional
% 1* - reorients via autoreoirent or center of mass (deprecated, replaced by
% cat12)
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
% v18: uses util.defs to move to MNI space rather than spatial.normalise
% v19: lets not speak about v19, except to say that Terry was most likely
% incorrect, but has not been yet proven so.
% v20: belatedly changes the segmentation and registration to cat12, rather
% than spm unified. Removes step 1, reorientation via centre of mass, thus
% requires manual approximate reorientation of PET to MRI
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
end
if ~isfield(options,'nclus')    
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

% the optional steps are done now, so we can abandon ImgStruct and use more
% interpretable variable names
for i_img=1:numel(in_img)
    eval( in_img{i_img} + "="  + '''' + ImgStruct.(in_img{i_img}).currentfname + '''' + ';'); 
end

%=====================================
% Process MRI
%=====================================  

%---------------------------------
% Segment using cat12
%---------------------------------
if any(options.steps==2)
    
    matlabbatch{1}.spm.tools.cat.estwrite.data = cellstr({spm_file(preop_anat_fname,'prefix','','ext','nii'); spm_file(postop_anat_fname,'prefix','','ext','nii')});
    matlabbatch{1}.spm.tools.cat.estwrite.data_wmh = {''};
    matlabbatch{1}.spm.tools.cat.estwrite.nproc = 6;
    matlabbatch{1}.spm.tools.cat.estwrite.useprior = '';
    matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {fullfile(spm('dir'),'tpm/TPM.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
    matlabbatch{1}.spm.tools.cat.estwrite.opts.biasacc = 0.5;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.optimal = [1 0.3];
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.setCOM = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1070;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.affmod = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.spm_kamap = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASmyostr = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 2;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.WMHC = 2;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = {fullfile(spm('dir'),'/toolbox/cat12/templates_MNI152NLin2009cAsym/Template_0_GS.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.bb = 12;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.SRP = 22;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.ignoreErrors = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.suit = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ownatlas = {''};
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.pp.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.pp.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.pp.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.las.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [1 1];
    matlabbatch{1}.spm.tools.cat.estwrite.output.rmat = 0;
    
    spm_jobman('run',matlabbatch); clear matlabbatch;
    
    for i_img=1:2
        system(['mv ',ImgStruct.(in_img{i_img}).outpath,'/mri/* ', ImgStruct.(in_img{i_img}).outpath]);
    end
    
    % rename deformation field to avoid name clash with longitudinal registration
    system(['mv ',spm_file(preop_anat_fname,'prefix','y_','ext','nii'),' ',spm_file(preop_anat_fname,'prefix','y_cat12_','ext','nii')]);
    system(['mv ',spm_file(postop_anat_fname,'prefix','y_','ext','nii'),' ',spm_file(postop_anat_fname,'prefix','y_cat12_','ext','nii')]);   
    system(['mv ',spm_file(preop_anat_fname,'prefix','iy_','ext','nii'),' ',spm_file(preop_anat_fname,'prefix','iy_cat12_','ext','nii')]);
    system(['mv ',spm_file(postop_anat_fname,'prefix','iy_','ext','nii'),' ',spm_file(postop_anat_fname,'prefix','iy_cat12_','ext','nii')]); 
end

%---------------------------------
% skullstrip T1, for better coregistration with PET
%---------------------------------
if any(options.steps==3) 
    job_imcalc( spm_file(preop_anat_fname,'prefix','','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','p1','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','p2','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','p3','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','brain_','ext','nii'),...
                'i1.*((i2+i3+i4)>0.5)');
    job_imcalc( spm_file(postop_anat_fname,'prefix','','ext','nii'),...
                spm_file(postop_anat_fname,'prefix','p1','ext','nii'),...
                spm_file(postop_anat_fname,'prefix','p2','ext','nii'),...
                spm_file(postop_anat_fname,'prefix','p3','ext','nii'),...
                spm_file(postop_anat_fname,'prefix','brain_','ext','nii'),...
                'i1.*((i2+i3+i4)>0.5)');        
end


% make mask 1
if any(options.steps==4) 
    job_imcalc( spm_file(preop_anat_fname,'prefix','p1','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','p2','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','p1p2_','ext','nii'),...
                'i1+i2');  
    job_imcalc( spm_file(postop_anat_fname,'prefix','p1','ext','nii'),...
                spm_file(postop_anat_fname,'prefix','p2','ext','nii'),...
                spm_file(postop_anat_fname,'prefix','p1p2_','ext','nii'),...
                'i1+i2');   
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

    % rename longitudinal registration deformation fields
    system(['mv ',spm_file(preop_anat_fname,'prefix','y_','ext','nii'),' ',spm_file(preop_anat_fname,'prefix','y_LR_','ext','nii')]);   % has to be y_, cant have yLR_ 
    system(['mv ',spm_file(postop_anat_fname,'prefix','y_','ext','nii'),' ',spm_file(postop_anat_fname,'prefix','y_LR_','ext','nii')]);

end

% Longitudinal Registration: deform
if any(options.steps==5.2) % register postop MRI to preop  

    disp(['Running step: Longitudinal Registration Deform']);

    images2LRpre = {    spm_file(postop_anat_fname,'prefix','','ext','nii');...
                        %spm_file(postop_anat_fname,'prefix','p1','ext','nii');...
                        %% if you want the GM and WM images in native space
                        %spm_file(postop_anat_fname,'prefix','p2','ext','nii');...
                        spm_file(postop_anat_fname,'prefix','p1p2_','ext','nii')}                                             
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
    resected_tissue_fname = spm_file(postop_anat_fname,'prefix',[options.RRRprefix,'p1p2loss_','thr_',options.LRprefix,'p1p2_'],'ext','nii');
elseif options.RRRprefix=="subthr"
    resected_tissue_fname = spm_file(postop_anat_fname,'prefix',[options.RRRprefix,'thr_','p1p2loss_',options.LRprefix,'p1p2_'],'ext','nii');
end

if any(options.steps==6)

    if options.RRRprefix=='thrsub'
        %threshold preop
        job_imcalc( spm_file(preop_anat_fname,'prefix','p1p2_','ext','nii'),...
                    spm_file(preop_anat_fname,'prefix','thr_p1p2_','ext','nii'),...
                    'i1>0.1');
        %threshold postop
        job_imcalc( spm_file(postop_anat_fname,'prefix',[options.LRprefix,'p1p2_'],'ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',['thr_',options.LRprefix,'p1p2_'],'ext','nii'),...
                    'i1>0.1'); 
        %subtract
        job_imcalc( spm_file(preop_anat_fname,'prefix','thr_p1p2_','ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',['thr_',options.LRprefix,'p1p2_'],'ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',['p1p2loss_','thr_',options.LRprefix,'p1p2_'],'ext','nii'),...
                    'i1-i2'); 
        %retain positive         
        job_imcalc( spm_file(postop_anat_fname,'prefix',['p1p2loss_','thr_',options.LRprefix,'p1p2_'],'ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',[options.RRRprefix,'p1p2loss_','thr_',options.LRprefix,'p1p2_'],'ext','nii'),...
                    'i1>0');    

    elseif options.RRRprefix=='subthr'
        %subtract


        job_imcalc( spm_file(preop_anat_fname,'prefix','p1p2_','ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',[options.LRprefix,'p1p2_'],'ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',['p1p2loss_',options.LRprefix,'p1p2_'],'ext','nii'),...
                    'i1-i2'); 

        %threshold 
        job_imcalc( spm_file(postop_anat_fname,'prefix',['p1p2loss_',options.LRprefix,'p1p2_'],'ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',['thr_','p1p2loss_',options.LRprefix,'p1p2_'],'ext','nii'),...
                    'i1>0.1'); 

        % retain positive 
        job_imcalc( spm_file(postop_anat_fname,'prefix',['thr_','p1p2loss_',options.LRprefix,'p1p2_'],'ext','nii'),...
                    spm_file(postop_anat_fname,'prefix',[options.RRRprefix,'thr_','p1p2loss_',options.LRprefix,'p1p2_'],'ext','nii'),...
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
                spm_file(preop_anat_fname,'prefix','p1','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','p2','ext','nii'),...
                spm_file(preop_anat_fname,'prefix','p3','ext','nii'),...
                spm_file(preop_pet_fname,'prefix','brain_cpre_','ext','nii'),...
                'i1.*((i2+i3+i4)>0.5)');
    end
end

%=====================================
% Move all images in preop space to MNI space
%=====================================

files_to_mni = {};  

if any(contains(in_img,'preop_anat_fname'));    files_to_mni{end+1} = spm_file(preop_anat_fname,'prefix','','ext','nii'); end
if any(contains(in_img,'preop_pet_fname'));     files_to_mni{end+1} = spm_file(preop_pet_fname,'prefix','cpre_','ext','nii'); files_to_mni{end+1} = spm_file(preop_pet_fname,'prefix','brain_cpre_','ext','nii'); end
if any(contains(in_img,'postop_anat_fname'));   files_to_mni{end+1} = spm_file(postop_anat_fname,'prefix',options.LRprefix,'ext','nii'); end
if any(contains(in_img,'postop_rsctman_fname'));files_to_mni{end+1} = spm_file(postop_rsctman_fname,'prefix',options.LRprefix,'ext','nii'); end
%files_to_mni{end+1} = spm_file(preop_anat_fname,'prefix','c1','ext','nii');


if any(options.steps==9.1)
    disp(['Running step 9.1: Normalise Write']);
    
    matlabbatch{1}.spm.util.defs.comp{1}.def = {spm_file(preop_anat_fname,'prefix','y_cat12_','ext','nii')}; 
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = files_to_mni';
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;   
    spm_jobman('run',matlabbatch); clear matlabbatch;

    
    % do nearest neighbour interpolation for the resection region image
    matlabbatch{1}.spm.util.defs.comp{1}.def = {spm_file(preop_anat_fname,'prefix','y_cat12_','ext','nii')}; 
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = resected_tissue_clean_fname';
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0; 
    spm_jobman('run',matlabbatch); clear matlabbatch;
    
end

%=====================================
% Move all relevant postop images directly to MNI space (postop T1 and manual resection cavity) 
% Note: this uses SPMs normalise function, and does not involve the
% longitudinal registration to account for postoperative tissue deformation
%=====================================

if any(options.steps==9.2)

    disp(['Running step 9.2: Normalise Write']);
    
    matlabbatch{1}.spm.util.defs.comp{1}.def = {spm_file(postop_anat_fname,'prefix','y_cat12_','ext','nii')}; 
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {spm_file(postop_anat_fname,'prefix','','ext','nii')};
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
    
    spm_jobman('run',matlabbatch); clear matlabbatch;

    if any(contains(in_img,'postop_rsctman_fname'))    
        
        matlabbatch{1}.spm.util.defs.comp{1}.def = {spm_file(postop_anat_fname,'prefix','y_cat12_','ext','nii')}; 
        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {spm_file(postop_rsctman_fname,'prefix','','ext','nii')};
        matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';
        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
        matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0; 

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

    % resample postop MRI to preop   WRONG
    % resample preop MRI to postop
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
