
code_path=''
slurm_path=$code_path/slurm_jobs
site="S4"
list2run='text file with subject names'

for subject_name in `cat $list2run`
do
	echo $subject_name
    
  preop_anat_fname=`ls $in_path/$subject_name/preop/anat/S*T1w.nii`
  postop_anat_fname=`ls $in_path/$subject_name/postop/anat/S*T1w.nii` # optional	
  preop_pet_fname=`ls $in_path/$subject_name/preop/pet/S*pet.nii` # optional

  rm ${slurm_path}/slurm_${subject_name}.job

  echo '#!/bin/bash' >> ${slurm_path}/slurm_${subject_name}.job
  echo "#SBATCH --job-name=${subject_name}" >> ${slurm_path}/slurm_${subject_name}.job
  echo "#SBATCH --account=wr47" >> ${slurm_path}/slurm_${subject_name}.job
  echo "#SBATCH --time=3:00:00" >> ${slurm_path}/slurm_${subject_name}.job
  echo "#SBATCH --ntasks=1" >> ${slurm_path}/slurm_${subject_name}.job
  echo "#SBATCH --mem-per-cpu=32000" >> ${slurm_path}/slurm_${subject_name}.job
  echo "#SBATCH --cpus-per-task=1" >> ${slurm_path}/slurm_${subject_name}.job
  echo "module load spm12" >> ${slurm_path}/slurm_${subject_name}.job

  # MRI preop + MRI postop + PET
  echo "matlab -nosplash -singleCompThread -r \"addpath(genpath('${code_path}')) ;  EPIC_CHOP('${preop_anat_fname}','${postop_anat_fname}','preop_pet_fname','${preop_pet_fname}','out_path','${out_path}','subject_name','${subject_name}') ; exit \" " >> ${slurm_path}/slurm_${subject_name}.job

  # MRI preop + MRI postop + manual
  #echo "matlab -nosplash -singleCompThread -r \"addpath(genpath('${code_path}')) ;  EPIC_CHOP('${preop_anat_fname}','${postop_anat_fname}','postop_rsctman_fname','${postop_rsctman_fname}','out_path','${out_path}','subject_name','${subject_name}') ; exit \" " >> ${slurm_path}/slurm_${subject_name}.job

  sbatch ${slurm_path}/slurm_${subject_name}.job

done
