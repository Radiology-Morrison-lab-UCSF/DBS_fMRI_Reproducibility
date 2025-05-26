import SimpleITK as sitk
import os
import sys

dir_script = os.path.dirname(os.path.abspath(__file__))


loc_atlas_gm = os.path.join(dir_script, "grey.nii.gz")
loc_atlas_wm = os.path.join(dir_script, "white.nii.gz")
loc_atlas_csf = os.path.join(dir_script, "csf.nii.gz")
loc_exclude = os.path.join(dir_script, "exclude.nii.gz")


def GetInput(subj:str, condition:str, repeat:str):

    loc =  os.path.join(source_dir,subj, "wart_mean_au" + subj + "_fmri_"+  condition+"_" + repeat+ ".nii")

    if not os.path.exists(loc):
        print("Could not find input:", loc)
        exit(1)

    return sitk.ReadImage(loc,sitk.sitkFloat32)


def ProcessSubj(subj):

    dir_out = os.path.join(dest_dir, subj)
    loc_out = os.path.join(dir_out, "lead-seg-otsu.nii.gz")

    if not os.path.exists(dir_out):
        os.mkdir(dir_out)

    mean = GetInput(subj, "off", "1")
    mask = ((sitk.ReadImage(loc_atlas_gm) + sitk.ReadImage(loc_atlas_wm)) > 0.5)

    mask = sitk.BinaryErode(sitk.BinaryErode(mask))

    mask = mask * (sitk.ReadImage(loc_exclude) == 0)

    thresholded = sitk.OtsuThreshold(mean)

    thresholded = thresholded * mask

    sitk.WriteImage(thresholded, loc_out)

    print("Done for ", subj, ". File at ", loc_out)


if len(sys.argv) < 3 or sys.argv[1].startswith("-h"):
    print("Runs otsu thresholding on an fmri image and crops the result")
    print("Use as:")
    print("python ./main.py <source_dir> <destination_dir>")
    print("source_dir should have one directory per subject")
    exit(1)


source_dir = sys.argv[1]
dest_dir = sys.argv[2]


if not os.path.exists(dest_dir):
    os.mkdir(dest_dir)

for subj in os.listdir(source_dir):
    if os.path.isdir(os.path.join(source_dir,subj)):
        try:
            ProcessSubj(subj)
        except:
            print("FAILED", subj)
