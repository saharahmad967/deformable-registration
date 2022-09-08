# Deformable Registration
Dynamic elasticity model based deformable image registration

# System Requirements
- Software: MATLAB 2019b or higher
- Memory: 16 GB RAM
- OS: Windows/Mac/Linux

# Demo
Example data are provided in the subfolders

- Run DEM_main.m in MATLAB
- Results (warped tissue segmentation map and displacement field) will be saved in the 'output' subfolder
- Code takes approximately 1 - 2 hours to finish, depending on the system

# Results Visualization
- Open moving/fixed/warped tissue segmentation maps in ITK-SNAP

# Instructions for Usage with New Data
1) Copy fixed tissue segmentation map (nifti format) to 'fixed' subfolder
2) Copy moving tissue segmentation map (nifti format), affine-aligned with fixed tissue map, to 'moving' subfolder
3) Run DEM_main.m 
4) Displacement field and warped tissue segmentation map will be saved in the 'output' subfolder
