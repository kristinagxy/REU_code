# REU_code
code for the project **Iterative Methods at Lower Precision** at Emory REU 2022

the modified codes for iterative methods at custom precision level from the IR Tools package -> [IRtools](https://github.com/kristinagxy/REU_code/tree/main/IRtools);
* chopped CGLS -> [IRcgls_chop](https://github.com/kristinagxy/REU_code/blob/main/IRtools/IRcgls_chop.m)
* CS -> [cs_o](https://github.com/kristinagxy/REU_code/blob/main/IRtools/cs_o.m)
* chopped CS -> [cs_original](https://github.com/kristinagxy/REU_code/blob/main/IRtools/cs_original.m)
* chopped CS with Tikhonov Regularization -> [cs_chop](https://github.com/kristinagxy/REU_code/blob/main/IRtools/cs_chop.m)

the modified basic operations for custom precision level -> [basic_operations](https://github.com/kristinagxy/REU_code/tree/main/basic_operations)

function for drawing error norms of the .mat files in a given directory, can choose an aspect of the problem to make comparisons ('precision','size','blurlevel','noise') -> [make_plots](https://github.com/kristinagxy/REU_code/blob/main/make_plots.m)

* Files must be saved with the format: problem name_precision level_size_blur level_noise leve_file type.mat
  (tomo_single_64_default_0.01_info.mat)

function for generating image from the solutions in a given directory -> [do_plots](https://github.com/kristinagxy/REU_code/blob/main/do_plots.m)

function for running CGLS without Tikhonov Regularization -> [run_cgls](https://github.com/kristinagxy/REU_code/blob/main/run_cgls.m)

function for running CGLS with Tikhonov Regularization -> [run_cgls_reg](https://github.com/kristinagxy/REU_code/blob/main/run_cgls_reg.m)

function for running CS with Tikhonov Regularization -> [run_cs](https://github.com/kristinagxy/REU_code/blob/main/run_cs.m)
