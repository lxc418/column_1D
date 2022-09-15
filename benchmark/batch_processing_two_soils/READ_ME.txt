SEQUENCE FOR .sh RUNING

1.run_simu (sbatch test.script)
(create subfolders;copy inp files to each subfolder;
 replace given variable;create inp file;start the simu)

2.run_check_end
(check the if the simu crashes by checking the end of .smy;
results store in check.csv; need to delete .csv before a new checking)

3.run_plot_figure (sbatch post_process.script)
(copy postprocess files to each subfolder;
 replace given variable;create postprocess file;start post)

4.(optional)run_copy_files
(copy the specified file to the specified folder)

5.(optional)run_print_salt_thickness (sbatch salt_thickness.script)
(run sl_salt_thickness.m;
store the solid salt thicknesses in M.xlsx)