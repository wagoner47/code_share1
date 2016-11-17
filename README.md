# code_share1
Code for first group code share

# To run on calvin with basic options
(Replace save_dir with the directory in which to save the output)

"python 3d_zbins_los_pert.py /calvin1/wagoner47/halo_catalogs/buzzard/v1.1/split_randoms/Buzzard_v1.1_halos_log10M200B_13.5_scatter0.01_randoms_f10_split2-1_seed0.fit /calvin1/wagoner47/halo_catalogs/buzzard/v1.1/split_randoms/Buzzard_v1.1_halos_log10M200B_13.5_scatter0.01_randoms_f10_split2-2_seed0.fit (save_dir) buzzard_cosmology.ini"

## To make a plot
Include "--plot" followed by the directory in which to save the plot.

## Other options I find helpful
"--mail_options" followed by ini file with username, password, etc. for mailing: Send an email when code finishes running

"--nohup" followed by file used for nohup output: Attach the nohup file to the notice email

## Remaining options
(Fill this in later)
