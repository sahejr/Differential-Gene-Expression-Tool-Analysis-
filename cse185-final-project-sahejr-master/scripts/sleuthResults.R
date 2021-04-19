library("sleuth")

# Set up the paths to our kallisto results
sample_id = c("SRX642051", "SRX642051_r", "SRX642055", "SRX642055_r")
kal_dirs = file.path(sample_id)

# Load metadata
s2c = read.table(file.path("/home/linux/ieng6/cs185s/sdrandha/project/exp_info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c = dplyr::mutate(s2c, path = kal_dirs)

# Create sleuth object
so = sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

# Fit each model and test
so = sleuth_fit(so, ~condition, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')

# Get output, write results to file
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
# Note, you may need to edit the output path below if your $HOME
# directory is not the same as your CSE185 course directory
write.table(sleuth_significant, "final_sleuth_results.tab", sep="\t", quote=FALSE)
