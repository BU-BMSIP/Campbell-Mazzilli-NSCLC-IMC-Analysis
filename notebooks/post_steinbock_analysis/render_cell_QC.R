library(rmarkdown)

output_dir <- "/restricted/projectnb/camplab/projects/20220504_Suzuki_LN/TumorIMC/analysis/distribution_reports/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

patient_ids <- as.character(11:19)
p20_dirs    <- c("p20_4_25_results", "p20_5_19_results")

for (pt in patient_ids) {
  message("Rendering Patient ", pt, "...")
  tryCatch({
    render(
      input       = "cell_distributions.Rmd",
      output_file = paste0("p", pt, "_cell_distributions.pdf"),
      output_dir  = output_dir,
      params      = list(patient_id = pt),
      envir       = new.env()
    )
  }, error = function(e) {
    message("Patient ", pt, " failed: ", e$message)
  })
}

for (d in p20_dirs) {
  message("Rendering ", d, "...")
  tryCatch({
    render(
      input       = "cell_distributions.Rmd",
      output_file = paste0(d, "_cell_distributions.pdf"),
      output_dir  = output_dir,
      params      = list(patient_id = d),
      envir       = new.env()
    )
  }, error = function(e) {
    message(d, " failed: ", e$message)
  })
}

message("Done. PDFs written to: ", output_dir)