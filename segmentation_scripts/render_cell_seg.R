# render_cell_seg.R

library(rmarkdown)

base_path <- "/restricted/projectnb/camplab/projects/20220504_Suzuki_LN/TumorIMC/steinbock"
out_dir   <- "/restricted/projectnb/camplab/projects/20220504_Suzuki_LN/TumorIMC/analysis/cell_seg_reports"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

patient_dirs <- list(
  p10  = file.path(base_path, "p10_results"),
  p11  = file.path(base_path, "p11_results"),
  p12  = file.path(base_path, "p12_results"),
  p13  = file.path(base_path, "p13_results"),
  p14  = file.path(base_path, "p14_results"),
  p15  = file.path(base_path, "p15_results"),
  p16  = file.path(base_path, "p16_results"),
  p17  = file.path(base_path, "p17_results"),
  p18  = file.path(base_path, "p18_results"),
  p19  = file.path(base_path, "p19_results"),
  p20a = file.path(base_path, "p20_4_25_results"),
  p20b = file.path(base_path, "p20_5_19_results")
)

for (pt_name in names(patient_dirs)) {
  cat("Rendering:", pt_name, "\n")
  
  rmarkdown::render(
    input       = "cell_segmentation.Rmd",
    params      = list(
      patient_id     = pt_name,
      steinbock_path = patient_dirs[[pt_name]]
    ),
    output_file = file.path(out_dir, paste0(pt_name, "_mask_visualization.pdf"))
  )
}