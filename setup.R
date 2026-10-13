# Run from the folder containing app.R, before uploading to GitHub.
options(repos = c(CRAN = "https://cloud.r-project.org"))
packages <- c("shiny", "shinydashboard", "echarts4r", "readr", "mlr3",
              "mlr3learners", "xgboost", "kernelshap", "shapviz", "ggplot2",
              "rsconnect")
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) install.packages(missing)
stopifnot(file.exists("app.R"))
parse("app.R")
rsconnect::writeManifest(appDir = ".", appFiles = "app.R")
message("manifest.json created. Test with shiny::runApp('.') and click RUN.")
