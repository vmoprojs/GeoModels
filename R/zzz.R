.onAttach <- function(lib, pkg) {
 
  pkg.info <- drop(read.dcf(file = system.file("DESCRIPTION", package = "GeoModels"),
                            fields = c("Title", "Version", "Date")))
  
  # Startup message
  packageStartupMessage(paste("--------------------------------------------------------------\n",
                              pkg.info["Title"], "\n",
                              "For an Introduction to GeoModels, go to https://vmoprojs.github.io/GeoModels-page/\n",
                              paste("GeoModels version ", pkg.info["Version"],
                                    " (built on ", pkg.info["Date"], ") is now loaded\n", sep = ""),
                              "--------------------------------------------------------------\n"))

}