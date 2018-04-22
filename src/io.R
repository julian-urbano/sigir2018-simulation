save.object <- function(obj, file) {
  dir <- dirname(file)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  save(file = file, obj)
}

load.object <- function(file) {
  env <- new.env()
  name <- load(file = file, envir = env)
  obj <- get(name, envir = env)
  return(obj)
}
