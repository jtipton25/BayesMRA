# schedule statszilla server jobs
library(here)
library(rstudioapi)

jobRunScript(here::here("scripts", "test_fit_resolution.R"))
jobRunScript(here::here("scripts", "test_fit.R"))
jobRunScript(here::here("scripts", "test_fit-integrated.R"))
jobRunScript(here::here("scripts", "heaton-fit.R"))
