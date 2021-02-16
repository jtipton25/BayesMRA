#' Temperature data from satellite measurements
#'
#' A dataset containing the temperature measurements and locations
#' from satellite data used in Heaton et al, A Case Study
#' Competition Among Methods for Analyzing Large Spatial Data,
#' Journal of Agricultural, Biological and Environmental Statistics,
#' 2018: http://link.springer.com/article/10.1007/s13253-018-00348-w
#'
#' @docType data
#'
#' @usage data("all_sat_temps")
#'
#' @format A data frame with 150000 rows and 4 variables:
#' \describe{
#' \item{Lon}{The longitude of the obseravtion}
#' \item{Lat}{The latitude of the obseravtion}
#' \item{MaskTemp}{The temperature observations with a mask that is used for training the model}
#' \item{TrueTemp}{The temperature observations without a mask that is used for evaluating the model}
#' }
#'
#'
#' @references Heaton et al, A Case Study
#' Competition Among Methods for Analyzing Large Spatial Data,
#' Journal of Agricultural, Biological and Environmental Statistics,
#'  2018:
#'
#' (\href{http://link.springer.com/article/10.1007/s13253-018-00348-w}{JABES})
#'
#' @source \href{https://github.com/finnlindgren/heatoncomparison}{GitHub}
#'
#' @examples
#' library(BayesMRA)
#' library(tidyverse)
#' data("all_sat_temps")
#' ggplot(all_sat_temps, aes(x = Lon, y = Lat, fill = TrueTemp)) +
#'  geom_raster()
"all_sat_temps"



#' Temperature data from simulation
#'
#' A dataset containing the temperature measurements and locations
#' from simulated data used in Heaton et al, A Case Study
#' Competition Among Methods for Analyzing Large Spatial Data,
#' Journal of Agricultural, Biological and Environmental Statistics,
#' 2018: http://link.springer.com/article/10.1007/s13253-018-00348-w
#'
#' @docType data
#'
#' @usage data("all_sim_data")
#'
#' @format A data frame with 150000 rows and 4 variables:
#' \describe{
#' \item{Lon}{The longitude of the obseravtion}
#' \item{Lat}{The latitude of the obseravtion}
#' \item{MaskTemp}{The simulated temperature observations with a mask that is used for training the model}
#' \item{TrueTemp}{The simulated temperature observations without a mask that is used for evaluating the model}
#' }
#'
#'
#' @references Heaton et al, A Case Study
#' Competition Among Methods for Analyzing Large Spatial Data,
#' Journal of Agricultural, Biological and Environmental Statistics,
#'  2018:
#'
#' (\href{http://link.springer.com/article/10.1007/s13253-018-00348-w}{JABES})
#'
#' @source \href{https://github.com/finnlindgren/heatoncomparison}{GitHub}
#'
#' @examples
#' library(BayesMRA)
#' library(tidyverse)
#' data("all_sim_data")
#' ggplot(all_sim_data, aes(x = Lon, y = Lat, fill = TrueTemp)) +
#'  geom_raster()
"all_sim_data"


#' Small temperature data from simulation for testing
#'
#' A dataset containing the temperature measurements and locations
#' from simulated data used in Heaton et al, A Case Study
#' Competition Among Methods for Analyzing Large Spatial Data,
#' Journal of Agricultural, Biological and Environmental Statistics,
#' 2018: http://link.springer.com/article/10.1007/s13253-018-00348-w
#'
#' @docType data
#'
#' @usage data("code_test")
#'
#' @format A data frame with 10000 rows and 4 variables:
#' \describe{
#' \item{Lon}{The longitude of the obseravtion}
#' \item{Lat}{The latitude of the obseravtion}
#' \item{MaskTemp}{The simulated temperature observations with a mask that is used for training the model}
#' \item{TrueTemp}{The simulated temperature observations without a mask that is used for evaluating the model}
#' }
#'
#'
#' @references Heaton et al, A Case Study
#' Competition Among Methods for Analyzing Large Spatial Data,
#' Journal of Agricultural, Biological and Environmental Statistics,
#'  2018:
#'
#' (\href{http://link.springer.com/article/10.1007/s13253-018-00348-w}{JABES})
#'
#' @source \href{https://github.com/finnlindgren/heatoncomparison}{GitHub}
#'
#' @examples
#' library(BayesMRA)
#' library(tidyverse)
#' data("code_test")
#' ggplot(code_test, aes(x = Lon, y = Lat, fill = TrueTemp)) +
#'  geom_raster()
"code_test"
