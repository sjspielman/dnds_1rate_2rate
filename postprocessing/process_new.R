require(readr)  # for read_csv()
require(dplyr)  # for mutate()
require(tidyr)  # for unnest()
require(purrr)  # for map(), reduce()

result.directory <- "/Users/sjspielman/Research/github/dnds_1rate_2rate/results/balancedtrees_results"

files <- dir(path = result.directory, pattern = "*FUBAR1.txt")
data <- data_frame(filename = files) %>%
mutate(file_contents = map(filename,          # read files into
           ~ read_csv(file.path(result.directory, .))) # a new data column
        )
data
