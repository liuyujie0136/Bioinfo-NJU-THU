## Use R package "conflicted" to detect conflicted functions in different R packages

install.packages("conflicted")

library(conflicted)


## Example

library(dplyr)
filter(mtcars, am & cyl == 8)

# ERROR: [conflicted] `filter` found in 2 packages.
# Either pick the one you want with `::` 
# * dplyr::filter
# * stats::filter
# Or declare a preference with `conflict_prefer()`
# * conflict_prefer("filter", "dplyr")
# * conflict_prefer("filter", "stats")
# Run `rlang::last_error()` to see where the error occurred.


## OR, add the following line in your ~/.Rprofile
require(conflicted)
