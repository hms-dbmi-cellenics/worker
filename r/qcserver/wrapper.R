# dummy object for placeholder
count_matrix <- 1

wrapper <- function(filter_name) {
    print(paste(count_matrix, "placeholder wrapper function, running filter", filter_name))
    count_matrix <<- count_matrix + 1
}
