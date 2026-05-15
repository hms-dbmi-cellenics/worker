---
applyTo: "r/**/*.R"
description: "R code formatting and style for Cellenics worker: 80-char lines, pipe |>, snake_case, roxygen2 docs, no library()"
---

# R Code Formatting & Style

This instruction applies to all R files in the `r/` directory.

## Before Making Changes

**Always auto-format first**, then fix any remaining issues:

```bash
# In R console or terminal:
styler::style_file("path/to/file.R")  # Auto-format
lintr::lint("path/to/file.R")         # Identify issues
```

Fix any linting issues that `styler` doesn't handle automatically (see "Common Linting Issues" below).

## Style Rules

### Line Length: < 80 characters

Hard limit. Break long lines by splitting on open/close brackets:

```r
# ✓ GOOD: Break on opening bracket
config <- list(
  name = input$experimentName,
  samples = input$sampleIds,
  organism = input$organism,
  input = list(type = input$input$type),
  sampleOptions = input$sampleOptions
)

# ✗ BAD: Line too long
config <- list(name = input$experimentName, samples = input$sampleIds, organism = input$organism)
```

```r
# ✓ GOOD: Break on commas in function calls
result <- dplyr::filter(
  data,
  condition_a == TRUE,
  condition_b < 100
)

# ✓ GOOD: Pipe with line breaks
data |>
  dplyr::filter(condition == TRUE) |>
  dplyr::select(col1, col2)
```

### Pipes: Use `|>` not `%>%`

Always use the built-in pipe `|>`:

```r
# ✓ GOOD
result <- data |>
  dplyr::filter(value > 10) |>
  dplyr::group_by(group) |>
  dplyr::summarize(mean_val = mean(value))

# ✗ BAD: Old magrittr pipe
result <- data %>%
  filter(value > 10) %>%
  group_by(group)
```

### Variable Names: snake_case

All variables must be lowercase with underscores:

```r
# ✓ GOOD
cell_ids <- 1:100
expression_matrix <- matrix(rnorm(1000), nrow = 100)
seurat_object <- CreateSeuratObject(counts = matrix_data)

# ✗ BAD: camelCase
cellIds <- 1:100
expressionMatrix <- matrix(rnorm(1000))
seuratObject <- CreateSeuratObject(counts = matrix_data)
```

### Strings: Use "" not ''

Always double quotes:

```r
# ✓ GOOD
message <- "Processing cluster data"
assay_name <- "RNA"
query <- 'SELECT * FROM table'

# ✗ BAD: Single quotes
message <- 'Processing cluster data'
assay_name <- 'RNA'
```

### Comments: Above Code, All Lowercase

Comments go **above** the code they describe. No inline comments. Start lowercase:

```r
# ✓ GOOD: Comment above, lowercase
# calculate weighted mean for each cluster
weighted_mean <- rowMeans(expression_matrix)

# subset to high-expressing genes
high_expr <- expression_matrix[rowMeans(expression_matrix) > 0.5, ]

# ✗ BAD: Inline comment
weighted_mean <- rowMeans(expression_matrix)  # calculate weighted mean

# ✗ BAD: Uppercase
# Calculate weighted mean for each cluster
weighted_mean <- rowMeans(expression_matrix)
```

### No `library()` or `require()`

Use explicit `package::function()` notation:

```r
# ✓ GOOD
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_matrix)
df_filtered <- dplyr::filter(df, value > 100)
ggplot2::ggplot(data, ggplot2::aes(x, y))

# ✗ BAD: library() call
library(Seurat)
library(dplyr)
seurat_obj <- CreateSeuratObject(counts = counts_matrix)
```

## Documentation: roxygen2

Functions must be documented with roxygen2 style comments. Use `#'` prefix:

```r
#' Add cluster information to Seurat object
#'
#' Integrates clustering results into the metadata and adds cluster
#' markers to the object for downstream analysis.
#'
#' @param seurat_object Seurat object to annotate
#' @param cluster_ids integer vector of cluster assignments
#' @param assay character name of assay to use (default: "RNA")
#'
#' @return Seurat object with updated cluster metadata
#' @export
#'
#' @examples
#' \dontrun{
#'   clustered_obj <- add_clusters(
#'     seurat_obj,
#'     cluster_ids = clusters$cluster
#'   )
#' }
add_clusters <- function(seurat_object, cluster_ids, assay = "RNA") {
  # implementation
}
```

**Guidelines**:
- Title: One-line summary (no period)
- Description: Full details (can span multiple lines)
- `@param`: Name and description of each parameter
- `@return`: What the function returns
- `@export`: If function is part of public API
- `@examples`: Optional; wrap in `\dontrun{}` if interactive/requires data

## Common Linting Issues

After running `lintr::lint()`, fix these manually (styler won't auto-fix):

| Issue | Fix |
|-------|-----|
| `object_name_linter`: Variable not snake_case | Rename to `snake_case` |
| `commented_code_linter`: Commented-out code | Remove or convert to proper comment |
| `assignment_linter`: Using `=` instead of `<-` | Use `<-` for assignment |
| `spaces_linter`: Missing spaces around operators | Add spaces: `a + b` not `a+b` |
| `T/F` usage | Use `TRUE`/`FALSE` (not `T`/`F`) |

Example:

```r
# ✗ DETECTED BY LINTER
cellSet <- c(1, 2, 3)  # variable name, assignment operator
isLarge=TRUE           # spacing

# ✓ AFTER FIX
cell_set <- c(1, 2, 3)
is_large <- TRUE
```

## Workflow

1. **Edit your code** with any changes
2. **Auto-format**: `styler::style_file("file.R")`
3. **Lint**: `lintr::lint("file.R")`
4. **Fix remaining issues** manually from linting output
5. **Run tests**: `make test-r-file FILE=test-myfile.R` (or full suite: `make test-r`)
6. **Commit** with signoff: `git commit -s -m "brief message"`

## References

- [roxygen2 documentation](https://roxygen2.r-lib.org/)
- [tidyverse R style guide](https://style.tidyverse.org/)
- [lintr configuration](https://lintr.r-lib.org/)
