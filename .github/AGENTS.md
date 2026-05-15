---
description: "Agent instructions for Cellenics worker project: mixed R/Python data analysis executor with Kubernetes deployment"
---

# Cellenics Worker - Agent Instructions

Cellenics worker is a dual-container data analysis executor for single-cell genomics. The **R component** performs statistical analysis on Seurat objects; the **Python component** manages task queuing, orchestration, and S3 integration.

For detailed architecture, patterns, and common pitfalls, see [CODEBASE_ANALYSIS.md](./CODEBASE_ANALYSIS.md).

## Quick Start

- **Build & test**: `make build test-r test-py`
- **Run locally**: `make run` (requires Docker, 10GB+ RAM, `GITHUB_API_TOKEN` env var)
- **R version**: 4.4.0+ with gFortran (macOS)
- **Python**: 3.8+

## Code Conventions

### R Code Style

**Formatting & Linting**:
- Lines **< 80 characters** (hard limit)
- Auto-format with: `styler::style_file("path/to/file.R")`
- Lint with: `lintr::lint("path/to/file.R")` and fix remaining issues
- Use built-in pipe `|>` (not magrittr `%>%`)
- All `library()` calls forbidden; use `package::function()` instead

**Naming & Comments**:
- Variables: `snake_case` (never camelCase)
- Comments: **above code**, all lowercase, no leading uppercase
  ```r
  # calculate the weighted mean for each cluster
  weighted_mean <- rowMeans(expression_matrix)
  ```
- Function documentation: roxygen2 style (`#' @param`, `#' @return`, `#' @export`)

**Strings & Syntax**:
- Strings: double quotes `""` (never `''`)
- Long function calls: break on opening/closing brackets
  ```r
  config <- list(
    name = input$experimentName,
    samples = input$sampleIds,
    organism = input$organism
  )
  ```

**Avoid**:
- `library()` or `require()` → use `package::function()`
- Inline comments → move above code
- Single quotes → always `""`

### Python Code Style

- Line length: **79 characters**
- Formatter: `black` (configured in Makefile)
- Import sorter: `isort`
- Linter: `flake8`
- Run: `make fmt` to format; `make check` to validate

### Commit Messages

- **One line**, brief, descriptive
- Always sign with `-s` flag: `git commit -s -m "message"`
- Example: `git commit -s -m "add clustering task handler"`

## R Component

See [r/README.md](./r/README.md) for R-specific setup and development.

**Key patterns**:
- Error handling: `generateErrorMessage(code, message)` with delimiter `:|:`
- Reproducibility: Set `ULTIMATE_SEED = 42` in tasks
- Seurat workflows: Auto-detect BPCells, handle Sketch assays
- Task execution: REST API listens on localhost (RestRserve)

**Testing**: `make test-r` or `make test-r-file FILE=test-expression.R`
- Framework: testthat (3.0.0+)
- Pattern: Mock Seurat objects, snapshot regression

## Python Component

See [python/README.md](./python/README.md) for Python-specific setup.

**Key patterns**:
- Task factory pattern with ABC base classes
- Config precedence: K8s labels → env vars → defaults
- HTTP retries: Exponential backoff (max 30s)
- S3 compression: >250KB uploaded to S3; ≤250KB via socket

**Testing**: `make test-py`
- Framework: pytest
- Pattern: `@mock.patch` for S3/SQS/Redis, fixture-based

## Deployment

Helm chart in [chart-infra/](./chart-infra/). Deployed to AWS Fargate via cluster Helm operator.

**Critical variables**:
- `CLUSTER_ENV`: deployment target (development/staging/production)
- `GITHUB_API_TOKEN`: required for Docker builds (silently fails if missing)
- Pod labels override env vars; env vars override defaults

See [ci-role-worker-cf.yaml](./ci-role-worker-cf.yaml) and Github Actions workflows.

## Common Development Issues

1. **Docker build fails silently**: Missing `GITHUB_API_TOKEN` env var
2. **Freezing/OOM**: Docker allocated <10GB RAM
3. **R version mismatch**: Check `renv.lock`, use 4.4.0+
4. **gFortran missing** (macOS): Download from [gFortran releases](https://github.com/fxcoudert/gfortran-for-macOS/releases)
5. **File staleness**: If reloading Seurat objects fails, check S3 timestamps (built-in retry loop handles transient failures)

For more, see [CODEBASE_ANALYSIS.md](./CODEBASE_ANALYSIS.md).

## Agents & Skills

- **Test Coverage Agent** (`.github/agents/test-coverage.agent.md`): Comparing branches, writing tests, ensuring coverage with best practices
