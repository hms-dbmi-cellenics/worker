---
applyTo: "python/**/*.py"
description: "Python code formatting for Cellenics worker: 79-char lines, black formatter, isort, flake8 linter"
---

# Python Code Formatting & Style

This instruction applies to all Python files in the `python/` directory.

## Formatting

**Line length**: 79 characters (enforced by Black)

**Formatter**: [Black](https://black.readthedocs.io/)
```bash
black --line-length 79 python/
```

**Import sorter**: [isort](https://pycqa.github.io/isort/)
```bash
isort python/
```

**Linter**: [flake8](https://flake8.pycqa.org/)
```bash
flake8 python/
```

**Quick workflow**:
```bash
make fmt   # Run black + isort
make check # Run flake8 validation
```

## Style Rules

### Naming Conventions

- **Functions & variables**: `snake_case`
- **Classes**: `PascalCase`
- **Constants**: `UPPER_SNAKE_CASE`
- **Private** (internal): prefix with `_`

```python
# ✓ GOOD
def process_task(task_data):
    config_value = get_config()
    
class TaskProcessor:
    RETRY_LIMIT = 3
    
    def _validate_input(self):
        pass

# ✗ BAD
def processTask(taskData):
    configValue = get_config()
    
class taskProcessor:  # Should be PascalCase
    retryLimit = 3    # Should be UPPER_SNAKE_CASE
```

### Imports

Use `isort` to organize imports. Typically:
1. Standard library
2. Third-party packages
3. Local modules

```python
# ✓ GOOD (isort will organize this)
import json
import logging
from pathlib import Path
from typing import Optional

import boto3
import pytest

from worker.tasks import TaskFactory
from worker.utils import logger
```

### Type Hints

Use type hints for clarity:

```python
# ✓ GOOD
def calculate_statistics(
    data: list[float],
    window_size: int = 10
) -> dict[str, float]:
    """Calculate rolling statistics."""
    return {"mean": 0.0}

# ✗ BAD: No type hints
def calculate_statistics(data, window_size=10):
    return {"mean": 0.0}
```

### Docstrings

Use docstrings for modules, classes, and functions:

```python
"""Module for task processing and queue management."""

class TaskProcessor:
    """Handle individual task execution and error reporting."""
    
    def process_task(self, task_id: str, timeout: int = 300) -> dict:
        """
        Process a single task from the queue.
        
        Args:
            task_id: Unique identifier for the task
            timeout: Maximum execution time in seconds
            
        Returns:
            Dictionary with result_status, output, and any error messages
            
        Raises:
            TimeoutError: If task exceeds timeout
            TaskValidationError: If task format is invalid
        """
        # implementation
```

### Comments

- **Above code**, not inline
- Start with lowercase
- Explain *why*, not *what*

```python
# ✓ GOOD
# exponential backoff helps reduce load on API during failures
retry_delay = initial_delay * (backoff_factor ** attempt)
perform_retry(request, delay=retry_delay)

# ✗ BAD: Inline, restating obvious
retry_delay = initial_delay * (backoff_factor ** attempt)  # multiply by backoff_factor

# ✗ BAD: Uppercase start
# Retry the request with exponential backoff
```

### String Formatting

Prefer f-strings:

```python
# ✓ GOOD
message = f"Processing {task_count} tasks for {experiment_id}"
log.info(f"Status: {status}, retries: {retry_count}")

# ✗ OLD: .format()
message = "Processing {} tasks for {}".format(task_count, experiment_id)

# ✗ OLD: % formatting
message = "Processing %s tasks for %s" % (task_count, experiment_id)
```

### Error Handling

Use specific exceptions and provide context:

```python
# ✓ GOOD
try:
    response = requests.post(api_url, json=payload, timeout=30)
    response.raise_for_status()
except requests.Timeout:
    logger.error(f"API request timed out for experiment {exp_id}")
    raise TaskExecutionError(f"API timeout: {exp_id}") from None
except requests.HTTPError as e:
    logger.error(f"API returned {e.response.status_code}: {e}")
    raise TaskExecutionError(f"API error: {e.response.text}") from e

# ✗ BAD: Too broad
try:
    response = requests.post(api_url, json=payload)
except Exception:
    print("error")  # Lost context
```

## Commit Messages

- **One line**, brief and descriptive
- Always sign with `-s` flag
- Example: `git commit -s -m "add retry logic for S3 uploads"`

## Testing

Tests go in `python/tests/` and follow pytest conventions:

```python
# test_task_processor.py
import pytest
from unittest import mock

from worker.tasks import TaskProcessor

def test_process_task_success():
    """Test successful task processing."""
    processor = TaskProcessor()
    result = processor.process_task("task123")
    assert result["status"] == "completed"

@mock.patch("worker.tasks.S3Client")
def test_process_task_with_s3_failure(mock_s3):
    """Test graceful handling of S3 errors."""
    mock_s3.return_value.upload.side_effect = Exception("S3 timeout")
    processor = TaskProcessor()
    
    with pytest.raises(TaskExecutionError):
        processor.process_task("task123")
```

Run tests with:
```bash
make test-py
```

## Configuration

Black and isort configuration is in `Makefile`:
- Black: `--line-length 79`

No additional config files needed; the Makefile handles setup.

## References

- [Black documentation](https://black.readthedocs.io/)
- [isort documentation](https://pycqa.github.io/isort/)
- [flake8 rules](https://flake8.pycqa.org/en/latest/user/error-codes.html)
- [PEP 8 style guide](https://www.python.org/dev/peps/pep-0008/)
