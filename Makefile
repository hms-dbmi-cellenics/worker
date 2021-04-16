#!make
#----------------------------------------
# Settings
#----------------------------------------
.DEFAULT_GOAL := help
#--------------------------------------------------
# Variables
#--------------------------------------------------
PYTHON_FILES?=$$(find python -name '*.py' | grep -v venv) # do not include virtual env files
# If unix name is Darwin, we are in MacOS => use regular docker-compose.yaml
# Otherwise assume Linux and add docker-compose.linux-dev.yaml overrides
ifneq ($(shell uname -s), Darwin)
  docker_files=-f docker-compose.yaml
else
  docker_files=-f docker-compose.yaml -f docker-compose.linux-dev.yaml
endif
#--------------------------------------------------
# Targets
#--------------------------------------------------
bootstrap: ## Installs requirements (python linter & formatter)
	@pip install flake8 black
fmt: ## Formats python files
	@echo "==> Formatting files..."
	@black $(PYTHON_FILES)
	@echo ""
check: ## Checks code for linting/construct errors
	@echo "==> Checking if files are well formatted..."
	flake8 $(PYTHON_FILES)
	@echo "    [✓]\n"
build: ## Builds the docker-compose environment
	@echo "==> Building docker image..."
	docker-compose $(docker_files) build
	@echo "    [✓]\n"
run: # Runs the docker environment
	@docker-compose $(docker_files) up
logs: # Shows live logs if the workers are running or logs from last running worker if they are not.
	@docker-compose $(docker_files) logs -f
kill: # Kills the currently running environment
	@docker-compose $(docker_files) kill
.PHONY: bootstrap fmt check build run clean help
clean: ## Cleans up temporary files
	@echo "==> Cleaning up ..."
	@find . -name "*.pyc" -exec rm -f {} \;
	@echo "    [✓]"
	@echo ""
help: ## Shows available targets
	@fgrep -h "## " $(MAKEFILE_LIST) | fgrep -v fgrep | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-13s\033[0m %s\n", $$1, $$2}'
