#!/bin/bash

# Usage: build.sh 

# Determine the folder of this script
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo BASE_DIR
docker build -t gp-docker \
	-f "$BASE_DIR/Dockerfile" 
