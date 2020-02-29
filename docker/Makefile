THIS_FILE := $(realpath $(lastword $(MAKEFILE_LIST)))
THIS_FILE_DIR := $(shell dirname $(THIS_FILE))
IMAGE = gauss:1.0
CONTAINER = gauss-container

build: 
	docker build --tag $(IMAGE)  \
	    -f Dockerfile \
        .

run:
	docker run --rm \
		-it \
		-p 8888:8888 \
	    -v $(THIS_FILE_DIR)/gp-tutorial:/tf/gp-tutorial \
	    --name $(CONTAINER) \
	    $(IMAGE)
