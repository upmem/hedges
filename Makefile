.PHONY: clean build

PWD= $(shell pwd)
DOCKER_IMG=bionic-hedges

build_docker:
	docker build -t ${DOCKER_IMG} dhub/bionic/python2.7-numpy/

test_docker_env:
	docker run -v ${PWD}:/work ${DOCKER_IMG} bash -c "mkdir -p /tmp/build && cd /tmp/build && cmake /work && make && PYTHONPATH=/tmp/build/src python -u /work/print_module_help_files.py"

clean:
	docker run -v ${PWD}:/work ${DOCKER_IMG} bash -c "rm -r /work/build/"

build:
	docker run -v ${PWD}:/work ${DOCKER_IMG} bash -c "mkdir -p /work/build && cd /work/build && cmake /work/ && make"

testprogramm:
	docker run -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work/ && PYTHONPATH=/work/build/src python -u test_program.py"

print_module_help_files:
	docker run -v ${PWD}:/work ${DOCKER_IMG} bash -c "cd /work/ && PYTHONPATH=/work/build/src python -u print_module_help_files.py"

.DEFAULT_GOAL := build