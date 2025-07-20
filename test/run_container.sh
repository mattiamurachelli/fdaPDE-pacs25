#/bin/bash

testing=$(pwd)/../

docker run --rm -v $testing:/root/fdaPDE-testing -ti aldoclemente/fdapde-docker:latest /bin/bash 

