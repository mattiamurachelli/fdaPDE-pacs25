#!/bin/bash

docker run --rm -d -p 8787:8787 -v $(pwd)/../../fdaPDE-pacs25:/home/user/fdaPDE-pacs25 --name rstudio -e PASSWORD=password aldoclemente/fdapde-docker:rstudio


# connect to http://localhost:8787
# username: user
# password: password
