FROM continuumio/miniconda3
WORKDIR /usr/src/app
RUN apt-get update && \
    apt-get -y install vim 

# Copy function code
ADD ../ /usr/CixTools

#RUN conda update -n base -c defaults conda
RUN conda init bash
WORKDIR /usr/CixTools
RUN conda create --name SignelIx_CIx --file requirements.txt -c conda-forge -c rdkit -c intel -c default 
RUN conda init bash
RUN conda env update --name SignelIx_CIx --file addtl_env.yml 

SHELL ["bash", "-lc"]

#docker build -t cixtools 
# docker run -w /usr/CixTools -i -t cixtools  /bin/bash
#conda activate SignelIx_CIx
#pip list --format=freeze > requirements.txt
