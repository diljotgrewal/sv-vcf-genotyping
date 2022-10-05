FROM quay.io/singlecellpipelinetest/miniconda3:4.8.2


ADD . /app


RUN cp -r /app /code && chmod -R 777 /code/*
ENV PATH="${PATH}:/code"

RUN apt-get --allow-releaseinfo-change update && apt install gcc zlib1g-dev g++ libbz2-dev build-essential -y
RUN cd /code/ && git clone https://github.com/fritzsedlazeck/SURVIVOR.git && cd SURVIVOR/Debug && pwd && ls && make
ENV PATH="${PATH}:/code/SURVIVOR/Debug"

RUN pip install click pandas sniffles pyvcf
