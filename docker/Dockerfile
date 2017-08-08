FROM ubuntu:14.04

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 libbz2-dev liblzma-dev libcurl4-openssl-dev \
    ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 git mercurial subversion \
    make gcc \
    postgresql postgresql-contrib python-psycopg2 libpq-dev libmysqlclient-dev 

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/archive/Anaconda2-4.4.0-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm ~/anaconda.sh

RUN /opt/conda/bin/pip install mysql-python


RUN wget --quiet https://bootstrap.pypa.io/ez_setup.py -O - | /opt/conda/bin/python 

RUN cd && \
 wget --quiet https://github.com/biopython/biopython/archive/biopython-169.tar.gz  && \
 tar zxvf biopython-169.tar.gz && \
 cd /root/biopython-biopython-169 && \
 /opt/conda/bin/python setup.py install && \
 cd 

RUN cd && \
 wget --quiet https://github.com/samtools/htslib/releases/download/1.5/htslib-1.5.tar.bz2 && \
 bunzip2 htslib-1.5.tar.bz2 && \
 tar xvf htslib-1.5.tar && \
 cd /root/htslib-1.5 && \
 ./configure && \
 make && \
 make install && \
 cd 


# The command that is twice is NOT an error.. 
RUN cd && \
 git clone https://github.com/kantale/MutationInfo && \
 cd /root/MutationInfo/ && \
 /opt/conda/bin/python setup.py install && \
 /opt/conda/bin/python setup.py install && \
 cd

ADD initialize_mi.py initialize_mi.py 

RUN /opt/conda/bin/python /initialize_mi.py  

RUN /opt/conda/bin/pip install Django==1.9

RUN cd /root/MutationInfo/web && \
  /opt/conda/bin/python manage.py migrate && \
  cd

ENTRYPOINT [ "/opt/conda/bin/python", "/root/MutationInfo/web/manage.py", "runserver", "0.0.0.0:4444" ]


