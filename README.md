# MeTarEnz

MeTarEnz (metagenomic targeted enzyme miner) is a multi-service software that enables the targeted screening of high-throughput metagenomic data with user-defined databases and bit-score cut-offs. Moreover, MeTarEnz contains an assembly module making it capable of providing its screening services while accepting various types of input data such as raw reads, assembled contigs, and translated protein sequences. As a complementary section of this software, specially designed for a lipase mining project, trained regression models for prediction of lipases’ temperature and pH optima are included as well. This tool is freely accessible in forms of the standalone toolkit, docker image, and python package. Figure 1 is a pictorial workflow of MeTarEnz’s services.

![alt text](https://github.com/mehdiforoozandeh/MeTarEnz/blob/master/Workflow.jpg?raw=true)

*Figure 1: Workflow of Services. This figure illustrates the general workflow of the MeTarEnz including its acceptable inputs, utilized software, user-defined parameters, etc.*

## Distributions

### Pre-built Docker Image

Firstly, to use MeTarEnz’s pre-built Docker image, you must have docker software installed (https://docs.docker.com/engine/install/). 
To download the MeTarEnz image, run the following command:

'docker pull mforooz/metarenz'

After pulling the image, you can run MeTarEnz through the command below:

'docker run -ti --entrypoint "" mforooz/metarenz python metarenz.py'

Or you can run a general test to check if it works fine by the following command:

'docker run -ti --entrypoint "" mforooz/metarenz python test.py'

For more detailed instructions on how to use docker images, containers, etc. please read the Docker documentation (https://docs.docker.com/).				

### Standalone Toolkit
MeTarEnz is also available in the form of pre-built executables for linux64. The Standalone distribution of this tool can be downloaded from cbb.ut.ac.ir/metarenz or https://github.com/mehdiforoozandeh/metarenz. 

### Python Package
All MeTarEnz’s codes, dependencies, models, etc. are available as a python package. This version enables users to modify and use this tool according to their needs by tweaking the source code, updating dependencies, and models. The python package of this tool can be downloaded from both cbb.ut.ac.ir/metarenz and https://github.com/mehdiforoozandeh/metarenz. 


## General Usage Guide	
Throughout this section, it is assumed that the MeTarEnz is executed from its main directory. If otherwise, it is recommended to use files’ full addresses. 

**MeTarEnz can also be executed interactively with -int or --interactive.**


'./metarenz --interactive' (from the same directory of the standalone version)
'python metarenz.py --interactive' (from the same directory of the python package version)
'''docker run -ti --entrypoint "" mforooz/metarenz python metarenz.py'''  --interactive(docker image)


In the usage guide below, for the non-interactive execution of MeTarEnz, arguments are separated with space. Figure 2 summarizes MeTarEnz’s different functions, inputs, and user-defined parameters. 

![alt text](https://github.com/mehdiforoozandeh/MeTarEnz/blob/master/Graphical%20Help.jpg?raw=true)

*Figure 2: Summary of Functions. This figure illustrates MeTarEnz’s different services, inputs, parameters, etc.*
