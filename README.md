<img src="https://raw.githubusercontent.com/skptic/nf-toxomix/master/assets/toxomix_logo.png" width="400">

<img src="https://raw.githubusercontent.com/skptic/nf-toxomix/master/assets/orn_logo.png" width="400">


# NF-toxomix
Pipeline for toxicology predictions based on transcriptomic profiles.

It acts a pilot workflow as part of the [OpenRiskNet](https://openrisknet.org/) project with the goal to incorperate genomic data into the OpenRiskNet infrastructure.

[![Build Status](https://travis-ci.org/skptic/NF-toxomix.svg?branch=master)](https://travis-ci.org/skptic/NF-toxomix)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](https://www.nextflow.io/)


### Introduction
NF-toxomix: Pipeline for toxicology predictions based on transcriptomic profiles

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.

The pipeline from Juma Bayjan at Maastricht University which in turn aimed to reproduce the article "A transcriptomics-based in vitro assay for predicting chemical genotoxicity in vivo", by C.Magkoufopoulou et. al.


### Documentation
The NF-toxomix pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Credits
This pipeline was written by Evan Floden ([skptic](https://github.com/skptic)) at [Center for Genomic Regulation (CRG)](http://www.crg.eu).
