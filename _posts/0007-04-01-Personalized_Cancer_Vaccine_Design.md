---
feature_text: |
  ## Precision Medicine
title: Personalized Cancer Vaccine Design
categories:
    - Module 7
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0007-04-01
---

## Using AWS to run Pvactools 

#### **Installation step**
_____________________

1. Configure and launch a ubuntu instance. Depending on how much data you plan on processing, the amount of memory and storage should be customized accordingly. For this small example, the configuration was: **t2.medium instance + 50G EBS storage**. 

2. Assuming a ubuntu instance with appropriate amount of storage and launched, login to designated remote computer through terminal and mount the storage to the data directory (`/data`).

    In order to fully utilize the attached storage, before installing any of the packages, please enter the data directory using the command:`cd /data`  

3. Installing pvactools package:

    * Check for Python version using python3 -V. If you launched an instance as specified, the screen should show:
    `Python 3.5.2`  
    * Next you will need to install pip3 using the following command:  
    `sudo apt-get install python3-setuptools`  
    `sudo easy_install3 pip`  
    * Now you are ready to install pvactools:  
    `sudo pip3 install pvactools`  
    * In order to confirm that pvactools has been successfully installed, use the command `pip3 list`, in the list of packages, you should be able to find `pvactools (1.0.3)`.  

4. Installing **IEDB binding prediction tools** locally is highly recommended for analysis of larger datasets or when the making predictions for multiple alleles, epitope lengths, or prediction algorithm. In order to use the local version of IEDB with pvactools, we will first need to install **Anaconda**.  
    * Inside the data mounted directory (`/data`), install anaconda by using the following commands:  
    `wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh`  
    `bash Anaconda3-5.0.1-Linux-x86_64.sh -b -p /data/anaconda`  
    `rm Anaconda3-5.0.1-Linux-x86_64.sh`  
    `echo 'export PATH="/data/anaconda/bin:$PATH"' >> ~/.bashrc`  
    `source ~/.bashrc`

    * In your data storage directory (with the attached volume mounted), download IEDB MHC Class I&II predictor tools using the commands:  
    `wget http://downloads.iedb.org/tools/mhci/LATEST/IEDB_MHC_I-2.17.tar.gz`  
    `wget http://downloads.iedb.org/tools/mhcii/LATEST/IEDB_MHC_II-2.17.3.tar.gz`  

    * IEDB also requires **tcsh**, install as following:  
    `sudo apt-get update`  
    `sudo apt-get install tcsh`  

    * Unpack the files downloaded from IEDB using following commands:  
    `tar -zxvf IEDB_MHC_I-2.17.tar.gz`  
    `tar -zxvf IEDB_MHC_II-2.17.3.tar.gz`  

    * Now you will need to setup a **Python2.7** environment:  
    `conda create -n python2.7 python=2.7`  

    * Enter the mch_i directory and activate python2.7 environment as following:  
    `cd mhc_i`  
    `source activate python2.7`  
    `./configure`  

    * Now enter the mch_ii directory (`cd ../mhc_ii`) and do the following:  
    `python configure.py`  

    * Moving forward, you will run pvactools in a Python3.5 environment (default), thus please deactivate the python2.7 environment:  
    `source deactivate python2.7`

**This concludes the necessary installations for running pVACtools.**

#### **Test Run Example**
__________________________  
To do a simple test run, you can download an example using the command:  
    `pvacseq download_example_data /data/`  
    In the current additional_input_file_list.yaml, you will need to change the directory of the files. If you followed the tutorial as specified, the directory of all four files would be: `/data/pvacseq_example_data/`  

Ensure that you are running in a Python3.5 environment, and run the following example command:  
    `pvacseq run pvacseq_example_data/input.vcf Test HLA-G*01:09,HLA-E*01:01,H2-IAb NetMHC PickPocket NNalign example_output/ -e 9,10 -i pvacseq_example_data/additional_input_file_list.yaml --tdna-vaf 20 --net-chop-method cterm --netmhc-stab --top-score-metric=lowest -d full --keep-tmp-files --iedb-install-directory /data/`  

Note that if the job has failed previously, you will need to delete the directories that have been created in the output file in order to rerun the command.  

If you would like to transfer your own data from your laptop/PC, you may use the scp command with -i and your private key for the instance.

For more specific details on usage options, please refer to the [pVACtools website](http://pvactools.readthedocs.io/en/latest/pvacseq/run.html).
