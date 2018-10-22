---
feature_text: |
  ## Precision Medicine
title: AWS Instance Setup
categories:
    - Module-01-Setup
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-02-02
---

Covered in this section: logging into AWC EC2 console, starting an instance from the course AMI, configuring it in the console (select instance AMI, instance type, instance details, storage volumes, tags, security group, and key pairs).

Basic intro to the instance (top, resources available, mount location of volumes, etc.).

#### Logging into an instance

In the previous section [AWS Intro](http://pmbio.org/module-01-setup/0001/02/01/AWS_Intro/) we reviewed fundamental concepts of cloud computing and some of the jargon and features specific to AWS. In this section we will learn how to launch an instance specifically for this course.

In order to launch your own instance you will either need to use your own personal AWS account (not recommended unless you are already familiar with and using AWS) OR you will need to be assigned an AWS account using the IAMS system. If neither of these is possible, the instructors will have to launch an instance for you and provide the login details.

Instance details:
* AMI: XXXXXXXXXXXX
* Based on OS: Ubuntu Server 18.04 LTS (HVM), SSD Volume Type
* Instance type: r5.2xlarge (8 CPUs, 64 GB RAM, Up to 10 Gigabit network speed, EBS storage only)
* Volumes: 2. A 500GB "root" ("/") volume, and a 2TB "workspace" volume 
* Security group and details: "pmbio ssh and http"
* Other notes: Protect from accidental termination.

NOTE: CHECK INSTANCE LIMITS FOR THE COURSE AWS REGION (EU LONDON) FOR THE INSTANCE TYPE WE CHOOSE. REQUESTED. CHECK BACK FOR APPROVAL
