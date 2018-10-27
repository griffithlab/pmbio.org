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

#### Launching an AWS instance for the course

In the previous section [AWS Intro](http://pmbio.org/module-01-setup/0001/02/01/AWS_Intro/) we reviewed fundamental concepts of cloud computing and some of the jargon and features specific to AWS. In this section we will learn how to launch an instance specifically for this course.

In order to launch your own instance you will either need to use your own personal AWS account (not recommended unless you are already familiar with and using AWS) OR you will need to be assigned an AWS account using the IAMS system. If neither of these is possible, the instructors will have to launch an instance for you and provide the login details.

Instance details:
* AMI: `pmbio ami v3` (`ami-08d777ebd086b116f`) (Publicly available in `London` Zone) 
* Based on OS: Ubuntu Server 18.04 LTS (HVM), SSD Volume Type
* Instance type: r5.2xlarge (8 CPUs, 64 GB RAM, Up to 10 Gigabit network speed, EBS storage only)
* Volumes: 2. A 500GB "root" ("/") volume, and a 2TB "workspace" volume 
* Security group and details: "pmbio ssh and http"
* Other notes: Protect from accidental termination.

#### Logging into your instance

For this course the instructors have launched instances for you (so that you will not be charged). You will be provided with the IP address and/or a shorthand name for your instance. Each student will have an independent instance. It is important that you always use your own instance and do not share with another student. You can use the 'who' command to see if there is more than one login active on your instance.

Example login commands (for Mac/linux laptops):
```bash
ssh -i pmbio.pem ubuntu@18.220.123.159

```

In this command, 'pmbio.pem' is the AWS key file used to authenticate SSH access to the machine, 'ubuntu' is the name of the user account on the machine (a common default user name for the Ubuntu linux operating system), and '18.220.123.159' is the IP address of the machine.  This last part will need to replace with your instance's unique IP address.


