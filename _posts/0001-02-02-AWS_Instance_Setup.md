---
feature_text: |
  ## Precision Medicine
title: AWS Instance Setup
categories:
    - Module-01-Setup
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-02-02
---

Covered in this section: logging into AWS EC2 console, starting an instance from the course AMI, configuring it in the console (select instance AMI, instance type, instance details, storage volumes, tags, security group, and key pairs).

Basic intro to the instance (top, resources available, mount location of volumes, etc.).

#### Launching an AWS instance for the course

In the previous section [AWS Intro](http://pmbio.org/module-01-setup/0001/02/01/AWS_Intro/) we reviewed fundamental concepts of cloud computing and some of the jargon and features specific to AWS. In this section we will learn how to launch an instance specifically for this course.

In order to launch your own instance you will either need to use your own personal AWS account (not recommended unless you are already familiar with and using AWS) OR you will need to be assigned an AWS account using the IAMS system. If neither of these is possible, the instructors will have to launch an instance for you and provide the login details.

Instance details:
* AMI: `pmbio ami v4` (`ami-05634ed111b3581c1`) (Publicly available in `London` Zone)
* Based on OS: Ubuntu Server 18.04 LTS (HVM), SSD Volume Type
* Instance type: r5.2xlarge (8 CPUs, 64 GB RAM, Up to 10 Gigabit network speed, EBS storage only)
* Volumes: 2. A 500GB "root" ("/") volume, and a 2TB "workspace" volume
* Security group and details: "pmbio ssh and http"
* Other notes: Protect from accidental termination.

#### Logging into your instance

For this course the instructors have launched instances for you (so that you will not be charged). You will be provided with the IP address and/or a shorthand name for your instance. Each student will have an independent instance. It is important that you always use your own instance and do not share with another student. You can use the 'who' command to see if there is more than one login active on your instance.

Example login commands (for Mac/linux laptops):
```bash
ssh -i pmbio.pem ubuntu@s1.pmbio.org

```

In this command, 'pmbio.pem' is the AWS key file used to authenticate SSH access to the machine, 'ubuntu' is the name of the user account on the machine (a common default user name for the Ubuntu linux operating system), and 's1.pmbio.org' is a domain name linked to the IP address of the machine.  This last part will need to replace with your instance's unique IP address or the shorthand domain name.

#### Use of a terminal multiplexer

Throughout this course we will be logged into an AWS terminal session via SSH connection from your laptop to the remote cloud instance. Some points on this:

* Your terminal connection relies on your laptop's WiFI connection (which could reset periodically)
* Your terminal connection also relies on your laptop being on. Run out of power or close the lid and it dies
* If your terminal connection dies and you haven't taken precautions, any command currently running will likely fail

There are various ways to get around this problem, but the most robust approach is to use a "terminal multiplexer". The purpose of these programs is to make your session persist regardless of your SSH connection.  If your connection gets killed, you can just log into the AWS instance again, attach the multiplexer session and it will be as if nothing happened.

If you are already familiar with `screen`, `tmux`, `byobu`, etc. you can use on of those. If you are not already familiar, we recommend `screen` for its simplicity.

To use it, simply type `screen` when you login. If you have already logged in previously and left a `screen` session running you may have to use `screen -r` or `screen -r -d`. When you want to leave your session, you can just close the window or use `screen -d` to detach your session.

Basic screen commands:
* `screen`. Start a new session.
* `screen -r`. Resume a session that was properly detachead.
* `screen -r -d`. Resume a session that was NOT properly detached.
* `screen -d`. Detach a session so that later you can resume.
* `<ctrl> <a> <ctrl> <c>`. Create an additional tab in your current session.
* `<ctrl> <space>`. Cycle through existing tabs.
* `<ctrl> <a> <ctrl> <a>`. Flip between two tabs.
* `<ctrl> <a> <">` View a list of tabs.
