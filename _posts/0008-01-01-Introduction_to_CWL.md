---
feature_text: |
  ## Precision Medicine
title: Introduction to the Common Workflow Language
categories:
    - Module 8
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-01-01
---

Often in data analysis a collection of programatic tools and resources are required to perform a specific task. This can take the form of bash scripts linking tools together however bash scripts offer no standardization and often need to be tweaked from one project to the next. The Common Workflow Language (CWL) is a specific for designing portable and scaleable workflows. It is open source and available on [github](https://github.com/common-workflow-language/common-workflow-language) under an apache 2.0 license. Using CWL instead of bash offers a number of advantages, these include:

1. Automated paralization of workflow sub tasks
2. Modular, can easily add and remove tools from a workflow
3. Portability, workflows will work with any environment with CWL installed
4. No need to monitor workflows, one sub-task failure won't cause the entire workflow to fail
5. Standardized format makes reading workflows easier

In this module we will use CWL with Docker to build an analysis pipeline to perform DNA alignment.
