---
feature_text: |
  ## Precision Medicine
title: RNAseq Fusion Filtering/Annotation/Review
categories:
    - Module-06-RNAseq
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0006-05-02
---

**DEV notes:**


## Collapsible markdown example

<details><summary>CLICK ME</summary>
<p>

#### yes, even hidden code blocks!

```python
print("hello world!")
```

</p>
</details>



---


# Introduction
# Setup
# Filtering Example
Pizzly generates outputs in fasta and json formats. The downloaded pizzly repository also contains a python script to convert .json output to a simple gene table:

```bash
cd ~/workspace/bin
 wget https://raw.githubusercontent.com/pmelsted/pizzly/master/scripts/flatten_json.py
python ~/workpace/bin/flatten_json.py <FILE>
```

# Basic Visualization in IGV
# Advanced Visualization in chimeraviz

```bash
R

source("https://bioconductor.org/biocLite.R")
biocLite("chimeraviz")
```
https://bioconductor.org/packages/release/bioc/html/chimeraviz.html
dependnecy: openssl error
sudo apt-get install libssl-dev
```
