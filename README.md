# Predicting the Frequencies of drug side effects

This repository contains two folders:

**Data collection**. Use the python code provided in this folder if you want to parse the original files from SIDER 4.1. Notice that this part of the code requires data from http://sideeffects.embl.de/ and ensure that each drug has an ATC code, drugs need to be further filtered using propietary data from the World Health Organization (WHO) https://www.whocc.no/atc_ddd_index/ (commercial).

**Data analysis**. Use the matlab code provided in this folder if you want to run our matrix decomposition algorithm. This is self-contained: it contains the .mat files of the matrices (frequencies and post-marketing associations) reproduce our study.

## Dependencies

### 1. Data Collection.

The data collection is illustrated in a jupyther notebook.

The dependencies required to run the algorithm are the following:
```python
import csv
import os
import pickle
import numpy, scipy.io
from collections import defaultdict
import numpy, scipy.io
from xml.etree.ElementTree import iterparse
from collections import defaultdict
from pprint import pprint
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles
from collections import Counter
import operator
```

### 2. Data Analysis

This folder contains the code to run our algorithm.

**Compatibility**: the code was implemented in Matlab R2018a. However, we used built-in functions
that had been available since previous versions:

```matlab
* perfcurve: introduced in R2009a.
* histogram: introduced in R2014b.
```

To run the code, follow the guidelines in */data analysis/readme.pdf*

#### Time to run the algorithm 

Running the algorithm takes roughly *32 seconds* in a 6-core Intel Xeon CPU E5-1650v3@3.5GHz with 32GB of DDR4. 

# Reference
Galeano, D., Li, S., Gerstein, M & Paccanaro, A. (2020). Predicting the frequencies to drug side effects. Nature Communications.

Galeano, D., & Paccanaro, A. (2019). Predicting the Frequency of Drug Side effects. bioRxiv, 594465.
