# Predicting the Frequencies of drug side effects

This repository contains two folders:

**Data collection**. If you want to parse the original files from SIDER 4.1. Notice that this part of the code requires data from http://sideeffects.embl.de/ and drugs need to be filtered using propietary data from https://www.whocc.no/atc_ddd_index/ (commercial).

**Data analysis**. Contains a step-by-step Matlab R2018a code to run our matrix decomposition algorithm. This also contains .mat files of the matrices (frequencies and post-marketing associations) built for our study.

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

Compatibility: the code was implemented in Matlab R2018a. However, we used built-in functions
that had been available since previous versions:

```matlab
* perfcurve: introduced in R2009a.
* histogram: introduced in R2014b.
```

To run the code, follow the guidelines in data analysis/readme.pdf

#### Time to run the algorithm 

Running the algorithm takes roughly *32 seconds* in a 6-core Intel Xeon CPU E5-1650v3@3.5GHz with 32GB of DDR4. 
