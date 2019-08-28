# Predicting the Frequencies of drug side effects

This repository contains the code to run our matrix decomposition algorithm to predict the frequencies of drug side effects.

The data collection (parse if SIDER 4.1 database) is contained in the data collection folder

## 1. Data Collection.

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

The code requires data from <http://sideeffects.embl.de/> which should be located in a /data/databases/ file.

## 2. Algorithm

Compatibility: the code was implemented in Matlab R2018a. However, we used built-in functions
that had been available since previous versions:

* perfcurve: introduced in R2009a.
* histogram: introduced in R2014b.

To run the code, follow the guidelines below.

1. Open the Matlab script DRIVER.m
2. The code is divided into a series of sections e.g. Initialization. It is recommended to run
section by section to follow the code. To run a given section, click on the section of the code
you want to run and then go to the Editor panel and click on Run Section.
(Alternatively, run the whole script but you will get all the outputs of each section
sequentially).
3. Sections are dependant of each other, so make sure to run each section before the next one.
Run sections in the following order:
  a. Initialization.
  b. Section 1. Load dataset.
  c. Section 2. Train the model.
  d. Section 3. Evaluation in HeldOut test set.
  e. Section 4. Evaluation in PostMarket test set.
  f. Section 5. Reproduce Figure 2 of the main paper.

### Time to run the algorithm 

Running the algorithm takes roughly 32 seconds in a 6-core Intel Xeon CPU E5-1650v3@3.5GHz with 32GB of DDR4. 
