# CATBOSS

CATBOSS: Cluster Analysis of Trajectories Based on Segment Splitting

# Pre-requisites:

In order to run CATBOSS, you will need the following dependencies:
1. The SIMPLEchangepoint Python package (https://github.com/DEShawResearch/SIMPLEchangepoint) [1]
2. The MATLAB wrapper for earth mover's distance calculation (https://github.com/garydoranjr/pyemd) [2,3]

# How to run:
1. Prepare the dataset to be analyzed in space-delimited format, with the first column being time, and the remaining columns being the variables of the dataset. You may also use sample data provided here for your convenience.
2. Run CATBOSS_detectChanges.py <input_file> --lam <lambda value> --alpha <alpha value>
(Optional: Import the change probability matrix and data files into your MATLAB workspace and run CATBOSS_slopeAnalysis.m to identify sloped segments)
3. Import the change probability matrix and data files into your MATLAB workspace and run CATBOSS_distanceMatrix.m
4. Import the distance matrix into your MATLAB workspace and run CATBOSS_cluster.m
5. Import the resulting cluster assignment into your MATLAB workspace and run CATBOSS_mapToPoints.m to get final assignments for each data point.

# References(BibTeX):
1. F@article {Fan7454,
	author = {Fan, Zhou and Dror, Ron O. and Mildorf, Thomas J. and Piana, Stefano and Shaw, David E.},
	title = {Identifying localized changes in large systems: Change-point detection for biomolecular simulations},
	volume = {112},
	number = {24},
	pages = {7454--7459},
	year = {2015},
	doi = {10.1073/pnas.1415846112},
	publisher = {National Academy of Sciences},
	issn = {0027-8424},
	URL = {https://www.pnas.org/content/112/24/7454},
	eprint = {https://www.pnas.org/content/112/24/7454.full.pdf},
	journal = {Proceedings of the National Academy of Sciences}
}

2. @Misc{,
  author =    {Gary Doran},
  title =     {{PyEMD}: Earth Mover's Distance for {Python}},
  year =      {2014--},
  url = "https://github.com/garydoranjr/pyemd",
  note = {[Online; accessed <today>]}
}
  
3. @INPROCEEDINGS{Pele-eccv2008,
author = {Ofir Pele and Michael Werman},
title = {A Linear Time Histogram Metric for Improved SIFT Matching},
booktitle = {ECCV},
year = {2008}
}

# License information

Copyright (c) 2021

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
