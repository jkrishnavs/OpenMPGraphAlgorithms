# 1.Introduction

Realization of popular Graph Algorithms using C and OpenMP.
The base algorithms are taken from GreenMarl[1] example algorithms.







# 2. Algorithm Details
i. communities
ii. conduct
iii. pagerank
iv. sssp
v. triangle_counting


# 3. Folder Structure
   inputGraphs: a set of real world and generated graphs in EDGE_LIST[1] format. We only handle EDGE_LIST format. 



# 4. Input Graphs
The input graphs are either downloaded from SNAP data set[2] or RMAT graphs generated using SNAP[3] or random graphs generated using GreenMarl graph gen program[1].

***** skipped due to size **** Please contact for set of input graphs.

SNAP data set graphs:
1. amazon.edge
2. higgs.edge
3. patent.edge
4. pokec.edge
5. pokecr.edge
6. rhiggs.edge
7. roadNet-CA.edge
8. rpatent.edge
9. web-BerkStan.edge
10. wiki-topcats.edge 

RMAT Graphs:
1. rmatv1e5e1e7.edge
2. rmatv7e5e5e6.edge
3. rrmatv1e5e1e7.edge
4. rrmatv7e5e5e6.edge 


Random Graphs:
1. rm1e6u8e6.edge
2. rm1e6u35e6.edge 



# 4.Licence
See GMLICENCE for GreenMarl Licence.

# 5. Compilation
   Run ./configure.sh before running make.
   
   make modes:
   i. make 
   ii. make intermediate
   iii: make debug


   To change chunk size
   make SET_CHUNKSIZE=<newchunkSize>
   

   Place Holder header:
   /include/energylib.h

   You can replace the function in energylib to track power
   and energy consumption of the kernel.




# 6.References
[1] https://github.com/stanford-ppl/Green-Marl  
[2] https://snap.stanford.edu/data/  
[3] SNAP: A General-Purpose Network Analysis and Graph-Mining Library

# 7. Contact
1. Jyothi Krishna V S, IIT Madras (jkrishna@cse.iitm.ac.in)
2. Rupesh Nasre, IIT Madras 


