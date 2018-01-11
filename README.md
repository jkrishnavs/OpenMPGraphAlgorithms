1.Introduction

Realization of popular Graph Algorithms using C and OpenMP.
The base algorithms are taken from GreenMarl[1] example algorithms.







2. Algorithm Details
i. communities
ii. conduct
iii. pagerank
iv. sssp
v. triangle_counting


3. Folder Structure
   inputGraphs: a set of real world and generated graphs in EDGE_LIST[1] format. We only handle EDGE_LIST format. 



4. Input Graphs
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



4.Licence
============== GreenMarl License Copy =============================
Copyright (c) 2011-2012 Stanford University, unless otherwise specified.
All rights reserved.

This software was developed by the Pervasive Parallelism Laboratory of
Stanford University, California, USA.

Permission to use, copy, modify, and distribute this software in source
or binary form for any purpose with or without fee is hereby granted,
provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright
         notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
         notice, this list of conditions and the following disclaimer in the
	       documentation and/or other materials provided with the distribution.

   3. Neither the name of Stanford University nor the names of its
         contributors may be used to endorse or promote products derived
	       from this software without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

5. Compilation
   i. make
   ii. make intermediate
   iii: make debug:
   

   Place Holder header:
   /include/energylib.h

   You can replace the function in energylib to track power
   and energy consumption of the kernel.




5.References
[1] https://github.com/stanford-ppl/Green-Marl
[2] https://snap.stanford.edu/data/
[3] SNAP: A General-Purpose Network Analysis and Graph-Mining Library


