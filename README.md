# VieCut 1.00 - Shared-memory Minimum Cuts

[![Build Status](https://travis-ci.com/alexnoe/VieCut.svg?branch=master)](https://travis-ci.com/alexnoe/VieCut)

This is the code repository to accompany our papers:

*Henzinger, M., Noe, A., Schulz, C. and Strash, D., 2018. Practical Minimum Cut Algorithms. arXiv preprint [arXiv:1708.06127.](https://arxiv.org/abs/1708.06127)*

*Henzinger, M., Noe, A. and Schulz, C., 2019. Shared-memory Exact Minimum Cuts. arXiv preprint [arXiv:1808.05458.](https://arxiv.org/abs/1808.05458)*

*Henzinger, M, Noe, A. and Schulz, C., 2019. Shared-Memory Branch-and-Reduce for Multiterminal Cuts. arXiv preprint [arXiv:1908.04141.](https://arxiv.org/abs/1908/04141)*

The papers can be freely accessed online in the arXiv.

If you use this code in the context of an academic publication, we ask that you cite the applicable papers:
```bibtex
@inproceedings{henzinger2018practical,
  title={Practical Minimum Cut Algorithms},
  author={Henzinger, Monika and Noe, Alexander and Schulz, Christian and Strash, Darren},
  booktitle={2018 Proceedings of the Twentieth Workshop on Algorithm Engineering and Experiments (ALENEX)},
  pages={48--61},
  year={2018},
  organization={SIAM}
}

@inproceedings{henzinger2019shared,
  author = {Henzinger, Monika and Noe, Alexander and Schulz, Christian},
  title = {{Shared-memory Exact Minimum Cuts}},
  booktitle={Proceedings of the 33rd International Parallel and Distributed Processing Symposium (IPDPS)},
  year = {2019}
}

@article{henzinger2019branch,
  title={Shared-Memory Branch-and-Reduce for Multiterminal Cuts},
  author={Henzinger, Monika and Noe, Alexander and Schulz, Christian},
  journal={arXiv preprint arXiv:1908.04141},
  year={2019}
}
```

## Introduction

The minimum cut problem for an undirected edge-weighted graph asks us to divide its set of nodes into two blocks while minimizing the weight sum of the cut edges.
 It is a fundamental graph problem with many applications in different fields, such as network reliability,
 where assuming equal failure probability edges, the smallest edge cut in the network has the highest chance to disconnect the network;
 in VLSI design, where the minimum cut can be used to minimize the number of connections between microprocessor blocks;
 and as a subproblem in branch-and-cut algorithms for the Traveling Salesman Problem (TSP) and other combinatorial problems.

In our work, we present fast shared-memory parallel algorithms for the minimum cut problem.
In the first paper, *Practical Minimum Cut Algorithm*, we present the fast shared-memory parallel heuristic algorithm `VieCut`.
While the algorithm can not guarantee solution quality, in practice it usually outputs cuts that are minimal or very close.

Based on this algorithm and the algorithm of Nagamochi et al.,
in our paper *Shared-memory Exact Minimum Cuts* we present a fast shared-memory exact algorithm for the minimum cut problem.
These algorithms significantly outperform the state of the art.

## Installation

### Prerequisites

In order to compile the code you need a version of the GCC that supports `c++-17`, such as `g++-7`, and `cmake`. 
As we use TCMalloc to allocate and deallocate memory, you need to have `google-perftools` installed.
If you haven't installed these dependencies, please do so via your package manager

```
sudo apt install gcc-7 g++-7 cmake libgoogle-perftools-dev
```

### Compiling

To compile the code use the following commands

```
  git submodule update --init --recursive
  mkdir build
  cd build
  cmake ..
  make
```

We also offer a compile script `compile.sh` which compiles the executables and runs tests.

All of our programs are compiled both for single threaded and shared-memory parallel use. 
The name of the parallel executable is indicated by appending it with `_parallel`. 
The executables can be found in subfolder `build`.


# Running the programs

## `mincut`

The main executable in our program is `mincut`. The shared-memory parallel executable using OpenMP is `mincut_parallel`
This executable can be used to compute the minimum cut of a given graph with different algorithms.
We use the [METIS graph format](http://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html) for all graphs.
Run any minimum cut algorithm using the following command

```
./build/mincut [options] /path/to/graph.metis <algorithm>

./build/mincut_parallel [options] /path/to/graph.metis <algorithm>
```

For <algorithm> use one of the following:

* `vc` - `VieCut` [HNSS'18]
* `noi` - Algorithm of Nagamochi et al. [NOI'94]
* `ks` - Algorithm of Karger and Stein [KS'96]
* `matula` - Approximation Algorithm of Matula [Matula'93]
* `pr` - Repeated application of Padberg-Rinaldi contraction rules [PR'91]
* `cactus` - Find _all_ minimum cuts and give the cactus that represents them. [ongoing work]

when parallelism is enabled, use one of the following:

* `inexact` - shared-memory parallel version of `VieCut` [HNSS'18]
* `exact` - exact shared-memory parallel minimum cut [HNS'19a]
* `cactus` - Find _all_ minimum cuts and give the cactus that represents them. [ongoing work]

#### (Optional) Program Options:

* `-q` - Priority queue implementation ('`bqueue`, `bstack`, `heap`, see [HNS'19a] for details)
* `-i` - Number of iterations (default: 1)
* `-l` - Disable limiting of values in priority queue (only relevant for `noi` and `exact`, see [HNS'19a])
* `-p` - [Only for `mincut_parallel`] Use `p` processors (multiple values possible)
* `-s` - Compute and save minimum cut. The cut will be written to disk in a file which contains one line per node, either `0` or `1` depending on which side of the cut the node is.
* `-o` - [`-s` needs to be enabled as well] Path of output file. If this is set, we print the minimum cut to file.
* `-b` - [Only for algorithm `cactus`, `-s` needs to be enabled as well] Find most balanced minimum cut and print its balance.

The following command

```
./build/mincut_parallel -q bqueue -i 3 -p 2 -p 12 /path/to/my/graph.metis exact
```

runs algorithm `exact` using the `BQueue` priority queue implementation for 3 iterations both with 2 and 12 processors.
For each of the runs we print running time and results, as well as a few informations about the graph and algorithm configuration.


### `multiterminal_cut`

The multiterminal cut of a graph G and a set of terminals T is to find the minimum cut of G that separates all terminals from each other. 
For |T|=2 this is equal to the minimum s-t-cut problem, for |T|>2, the problem is NP-hard. 
We solve the problem using a branch-and-reduce approach which finds a hard kernel by applying reduction rules and then branching on edges adjacent to a terminal. For a more detailed description, we refer the reader to [HNS'19b]. The algorithm is shared-memory parallel and uses OpenMP.

#### Usage:

```
./build/multiterminal_cut /path/to/graph.metis [options] 
```

#### (Optional) Program Options:
* `-f` - Path to partition file. This file has one line for each vertex. A value of 0 to |T|-1 indicates which terminal a vertex belongs to, otherwise a value of |T| indicates that the vertex does not belong to a terminal.
* `-t` - Add vertex `t` as a terminal.
* `-k` - Find multiterminal cut between `k` vertices with highest vertex degree.
* `-r` - Find multiterminal cut between `r` random vertices.
* `-b` - Run BFS around each terminal and add up to `b` vertices discovered first to each terminal.
* `-p` - Number of threads (default: OMP_NUM_THREADS, which defaults to the number of hardware threads).
* `-c` - Disable kernelization variants [values in 0-4] (default: 0 - all enabled). 


The following command

```
./build/multiterminal_cut /path/to/my/graph.metis -f /path/to/my/partition_file -p 12 -c 1
```

finds the multiterminal cut as given by the graph and partition file using 12 threads and all except for the high-connectivity kernelization (see [HNS'19b] for further details) rules enabled.





## Other Executables

### `kcore`

As most real-world graphs contain vertices with degree 1 and multiple connected components, finding the minimum cut is
as easy as finding the minimum degree or checking whether the graph has multiple connected components.
In order to create harder instances for [HNS'19a] and [HNSS'18] we use the cores decomposition of the graph.
The k-core of a graph is the largest subgraph of the graph, in which every node has at least degree k in the k-core.
We use the executable `kcore` to find k-cores of a graph where the minimum cut is not equal to the minimum degree.
If the minimum cut is not equal to the minimum degree, the k-core graph is written both in METIS and in DIMACS format.

#### Usage:

```
./build/kcore <options> /path/to/graph.metis
```

with following options:

* `-l` - search for the lowest value of k where the minimum cut of the k-core is not equal to the minimum degree
* `-c` - disable testing for the minimum cut, just compute cores decomposition and write k-core graphs to disk
* `-k` - compute k-core for k

For example, to compute the 5- and 10-cores of a graph without minimum cut testing, use the following command

```
./build/kcore -c -k 5 -k 10 /path/to/graph.metis
```

### `mincut_contract`

The executable `mincut_contract` runs a version of `mincut` that begins with contracting random edges.
This idea is taken from the algorithm or Karger and Stein [KS96], which contracts random edges and recurses on the contracted graph.
Mincut contract begins by contracting random edges until a user-defined number of vertices is left.
Afterwards we run a minimum cut algorithm on the contracted graph.
Usage is similar to mincut and can be combined with any minimum cut algorithm.

```
./build/mincut_contract [options] /path/to/graph.metis <algorithm>
```

#### Program Options:

* `-q` - Priority queue implementation ('`bqueue`, `bstack`, `heap`, see [HNS'19a] for details)
* `-i` - Number of iterations (default: 1)
* `-l` - Disable limiting of values in priority queue (only relevant for `noi` and `exact`, see [HNS'19a])
* `-p` - Use `p` processors (multiple values possible)
* `-c` - Contraction factor: we contract until only n*(1-c) vertices are left.
* `-s` - Compute and save minimum cut.

### `mincut_recursive`

The executable `mincut_recusive` runs an exact minimum cut algorithm (`exact` in parallel, `noi` otherwise)
on a graph and creates a graph for the largest SCC of the larger block of the cut.
We then run the minimum cut algorithm on this graph and repeat this process until the minimum cut is equal to the minimum degree.
This can be used to create subgraphs of a graph that have different minimum cuts.

```
./build/mincut_recursive [options] /path/to/graph.metis
```

#### Program Options:

* `-o` - Write all graphs to disk (DIMACS and METIS format) where the minimum cut is larger than the minimum cut of the previous graph.

## References

[BZ'03] - *Batagelj, V. and Zaversnik, M., 2003. An O(m) algorithm for cores decomposition of networks.*

[HNS'19a] - *Henzinger, M., Noe, A. and Schulz, C., 2019. Shared-memory Exact Minimum Cuts.*

[HNS'19b] - *Henzinger, M., Noe, A. and Schulz, C., 2019. Shared-memory Branch-and-reduce for Multiterminal Cuts.*

[HNSS'18] - *Henzinger, M., Noe, A., Schulz, C. and Strash, D., 2018. Practical Minimum Cut Algorithms.*

[KS'96] - *Karger, D. and Stein, C., 1996. A new approach to the minimum cut problem.*

[Matula'93] - *Matula, D., 1993. A linear time (2 + ε)-approximation algorithm for edge connectivity*

[NOI'94] - *Nagamochi, H., Ono, T. and Ibaraki, T., 1994. Implementing an efficient minimum capacity cut algorithm.*

[PR'91] - *Padberg, M. and Rinaldi G., 1991. An efficient algorithm for the minimum capacity cut problem.*

[SW'97] - *Stoer, M. and Wagner, F., 1997. A simple min-cut algorithm.*



