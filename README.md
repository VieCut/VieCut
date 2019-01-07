# VieCut 1.00 - Shared-memory Minimum Cuts

This is the code repository to accompany our papers:

*Henzinger, M., Noe, A., Schulz, C. and Strash, D., 2018. Practical Minimum Cut Algorithms. arXiv preprint [arXiv:1708.06127.](https://arxiv.org/abs/1708.06127)*
*Henzinger, M., Noe, A. and Schulz, C., 2018. Shared-memory Exact Minimum Cuts. arXiv preprint [arXiv:1808.05458.](https://arxiv.org/abs/1808.05458)*

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

@article{henzinger2018shared,
  title={Shared-memory Exact Minimum Cuts},
  author={Henzinger, Monika and Noe, Alexander and Schulz, Christian},
  journal={arXiv preprint arXiv:1808.05458},
  year={2018}
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
If you haven't installed these dependencies, please do so via your package manager

```
sudo apt install gcc-7 g++-7 cmake
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

All of our programs are compiled both for single threaded and shared-memory parallel use. The name of the parallel executable is indicated by appending it with `_parallel`.
When compiling with `./compile.sh -DSAVE_CUTS=ON` the algorithms do not only store the value of the minimum cut but also write a cut file to disk.
This cut file has one line per vertex, either '0' or '1', depending on which side of the cut the vertex is.
The executables can be found in subfolder `/build`.


# Running the programs

## `mincut`

The main executable in our program is `mincut`.
This executable can be used to compute the minimum cut of a given graph with different algorithms.
We use the [METIS graph format](http://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html) for all graphs.
Run any minimum cut algorithm using the following command

```
./build/mincut [options] /path/to/graph.metis <algorithm>
```

For <algorithm> use one of the following:

* `vc` - `VieCut` [HNSS'18]
* `noi` - Algorithm of Nagamochi et al. [NOI'94]
* `ks` - Algorithm of Karger and Stein [KS'96]
* `sw` - Algorithm of Stoer and Wagner [SW'97]
* `pr` - Repeated application of Padberg-Rinaldi contraction rules [PR'91]

when parallelism is enabled, use one of the following:

* `vc` - shared-memory parallel `VieCut` [HNSS'18]
* `exact` - exact shared-memory parallel minimum cut [HNS'18]

#### (Optional) Program Options:

* `-q` - Priority queue implementation ('`bqueue`, `bstack`, `heap`, see [HNS'18] for details)
* `-i` - Number of iterations (default: 1)
* `-l` - Disable limiting of values in priority queue (only relevant for `noi` and `exact`, see [HNS'18])
* `-p` - Use `p` processors (multiple values possible)

The following command

```
./build/mincut -q bqueue -i 3 -p 2 -p 12 /path/to/my/graph.metis exact
```

runs algorithm `exact` using the `BQueue` priority queue implementation for 3 iterations both with 2 and 12 processors.
For each of the runs we print running time and results, as well as a few informations about the graph and algorithm configuration.

## Other Executables

### `kcore`

As most real-world graphs contain vertices with degree 1 and multiple connected components, finding the minimum cut is
as easy as finding the minimum degree or checking whether the graph has multiple connected components.
In order to create harder instances for [HNS'18] and [HNSS'18] we use the cores decomposition of the graph.
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

* `-q` - Priority queue implementation ('`bqueue`, `bstack`, `heap`, see [HNS'18] for details)
* `-i` - Number of iterations (default: 1)
* `-l` - Disable limiting of values in priority queue (only relevant for `noi` and `exact`, see [HNS'18])
* `-p` - Use `p` processors (multiple values possible)
* `-c` - Contraction factor: we contract until only n*(1-c) vertices are left.

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

[HNS'18] - *Henzinger, M., Noe, A. and Schulz, C., 2018. Shared-memory Exact Minimum Cuts.*

[HNSS'18] - *Henzinger, M., Noe, A., Schulz, C. and Strash, D., 2018. Practical Minimum Cut Algorithms.*

[KS'96] - *Karger, D. and Stein, C., 1996. A new approach to the minimum cut problem.*

[NOI'94] - *Nagamochi, H., Ono, T. and Ibaraki, T., 1994. Implementing an efficient minimum capacity cut algorithm.*

[PR'91] - *Padberg, M. and Rinaldi G., 1991. An efficient algorithm for the minimum capacity cut problem.*

[SW'97] - *Stoer, M. and Wagner, F., 1997. A simple min-cut algorithm.*



