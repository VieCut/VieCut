/******************************************************************************
 * multiterminal_cut_test.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <type_traits>

#include "algorithms/multicut/multiterminal_cut.h"
#include "gperftools/malloc_extension.h"
#include "gtest/gtest.h"
#include "io/graph_io.h"

class MultiterminalCutTest : public ::testing::Test {
 public:
    static void SetUpTestCase() {
        MPI_Init(NULL, NULL);
    }

    static void TearDownTestCase() {
        MPI_Finalize();
    }
};

TEST_F(MultiterminalCutTest, FourClusters) {
    std::vector<size_t> sizes = { 5, 10, 50, 100 };
    for (size_t cluster_size : sizes) {
        std::shared_ptr<graph_access> G = std::make_shared<graph_access>();

        G->start_construction(4 * cluster_size,
                              2 * cluster_size * (cluster_size - 1) * 4 + 12);

        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                if (i != j) {
                    G->new_edge(i * cluster_size, j * cluster_size);
                }
            }

            for (size_t j = 0; j < cluster_size; ++j) {
                for (size_t k = 0; k < cluster_size; ++k) {
                    if (j != k) {
                        NodeID base = cluster_size * i;
                        G->new_edge(base + j, base + k);
                    }
                }
            }
        }

        G->finish_construction();
        auto mG = mutable_graph::from_graph_access(G);
        std::vector<NodeID> terminals;
        std::random_device rd;
        std::mt19937 eng(rd());
        std::uniform_int_distribution<> distribution(0, cluster_size - 1);
        for (size_t i = 0; i < 4; ++i) {
            terminals.emplace_back(i * cluster_size + distribution(eng));
        }

        multiterminal_cut mct;
        FlowType f = mct.multicut(mG, terminals);

        ASSERT_EQ(f, (FlowType)6);
    }
}

TEST_F(MultiterminalCutTest, FourClustersWeighted) {
    std::vector<size_t> sizes = { 20, 50, 100 };
    for (size_t cluster_size : sizes) {
        std::shared_ptr<graph_access> G = std::make_shared<graph_access>();

        G->start_construction(4 * cluster_size,
                              2 * cluster_size * (cluster_size - 1) * 4 + 12);

        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                if (i != j) {
                    EdgeID e = G->new_edge(i * cluster_size, j * cluster_size);
                    G->setEdgeWeight(e, i + j);
                }
            }

            for (size_t j = 0; j < cluster_size; ++j) {
                for (size_t k = 0; k < cluster_size; ++k) {
                    if (j != k) {
                        NodeID base = cluster_size * i;
                        G->new_edge(base + j, base + k);
                    }
                }
            }
        }

        G->finish_construction();
        auto mG = mutable_graph::from_graph_access(G);

        multiterminal_cut mct;

        std::vector<NodeID> terminals;

        std::random_device rd;
        std::mt19937 eng(rd());
        std::uniform_int_distribution<> distribution(0, cluster_size - 1);
        for (size_t i = 0; i < 4; ++i) {
            terminals.emplace_back(i * cluster_size + distribution(eng));
        }

        FlowType f = mct.multicut(mG, terminals);

        ASSERT_EQ(f, (FlowType)18);
    }
}

TEST_F(MultiterminalCutTest, FourClustersMinCutUnequal) {
    std::vector<size_t> sizes = { 1, 4, 10, 50, 100 };
    for (size_t cluster_size : sizes) {
        std::shared_ptr<graph_access> G = std::make_shared<graph_access>();

        G->start_construction(8 * cluster_size,
                              2 * cluster_size * (cluster_size - 1) * 8 + 40);

        for (size_t i = 0; i < 8; ++i) {
            for (size_t j = 0; j < 8; ++j) {
                if (i / 2 == j / 2 && i != j) {
                    EdgeID e = G->new_edge(i * cluster_size, j * cluster_size);
                    G->setEdgeWeight(e, 3);
                }

                if (i != j && i % 2 == 1 && j % 2 == 1 && i + j != 8) {
                    EdgeID e = G->new_edge(i * cluster_size, j * cluster_size);
                    G->setEdgeWeight(e, 2);
                }
            }

            for (size_t j = 0; j < cluster_size; ++j) {
                for (size_t k = 0; k < cluster_size; ++k) {
                    if (j != k) {
                        NodeID base = cluster_size * i;
                        G->new_edge(base + j, base + k);
                    }
                }
            }
        }

        G->finish_construction();
        auto mG = mutable_graph::from_graph_access(G);

        multiterminal_cut mct;

        std::vector<NodeID> terminals;

        std::random_device rd;
        std::mt19937 eng(rd());
        std::uniform_int_distribution<> distribution(0, cluster_size - 1);
        for (size_t i = 0; i < 4; ++i) {
            terminals.emplace_back((2 * i) * cluster_size + distribution(eng));
        }

        FlowType f = mct.multicut(mG, terminals);

        ASSERT_EQ(f, (FlowType)8);
    }
}

TEST_F(MultiterminalCutTest, TotallyDisconnected) {
    std::vector<size_t> sizes = { 1, 5, 10, 50, 100 };
    for (size_t cluster_size : sizes) {
        std::shared_ptr<graph_access> G = std::make_shared<graph_access>();

        G->start_construction(4 * cluster_size,
                              2 * cluster_size * (cluster_size - 1) * 4);

        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < cluster_size; ++j) {
                for (size_t k = 0; k < cluster_size; ++k) {
                    if (j != k) {
                        NodeID base = cluster_size * i;
                        G->new_edge(base + j, base + k);
                    }
                }
            }
        }

        G->finish_construction();
        auto mG = mutable_graph::from_graph_access(G);

        std::vector<NodeID> terminals;
        std::random_device rd;
        std::mt19937 eng(rd());
        std::uniform_int_distribution<> distribution(0, cluster_size - 1);
        for (size_t i = 0; i < 4; ++i) {
            terminals.emplace_back(i * cluster_size + distribution(eng));
        }
        multiterminal_cut mct;

        FlowType f = mct.multicut(mG, terminals);

        ASSERT_EQ(f, (FlowType)0);
    }
}

TEST_F(MultiterminalCutTest, PartiallyDisconnected) {
    std::vector<size_t> sizes = { 1, 5, 10, 50, 100 };
    for (size_t cluster_size : sizes) {
        std::shared_ptr<graph_access> G = std::make_shared<graph_access>();

        G->start_construction(4 * cluster_size,
                              2 * cluster_size * (cluster_size - 1) * 4 + 4);

        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                if (i / 2 == j / 2) {
                    G->new_edge(i * cluster_size, j * cluster_size);
                }
            }

            for (size_t j = 0; j < cluster_size; ++j) {
                for (size_t k = 0; k < cluster_size; ++k) {
                    if (j != k) {
                        NodeID base = cluster_size * i;
                        G->new_edge(base + j, base + k);
                    }
                }
            }
        }

        G->finish_construction();
        auto mG = mutable_graph::from_graph_access(G);

        std::vector<NodeID> terminals;
        std::random_device rd;
        std::mt19937 eng(rd());
        std::uniform_int_distribution<> distribution(0, cluster_size - 1);
        for (size_t i = 0; i < 4; ++i) {
            terminals.emplace_back(i * cluster_size + distribution(eng));
        }
        multiterminal_cut mct;

        FlowType f = mct.multicut(mG, terminals);

        ASSERT_EQ(f, (FlowType)2);
    }
}
