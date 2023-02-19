/******************************************************************************
 * cactus_mincut.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/global_mincut/cactus/most_balanced_minimum_cut.h"
#include "algorithms/global_mincut/cactus/recursive_cactus.h"
#include "algorithms/global_mincut/minimum_cut.h"
#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "algorithms/global_mincut/viecut.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "io/graph_io.h"
#include "tools/string.h"
#include "tools/timer.h"

#ifdef PARALLEL
#include "parallel/coarsening/contract_graph.h"
#else
#include "coarsening/contract_graph.h"
#endif

template <class GraphPtr>
class cactus_mincut : public minimum_cut {
 public:
    typedef GraphPtr GraphPtrType;
    cactus_mincut() { }
    virtual ~cactus_mincut() { }
    static constexpr bool debug = false;

    bool timing = configuration::getConfig()->verbose;

    EdgeWeight perform_minimum_cut(GraphPtr G) {
        // compatibility with min cut interface
        return std::get<0>(findAllMincuts(G));
    }

    std::tuple<EdgeWeight, mutableGraphPtr,
               std::vector<std::pair<NodeID, EdgeID> > >
    findAllMincuts(GraphPtr G, EdgeWeight known_mincut = UNDEFINED_NODE) {
        std::vector<GraphPtr> v = { G };
        return findAllMincuts(v, known_mincut);
    }

    std::tuple<EdgeWeight, mutableGraphPtr,
               std::vector<std::pair<NodeID, EdgeID> > >
    findAllMincuts(std::vector<GraphPtr> graphs,
                   EdgeWeight known_mincut = UNDEFINED_NODE) {
        if (graphs.size() == 0 || !graphs.back()) {
            mutableGraphPtr empty;
            return std::make_tuple(
                -1, empty, std::vector<std::pair<NodeID, EdgeID> > { });
        }
        recursive_cactus<GraphPtr> rc;
        EdgeWeight mincut = graphs.back()->getMinDegree();
        timer t;
        if (known_mincut == UNDEFINED_NODE) {
            viecut<GraphPtr> vc;
            mincut = vc.perform_minimum_cut(graphs.back());
        } else {
            mincut = known_mincut;
        }
        noi_minimum_cut<GraphPtr> noi;
        std::vector<std::vector<std::pair<NodeID, NodeID> > > guaranteed_edges;
        std::vector<size_t> ge_ids;

        minimum_cut_helpers<GraphPtr>::setInitialCutValues(graphs);

        NodeID previous_size = UNDEFINED_NODE;
        while (graphs.back()->number_of_nodes() * 1.01 < previous_size) {
            previous_size = graphs.back()->number_of_nodes();
            auto current_graph = graphs.back();
            EdgeWeight current_mincut = mincut;
            ge_ids.emplace_back(graphs.size() - 1);
            guaranteed_edges.emplace_back();
            LOGC(timing) << "t " << t.elapsed()
                         << " n " << current_graph->number_of_nodes()
                         << " m " << current_graph->number_of_edges()
                         << " cut " << mincut;

            std::vector<std::pair<NodeID, NodeID> > contractable;

            timer time;
            auto uf = noi.modified_capforest(current_graph, mincut + 1);

            for (NodeID n : current_graph->nodes()) {
                EdgeID e = current_graph->get_first_edge(n);
                if (current_graph->get_first_invalid_edge(n) - e == 1) {
                    if ((current_graph->getEdgeWeight(n, e) == mincut)
                        && uf.n() > 1) {
                        NodeID t = current_graph->getEdgeTarget(n, e);
                        uf.Union(n, t);
                        guaranteed_edges.back().emplace_back(n, t);
                    }
                }
            }

            LOGC(timing) << "t-noi " << t.elapsed() << " contract "
                         << graphs.back()->number_of_nodes()
                         << " to " << uf.n();

            if (uf.n() < current_graph->number_of_nodes()) {
                auto newg =
                    contraction::fromUnionFind(current_graph, &uf, true);
                graphs.emplace_back(newg);
                mincut = minimum_cut_helpers<GraphPtr>::updateCut(
                    graphs, mincut);
            }

            union_find uf12 = tests::prTests12(
                graphs.back(), mincut + 1, true);
            LOGC(timing) << "t12 " << t.elapsed() << " contract "
                         << graphs.back()->number_of_nodes()
                         << " to " << uf12.n();
            if (uf12.n() < graphs.back()->number_of_nodes()) {
                auto g12 = contraction::fromUnionFind(
                    graphs.back(), &uf12, true);
                graphs.push_back(g12);
                mincut = minimum_cut_helpers<GraphPtr>::updateCut(
                    graphs, mincut);
            }

            union_find uf34 = tests::prTests34(
                graphs.back(), mincut + 1, true);
            LOGC(timing) << "t34 " << t.elapsed() << " contract "
                         << graphs.back()->number_of_nodes()
                         << " to " << uf34.n();
            if (uf34.n() < graphs.back()->number_of_nodes()) {
                auto g34 = contraction::fromUnionFind(
                    graphs.back(), &uf34, true);
                graphs.push_back(g34);
                mincut = minimum_cut_helpers<GraphPtr>::updateCut(
                    graphs, mincut);
            }

            if (current_mincut > mincut) {
                // mincut has improved
                // so all the edges that were in guaranteed edges before
                // are now _not_ cactus edges as there is a lighter cut
                guaranteed_edges.clear();
                ge_ids.clear();
            }
        }

        if (graphs.back()->number_of_nodes() > 1)
            mincut = noi.perform_minimum_cut(graphs.back());

        rc.setMincut(mincut);
        auto out_graph = rc.flowMincut(graphs);  // This is the cactus graph!

        minimum_cut_helpers<GraphPtr>::setVertexLocations(
            out_graph, graphs, ge_ids, guaranteed_edges, mincut);
        LOGC(timing)
            << "unpacked - n " << out_graph->n() << " m " << out_graph->m();

        std::vector<std::pair<NodeID, EdgeID> > mb_edges;
        if (configuration::getConfig()->find_most_balanced_cut) {
            most_balanced_minimum_cut<GraphPtr> mbmc;
            mb_edges = mbmc.findCutFromCactus(out_graph, mincut, graphs[0]);
        }

	if (configuration::getConfig()->cactus_filename != "") {
	    _write_cactus_graphml_file(out_graph, configuration::getConfig()->cactus_filename);
	}

        return std::make_tuple(mincut, out_graph, mb_edges);
    }
};


void _write_cactus_graphml_file(mutableGraphPtr cactus, const std::string& cactus_filename) {
    std::ostringstream oss;

    oss << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
	<< "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n"
	<< "         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
	<< "         xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n";

    oss << "  <key id=\"containedVertices\" for=\"node\" attr.name=\"containedVertices\" attr.type=\"string\"/>\n"
	<< "  <key id=\"weight\" for=\"edge\" attr.name=\"weight\" attr.type=\"double\"/>\n";

    oss << "  <graph id=\"cactus\" edgedefault=\"undirected\">\n";

    // Generate XML for nodes.
    auto nodes = cactus->nodes();

    for (auto node_it = nodes.begin(); node_it != nodes.end(); ++node_it) {
        oss << "    <node id=\"" << *node_it << "\">\n"
	    << "      <data key=\"containedVertices\">";

	auto vertices = cactus->containedVertices(*node_it); // needed as containedVertices() creates copy, not reference

        for (auto vertex_it = vertices.begin(); vertex_it != vertices.end(); ) {
            oss << *vertex_it;
	    if (++vertex_it != vertices.end())
		oss << ",";
        }

        oss << "</data>\n";
	oss << "    </node>\n";
    }

    // Generate XML for edges.
    for (NodeID n : cactus->nodes()) {
        for (EdgeID e : cactus->edges_of(n)) {
	    if (n > cactus->getEdgeTarget(n, e))
		continue;

	    std::ostringstream edge_id;
	    edge_id << n << "-" << cactus->getEdgeTarget(n, e) << "-" << e;

	    oss << "    <edge id=\"" << edge_id.str() << "\" source=\"" << n << "\" target=\"" << cactus->getEdgeTarget(n, e) << "\">\n"
		<< "      <data key=\"weight\">" << cactus->getEdgeWeight(n, e) << "</data>\n"
		<< "    </edge>\n";
        }
    }

    oss << "  </graph>\n";
    oss << "</graphml>\n";

    std::ofstream cactus_file(cactus_filename);
    cactus_file << oss.str();
    cactus_file.close();
}
