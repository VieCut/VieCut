/******************************************************************************
 * excess_scaling.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once

#include "coarsening/unit_flow.h"
#include "data_structure/flow_graph.h"
#include "data_structure/graph_access.h"
#include "definitions.h"
#include "tlx/logger.hpp"

#include <cstdint>
#include <cstdlib>
#include <functional>
#include <queue>

class excess_scaling
{

    static constexpr bool debug = false;

public:
    excess_scaling() { }
    virtual ~excess_scaling() { }

    void init(flow_graph& fg, std::vector<EdgeWeight>& delta_src,
              float tau, FlowType U, NodeID max_height) {
        m_fg = &fg;
        m_delta_src = delta_src;
        m_tau = tau;
        m_U = U;
        m_max_height = max_height;
        m_F = 0.f;

        for (size_t n = 0; n < fg.number_of_nodes(); ++n) {
            float div = m_delta_src[n] / (float)(2 * fg.getCapacity(n));
            m_F = std::max(m_F, div);
        }
        m_mu = std::ceil(m_F);
        m_j = 0;
        m_flows.resize(fg.number_of_edges(), 0);
        m_fj.resize(fg.number_of_nodes(), 0);
    }

    void run() {
        LOG << "START: graph with n=" << m_fg->number_of_nodes()
            << " m=" << m_fg->number_of_edges();

        std::vector<EdgeWeight> flowsrc(m_delta_src.size());
        for (size_t i = 0; i < m_delta_src.size(); ++i) {
            // ceil(delta div mu)
            flowsrc[i] = (m_delta_src[i] + m_mu - 1) / m_mu;
        }
        while (m_mu > 4) {

            size_t discard = 0;
            unit_flow uf;

            uf.init(*m_fg, flowsrc, m_U, m_max_height, 2);

            uf.run();
            size_t ctr = 0;

            for (NodeID n : m_fg->nodes()) {

                for (EdgeID e : m_fg->edges_of(n)) {
                    EdgeID rev = m_fg->getReverseEdge(n, e);
                    NodeID tgt = m_fg->getEdgeTarget(n, e);

                    m_flows[ctr++] += m_fg->getEdgeFlow(n, e) * m_mu;
                    m_fg->setEdgeFlow(n, e, 0);
                    m_fg->setEdgeFlow(tgt, rev, 0);
                }

                if (uf.excess(n) > 0) {
                    discard += uf.excess(n);
                    m_delta_src[n] = m_delta_src[n] - uf.excess(n);
                    flowsrc[n] = m_delta_src[n] / m_mu;
                }
            }
            m_mu /= 2;
        }
    }

private:
private:
    flow_graph* m_fg;
    std::vector<EdgeWeight> m_delta_src;
    std::vector<FlowType> m_flows;
    std::vector<FlowType> m_fj;
    float m_tau;
    FlowType m_U;
    NodeID m_max_height;
    float m_F;
    FlowType m_mu;
    size_t m_j;
};
