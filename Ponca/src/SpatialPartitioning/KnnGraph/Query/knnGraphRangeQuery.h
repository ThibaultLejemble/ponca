/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../../query.h"
#include "../Iterator/knnGraphRangeIterator.h"

#include <vector>
#include <stack>

namespace Ponca {
template <typename Traits> class KnnGraphBase;

template <typename Traits>
class KnnGraphRangeQuery : public RangeIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>
{
protected:
    using QueryType = RangeIndexQuery<typename Traits::IndexType, typename Traits::DataPoint::Scalar>;
    friend class KnnGraphRangeIterator<Traits>; // This type must be equal to KnnGraphRangeQuery::Iterator

public:
    using DataPoint  = typename Traits::DataPoint;
    using IndexType  = typename Traits::IndexType;
    using Scalar     = typename DataPoint::Scalar;
    using VectorType = typename DataPoint::VectorType;
    using Iterator   = KnnGraphRangeIterator<Traits>;

public:
    inline KnnGraphRangeQuery(const KnnGraphBase<Traits>* graph, Scalar radius, int index):
            QueryType(radius, index),
            m_graph(graph),
            m_flag(graph->size()),
            m_stack() {}

public:
    inline Iterator begin(){
        KnnGraphRangeIterator it(this);
        this->initialize(it);
        this->advance(it);
        return it;
    }

    inline Iterator end(){
        return KnnGraphRangeIterator(this, m_graph->size());
    }

protected:
    inline void initialize(Iterator& iterator){
        m_flag.resize(m_graph->size());
        std::fill(m_flag.begin(), m_flag.end(), false);

        PONCA_DEBUG_ASSERT(m_stack.empty());
        m_stack.push(QueryType::input());
        m_flag[QueryType::input()] = true;

        iterator.m_index = -1;
    }

    inline void advance(Iterator& iterator){
        const auto& points  = m_graph->m_kdTreePoints;
        const auto& point   = points[QueryType::input()].pos();

        if(! (iterator != end())) return;

        if(m_stack.empty())
        {
            iterator = end();
        }
        else
        {
            int idx_current = m_stack.top();
            m_stack.pop();

            PONCA_DEBUG_ASSERT((point - points[idx_current]).squaredNorm() < m_squared_radius);

            iterator.m_index = idx_current;

            for(int idx_nei : m_graph->k_nearest_neighbors(idx_current))
            {
                PONCA_DEBUG_ASSERT(idx_nei>0);
                Scalar d  = (point - points[idx_nei].pos()).squaredNorm();
                Scalar th = QueryType::descentDistanceThreshold();
                if(!m_flag[idx_nei] && (point - points[idx_nei].pos()).squaredNorm() < QueryType::descentDistanceThreshold())
                {
                    m_flag[idx_nei] = true;
                    m_stack.push(idx_nei);
                }
            }
            if (iterator.m_index == QueryType::input()) advance(iterator); // query is not included in returned set
        }
    }

protected:
    const KnnGraphBase<Traits>*   m_graph {nullptr};
    std::vector<bool> m_flag;  ///< hold ids status (ids range from 0 to point cloud size)
    std::stack<int>   m_stack; ///< hold ids (ids range from 0 to point cloud size)
};

} // namespace Ponca
