/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template<typename Traits>
auto KdTreeKNearestIndexQuery<Traits>::begin()
    -> KdTreeKNearestIterator<IndexType, DataPoint>
{
    QueryAccelType::reset();
    QueryType::reset();
    this->search();
    return KdTreeKNearestIterator<IndexType, DataPoint>(QueryType::m_queue.begin());
}

template<typename Traits>
auto KdTreeKNearestIndexQuery<Traits>::end()
    -> KdTreeKNearestIterator<IndexType, DataPoint>
{
    return KdTreeKNearestIterator<IndexType, DataPoint>(QueryType::m_queue.end());
}

template<typename Traits>
void KdTreeKNearestIndexQuery<Traits>::search()
{
    KdTreeQuery<Traits>::search_internal(QueryAccelType::m_kdtree->point_data()[QueryType::input()].pos(),
                                        [](IndexType, IndexType){},
                                        [this](){return QueryType::descentDistanceThreshold();},
                                        [this](IndexType idx){return QueryType::input() == idx;},
                                        [this](IndexType idx, IndexType, Scalar d){QueryType::m_queue.push({idx, d}); return false;}
                                        );
}
