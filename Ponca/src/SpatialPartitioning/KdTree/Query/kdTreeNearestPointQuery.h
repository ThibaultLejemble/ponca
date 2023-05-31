/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "../kdTreeQuery.h"
#include "../../query.h"
#include "../Iterator/kdTreeNearestIterator.h"

namespace Ponca {

template <class DataPoint, class Adapter>
class KdTreeNearestPointQuery : public KdTreeQuery<DataPoint, Adapter>,
    public NearestPointQuery<typename Adapter::IndexType, DataPoint>
{
public:
    using IndexType       = typename Adapter::IndexType;
    using Scalar          = typename DataPoint::Scalar;
    using VectorType      = typename DataPoint::VectorType;
    using QueryType       = NearestPointQuery<IndexType, DataPoint>;
    using QueryAccelType  = KdTreeQuery<DataPoint, Adapter>;

    KdTreeNearestPointQuery(const KdTree<DataPoint, Adapter>* kdtree, const VectorType& point) :
        KdTreeQuery<DataPoint, Adapter>(kdtree), NearestPointQuery<IndexType, DataPoint>(point)
    {
    }

public:
    KdTreeNearestIterator<IndexType> begin();
    KdTreeNearestIterator<IndexType> end();

protected:
    void search();
};

#include "./kdTreeNearestPointQuery.hpp"
} // namespace ponca