/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

namespace Ponca {

class KnnGraphRangeQuery;

class KnnGraphRangeIterator
{
protected:
    friend class KnnGraphRangeQuery;

public:
    KnnGraphRangeIterator();
    KnnGraphRangeIterator(KnnGraphRangeQuery* query);
    KnnGraphRangeIterator(KnnGraphRangeQuery* query, int index);

public:
    bool operator != (const KnnGraphRangeIterator& other) const;
    void operator ++ ();
    int  operator *  () const;

protected:
    KnnGraphRangeQuery* m_query;
    int m_index;
};

} // namespace Ponca
