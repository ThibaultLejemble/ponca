#include "mlsProjection.h"

namespace Ponca {

template<class FitT, class FitFinalT>
template<class NeighborQueryT>
typename QuasiOrthogonalMLSProjection<FitT,FitFinalT>::VectorType
QuasiOrthogonalMLSProjection<FitT,FitFinalT>::project(const VectorType& p, const NeighborQueryT& nei)
{
    auto q = p; // moving point

    FitT fit;
    m_fitFinal = FitFinalT();

    fit.setWeightFunc(WFunctor(m_scale));
    m_fitFinal.setWeightFunc(WFunctorFinal(m_scale));

    const auto epsilon = m_convRatio * m_scale;
    const auto epsilon2 = epsilon * epsilon;

    // Three stopping conditions:
    //   1. reach convergence (displacement is below epsilon)
    //   2. reach step max - 1 (save one final step)
    //   3. fitting is not stable (see return statement below)
    m_converged = false;
    for(m_step = 0; m_step < m_stepMax-1 and not m_converged; ++m_step)
    {
        fit.init(q);

        auto status = Ponca::FIT_RESULT::NEED_OTHER_PASS;
        while(status == Ponca::FIT_RESULT::NEED_OTHER_PASS)
        {
            for(const auto& n : nei(q, m_scale))
            {
                fit.addNeighbor(n);
            }
            status = fit.finalize();
        }

        if(status != Ponca::FIT_RESULT::STABLE)
        {
            return q;
        }

        const VectorType proj = fit.project(p); // project the initial point
        const Scalar dist2 = (proj - q).squaredNorm();
        m_converged = dist2 < epsilon2;
        q = proj;
    }

    // final step
    {
        ++m_step;
        m_fitFinal.init(q);

        auto status = Ponca::FIT_RESULT::NEED_OTHER_PASS;
        while(status == Ponca::FIT_RESULT::NEED_OTHER_PASS)
        {
            for(const auto& n : nei(q, m_scale))
            {
                m_fitFinal.addNeighbor(n);
            }
            status = m_fitFinal.finalize();
        }
        return status == Ponca::FIT_RESULT::STABLE ? m_fitFinal.project(p) : q;
    }
}


template<class FitT, class FitFinalT>
template<class NeighborQueryT>
typename QuasiOrthogonalMLSProjection<FitT,FitFinalT>::VectorType
QuasiOrthogonalMLSProjection<FitT,FitFinalT>::project(const VectorType& p, const VectorType& n, const NeighborQueryT& nei)
{
    auto q = p; // moving point
    auto m = n;

    FitT fit;
    // m_fitFinal = FitFinalT();

    fit.setWeightFunc(WFunctor(m_scale));
    m_fitFinal.setWeightFunc(WFunctorFinal(m_scale));

    const auto epsilon = m_convRatio * m_scale;
    const auto epsilon2 = epsilon * epsilon;

    // Three stopping conditions:
    //   1. reach convergence (displacement is below epsilon)
    //   2. reach step max - 1 (save one final step)
    //   3. fitting is not stable (see return statement below)
    m_converged = false;
    for(m_step = 0; m_step < m_stepMax-1 and not m_converged; ++m_step)
    {
        fit.init(q, m);

        auto status = Ponca::FIT_RESULT::NEED_OTHER_PASS;
        while(status == Ponca::FIT_RESULT::NEED_OTHER_PASS)
        {
            for(const auto& n : nei(q, m_scale))
            {
                fit.addNeighbor(n);
            }
            status = fit.finalize();
        }

        if(status != Ponca::FIT_RESULT::STABLE)
        {
            return q;
        }

        const VectorType proj = fit.project(p); // project the initial point
        const Scalar dist2 = (proj - q).squaredNorm();
        m_converged = dist2 < epsilon2;
        q = proj;
        m = fit.primitiveGradient(q).normalized();
    }

    // final step
    {
        ++m_step;
        m_fitFinal.init(q, m);

        auto status = Ponca::FIT_RESULT::NEED_OTHER_PASS;
        while(status == Ponca::FIT_RESULT::NEED_OTHER_PASS)
        {
            for(const auto& n : nei(q, m_scale))
            {
                m_fitFinal.addNeighbor(n);
            }
            status = m_fitFinal.finalize();
        }
        return status == Ponca::FIT_RESULT::STABLE ? m_fitFinal.project(p) : q;
    }
}

} // namespace Ponca