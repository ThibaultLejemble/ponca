#pragma once

#include <Ponca/Fitting>

namespace Ponca {

/*!
    \brief Quasi-Orthogonal Moving Least Squares projection
    
    See Section 4.2 of
    On Normals and Projection Operators for Surfaces Defined by Point Sets
    Marc Alexa and Anders Adamson
    Eurographics Symposium on Point-Based Graphics (2004)

    \note Also called 'almost' orthogonal projection
    \tparam FitT is the main Concept::FittingProcedureConcept used during the first N-1 steps
    \tparam FitFinalT is another Concept::FittingProcedureConcept used for the N-th fitting step (default is FitT)
 */
template<class FitT, class FitFinalT = FitT>
class QuasiOrthogonalMLSProjection
{
public:
    using Scalar        = typename FitT::Scalar;
    using VectorType    = typename FitT::VectorType;
    using WFunctor      = typename FitT::WFunctor;
    using WFunctorFinal = typename FitFinalT::WFunctor;

public:
    /*!
        \brief Perform the fitting to nei and return the projection of p
        \tparam NeighborQueryT must provide an operator()(const VectorType& q, Scalar t)
        that returns a range over the neighbors of q at a distance t (of type FitT::DataPoint)
    */
    template<class NeighborQueryT>
    VectorType project(const VectorType& p, const NeighborQueryT& nei);

    template<class NeighborQueryT>
    VectorType project(const VectorType& p, const VectorType& n, const NeighborQueryT& nei);

public:
    void setScale(Scalar scale) {m_scale = scale;}
    void setStepMax(int step) {m_stepMax = step;}
    void setConvRatio(Scalar convRatio) {m_convRatio = convRatio;}

protected:
    int    m_stepMax{20};       /*!< \brief Maximal number of iterations */
    Scalar m_convRatio{0.001};  /*!< \brief Distance threshold for convergence (ratio of scale) */
    Scalar m_scale{0.0};        /*!< \brief Weighting kernel support size */

public:
    FitFinalT m_fitFinal;       /*!< \brief The final fit can be used after a call to project */
    int m_step{0};
    bool m_converged{false};
};

} // namespace Ponca

#include "mlsProjection.hpp"