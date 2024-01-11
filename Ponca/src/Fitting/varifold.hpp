#include "varifold.h"

namespace Ponca
{

namespace v0 
{

template<class DataPoint, class _WFunctor, typename T>
void Varifold<DataPoint, _WFunctor, T>::init(
    const VectorType& x_l0,
    const VectorType& n_l0)
{
    Base::init(x_l0);

    m_n_l0 = n_l0;
    m_P_l0 = MatrixType::Identity() - n_l0 * n_l0.transpose();

    for(auto k = 0; k < 3; ++k)
        m_B_ijk[k] = MatrixType::Zero();
    m_deno = Scalar(0);

    m_k1 = Scalar(0);
    m_dir1 = VectorType::Zero();
    m_k2 = Scalar(0);
    m_dir2 = VectorType::Zero();
}

template<class DataPoint, class _WFunctor, typename T>
bool Varifold<DataPoint, _WFunctor, T>::addLocalNeighbor(
    Scalar w, const VectorType &localQ, const DataPoint &attributes)
{
    if(not Base::addLocalNeighbor(w, localQ, attributes))
        return false;

    if(localQ.squaredNorm() < Eigen::NumTraits<Scalar>::dummy_precision())
        return false;
    
    const Scalar t = Base::m_w.evalScale();

    // x_l0 = 0 in centered basis 
    const VectorType& x_l = localQ;
    const VectorType x_l0 = VectorType::Zero();
    const VectorType n_l = attributes.normal();
    const MatrixType P_l = MatrixType::Identity() - n_l * n_l.transpose();
    constexpr Scalar m_l = 1.0;
    constexpr Scalar d = 2.0;
    constexpr Scalar n = 3.0;

    const Scalar rho_prime = rho_der((x_l - x_l0).norm() / t);
    const Scalar xi = - t * rho_prime / 3;

    const Scalar coeff = d/n 
        * m_l
        * rho_prime
        / (x_l0 - x_l).norm()
        * 1.0/2.0;

    m_deno += m_l * xi;

    const VectorType u = P_l * (x_l0 - x_l);
    const MatrixType Pdiff = P_l - m_P_l0;

    for(auto i = 0; i < 3; ++i)
    {
        VectorType e_i = VectorType::Zero();
        e_i[i] = 1;
        for(auto j = 0; j < 3; ++j)
        {
            VectorType e_j = VectorType::Zero();
            e_j[j] = 1;
            for(auto k = 0; k < 3; ++k)
            {
                VectorType e_k = VectorType::Zero();
                e_k[k] = 1;
                {
                    const VectorType v1 = Pdiff(j,k) * e_i;
                    const VectorType v2 = Pdiff(i,k) * e_j;
                    const VectorType v3 = Pdiff(i,j) * e_k;
                    const VectorType v = v1 + v2 - v3;

                    m_B_ijk[k](i,j) += coeff * u.dot(v);
                }
            }
        }
    }
    return true;
}

template<class DataPoint, class _WFunctor, typename T>
FIT_RESULT Varifold<DataPoint, _WFunctor, T>::finalize()
{
    if(Base::finalize() != STABLE)
        return Base::m_eCurrentState; 

    constexpr Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();

    if(std::abs(m_deno) < epsilon)
        return UNSTABLE;

    for(auto k = 0; k < 3; ++k) {
        m_B_ijk[k] /= m_deno;
    }

    MatrixType B_ij = MatrixType::Zero();
    for(auto i = 0; i < 3; ++i)
    for(auto j = 0; j < 3; ++j)
    for(auto k = 0; k < 3; ++k)
    {
        B_ij(i,j) += m_B_ijk[k](i,j) * m_n_l0[k];
    }

    const Mat32 Q = tangentPlane();
    const Mat22 B = - Q.transpose() * B_ij * Q; // minus for sign convention

    Eigen::SelfAdjointEigenSolver<Mat22> solver;
    solver.computeDirect(B);

    if(solver.info() != Eigen::Success)
        return UNSTABLE;

    m_k1 = solver.eigenvalues()[0];
    m_dir1 = Q * solver.eigenvectors().col(0);
    m_k2 = solver.eigenvalues()[1];
    m_dir2 = Q * solver.eigenvectors().col(1);

    // if(std::abs(m_k1) < std::abs(m_k2))
    // {
    //     std::swap(m_k1, m_k1);
    //     std::swap(m_dir1, m_dir2);
    // }

    return STABLE;
}

template<class DataPoint, class _WFunctor, typename T>
typename Varifold<DataPoint, _WFunctor, T>::Mat32 Varifold<DataPoint, _WFunctor, T>::tangentPlane() const
{
    Mat32 B;
    int i0=-1, i1=-1, i2=-1;
    const VectorType& n = m_n_l0;
    n.array().abs().minCoeff(&i0); // i0: dimension where n extends the least
    i1 = (i0+1)%3;
    i2 = (i0+2)%3;

    B.col(0)[i0] = 0;
    B.col(0)[i1] = n[i2];
    B.col(0)[i2] = -n[i1];

    B.col(0).normalize();
    B.col(1) = B.col(0).cross(n);
    return B;
}

template<class DataPoint, class _WFunctor, typename T>
typename Varifold<DataPoint, _WFunctor, T>::Scalar Varifold<DataPoint, _WFunctor, T>::rho_der(Scalar _x)
{
    const Scalar g = Scalar(1) / (Scalar(1) - _x*_x);
    return - Scalar(2) * _x  * g * g * std::exp(-g);
}

template<class DataPoint, class _WFunctor, typename T>
typename Varifold<DataPoint, _WFunctor, T>::VectorType Varifold<DataPoint, _WFunctor, T>::project(const VectorType& p) const
{
    return Base::m_w.basisCenter() + (p - (p - Base::m_w.basisCenter()).dot(m_n_l0) * m_n_l0);
}

template<class DataPoint, class _WFunctor, typename T>
typename Varifold<DataPoint, _WFunctor, T>::VectorType Varifold<DataPoint, _WFunctor, T>::primitiveGradient(const VectorType&) const
{
    return m_n_l0;
}

} //namespace v0

// ============================================================================

namespace v1
{

template<class DataPoint, class _WFunctor, typename T>
void Varifold<DataPoint, _WFunctor, T>::init(
    const VectorType& x_l0,
    const VectorType& n_l0)
{
    Base::init(x_l0);

    m_n_l0 = n_l0;
    m_P_l0 = MatrixType::Identity() - n_l0 * n_l0.transpose();

    m_nume = MatrixType::Zero();
    m_deno = Scalar(0);

    m_k1 = Scalar(0);
    m_dir1 = VectorType::Zero();
    m_k2 = Scalar(0);
    m_dir2 = VectorType::Zero();
}

template<class DataPoint, class _WFunctor, typename T>
bool Varifold<DataPoint, _WFunctor, T>::addLocalNeighbor(
    Scalar w, const VectorType &localQ, const DataPoint &attributes)
{
    if(not Base::addLocalNeighbor(w, localQ, attributes))
        return false;

    if(localQ.squaredNorm() < Eigen::NumTraits<Scalar>::dummy_precision())
        return false;

    const VectorType n_l = attributes.normal();
    const MatrixType P_l = MatrixType::Identity() - n_l * n_l.transpose();
    const Scalar d_l = localQ.norm();
    const VectorType delta_x = localQ / d_l; // = localQ.normalized()
    const VectorType u = P_l * delta_x;
    const MatrixType DeltaP = P_l - m_P_l0;
    const VectorType Pl_nl0 = P_l * m_n_l0;

    const MatrixType a1 = u * Pl_nl0.transpose();
    const MatrixType a2 = u.dot(m_n_l0) * DeltaP;

    w *= -1; // see WARNING in VarifoldWeightKernel::f() (nor necessary because of the ratio nume/deno)

    m_nume += w * (a1 + a1.transpose() - a2);
    m_deno += w * d_l;
    
    return true;
}

template<class DataPoint, class _WFunctor, typename T>
FIT_RESULT Varifold<DataPoint, _WFunctor, T>::finalize()
{
    if(Base::finalize() != STABLE)
        return Base::m_eCurrentState;

    constexpr Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();

    if(std::abs(m_deno) < epsilon)
        return UNSTABLE;

    const MatrixType A = m_nume / m_deno * Base::m_w.evalScale(); // d/2 = 1

    const Mat32 Q = tangentPlane();
    const Mat22 B = - Q.transpose() * A * Q; // minus for sign convention

    Eigen::SelfAdjointEigenSolver<Mat22> solver;
    solver.computeDirect(B);

    if(solver.info() != Eigen::Success)
        return UNSTABLE;

    m_k1 = solver.eigenvalues()[0];
    m_dir1 = Q * solver.eigenvectors().col(0);
    m_k2 = solver.eigenvalues()[1];
    m_dir2 = Q * solver.eigenvectors().col(1);

    // if(std::abs(m_k1) < std::abs(m_k2))
    // {
    //     std::swap(m_k1, m_k1);
    //     std::swap(m_dir1, m_dir2);
    // }

    return STABLE;
}

template<class DataPoint, class _WFunctor, typename T>
typename Varifold<DataPoint, _WFunctor, T>::Mat32 Varifold<DataPoint, _WFunctor, T>::tangentPlane() const
{
    Mat32 B;
    int i0=-1, i1=-1, i2=-1;
    const VectorType& n = m_n_l0;
    n.array().abs().minCoeff(&i0); // i0: dimension where n extends the least
    i1 = (i0+1)%3;
    i2 = (i0+2)%3;

    B.col(0)[i0] = 0;
    B.col(0)[i1] = n[i2];
    B.col(0)[i2] = -n[i1];

    B.col(0).normalize();
    B.col(1) = B.col(0).cross(n);
    return B;
}

template<class DataPoint, class _WFunctor, typename T>
typename Varifold<DataPoint, _WFunctor, T>::VectorType Varifold<DataPoint, _WFunctor, T>::project(const VectorType& p) const
{
    return Base::m_w.basisCenter() + (p - (p - Base::m_w.basisCenter()).dot(m_n_l0) * m_n_l0);
}

template<class DataPoint, class _WFunctor, typename T>
typename Varifold<DataPoint, _WFunctor, T>::VectorType Varifold<DataPoint, _WFunctor, T>::primitiveGradient(const VectorType&) const
{
    return m_n_l0;
}

} //namespace v1
} //namespace Ponca