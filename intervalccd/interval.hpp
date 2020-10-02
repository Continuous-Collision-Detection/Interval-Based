// An interval object.
#pragma once

#include <string>

#include <Eigen/Core>
#include <boost/numeric/interval.hpp>

namespace intervalccd {

namespace interval_options {
    typedef boost::numeric::interval_lib::checking_base<double> CheckingPolicy;
} // namespace interval_options

typedef boost::numeric::interval<
    double,
    boost::numeric::interval_lib::policies<
        boost::numeric::interval_lib::save_state<
            boost::numeric::interval_lib::rounded_arith_std<double>>,
        interval_options::CheckingPolicy>>
    Interval;

} // namespace intervalccd

namespace Eigen {
typedef Matrix<intervalccd::Interval, 2, 1> Vector2I;
typedef Matrix<intervalccd::Interval, 3, 1> Vector3I;
typedef Matrix<intervalccd::Interval, Dynamic, 1> VectorXI;
typedef Matrix<intervalccd::Interval, Dynamic, 1, ColMajor, 3, 1> VectorX3I;
typedef Matrix<double, 3, 1> Vector3d;
typedef Matrix<double, Dynamic, 1, ColMajor, 3, 1> VectorX3d;

template <typename T, int dim, int max_dim = dim>
using Vector = Matrix<T, dim, 1, ColMajor, max_dim, 1>;

} // namespace Eigen
