#include <Fastor/Fastor.h>
#include <cmath>
#include <iostream>

using namespace Fastor;

class SKS {
public:
    SKS(double sep, double mass1, double mass2, double spin1, double spin2) 
        : b(sep), 
          M{mass1, mass2},
          a{spin1,spin2},
          Omega(std::sqrt((M[0] + M[1])/(b * b * b))) {}
    
    enum class bh_index {bh_1, bh_2};

    Tensor<double, 4, 4> metric(double t, double x, double y, double z) {
        const Tensor<double, 4, 4> eta = {{-1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
        
        const auto H_hat_1 = H_hat<bh_index::bh_1>(t, x, y, z);
        const auto H_hat_2 = H_hat<bh_index::bh_2>(t, x, y, z);

        const auto l_hat_1 = l_hat<bh_index::bh_1>(t, x, y, z);
        const auto l_hat_2 = l_hat<bh_index::bh_2>(t, x, y, z);

        return eta + 2 * H_hat_1 * l_hat_1 + 2 * H_hat_2 * l_hat_2;
    }

private:
    double b;
    Tensor<double,2> M;
    Tensor<double,2> a;
    double Omega;

    template<bh_index i>
    double x_K(double t) {
        if constexpr (i == bh_index::bh_1)
            return b/2 * std::cos(Omega * t);
        else
            return -b/2 * std::cos(Omega * t);
    }

    template<bh_index i>
    double y_K(double t) {
        if constexpr (i == bh_index::bh_1)
            return b/2 * std::sin(Omega * t);
        else
            return -b/2 * std::sin(Omega * t);
    }

    template<bh_index n>
    double rKS(double x, double y, double z) {
        constexpr const auto i = static_cast<int>(n);
        
        const double rho = std::sqrt(x*x + y*y + z*z);
        
        const auto p1 = rho * rho - a[i]*a[i];
        const auto p2 = 0.5 * p1 + std::sqrt(0.25 * p1 * p1 + a[i]*a[i] * z * z);

        return p2;
    }

    template<bh_index n>
    double H(double x, double y, double z) {
        constexpr const auto i = static_cast<int>(n);

        const auto rks = rKS<n>(x, y, z);
        const auto cos_theta = z/rks;

        return (M[i] * rks) / (rks * rks + a[i] * a[i] * cos_theta * cos_theta);
    }

    template<bh_index n>
    Tensor<double, 4> l(double x, double y, double z) {
        constexpr const auto i = static_cast<int>(n);

        const auto rks = rKS<n>(x, y, z);
        const auto den = rks * rks + a[i] * a[i];

        return Tensor<double,4> {
            1.0,
            (rks * x + a[i] * y)/den,
            (rks * y - a[i] * x)/den,
            z/rks
        };
    }

    template<bh_index n>
    Tensor<double, 3> v_i(double t) {
        if constexpr (n == bh_index::bh_1) {
            return Tensor<double, 3> {
                -b/2 * std::sin(Omega * t) * Omega,
                b/2 * std::cos(Omega * t) * Omega,
                0.0
            };
        }
        else {
            return Tensor<double, 3> {
                b/2 * std::sin(Omega * t) * Omega,
                -b/2 * std::cos(Omega * t) * Omega,
                0.0
            };
        }
    }

    template<bh_index n>
    Tensor<double, 3> v_hat_i(double t) {
        const auto vi = v_i<n>(t);
        const auto norm_vi = norm(vi);
        return vi/norm_vi;
    }

    template<bh_index n>
    double gamma(double t) {
        const auto v = v_i<n>(t);
        return std::sqrt(1 - (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
    }

    template<bh_index n>
    Tensor<double, 4> Lambda_circ(double th, double xh, double yh, double zh) {
        const auto g = gamma<n>(th);
        const auto v = v_hat_i<n>(th);

        return Tensor<double, 4> {
            g * th - g * v[0] * xh - g * v[1] * yh,
            -x_K<n>(th) + (1 + (g - 1)*v[0]*v[0])*xh + ((g - 1) * v[0]*v[1])*yh,
            -y_K<n>(th) + ((g - 1)*v[1]*v[0])*xh + (1 + (g - 1)*v[1]*v[1])*yh,
            zh,
        };
    }

    template<bh_index n>
    Tensor<double, 4,4> Lambda(double t) {
        const auto g = gamma<n>(t);
        const auto v = v_i<n>(t);
        const auto vh = v_hat_i<n>(t);

        Tensor<double, 4, 4> result{
            {g      , -g*v[0]            , -g*v[1]            , -g*v[2]            },
            {-g*v[0], 1+(g-1)*vh[0]*vh[0], (g-1)*vh[0]*vh[1]  , (g-1)*vh[0]*vh[2]  },
            {-g*v[1], (g-1)*vh[1]*vh[0]  , 1+(g-1)*vh[1]*vh[1], (g-1)*vh[1]*vh[2]  },
            {-g*v[2], (g-1)*vh[2]*vh[0]  , (g-1)*vh[2]*vh[1]  , 1+(g-1)*vh[2]*vh[2]},
        };

        return result;
    }

    template<bh_index n>
    Tensor<double, 3> l_hat(double th, double xh, double yh, double zh) {
        enum {mu, nu};
        const auto L = Lambda<n>(th);
        const auto L_circ = Lambda_circ<n>(th, xh, yh, zh);
        const auto old_l = l<n>(L_circ[1], L_circ[2], zh);
        return einsum<Index<nu,mu>, Index<nu>>(L, old_l);
    }

    template<bh_index n>
    double H_hat(double th, double xh, double yh, double zh) {
        const auto L_circ = Lambda_circ<n>(th, xh, yh, zh);
        return H<n>(L_circ[1], L_circ[2], zh);
    }
};

int main() {
    SKS g(1, 1, 1, 0.5, 0.5);
    g.metric(0, 0, 0, 0);
    return 0;
}
