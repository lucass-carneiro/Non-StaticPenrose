#include "KerrSchild_Kerr.hpp"
#include "aux_functions.hpp"

using grlensing::KerrSchild_Kerr;
using grlensing::metric_server;
using ksk_aux::d_r_KS_dx;
using ksk_aux::d_r_KS_dy;
using ksk_aux::d_r_KS_dz;
using ksk_aux::Power;
using ksk_aux::r_KS;
using ksk_aux::Sqrt;

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::grad_lapse(double, double x, double y,
                                                                      double z)
    -> metric_server::spatial_vector {

  auto r = r_KS(a, x, y, z);
  auto dr_dx = d_r_KS_dx(a, x, y, z);
  auto dr_dy = d_r_KS_dy(a, x, y, z);
  auto dr_dz = d_r_KS_dz(a, x, y, z);

  metric_server::spatial_vector grad_lapse{};

  grad_lapse[0]
      = (M * Power(r, 2)
         * (4 * M * x * Power(r, 4) * (Power(a, 2) + Power(r, 2))
            - (3 * Power(a, 6) * Power(z, 2)
               + r
                     * (4 * Power(a, 4) * M * Power(z, 2)
                        + r
                              * (6 * Power(a, 4) * Power(z, 2)
                                 + r
                                       * (8 * Power(a, 2) * M * Power(z, 2)
                                          + r
                                                * (-Power(a, 4) + 3 * Power(a, 2) * Power(z, 2)
                                                   + 4 * M
                                                         * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                         * r
                                                   - 2 * Power(a, 2) * Power(r, 2)
                                                   - Power(r, 4))))))
                  * dr_dx))
        / ((Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
            + Power(a, 2) * Power(z, 2) * Power(r, 2)
            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
            + Power(a, 2) * Power(r, 4) + Power(r, 6))
           * Sqrt(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                  + Power(a, 2) * Power(z, 2) * Power(r, 2)
                  + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                  + Power(a, 2) * Power(r, 4) + Power(r, 6))
           * Sqrt(Power(a, 4) * Power(z, 2)
                  + r
                        * (2 * Power(a, 2) * M * Power(z, 2)
                           + r
                                 * (Power(a, 2) * Power(z, 2)
                                    + r
                                          * (2 * M
                                                 * (-Power(a, 2) + Power(x, 2) + Power(y, 2)
                                                    + Power(z, 2))
                                             + r * (Power(a, 2) - 2 * M * r + Power(r, 2)))))));

  grad_lapse[1]
      = (M * Power(r, 2)
         * (4 * M * y * Power(r, 4) * (Power(a, 2) + Power(r, 2))
            - (3 * Power(a, 6) * Power(z, 2)
               + r
                     * (4 * Power(a, 4) * M * Power(z, 2)
                        + r
                              * (6 * Power(a, 4) * Power(z, 2)
                                 + r
                                       * (8 * Power(a, 2) * M * Power(z, 2)
                                          + r
                                                * (-Power(a, 4) + 3 * Power(a, 2) * Power(z, 2)
                                                   + 4 * M
                                                         * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                         * r
                                                   - 2 * Power(a, 2) * Power(r, 2)
                                                   - Power(r, 4))))))
                  * dr_dy))
        / ((Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
            + Power(a, 2) * Power(z, 2) * Power(r, 2)
            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
            + Power(a, 2) * Power(r, 4) + Power(r, 6))
           * Sqrt(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                  + Power(a, 2) * Power(z, 2) * Power(r, 2)
                  + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                  + Power(a, 2) * Power(r, 4) + Power(r, 6))
           * Sqrt(Power(a, 4) * Power(z, 2)
                  + r
                        * (2 * Power(a, 2) * M * Power(z, 2)
                           + r
                                 * (Power(a, 2) * Power(z, 2)
                                    + r
                                          * (2 * M
                                                 * (-Power(a, 2) + Power(x, 2) + Power(y, 2)
                                                    + Power(z, 2))
                                             + r * (Power(a, 2) - 2 * M * r + Power(r, 2)))))));

  grad_lapse[2]
      = (M * Power(r, 2)
         * (2 * z * r * (Power(a, 2) + 2 * M * r) * Power(Power(a, 2) + Power(r, 2), 2)
            - (3 * Power(a, 6) * Power(z, 2)
               + r
                     * (4 * Power(a, 4) * M * Power(z, 2)
                        + r
                              * (6 * Power(a, 4) * Power(z, 2)
                                 + r
                                       * (8 * Power(a, 2) * M * Power(z, 2)
                                          + r
                                                * (-Power(a, 4) + 3 * Power(a, 2) * Power(z, 2)
                                                   + 4 * M
                                                         * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                         * r
                                                   - 2 * Power(a, 2) * Power(r, 2)
                                                   - Power(r, 4))))))
                  * dr_dz))
        / ((Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
            + Power(a, 2) * Power(z, 2) * Power(r, 2)
            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
            + Power(a, 2) * Power(r, 4) + Power(r, 6))
           * Sqrt(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                  + Power(a, 2) * Power(z, 2) * Power(r, 2)
                  + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                  + Power(a, 2) * Power(r, 4) + Power(r, 6))
           * Sqrt(Power(a, 4) * Power(z, 2)
                  + r
                        * (2 * Power(a, 2) * M * Power(z, 2)
                           + r
                                 * (Power(a, 2) * Power(z, 2)
                                    + r
                                          * (2 * M
                                                 * (-Power(a, 2) + Power(x, 2) + Power(y, 2)
                                                    + Power(z, 2))
                                             + r * (Power(a, 2) - 2 * M * r + Power(r, 2)))))));

  return grad_lapse;
}