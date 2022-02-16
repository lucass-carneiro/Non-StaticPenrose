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

  auto r = [&](double x, double y, double z) { return r_KS(a, x, y, z); };
  auto dr_dx = [&](double x, double y, double z) { return d_r_KS_dx(a, x, y, z); };
  auto dr_dy = [&](double x, double y, double z) { return d_r_KS_dy(a, x, y, z); };
  auto dr_dz = [&](double x, double y, double z) { return d_r_KS_dz(a, x, y, z); };

  return metric_server::spatial_vector{
      (M * Power(r(x, y, z), 2)
       * (4 * M * x * Power(r(x, y, z), 4) * (Power(a, 2) + Power(r(x, y, z), 2))
          - (3 * Power(a, 6) * Power(z, 2)
             + r(x, y, z)
                   * (4 * Power(a, 4) * M * Power(z, 2)
                      + r(x, y, z)
                            * (6 * Power(a, 4) * Power(z, 2)
                               + r(x, y, z)
                                     * (8 * Power(a, 2) * M * Power(z, 2)
                                        + r(x, y, z)
                                              * (-Power(a, 4) + 3 * Power(a, 2) * Power(z, 2)
                                                 + 4 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                       * r(x, y, z)
                                                 - 2 * Power(a, 2) * Power(r(x, y, z), 2)
                                                 - Power(r(x, y, z), 4))))))
                * dr_dx(x, y, z)))
          / (Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                       + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                       + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6),
                   2)
             * Sqrt(
                 1
                 - (2 * M * Power(r(x, y, z), 3) * (Power(a, 2) + Power(r(x, y, z), 2)))
                       / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                          + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                          + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6)))),
      (M * Power(r(x, y, z), 2)
       * (4 * M * y * Power(r(x, y, z), 4) * (Power(a, 2) + Power(r(x, y, z), 2))
          - (3 * Power(a, 6) * Power(z, 2)
             + r(x, y, z)
                   * (4 * Power(a, 4) * M * Power(z, 2)
                      + r(x, y, z)
                            * (6 * Power(a, 4) * Power(z, 2)
                               + r(x, y, z)
                                     * (8 * Power(a, 2) * M * Power(z, 2)
                                        + r(x, y, z)
                                              * (-Power(a, 4) + 3 * Power(a, 2) * Power(z, 2)
                                                 + 4 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                       * r(x, y, z)
                                                 - 2 * Power(a, 2) * Power(r(x, y, z), 2)
                                                 - Power(r(x, y, z), 4))))))
                * dr_dy(x, y, z)))
          / (Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                       + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                       + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6),
                   2)
             * Sqrt(
                 1
                 - (2 * M * Power(r(x, y, z), 3) * (Power(a, 2) + Power(r(x, y, z), 2)))
                       / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                          + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                          + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6)))),
      (M * Power(r(x, y, z), 2)
       * (2 * z * r(x, y, z) * (Power(a, 2) + 2 * M * r(x, y, z))
              * Power(Power(a, 2) + Power(r(x, y, z), 2), 2)
          - (3 * Power(a, 6) * Power(z, 2)
             + r(x, y, z)
                   * (4 * Power(a, 4) * M * Power(z, 2)
                      + r(x, y, z)
                            * (6 * Power(a, 4) * Power(z, 2)
                               + r(x, y, z)
                                     * (8 * Power(a, 2) * M * Power(z, 2)
                                        + r(x, y, z)
                                              * (-Power(a, 4) + 3 * Power(a, 2) * Power(z, 2)
                                                 + 4 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                       * r(x, y, z)
                                                 - 2 * Power(a, 2) * Power(r(x, y, z), 2)
                                                 - Power(r(x, y, z), 4))))))
                * dr_dz(x, y, z)))
          / (Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                       + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                       + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6),
                   2)
             * Sqrt(
                 1
                 - (2 * M * Power(r(x, y, z), 3) * (Power(a, 2) + Power(r(x, y, z), 2)))
                       / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                          + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                          + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))))};
}