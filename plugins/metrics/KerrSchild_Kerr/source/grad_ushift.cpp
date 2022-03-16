#include "KerrSchild_Kerr.hpp"
#include "aux_functions.hpp"

using grlensing::KerrSchild_Kerr;
using grlensing::metric_server;
using ksk_aux::d_r_KS_dx;
using ksk_aux::d_r_KS_dy;
using ksk_aux::d_r_KS_dz;
using ksk_aux::Power;
using ksk_aux::r_KS;

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::grad_ushift(double, double x, double y,
                                                                       double z)
    -> metric_server::spatial_matrix {

  auto r = [&](double x, double y, double z) { return r_KS(a, x, y, z); };
  auto dr_dx = [&](double x, double y, double z) { return d_r_KS_dx(a, x, y, z); };
  auto dr_dy = [&](double x, double y, double z) { return d_r_KS_dy(a, x, y, z); };
  auto dr_dz = [&](double x, double y, double z) { return d_r_KS_dz(a, x, y, z); };

  metric_server::spatial_matrix matrix{};

  matrix[0][0]
      = (2 * M * Power(r(x, y, z), 2)
         * (Power(r(x, y, z), 2)
                * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                   + a * (-4 * M * x * y + a * Power(z, 2)) * Power(r(x, y, z), 2)
                   + 2 * M * (-Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                   + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))
            + (3 * Power(a, 5) * y * Power(z, 2)
               + 4 * Power(a, 3) * (a * x + M * y) * Power(z, 2) * r(x, y, z)
               + Power(a, 2) * (6 * M * x + a * y) * Power(z, 2) * Power(r(x, y, z), 2)
               + 2 * Power(a, 2) * x * Power(z, 2) * Power(r(x, y, z), 3)
               + (-(Power(a, 3) * y) + 2 * M * x * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                     * Power(r(x, y, z), 4)
               - 3 * a * y * Power(r(x, y, z), 6) - 2 * x * Power(r(x, y, z), 7))
                  * dr_dx(x, y, z)))
        / Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                    + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                    + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                    + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6),
                2);

  matrix[0][1]
      = (2 * M * Power(r(x, y, z), 2)
         * (r(x, y, z)
                * (Power(a, 5) * Power(z, 2) + 2 * Power(a, 3) * M * Power(z, 2) * r(x, y, z)
                   + Power(a, 3) * Power(z, 2) * Power(r(x, y, z), 2)
                   + 2 * a * M * (Power(x, 2) - Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                   + (Power(a, 3) - 4 * M * x * y) * Power(r(x, y, z), 4)
                   + a * Power(r(x, y, z), 6))
            + (3 * Power(a, 5) * y * Power(z, 2)
               + 4 * Power(a, 3) * (a * x + M * y) * Power(z, 2) * r(x, y, z)
               + Power(a, 2) * (6 * M * x + a * y) * Power(z, 2) * Power(r(x, y, z), 2)
               + 2 * Power(a, 2) * x * Power(z, 2) * Power(r(x, y, z), 3)
               + (-(Power(a, 3) * y) + 2 * M * x * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                     * Power(r(x, y, z), 4)
               - 3 * a * y * Power(r(x, y, z), 6) - 2 * x * Power(r(x, y, z), 7))
                  * dr_dy(x, y, z)))
        / Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                    + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                    + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                    + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6),
                2);

  matrix[0][2]
      = (2 * M * Power(r(x, y, z), 2)
         * (-2 * z * r(x, y, z) * (Power(a, 2) + 2 * M * r(x, y, z)) * (a * y + x * r(x, y, z))
                * (Power(a, 2) + Power(r(x, y, z), 2))
            + (3 * Power(a, 5) * y * Power(z, 2)
               + 4 * Power(a, 3) * (a * x + M * y) * Power(z, 2) * r(x, y, z)
               + Power(a, 2) * (6 * M * x + a * y) * Power(z, 2) * Power(r(x, y, z), 2)
               + 2 * Power(a, 2) * x * Power(z, 2) * Power(r(x, y, z), 3)
               + (-(Power(a, 3) * y) + 2 * M * x * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                     * Power(r(x, y, z), 4)
               - 3 * a * y * Power(r(x, y, z), 6) - 2 * x * Power(r(x, y, z), 7))
                  * dr_dz(x, y, z)))
        / Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                    + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                    + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                    + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6),
                2);

  matrix[1][0]
      = (2 * M * Power(r(x, y, z), 2)
         * (r(x, y, z)
                * (-(Power(a, 5) * Power(z, 2)) - 2 * Power(a, 3) * M * Power(z, 2) * r(x, y, z)
                   - Power(a, 3) * Power(z, 2) * Power(r(x, y, z), 2)
                   - 2 * a * M * (-Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                   - (Power(a, 3) + 4 * M * x * y) * Power(r(x, y, z), 4)
                   - a * Power(r(x, y, z), 6))
            + (-3 * Power(a, 5) * x * Power(z, 2)
               + 4 * Power(a, 3) * (-(M * x) + a * y) * Power(z, 2) * r(x, y, z)
               + Power(a, 2) * (-(a * x) + 6 * M * y) * Power(z, 2) * Power(r(x, y, z), 2)
               + 2 * Power(a, 2) * y * Power(z, 2) * Power(r(x, y, z), 3)
               + (Power(a, 3) * x + 2 * M * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                     * Power(r(x, y, z), 4)
               + 3 * a * x * Power(r(x, y, z), 6) - 2 * y * Power(r(x, y, z), 7))
                  * dr_dx(x, y, z)))
        / Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                    + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                    + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                    + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6),
                2);

  matrix[1][0]
      = (2 * M * Power(r(x, y, z), 2)
         * (Power(r(x, y, z), 2)
                * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                   + a * (4 * M * x * y + a * Power(z, 2)) * Power(r(x, y, z), 2)
                   + 2 * M * (Power(x, 2) - Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                   + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))
            + (-3 * Power(a, 5) * x * Power(z, 2)
               + 4 * Power(a, 3) * (-(M * x) + a * y) * Power(z, 2) * r(x, y, z)
               + Power(a, 2) * (-(a * x) + 6 * M * y) * Power(z, 2) * Power(r(x, y, z), 2)
               + 2 * Power(a, 2) * y * Power(z, 2) * Power(r(x, y, z), 3)
               + (Power(a, 3) * x + 2 * M * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                     * Power(r(x, y, z), 4)
               + 3 * a * x * Power(r(x, y, z), 6) - 2 * y * Power(r(x, y, z), 7))
                  * dr_dy(x, y, z)))
        / Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                    + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                    + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                    + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6),
                2);

  matrix[1][2]
      = (2 * M * Power(r(x, y, z), 2)
         * (2 * z * r(x, y, z) * (Power(a, 2) + 2 * M * r(x, y, z)) * (a * x - y * r(x, y, z))
                * (Power(a, 2) + Power(r(x, y, z), 2))
            + (-3 * Power(a, 5) * x * Power(z, 2)
               + 4 * Power(a, 3) * (-(M * x) + a * y) * Power(z, 2) * r(x, y, z)
               + Power(a, 2) * (-(a * x) + 6 * M * y) * Power(z, 2) * Power(r(x, y, z), 2)
               + 2 * Power(a, 2) * y * Power(z, 2) * Power(r(x, y, z), 3)
               + (Power(a, 3) * x + 2 * M * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                     * Power(r(x, y, z), 4)
               + 3 * a * x * Power(r(x, y, z), 6) - 2 * y * Power(r(x, y, z), 7))
                  * dr_dz(x, y, z)))
        / Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                    + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                    + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                    + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6),
                2);

  matrix[2][0]
      = (4 * M * z * r(x, y, z)
         * (-2 * M * x * Power(r(x, y, z), 4) * (Power(a, 2) + Power(r(x, y, z), 2))
            + (Power(a, 6) * Power(z, 2)
               + r(x, y, z)
                     * (Power(a, 4) * M * Power(z, 2)
                        + r(x, y, z)
                              * (2 * Power(a, 4) * Power(z, 2)
                                 + r(x, y, z)
                                       * (-(Power(a, 2) * M
                                            * (Power(x, 2) + Power(y, 2) - 2 * Power(z, 2)))
                                          + r(x, y, z)
                                                * (-Power(a, 4) + Power(a, 2) * Power(z, 2)
                                                   + M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                         * r(x, y, z)
                                                   - 2 * Power(a, 2) * Power(r(x, y, z), 2)
                                                   - Power(r(x, y, z), 4))))))
                  * dr_dx(x, y, z)))
        / Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                    + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                    + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                    + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6),
                2);
  matrix[2][1]
      = (4 * M * z * r(x, y, z)
         * (-2 * M * y * Power(r(x, y, z), 4) * (Power(a, 2) + Power(r(x, y, z), 2))
            + (Power(a, 6) * Power(z, 2)
               + r(x, y, z)
                     * (Power(a, 4) * M * Power(z, 2)
                        + r(x, y, z)
                              * (2 * Power(a, 4) * Power(z, 2)
                                 + r(x, y, z)
                                       * (-(Power(a, 2) * M
                                            * (Power(x, 2) + Power(y, 2) - 2 * Power(z, 2)))
                                          + r(x, y, z)
                                                * (-Power(a, 4) + Power(a, 2) * Power(z, 2)
                                                   + M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                         * r(x, y, z)
                                                   - 2 * Power(a, 2) * Power(r(x, y, z), 2)
                                                   - Power(r(x, y, z), 4))))))
                  * dr_dy(x, y, z)))
        / Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                    + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                    + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                    + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6),
                2);
  matrix[2][2]
      = (2 * M * r(x, y, z)
         * (r(x, y, z) * (Power(a, 2) + Power(r(x, y, z), 2))
                * (-(Power(a, 4) * Power(z, 2)) - 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                   - Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                   + 2 * M * (Power(x, 2) + Power(y, 2) - Power(z, 2)) * Power(r(x, y, z), 3)
                   + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))
            + 2 * z
                  * (Power(a, 6) * Power(z, 2)
                     + r(x, y, z)
                           * (Power(a, 4) * M * Power(z, 2)
                              + r(x, y, z)
                                    * (2 * Power(a, 4) * Power(z, 2)
                                       + r(x, y, z)
                                             * (-(Power(a, 2) * M
                                                  * (Power(x, 2) + Power(y, 2) - 2 * Power(z, 2)))
                                                + r(x, y, z)
                                                      * (-Power(a, 4) + Power(a, 2) * Power(z, 2)
                                                         + M
                                                               * (Power(x, 2) + Power(y, 2)
                                                                  + Power(z, 2))
                                                               * r(x, y, z)
                                                         - 2 * Power(a, 2) * Power(r(x, y, z), 2)
                                                         - Power(r(x, y, z), 4))))))
                  * dr_dz(x, y, z)))
        / Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                    + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                    + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                    + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6),
                2);

  return matrix;
}