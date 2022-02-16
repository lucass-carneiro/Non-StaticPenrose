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

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::ul_extrinsic(double, double x, double y,
                                                                        double z)
    -> metric_server::spatial_matrix {
  auto r = [&](double x, double y, double z) { return r_KS(a, x, y, z); };
  auto dr_dx = [&](double x, double y, double z) { return d_r_KS_dx(a, x, y, z); };
  auto dr_dy = [&](double x, double y, double z) { return d_r_KS_dy(a, x, y, z); };
  auto dr_dz = [&](double x, double y, double z) { return d_r_KS_dz(a, x, y, z); };

  metric_server::spatial_matrix mixedK;
  mixedK[0][0]
      = (2 * M * Power(r(x, y, z), 2)
         * (Power(r(x, y, z), 18) + 3 * Power(a, 13) * y * Power(z, 4) * dr_dx(x, y, z)
            + 2 * Power(a, 11) * (2 * a * x + 3 * M * y) * Power(z, 4) * r(x, y, z) * dr_dx(x, y, z)
            - 2 * x * Power(r(x, y, z), 17) * dr_dx(x, y, z)
            + 2 * Power(r(x, y, z), 15)
                  * (M * (Power(y, 2) + Power(z, 2)) - 3 * Power(a, 2) * x * dr_dx(x, y, z))
            + a * Power(r(x, y, z), 16) * (4 * a - 3 * y * dr_dx(x, y, z))
            + Power(r(x, y, z), 14)
                  * (2 * Power(a, 2) * (3 * Power(a, 2) + Power(z, 2))
                     - M * Power(x, 2) * (z * dr_dz(x, y, z) + y * dr_dy(x, y, z))
                     - (M * Power(x, 3) + 10 * Power(a, 3) * y) * dr_dx(x, y, z))
            + Power(a, 9) * Power(z, 4) * Power(r(x, y, z), 2)
                  * (Power(a, 3)
                     + 2 * (4 * a * M * x + 5 * Power(a, 2) * y + Power(M, 2) * y) * dr_dx(x, y, z))
            + Power(a, 8) * Power(z, 2) * Power(r(x, y, z), 3)
                  * (a * M * z * (2 * a * z + x * y * dr_dz(x, y, z))
                     + (3 * a * M * y * (Power(x, 2) + Power(y, 2))
                        + 2 * (7 * Power(a, 2) * x + Power(M, 2) * x + 10 * a * M * y)
                              * Power(z, 2))
                           * dr_dx(x, y, z))
            + Power(a, 4) * Power(r(x, y, z), 7)
                  * (M * z
                         * (a * x * y * (Power(a, 2) - 5 * Power(z, 2))
                            + 2 * M
                                  * (Power(x, 4) - Power(x, 2) * Power(y, 2) - 2 * Power(y, 4)
                                     + (Power(x, 2) - 4 * Power(y, 2)) * Power(z, 2)))
                         * dr_dz(x, y, z)
                     + 2 * M
                           * (a
                                  * (-2 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                     + a * Power(z, 2)
                                           * (Power(a, 2) - 2 * Power(x, 2) + 5 * Power(y, 2)
                                              + 6 * Power(z, 2)))
                              - M * y
                                    * (Power(x, 4) + 2 * Power(y, 2) * Power(z, 2)
                                       + Power(x, 2) * (Power(y, 2) - 3 * Power(z, 2)))
                                    * dr_dy(x, y, z))
                     + (M * y * (-Power(a, 3) + 2 * M * x * y) * (Power(x, 2) + Power(y, 2))
                        + (12 * Power(a, 4) * x + 10 * Power(M, 2) * Power(x, 3)
                           + 4 * Power(a, 3) * M * y + 5 * a * M * y * (Power(x, 2) + Power(y, 2)))
                              * Power(z, 2)
                        + 2 * (5 * Power(a, 2) * x + 6 * Power(M, 2) * x + 6 * a * M * y)
                              * Power(z, 4))
                           * dr_dx(x, y, z))
            + a * Power(r(x, y, z), 13)
                  * (M
                         * (-2 * a * Power(x, 2) + 8 * a * (Power(y, 2) + Power(z, 2))
                            - 3 * x * y * z * dr_dz(x, y, z)
                            + x * (Power(x, 2) - 3 * Power(y, 2)) * dr_dy(x, y, z))
                     - 2
                           * (3 * Power(a, 3) * x
                              + M * y * (3 * Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dx(x, y, z))
            + Power(a, 2) * Power(r(x, y, z), 12)
                  * (4 * a * (Power(a, 3) - M * x * y + 2 * a * Power(z, 2))
                     - M * (Power(x, 2) + 2 * Power(y, 2)) * z * dr_dz(x, y, z)
                     + M * y * (3 * Power(x, 2) - 2 * Power(y, 2)) * dr_dy(x, y, z)
                     + (-2 * a * y * (6 * Power(a, 2) + Power(z, 2))
                        + M * x * (Power(x, 2) - 4 * Power(y, 2) + 4 * Power(z, 2)))
                           * dr_dx(x, y, z))
            + Power(r(x, y, z), 11)
                  * (2 * a * M
                         * (a * Power(z, 2) * (Power(y, 2) + Power(z, 2))
                            + 2 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                            + Power(a, 3) * (-2 * Power(x, 2) + 5 * Power(y, 2) + 6 * Power(z, 2)))
                     + M * x
                           * (-(z
                                * (5 * Power(a, 3) * y
                                   + 2 * M * x * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                * dr_dz(x, y, z))
                              - 2 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                    * dr_dy(x, y, z))
                     + (-2 * Power(a, 6) * x + 4 * Power(a, 4) * x * Power(z, 2)
                        + 2 * Power(M, 2) * x * (Power(y, 2) + Power(z, 2))
                              * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                        - Power(a, 3) * M * y * (7 * (Power(x, 2) + Power(y, 2)) + 4 * Power(z, 2)))
                           * dr_dx(x, y, z))
            + Power(a, 6) * Power(z, 2) * Power(r(x, y, z), 5)
                  * (-2 * a * M
                         * (2 * M * x * y + a * (Power(x, 2) - 2 * Power(y, 2) - 4 * Power(z, 2)))
                     + M * (-(a * x * y) + 2 * M * (Power(x, 2) - 2 * Power(y, 2))) * z
                           * dr_dz(x, y, z)
                     - M * x * (2 * M * x * y + a * (Power(x, 2) - 3 * Power(y, 2)))
                           * dr_dy(x, y, z)
                     + 2
                           * (2 * Power(a, 4) * x + Power(a, 3) * M * y
                              + 9 * Power(a, 2) * x * Power(z, 2)
                              + Power(M, 2) * x * (Power(x, 2) + 2 * Power(y, 2) + 4 * Power(z, 2))
                              + a * M * y * (5 * Power(x, 2) + 3 * Power(y, 2) + 12 * Power(z, 2)))
                           * dr_dx(x, y, z))
            + Power(a, 7) * Power(z, 2) * Power(r(x, y, z), 4)
                  * (M * (2 * M * x * y + a * (Power(x, 2) - 2 * Power(y, 2))) * z * dr_dz(x, y, z)
                     + a
                           * (4 * a * (-(M * x * y) + a * Power(z, 2))
                              - M * Power(x, 2) * y * dr_dy(x, y, z))
                     + (2 * Power(a, 4) * y + 12 * Power(a, 2) * y * Power(z, 2)
                        + 2 * Power(M, 2) * y * (Power(x, 2) + Power(y, 2) + 4 * Power(z, 2))
                        + a * M * x * (5 * Power(x, 2) + 6 * Power(y, 2) + 28 * Power(z, 2)))
                           * dr_dx(x, y, z))
            + Power(a, 3) * Power(r(x, y, z), 8)
                  * (4 * a
                         * (-(a * M * x * y * (Power(a, 2) + Power(z, 2)))
                            + Power(a, 2) * Power(z, 2) * (2 * Power(a, 2) + Power(z, 2))
                            - Power(M, 2) * (x - y) * (x + y)
                                  * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                     - M * z
                           * (-(Power(a, 3) * (Power(x, 2) - 2 * Power(y, 2)))
                              + a * (Power(x, 2) + 2 * Power(y, 2)) * Power(z, 2)
                              + 2 * M * x * y * (2 * (Power(x, 2) + Power(y, 2)) + 5 * Power(z, 2)))
                           * dr_dz(x, y, z)
                     + M
                           * (x
                                  * (-(Power(a, 3) * x * y)
                                     - 2 * M * (Power(x, 2) - 3 * Power(y, 2))
                                           * (Power(x, 2) + Power(y, 2)))
                              + a * y * (3 * Power(x, 2) - 2 * Power(y, 2)) * Power(z, 2))
                           * dr_dy(x, y, z)
                     + (-(Power(a, 6) * y) + Power(a, 2) * y * Power(z, 4)
                        + Power(a, 3) * M * x * (Power(x, 2) + 2 * Power(y, 2) + 12 * Power(z, 2))
                        + a * M * x * Power(z, 2)
                              * (13 * Power(x, 2) + 8 * Power(y, 2) + 20 * Power(z, 2))
                        + 2 * Power(M, 2) * y
                              * (3 * Power(x, 4) + 2 * Power(x, 2) * Power(y, 2) - Power(y, 4)
                                 + 3 * (Power(x, 2) + Power(y, 2)) * Power(z, 2) + 4 * Power(z, 4)))
                           * dr_dx(x, y, z))
            + a * Power(r(x, y, z), 10)
                  * (Power(a, 7) - 8 * Power(a, 4) * M * x * y + 12 * Power(a, 5) * Power(z, 2)
                     + Power(a, 3) * Power(z, 4)
                     - 4 * a * Power(M, 2) * (x - y) * (x + y)
                           * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                     - M * z
                           * (-(Power(a, 3) * (Power(x, 2) - 4 * Power(y, 2)))
                              + a * Power(x, 2) * Power(z, 2)
                              + 6 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dz(x, y, z)
                     + M
                           * (Power(a, 3) * y * (3 * Power(x, 2) - 2 * Power(y, 2))
                              - a * Power(x, 2) * y * Power(z, 2)
                              + 2 * M * x * (Power(x, 2) - 3 * Power(y, 2))
                                    * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dy(x, y, z)
                     + (-6 * Power(a, 6) * y - 4 * Power(a, 4) * y * Power(z, 2)
                        + 2 * Power(M, 2) * y * (-3 * Power(x, 2) + Power(y, 2) + Power(z, 2))
                              * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                        + Power(a, 3) * M * x
                              * (3 * Power(x, 2) - 2 * Power(y, 2) + 12 * Power(z, 2))
                        + a * M * x * Power(z, 2)
                              * (3 * Power(x, 2) + 4 * (Power(y, 2) + Power(z, 2))))
                           * dr_dx(x, y, z))
            + Power(a, 2) * Power(r(x, y, z), 9)
                  * (2 * a * M
                         * (2 * M * x * y * Power(z, 2)
                            + Power(a, 3) * (-Power(x, 2) + 2 * Power(y, 2) + 4 * Power(z, 2))
                            + a * Power(z, 2) * (-Power(x, 2) + 4 * (Power(y, 2) + Power(z, 2))))
                     - M * z
                           * (y * (Power(a, 3) * x + 4 * M * y * (Power(x, 2) + Power(y, 2)))
                              + (3 * a * x * y + 2 * M * (Power(x, 2) + 2 * Power(y, 2)))
                                    * Power(z, 2))
                           * dr_dz(x, y, z)
                     + M
                           * (-(Power(a, 3) * (Power(x, 3) - 3 * x * Power(y, 2)))
                              + a * x * (Power(x, 2) - 3 * Power(y, 2)) * Power(z, 2)
                              + 2 * M * y
                                    * (4 * Power(x, 4)
                                       - 2 * Power(y, 2) * (Power(y, 2) + Power(z, 2))
                                       + Power(x, 2) * (2 * Power(y, 2) + 3 * Power(z, 2))))
                           * dr_dy(x, y, z)
                     + 2
                           * (Power(a, 2) * x * Power(z, 2) * (6 * Power(a, 2) + Power(z, 2))
                              + a * M * y
                                    * (-(Power(a, 2) * (Power(x, 2) + 3 * Power(y, 2)))
                                       + Power(z, 2) * (-Power(x, 2) + Power(y, 2) + Power(z, 2)))
                              + Power(M, 2) * x
                                    * (2 * Power(x, 4) - 4 * Power(y, 4) + 4 * Power(z, 4)
                                       + Power(x, 2) * (-2 * Power(y, 2) + 5 * Power(z, 2))))
                           * dr_dx(x, y, z))
            + Power(a, 5) * z * Power(r(x, y, z), 6)
                  * (M
                         * (a * (Power(x, 2) - 4 * Power(y, 2)) * Power(z, 2)
                            + 2 * M * x * y * (Power(x, 2) + Power(y, 2) - Power(z, 2)))
                         * dr_dz(x, y, z)
                     + z
                           * (2 * a
                                  * (Power(a, 4) - 4 * a * M * x * y
                                     + 2 * Power(M, 2) * (-Power(x, 2) + Power(y, 2))
                                     + 3 * Power(a, 2) * Power(z, 2))
                              + M
                                    * (-2 * M * Power(x, 3) + 3 * a * Power(x, 2) * y
                                       + 6 * M * x * Power(y, 2) - 2 * a * Power(y, 3))
                                    * dr_dy(x, y, z)
                              + (4 * Power(a, 3) * M * x + 4 * Power(a, 4) * y
                                 + 6 * Power(a, 2) * y * Power(z, 2)
                                 + 4 * Power(M, 2) * y
                                       * (3 * Power(x, 2) + Power(y, 2) + 3 * Power(z, 2))
                                 + a * M * x
                                       * (15 * Power(x, 2) + 10 * Power(y, 2) + 36 * Power(z, 2)))
                                    * dr_dx(x, y, z)))))
        / (Power(Power(a, 2) + Power(r(x, y, z), 2), 3)
           * (Power(a, 2) * Power(z, 2) + Power(r(x, y, z), 4))
           * Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
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
                        + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))));
  mixedK[0][1]
      = (M * Power(r(x, y, z), 2)
         * (3 * Power(a, 13) * Power(z, 4) * (y * dr_dy(x, y, z) - x * dr_dx(x, y, z))
            + 3 * a * Power(r(x, y, z), 16) * (-(y * dr_dy(x, y, z)) + x * dr_dx(x, y, z))
            - 2 * Power(r(x, y, z), 17) * (x * dr_dy(x, y, z) + y * dr_dx(x, y, z))
            + 4 * Power(a, 11) * Power(z, 4) * r(x, y, z)
                  * ((a * x + M * y) * dr_dy(x, y, z) + (-2 * M * x + a * y) * dr_dx(x, y, z))
            + 2 * Power(a, 9) * Power(z, 4) * Power(r(x, y, z), 2)
                  * (a * (3 * M * x + 5 * a * y) * dr_dy(x, y, z)
                     + (-5 * Power(a, 2) * x - 2 * Power(M, 2) * x + 5 * a * M * y)
                           * dr_dx(x, y, z))
            + 2 * Power(a, 8) * Power(z, 2) * Power(r(x, y, z), 3)
                  * (a * z
                         * (-(M * Power(x, 2) * dr_dz(x, y, z))
                            + (7 * a * x + 6 * M * y) * z * dr_dy(x, y, z))
                     + (-3 * a * M * x * (Power(x, 2) + Power(y, 2))
                        + (-14 * a * M * x + 7 * Power(a, 2) * y + 2 * Power(M, 2) * y)
                              * Power(z, 2))
                           * dr_dx(x, y, z))
            + 2 * Power(a, 4) * Power(r(x, y, z), 7)
                  * (2 * a * M
                         * (2 * M * Power(x, 2) * (Power(x, 2) + Power(y, 2))
                            + (-7 * a * x * y + 2 * M * (2 * Power(x, 2) + Power(y, 2)))
                                  * Power(z, 2))
                     + M * z
                           * (-(Power(a, 3) * Power(x, 2))
                              + a * (Power(x, 2) - 4 * Power(y, 2)) * Power(z, 2)
                              + 2 * M * x * y * (3 * (Power(x, 2) + Power(y, 2)) + 5 * Power(z, 2)))
                           * dr_dz(x, y, z)
                     + (2 * Power(M, 2) * Power(x, 3) * (Power(x, 2) + Power(y, 2))
                        + 2
                              * (3 * Power(a, 4) * x + 3 * Power(a, 3) * M * y
                                 + 5 * Power(M, 2) * x * Power(y, 2)
                                 - a * M * y * (Power(x, 2) + Power(y, 2)))
                              * Power(z, 2)
                        + a * (5 * a * x + 2 * M * y) * Power(z, 4))
                           * dr_dy(x, y, z)
                     + (M * x * (Power(a, 3) - 2 * M * x * y) * (Power(x, 2) + Power(y, 2))
                        + (2 * Power(a, 3) * M * x + 6 * Power(a, 4) * y
                           - 7 * a * M * x * (Power(x, 2) + Power(y, 2))
                           + 2 * Power(M, 2) * y * (7 * Power(x, 2) + 2 * Power(y, 2)))
                              * Power(z, 2)
                        + (-10 * a * M * x + 5 * Power(a, 2) * y + 12 * Power(M, 2) * y)
                              * Power(z, 4))
                           * dr_dx(x, y, z))
            + 2 * Power(a, 2) * Power(r(x, y, z), 12)
                  * (2 * a * M * (x - y) * (x + y) + M * x * y * z * dr_dz(x, y, z)
                     + (-(a * y * (6 * Power(a, 2) + Power(z, 2)))
                        + M * x * (2 * Power(x, 2) + 7 * Power(y, 2) + 6 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     + (6 * Power(a, 3) * x + a * x * Power(z, 2)
                        + M * y * (Power(x, 2) - 4 * Power(y, 2) - 2 * Power(z, 2)))
                           * dr_dx(x, y, z))
            + 2 * Power(r(x, y, z), 14)
                  * (-(M * x * y * z * dr_dz(x, y, z))
                     + (-5 * Power(a, 3) * y + M * x * (Power(x, 2) + Power(z, 2))) * dr_dy(x, y, z)
                     - (-5 * Power(a, 3) * x
                        + M * y * (2 * Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dx(x, y, z))
            + 2 * a * Power(r(x, y, z), 13)
                  * (-10 * a * M * x * y + M * (Power(x, 2) - 2 * Power(y, 2)) * z * dr_dz(x, y, z)
                     + (-3 * Power(a, 3) * x + 2 * M * (x - y) * y * (x + y)) * dr_dy(x, y, z)
                     + (-3 * Power(a, 3) * y
                        + M * x * (3 * Power(x, 2) - Power(y, 2) + 2 * Power(z, 2)))
                           * dr_dx(x, y, z))
            + 2 * Power(r(x, y, z), 11)
                  * (2 * a * M * y
                         * (-7 * Power(a, 3) * x - a * x * Power(z, 2)
                            + 2 * M * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                     - M * z
                           * (-(Power(a, 3) * (Power(x, 2) - 4 * Power(y, 2)))
                              + 2 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dz(x, y, z)
                     - (Power(a, 6) * x - 2 * Power(a, 4) * x * Power(z, 2)
                        + 2 * Power(a, 3) * M * y * (Power(x, 2) + Power(y, 2) - Power(z, 2))
                        + 2 * Power(M, 2) * x * Power(y, 2)
                              * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dy(x, y, z)
                     + (-(Power(a, 6) * y) + 2 * Power(a, 4) * y * Power(z, 2)
                        + 2 * Power(M, 2) * y * (Power(y, 2) + Power(z, 2))
                              * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                        + Power(a, 3) * M * x * (5 * (Power(x, 2) + Power(y, 2)) + 6 * Power(z, 2)))
                           * dr_dx(x, y, z))
            + 2 * Power(a, 7) * Power(z, 2) * Power(r(x, y, z), 4)
                  * (M * x * (-2 * M * x + 3 * a * y) * z * dr_dz(x, y, z)
                     + a
                           * (2 * a * M * (x - y) * (x + y)
                              + (a * y * (Power(a, 2) + 6 * Power(z, 2))
                                 + M * x * (2 * Power(x, 2) + Power(y, 2) + 10 * Power(z, 2)))
                                    * dr_dy(x, y, z))
                     - (Power(a, 4) * x + 6 * Power(a, 2) * x * Power(z, 2)
                        + 2 * Power(M, 2) * x * (Power(x, 2) + Power(y, 2) + 4 * Power(z, 2))
                        - a * M * y * (3 * Power(x, 2) + 4 * Power(y, 2) + 18 * Power(z, 2)))
                           * dr_dx(x, y, z))
            + 2 * Power(a, 6) * Power(z, 2) * Power(r(x, y, z), 5)
                  * (2 * a * M * x * (2 * M * x - 3 * a * y)
                     + M * (6 * M * x * y - a * (Power(x, 2) + 2 * Power(y, 2))) * z
                           * dr_dz(x, y, z)
                     + (2
                            * (Power(a, 4) * x + Power(M, 2) * Power(x, 3) + Power(a, 3) * M * y
                               - 2 * a * M * Power(x, 2) * y)
                        + 3 * a * (3 * a * x + 2 * M * y) * Power(z, 2))
                           * dr_dy(x, y, z)
                     + (2 * Power(a, 4) * y + 9 * Power(a, 2) * y * Power(z, 2)
                        + 2 * Power(M, 2) * y * (Power(y, 2) + 4 * Power(z, 2))
                        - a * M * x * (9 * Power(x, 2) + 5 * Power(y, 2) + 18 * Power(z, 2)))
                           * dr_dx(x, y, z))
            + Power(a, 3) * Power(r(x, y, z), 8)
                  * (4 * a * M
                         * (Power(a, 3) * (x - y) * (x + y) + a * (x - y) * (x + y) * Power(z, 2)
                            - 4 * M * x * y * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                     + 2 * M * z
                           * (y * (3 * Power(a, 3) * x - 4 * M * y * (Power(x, 2) + Power(y, 2)))
                              + (a * x * y + 2 * M * (Power(x, 2) - 4 * Power(y, 2))) * Power(z, 2))
                           * dr_dz(x, y, z)
                     + (-(Power(a, 6) * y) + Power(a, 2) * y * Power(z, 4)
                        - 8 * Power(M, 2) * y * (Power(x, 2) + Power(y, 2))
                              * (2 * Power(x, 2) + Power(z, 2))
                        + 2 * a * M * x * Power(z, 2)
                              * (2 * Power(x, 2) + 7 * Power(y, 2) + 6 * Power(z, 2))
                        + 2 * Power(a, 3) * M * x
                              * (2 * Power(x, 2) + Power(y, 2) + 10 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     + (Power(a, 6) * x - Power(a, 2) * x * Power(z, 4)
                        - 2 * Power(a, 3) * M * y * (Power(x, 2) - 2 * Power(z, 2))
                        + 2 * a * M * y * Power(z, 2)
                              * (13 * Power(x, 2) + 8 * Power(y, 2) + 14 * Power(z, 2))
                        - 4 * Power(M, 2) * x
                              * (2 * Power(x, 4) - 2 * Power(y, 4)
                                 + 5 * (Power(x, 2) + Power(y, 2)) * Power(z, 2) + 4 * Power(z, 4)))
                           * dr_dx(x, y, z))
            + 2 * a * Power(r(x, y, z), 10)
                  * (4 * a * M
                         * (Power(a, 3) * (x - y) * (x + y)
                            - 2 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                     + M * z
                           * (a * x * y * (5 * Power(a, 2) - Power(z, 2))
                              + 2 * M * (Power(x, 2) - 2 * Power(y, 2))
                                    * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dz(x, y, z)
                     + (-3 * Power(a, 6) * y - 2 * Power(a, 4) * y * Power(z, 2)
                        + a * M * x * Power(z, 2) * (Power(x, 2) + Power(z, 2))
                        + 4 * Power(M, 2) * (x - y) * y * (x + y)
                              * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                        + Power(a, 3) * M * x
                              * (3 * Power(x, 2) + 8 * Power(y, 2) + 12 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     + (3 * Power(a, 6) * x
                        + Power(a, 3) * M * y * (2 * Power(x, 2) - 3 * Power(y, 2))
                        + 2 * Power(a, 4) * x * Power(z, 2)
                        - 2 * Power(M, 2) * x * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                              * (4 * Power(y, 2) + Power(z, 2))
                        + a * M * y * Power(z, 2)
                              * (2 * Power(x, 2) + 3 * (Power(y, 2) + Power(z, 2))))
                           * dr_dx(x, y, z))
            + 2 * Power(a, 2) * Power(r(x, y, z), 9)
                  * (2 * a * M
                         * (-(a * x * y * (3 * Power(a, 2) + 5 * Power(z, 2)))
                            + 2 * M
                                  * (Power(Power(x, 2) + Power(y, 2), 2)
                                     + (Power(x, 2) + 2 * Power(y, 2)) * Power(z, 2)))
                     + M * z
                           * (-(Power(a, 3) * (Power(x, 2) + 2 * Power(y, 2)))
                              + a * (Power(x, 2) - 2 * Power(y, 2)) * Power(z, 2)
                              + 2 * M * x * y * (2 * (Power(x, 2) + Power(y, 2)) + Power(z, 2)))
                           * dr_dz(x, y, z)
                     + (-2 * M * x
                            * (2 * Power(a, 3) * x * y
                               + M * (Power(x, 2) - 5 * Power(y, 2)) * (Power(x, 2) + Power(y, 2)))
                        + 2
                              * (3 * Power(a, 4) * x + 3 * Power(a, 3) * M * y
                                 + a * M * (x - y) * y * (x + y)
                                 - Power(M, 2) * x * (Power(x, 2) - 4 * Power(y, 2)))
                              * Power(z, 2)
                        + Power(a, 2) * x * Power(z, 4))
                           * dr_dy(x, y, z)
                     + (6 * Power(a, 4) * y * Power(z, 2) + Power(a, 2) * y * Power(z, 4)
                        - a * M * x * Power(z, 2)
                              * (Power(x, 2) + 5 * Power(y, 2) + 2 * Power(z, 2))
                        + Power(a, 3) * M * x
                              * (3 * Power(x, 2) + 7 * Power(y, 2) + 6 * Power(z, 2))
                        + 2 * Power(M, 2) * y
                              * (5 * Power(x, 4) - Power(y, 4) + 3 * Power(y, 2) * Power(z, 2)
                                 + 4 * Power(z, 4)
                                 + 4 * Power(x, 2) * (Power(y, 2) + 2 * Power(z, 2))))
                           * dr_dx(x, y, z))
            - 2 * Power(r(x, y, z), 15)
                  * (2 * M * x * y + 3 * Power(a, 2) * (x * dr_dy(x, y, z) + y * dr_dx(x, y, z)))
            - 2 * Power(a, 5) * z * Power(r(x, y, z), 6)
                  * (M
                         * (-5 * a * x * y * Power(z, 2)
                            + 2 * M
                                  * (Power(x, 4) + 2 * Power(y, 2) * Power(z, 2)
                                     + Power(x, 2) * (Power(y, 2) + Power(z, 2))))
                         * dr_dz(x, y, z)
                     + z
                           * (4 * a * M * (2 * M * x * y + a * (-Power(x, 2) + Power(y, 2)))
                              - (3 * Power(a, 3) * M * x + 3 * a * M * Power(x, 3)
                                 + 2 * Power(a, 4) * y - 8 * Power(M, 2) * Power(x, 2) * y
                                 + 8 * a * M * x * Power(y, 2)
                                 + 3 * a * (4 * M * x + a * y) * Power(z, 2))
                                    * dr_dy(x, y, z)
                              + (2 * Power(a, 4) * x - Power(a, 3) * M * y
                                 + 3 * Power(a, 2) * x * Power(z, 2)
                                 + 2 * Power(M, 2) * x
                                       * (5 * Power(x, 2) + Power(y, 2) + 6 * Power(z, 2))
                                 - a * M * y
                                       * (14 * Power(x, 2) + 9 * Power(y, 2) + 24 * Power(z, 2)))
                                    * dr_dx(x, y, z)))))
        / (Power(Power(a, 2) + Power(r(x, y, z), 2), 3)
           * (Power(a, 2) * Power(z, 2) + Power(r(x, y, z), 4))
           * Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
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
                        + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))));
  mixedK[0][2]
      = (M * r(x, y, z)
         * (-3 * a * y * Power(r(x, y, z), 15) * dr_dz(x, y, z)
            + 2 * Power(a, 12) * Power(z, 5) * dr_dx(x, y, z)
            - 2 * Power(r(x, y, z), 16) * (x * dr_dz(x, y, z) + z * dr_dx(x, y, z))
            + 3 * Power(a, 10) * Power(z, 4) * r(x, y, z)
                  * (a * y * dr_dz(x, y, z) + 2 * M * z * dr_dx(x, y, z))
            + Power(r(x, y, z), 13)
                  * ((-7 * Power(a, 3) * y + 2 * M * x * (Power(x, 2) + Power(y, 2)))
                         * dr_dz(x, y, z)
                     - 2 * x * z * (Power(a, 2) + M * y * dr_dy(x, y, z))
                     - 2 * M * z * (2 * Power(x, 2) + Power(y, 2) + Power(z, 2)) * dr_dx(x, y, z))
            + Power(a, 8) * Power(z, 3) * Power(r(x, y, z), 3)
                  * (-2 * a * (a * x + 2 * M * y) + (8 * M * x + 7 * a * y) * z * dr_dz(x, y, z)
                     + 4 * M * (Power(x, 2) + Power(y, 2) + 6 * Power(z, 2)) * dr_dx(x, y, z))
            - Power(a, 3) * z * Power(r(x, y, z), 7)
                  * (2 * a
                         * (Power(a, 4) * x + 2 * Power(a, 3) * M * y
                            + 3 * Power(a, 2) * x * Power(z, 2) + 2 * a * M * y * Power(z, 2)
                            + 4 * Power(M, 2) * x * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                     - z
                           * (8 * Power(a, 3) * M * x + 2 * Power(a, 4) * y
                              + Power(a, 2) * y * Power(z, 2)
                              + 4 * a * M * x * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2))
                              - 8 * Power(M, 2) * y * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                           * dr_dz(x, y, z)
                     + 4 * M
                           * (M * Power(x, 2) * (Power(x, 2) + Power(y, 2))
                              + y * (-(a * x) + 2 * M * y) * Power(z, 2))
                           * dr_dy(x, y, z)
                     - 4 * M
                           * (-((Power(a, 3) - M * x * y) * (Power(x, 2) + Power(y, 2)))
                              + (-2 * Power(a, 3) + 4 * a * Power(x, 2) - 2 * M * x * y
                                 + 3 * a * Power(y, 2))
                                    * Power(z, 2)
                              + 6 * a * Power(z, 4))
                           * dr_dx(x, y, z))
            - 2 * z * Power(r(x, y, z), 10)
                  * (a
                         * (10 * Power(a, 3) * M * x + 3 * Power(a, 4) * y
                            + 2 * a * M * x * Power(z, 2)
                            - 4 * Power(M, 2) * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                     + 2 * z
                           * (-(Power(a, 4) * x) + Power(a, 3) * M * y
                              + Power(M, 2) * x * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dz(x, y, z)
                     + 2 * M * y
                           * (Power(a, 3) * y + M * x * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dy(x, y, z)
                     + 2
                           * (2 * Power(a, 6) + Power(a, 3) * M * x * y
                              - Power(M, 2) * (Power(y, 2) + Power(z, 2))
                                    * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dx(x, y, z))
            + 2 * Power(a, 6) * Power(z, 3) * Power(r(x, y, z), 4)
                  * ((5 * Power(a, 2) * x + 2 * Power(M, 2) * x + 2 * a * M * y) * z
                         * dr_dz(x, y, z)
                     - a * (a * (4 * M * x + 3 * a * y) + M * Power(x, 2) * dr_dy(x, y, z))
                     + (a * M * x * y + 6 * Power(a, 2) * Power(z, 2)
                        + 2 * Power(M, 2) * (Power(x, 2) + Power(y, 2) + 4 * Power(z, 2)))
                           * dr_dx(x, y, z))
            + 2 * Power(a, 2) * z * Power(r(x, y, z), 8)
                  * (-(a
                       * (4 * Power(a, 3) * M * x + 3 * Power(a, 4) * y
                          + 8 * a * M * x * Power(z, 2) + Power(a, 2) * y * Power(z, 2)
                          - 4 * Power(M, 2) * y * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2))))
                     + z
                           * (2 * Power(a, 3) * (2 * a * x + M * y)
                              + (Power(a, 2) * x - 2 * Power(M, 2) * x - 2 * a * M * y)
                                    * Power(z, 2))
                           * dr_dz(x, y, z)
                     + M
                           * (x * (-(Power(a, 3) * x) + 6 * M * y * (Power(x, 2) + Power(y, 2)))
                              + (4 * M * x * y + a * (Power(x, 2) - 2 * Power(y, 2))) * Power(z, 2))
                           * dr_dy(x, y, z)
                     + (-Power(a, 6) + Power(a, 3) * M * x * y - 3 * a * M * x * y * Power(z, 2)
                        + Power(a, 2) * Power(z, 4)
                        + 2 * Power(M, 2)
                              * (2 * Power(x, 4) - Power(y, 4) + 3 * Power(y, 2) * Power(z, 2)
                                 + 4 * Power(z, 4) + Power(x, 2) * (Power(y, 2) + 5 * Power(z, 2))))
                           * dr_dx(x, y, z))
            - 2 * a * Power(r(x, y, z), 12)
                  * ((Power(a, 3) * x + 2 * M * y * Power(z, 2)) * dr_dz(x, y, z)
                     + z
                           * (a * (8 * M * x + a * y)
                              - M * (Power(x, 2) - 2 * Power(y, 2)) * dr_dy(x, y, z)
                              + 3 * (2 * Power(a, 3) + M * x * y) * dr_dx(x, y, z)))
            - 4 * Power(r(x, y, z), 14)
                  * (M * x * z + Power(a, 2) * (x * dr_dz(x, y, z) + 2 * z * dr_dx(x, y, z)))
            + 2 * Power(a, 8) * Power(z, 3) * Power(r(x, y, z), 2)
                  * (-(Power(a, 3) * y)
                     + 2 * z
                           * (a * (a * x + M * y) * dr_dz(x, y, z)
                              + (2 * Power(a, 2) + Power(M, 2)) * z * dr_dx(x, y, z)))
            + Power(a, 2) * Power(r(x, y, z), 11)
                  * (-2 * a * (3 * a * x + 2 * M * y) * z
                     + (-5 * Power(a, 3) * y - 2 * a * y * Power(z, 2)
                        + 4 * M * x * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                           * dr_dz(x, y, z)
                     + 4 * M * z
                           * (x * y * dr_dy(x, y, z)
                              - (2 * Power(x, 2) + 3 * Power(y, 2) + 2 * Power(z, 2))
                                    * dr_dx(x, y, z)))
            + Power(a, 5) * Power(z, 2) * Power(r(x, y, z), 5)
                  * ((2 * Power(a, 4) * y + 5 * Power(a, 2) * y * Power(z, 2)
                      - 8 * Power(M, 2) * y * Power(z, 2)
                      + 2 * a * M * x * (Power(x, 2) + Power(y, 2) + 8 * Power(z, 2)))
                         * dr_dz(x, y, z)
                     + 2 * z
                           * (-(a * (3 * Power(a, 2) * x + 4 * Power(M, 2) * x + 4 * a * M * y))
                              + M * x * (-2 * M * x + 3 * a * y) * dr_dy(x, y, z)
                              + M
                                    * (-Power(a, 3) + 8 * a * Power(x, 2) + 2 * M * x * y
                                       + 5 * a * Power(y, 2) + 18 * a * Power(z, 2))
                                    * dr_dx(x, y, z)))
            + 2 * Power(a, 4) * z * Power(r(x, y, z), 6)
                  * (-(Power(a, 5) * y)
                     + a * (-10 * a * M * x - 3 * Power(a, 2) * y + 4 * Power(M, 2) * y)
                           * Power(z, 2)
                     + 2 * z
                           * ((Power(a, 4) * x + Power(a, 3) * M * y
                               + 2 * Power(a, 2) * x * Power(z, 2) - a * M * y * Power(z, 2)
                               + Power(M, 2) * x * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                  * dr_dz(x, y, z)
                              + M * y * (3 * M * x - a * y) * z * dr_dy(x, y, z)
                              + z
                                    * (M * (5 * M * Power(x, 2) - a * x * y + 2 * M * Power(y, 2))
                                       + 2 * (Power(a, 2) + 3 * Power(M, 2)) * Power(z, 2))
                                    * dr_dx(x, y, z)))
            - a * Power(r(x, y, z), 9)
                  * ((Power(a, 6) * y + 2 * Power(a, 4) * y * Power(z, 2)
                      - 2 * a * M * x * (Power(x, 2) + Power(y, 2)) * Power(z, 2)
                      + 8 * Power(M, 2) * y * Power(z, 2)
                            * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                      - 2 * Power(a, 3) * M * x * (Power(x, 2) + Power(y, 2) + 8 * Power(z, 2)))
                         * dr_dz(x, y, z)
                     + 2 * z
                           * (3 * Power(a, 5) * x + 4 * Power(a, 4) * M * y
                              + Power(a, 3) * x * Power(z, 2)
                              + 4 * a * Power(M, 2) * x * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                              + M
                                    * (a * x * y * (-3 * Power(a, 2) + Power(z, 2))
                                       - 2 * M * (Power(x, 2) - 2 * Power(y, 2))
                                             * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                    * dr_dy(x, y, z)
                              + M
                                    * (6 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                       + Power(a, 3)
                                             * (4 * Power(x, 2) + 7 * Power(y, 2) + 6 * Power(z, 2))
                                       - a * Power(z, 2)
                                             * (2 * Power(x, 2) + 3 * (Power(y, 2) + Power(z, 2))))
                                    * dr_dx(x, y, z)))))
        / (Power(Power(a, 2) + Power(r(x, y, z), 2), 2)
           * (Power(a, 2) * Power(z, 2) + Power(r(x, y, z), 4))
           * Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
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
                        + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))));
  mixedK[1][0]
      = (M * Power(r(x, y, z), 2)
         * (3 * Power(a, 13) * Power(z, 4) * (y * dr_dy(x, y, z) - x * dr_dx(x, y, z))
            + 3 * a * Power(r(x, y, z), 16) * (-(y * dr_dy(x, y, z)) + x * dr_dx(x, y, z))
            - 2 * Power(r(x, y, z), 17) * (x * dr_dy(x, y, z) + y * dr_dx(x, y, z))
            + 4 * Power(a, 11) * Power(z, 4) * r(x, y, z)
                  * ((a * x + 2 * M * y) * dr_dy(x, y, z) + (-(M * x) + a * y) * dr_dx(x, y, z))
            + 2 * Power(a, 9) * Power(z, 4) * Power(r(x, y, z), 2)
                  * ((5 * a * M * x + 5 * Power(a, 2) * y + 2 * Power(M, 2) * y) * dr_dy(x, y, z)
                     + a * (-5 * a * x + 3 * M * y) * dr_dx(x, y, z))
            - 2 * a * Power(r(x, y, z), 13)
                  * (10 * a * M * x * y + M * (-2 * Power(x, 2) + Power(y, 2)) * z * dr_dz(x, y, z)
                     + (3 * Power(a, 3) * x
                        + M * y * (-Power(x, 2) + 3 * Power(y, 2) + 2 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     + (3 * Power(a, 3) * y + 2 * M * x * (-Power(x, 2) + Power(y, 2)))
                           * dr_dx(x, y, z))
            + 2 * Power(a, 8) * Power(z, 2) * Power(r(x, y, z), 3)
                  * (a * M * Power(y, 2) * z * dr_dz(x, y, z)
                     + (3 * a * M * y * (Power(x, 2) + Power(y, 2))
                        + (7 * Power(a, 2) * x + 2 * Power(M, 2) * x + 14 * a * M * y)
                              * Power(z, 2))
                           * dr_dy(x, y, z)
                     + a * (-6 * M * x + 7 * a * y) * Power(z, 2) * dr_dx(x, y, z))
            + 2 * Power(a, 6) * Power(z, 2) * Power(r(x, y, z), 5)
                  * (-2 * a * M * y * (3 * a * x + 2 * M * y)
                     + M * (6 * M * x * y + a * (2 * Power(x, 2) + Power(y, 2))) * z
                           * dr_dz(x, y, z)
                     + (2 * Power(a, 4) * x + 9 * Power(a, 2) * x * Power(z, 2)
                        + 2 * Power(M, 2) * x * (Power(x, 2) + 4 * Power(z, 2))
                        + a * M * y * (5 * Power(x, 2) + 9 * Power(y, 2) + 18 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     + (2
                            * (-(Power(a, 3) * M * x) + Power(a, 4) * y
                               + 2 * a * M * x * Power(y, 2) + Power(M, 2) * Power(y, 3))
                        + 3 * a * (-2 * M * x + 3 * a * y) * Power(z, 2))
                           * dr_dx(x, y, z))
            + 2 * Power(a, 2) * Power(r(x, y, z), 12)
                  * (2 * a * M * (x - y) * (x + y) + M * x * y * z * dr_dz(x, y, z)
                     - (4 * M * Power(x, 3) + 6 * Power(a, 3) * y - M * x * Power(y, 2)
                        + (2 * M * x + a * y) * Power(z, 2))
                           * dr_dy(x, y, z)
                     + (6 * Power(a, 3) * x + 7 * M * Power(x, 2) * y + 2 * M * Power(y, 3)
                        + (a * x + 6 * M * y) * Power(z, 2))
                           * dr_dx(x, y, z))
            - 2 * Power(r(x, y, z), 14)
                  * (M * x * y * z * dr_dz(x, y, z)
                     + (5 * Power(a, 3) * y + M * x * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2)))
                           * dr_dy(x, y, z)
                     - (5 * Power(a, 3) * x + M * y * (Power(y, 2) + Power(z, 2))) * dr_dx(x, y, z))
            - 2 * Power(r(x, y, z), 11)
                  * (2 * a * M * x
                         * (a * y * (7 * Power(a, 2) + Power(z, 2))
                            + 2 * M * x * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                     + M * z
                           * (Power(a, 3) * (-4 * Power(x, 2) + Power(y, 2))
                              + 2 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dz(x, y, z)
                     + (Power(a, 6) * x - 2 * Power(a, 4) * x * Power(z, 2)
                        - 2 * Power(M, 2) * x * (Power(x, 2) + Power(z, 2))
                              * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                        + Power(a, 3) * M * y * (5 * (Power(x, 2) + Power(y, 2)) + 6 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     + (Power(a, 6) * y - 2 * Power(a, 4) * y * Power(z, 2)
                        - 2 * Power(a, 3) * M * x * (Power(x, 2) + Power(y, 2) - Power(z, 2))
                        + 2 * Power(M, 2) * Power(x, 2) * y
                              * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dx(x, y, z))
            + Power(a, 3) * Power(r(x, y, z), 8)
                  * (4 * a * M
                         * (Power(a, 3) * (x - y) * (x + y) + a * (x - y) * (x + y) * Power(z, 2)
                            - 4 * M * x * y * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                     + 2 * M * z
                           * (x * (3 * Power(a, 3) * y + 4 * M * x * (Power(x, 2) + Power(y, 2)))
                              + (8 * M * Power(x, 2) + a * x * y - 2 * M * Power(y, 2))
                                    * Power(z, 2))
                           * dr_dz(x, y, z)
                     + (-(Power(a, 6) * y) + Power(a, 2) * y * Power(z, 4)
                        - 2 * Power(a, 3) * M * x * (Power(y, 2) - 2 * Power(z, 2))
                        + 2 * a * M * x * Power(z, 2)
                              * (8 * Power(x, 2) + 13 * Power(y, 2) + 14 * Power(z, 2))
                        + 4 * Power(M, 2) * y
                              * (-2 * Power(x, 4) + 2 * Power(y, 4)
                                 + 5 * (Power(x, 2) + Power(y, 2)) * Power(z, 2) + 4 * Power(z, 4)))
                           * dr_dy(x, y, z)
                     + (Power(a, 6) * x - Power(a, 2) * x * Power(z, 4)
                        + 8 * Power(M, 2) * x * (Power(x, 2) + Power(y, 2))
                              * (2 * Power(y, 2) + Power(z, 2))
                        + 2 * a * M * y * Power(z, 2)
                              * (7 * Power(x, 2) + 2 * Power(y, 2) + 6 * Power(z, 2))
                        + 2 * Power(a, 3) * M * y
                              * (Power(x, 2) + 2 * Power(y, 2) + 10 * Power(z, 2)))
                           * dr_dx(x, y, z))
            + 2 * a * Power(r(x, y, z), 10)
                  * (4 * a * M
                         * (Power(a, 3) * (x - y) * (x + y)
                            - 2 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                     + M * z
                           * (a * x * y * (5 * Power(a, 2) - Power(z, 2))
                              + 2 * M * (2 * Power(x, 2) - Power(y, 2))
                                    * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dz(x, y, z)
                     + (-3 * Power(a, 6) * y
                        + Power(a, 3) * M * x * (-3 * Power(x, 2) + 2 * Power(y, 2))
                        - 2 * Power(a, 4) * y * Power(z, 2)
                        + 2 * Power(M, 2) * y * (4 * Power(x, 2) + Power(z, 2))
                              * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                        + a * M * x * Power(z, 2)
                              * (3 * Power(x, 2) + 2 * Power(y, 2) + 3 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     + (3 * Power(a, 6) * x + 2 * Power(a, 4) * x * Power(z, 2)
                        + a * M * y * Power(z, 2) * (Power(y, 2) + Power(z, 2))
                        + 4 * Power(M, 2) * x * (x - y) * (x + y)
                              * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                        + Power(a, 3) * M * y
                              * (8 * Power(x, 2) + 3 * Power(y, 2) + 12 * Power(z, 2)))
                           * dr_dx(x, y, z))
            - 2 * Power(a, 2) * Power(r(x, y, z), 9)
                  * (2 * a * M
                         * (a * x * y * (3 * Power(a, 2) + 5 * Power(z, 2))
                            + 2 * M
                                  * (Power(Power(x, 2) + Power(y, 2), 2)
                                     + (2 * Power(x, 2) + Power(y, 2)) * Power(z, 2)))
                     - M * z
                           * (Power(a, 3) * (2 * Power(x, 2) + Power(y, 2))
                              + a * (2 * Power(x, 2) - Power(y, 2)) * Power(z, 2)
                              + 2 * M * x * y * (2 * (Power(x, 2) + Power(y, 2)) + Power(z, 2)))
                           * dr_dz(x, y, z)
                     + (-(Power(a, 2) * x * Power(z, 2) * (6 * Power(a, 2) + Power(z, 2)))
                        + 2 * Power(M, 2) * x
                              * (Power(x, 4) - 5 * Power(y, 4) - 8 * Power(y, 2) * Power(z, 2)
                                 - 4 * Power(z, 4)
                                 - Power(x, 2) * (4 * Power(y, 2) + 3 * Power(z, 2)))
                        + a * M * y
                              * (-(Power(z, 2) * (5 * Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                                 + Power(a, 2)
                                       * (7 * Power(x, 2) + 3 * Power(y, 2) + 6 * Power(z, 2))))
                           * dr_dy(x, y, z)
                     - (4 * Power(a, 3) * M * x * Power(y, 2)
                        - 2 * a * M * x * (3 * Power(a, 2) - Power(x, 2) + Power(y, 2))
                              * Power(z, 2)
                        + Power(a, 2) * y * Power(z, 2) * (6 * Power(a, 2) + Power(z, 2))
                        - 2 * Power(M, 2) * y
                              * (-5 * Power(x, 4) - 4 * Power(x, 2) * (Power(y, 2) + Power(z, 2))
                                 + Power(y, 2) * (Power(y, 2) + Power(z, 2))))
                           * dr_dx(x, y, z))
            + 2 * Power(a, 4) * Power(r(x, y, z), 7)
                  * (2 * a * M
                         * (-7 * a * x * y * Power(z, 2)
                            - 2 * M
                                  * (Power(y, 4) + 2 * Power(y, 2) * Power(z, 2)
                                     + Power(x, 2) * (Power(y, 2) + Power(z, 2))))
                     + M * z
                           * (y * (Power(a, 3) * y + 6 * M * x * (Power(x, 2) + Power(y, 2)))
                              + (4 * a * Power(x, 2) + 10 * M * x * y - a * Power(y, 2))
                                    * Power(z, 2))
                           * dr_dz(x, y, z)
                     + (-(M * y * (Power(a, 3) + 2 * M * x * y) * (Power(x, 2) + Power(y, 2)))
                        + (6 * Power(a, 4) * x - 2 * Power(a, 3) * M * y
                           + 7 * a * M * y * (Power(x, 2) + Power(y, 2))
                           + 2 * Power(M, 2) * x * (2 * Power(x, 2) + 7 * Power(y, 2)))
                              * Power(z, 2)
                        + (5 * Power(a, 2) * x + 12 * Power(M, 2) * x + 10 * a * M * y)
                              * Power(z, 4))
                           * dr_dy(x, y, z)
                     + (2 * a * M * x * Power(z, 2)
                            * (-3 * Power(a, 2) + Power(x, 2) + Power(y, 2) - Power(z, 2))
                        + Power(a, 2) * y * Power(z, 2) * (6 * Power(a, 2) + 5 * Power(z, 2))
                        + 2 * Power(M, 2)
                              * (Power(y, 5) + Power(x, 2) * (Power(y, 3) + 5 * y * Power(z, 2))))
                           * dr_dx(x, y, z))
            - 2 * Power(r(x, y, z), 15)
                  * (2 * M * x * y + 3 * Power(a, 2) * (x * dr_dy(x, y, z) + y * dr_dx(x, y, z)))
            + 2 * Power(a, 7) * Power(z, 2) * Power(r(x, y, z), 4)
                  * (M * y * (3 * a * x + 2 * M * y) * z * dr_dz(x, y, z)
                     + (Power(a, 4) * y + 6 * Power(a, 2) * y * Power(z, 2)
                        + 2 * Power(M, 2) * y * (Power(x, 2) + Power(y, 2) + 4 * Power(z, 2))
                        + a * M * x * (4 * Power(x, 2) + 3 * Power(y, 2) + 18 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     + a
                           * (2 * a * M * (x - y) * (x + y)
                              + (-(Power(a, 3) * x) - 6 * a * x * Power(z, 2)
                                 + M * y * (Power(x, 2) + 2 * Power(y, 2) + 10 * Power(z, 2)))
                                    * dr_dx(x, y, z)))
            + 2 * Power(a, 5) * z * Power(r(x, y, z), 6)
                  * (M
                         * (2 * M * Power(y, 2) * (Power(x, 2) + Power(y, 2))
                            + (5 * a * x * y + 2 * M * (2 * Power(x, 2) + Power(y, 2)))
                                  * Power(z, 2))
                         * dr_dz(x, y, z)
                     + z
                           * (4 * a * M * (-2 * M * x * y + a * (x - y) * (x + y))
                              + (Power(a, 3) * M * x + 2 * Power(a, 4) * y
                                 + 3 * Power(a, 2) * y * Power(z, 2)
                                 + 2 * Power(M, 2) * y
                                       * (Power(x, 2) + 5 * Power(y, 2) + 6 * Power(z, 2))
                                 + a * M * x
                                       * (9 * Power(x, 2) + 14 * Power(y, 2) + 24 * Power(z, 2)))
                                    * dr_dy(x, y, z)
                              + (-2 * Power(a, 4) * x + 3 * Power(a, 3) * M * y
                                 + 8 * Power(M, 2) * x * Power(y, 2)
                                 - 3 * Power(a, 2) * x * Power(z, 2)
                                 + a * M * y
                                       * (8 * Power(x, 2) + 3 * Power(y, 2) + 12 * Power(z, 2)))
                                    * dr_dx(x, y, z)))))
        / (Power(Power(a, 2) + Power(r(x, y, z), 2), 3)
           * (Power(a, 2) * Power(z, 2) + Power(r(x, y, z), 4))
           * Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
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
                        + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))));
  mixedK[1][1]
      = (2 * M * Power(r(x, y, z), 2)
         * (Power(r(x, y, z), 18) - 3 * Power(a, 13) * x * Power(z, 4) * dr_dy(x, y, z)
            + 2 * Power(a, 11) * (-3 * M * x + 2 * a * y) * Power(z, 4) * r(x, y, z)
                  * dr_dy(x, y, z)
            - 2 * y * Power(r(x, y, z), 17) * dr_dy(x, y, z)
            + a * Power(r(x, y, z), 16) * (4 * a + 3 * x * dr_dy(x, y, z))
            + 2 * Power(r(x, y, z), 15)
                  * (M * (Power(x, 2) + Power(z, 2)) - 3 * Power(a, 2) * y * dr_dy(x, y, z))
            + Power(a, 9) * Power(z, 4) * Power(r(x, y, z), 2)
                  * (Power(a, 3)
                     - 2 * (5 * Power(a, 2) * x + Power(M, 2) * x - 4 * a * M * y) * dr_dy(x, y, z))
            + Power(a, 8) * Power(z, 2) * Power(r(x, y, z), 3)
                  * (a * M * z * (2 * a * z - x * y * dr_dz(x, y, z))
                     + (-3 * a * M * x * (Power(x, 2) + Power(y, 2))
                        + 2 * (-10 * a * M * x + 7 * Power(a, 2) * y + Power(M, 2) * y)
                              * Power(z, 2))
                           * dr_dy(x, y, z))
            - Power(a, 7) * Power(z, 2) * Power(r(x, y, z), 4)
                  * (-4 * Power(a, 2) * (M * x * y + a * Power(z, 2))
                     + M * (2 * a * Power(x, 2) + 2 * M * x * y - a * Power(y, 2)) * z
                           * dr_dz(x, y, z)
                     + (2 * Power(a, 4) * x + 12 * Power(a, 2) * x * Power(z, 2)
                        + 2 * Power(M, 2) * x * (Power(x, 2) + Power(y, 2) + 4 * Power(z, 2))
                        - a * M * y * (6 * Power(x, 2) + 5 * Power(y, 2) + 28 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     + a * M * x * Power(y, 2) * dr_dx(x, y, z))
            + a * Power(r(x, y, z), 13)
                  * (2 * a * M * (4 * Power(x, 2) - Power(y, 2) + 4 * Power(z, 2))
                     + 3 * M * x * y * z * dr_dz(x, y, z)
                     + 2
                           * (-3 * Power(a, 3) * y
                              + M * x * (Power(x, 2) + 3 * Power(y, 2) + Power(z, 2)))
                           * dr_dy(x, y, z)
                     - M * y * (-3 * Power(x, 2) + Power(y, 2)) * dr_dx(x, y, z))
            + Power(a, 2) * Power(r(x, y, z), 12)
                  * (4 * a * (Power(a, 3) + M * x * y + 2 * a * Power(z, 2))
                     - M * (2 * Power(x, 2) + Power(y, 2)) * z * dr_dz(x, y, z)
                     + (12 * Power(a, 3) * x + 2 * a * x * Power(z, 2)
                        + M * y * (-4 * Power(x, 2) + Power(y, 2) + 4 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     + M * x * (-2 * Power(x, 2) + 3 * Power(y, 2)) * dr_dx(x, y, z))
            + Power(a, 6) * Power(z, 2) * Power(r(x, y, z), 5)
                  * (2 * a * M
                         * (2 * M * x * y + a * (2 * Power(x, 2) - Power(y, 2) + 4 * Power(z, 2)))
                     + M * (a * x * y + 2 * M * (-2 * Power(x, 2) + Power(y, 2))) * z
                           * dr_dz(x, y, z)
                     + 2
                           * (-(Power(a, 3) * M * x) + 2 * Power(a, 4) * y
                              + 9 * Power(a, 2) * y * Power(z, 2)
                              + Power(M, 2) * y * (2 * Power(x, 2) + Power(y, 2) + 4 * Power(z, 2))
                              - a * M * x * (3 * Power(x, 2) + 5 * Power(y, 2) + 12 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     + M * y * (-2 * M * x * y + a * (-3 * Power(x, 2) + Power(y, 2)))
                           * dr_dx(x, y, z))
            + Power(r(x, y, z), 11)
                  * (2 * a * M
                         * (a * Power(z, 2) * (Power(x, 2) + Power(z, 2))
                            - 2 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                            + Power(a, 3) * (5 * Power(x, 2) - 2 * Power(y, 2) + 6 * Power(z, 2)))
                     - M * y * z
                           * (-5 * Power(a, 3) * x
                              + 2 * M * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dz(x, y, z)
                     + (-2 * Power(a, 6) * y + 4 * Power(a, 4) * y * Power(z, 2)
                        + 2 * Power(M, 2) * y * (Power(x, 2) + Power(z, 2))
                              * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                        + Power(a, 3) * M * x * (7 * (Power(x, 2) + Power(y, 2)) + 4 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     - 2 * Power(M, 2) * x * Power(y, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                           * dr_dx(x, y, z))
            + Power(a, 3) * Power(r(x, y, z), 8)
                  * (4 * a
                         * (a * M * x * y * (Power(a, 2) + Power(z, 2))
                            + Power(a, 2) * Power(z, 2) * (2 * Power(a, 2) + Power(z, 2))
                            + Power(M, 2) * (x - y) * (x + y)
                                  * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                     + M * z
                           * (Power(a, 3) * (-2 * Power(x, 2) + Power(y, 2))
                              - a * (2 * Power(x, 2) + Power(y, 2)) * Power(z, 2)
                              + 2 * M * x * y * (2 * (Power(x, 2) + Power(y, 2)) + 5 * Power(z, 2)))
                           * dr_dz(x, y, z)
                     + (Power(a, 6) * x - Power(a, 2) * x * Power(z, 4)
                        + Power(a, 3) * M * y * (2 * Power(x, 2) + Power(y, 2) + 12 * Power(z, 2))
                        + a * M * y * Power(z, 2)
                              * (8 * Power(x, 2) + 13 * Power(y, 2) + 20 * Power(z, 2))
                        + 2 * Power(M, 2) * x
                              * (Power(x, 4) - 2 * Power(x, 2) * Power(y, 2) - 3 * Power(y, 4)
                                 - 3 * (Power(x, 2) + Power(y, 2)) * Power(z, 2) - 4 * Power(z, 4)))
                           * dr_dy(x, y, z)
                     - M
                           * (Power(a, 3) * x * Power(y, 2)
                              + M
                                    * (6 * Power(x, 4) * y + 4 * Power(x, 2) * Power(y, 3)
                                       - 2 * Power(y, 5))
                              + a * x * (2 * Power(x, 2) - 3 * Power(y, 2)) * Power(z, 2))
                           * dr_dx(x, y, z))
            + a * Power(r(x, y, z), 10)
                  * (Power(a, 7) + 8 * Power(a, 4) * M * x * y + 12 * Power(a, 5) * Power(z, 2)
                     + Power(a, 3) * Power(z, 4)
                     + 4 * a * Power(M, 2) * (x - y) * (x + y)
                           * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                     + M * z
                           * (Power(a, 3) * (-4 * Power(x, 2) + Power(y, 2))
                              - a * Power(y, 2) * Power(z, 2)
                              + 6 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dz(x, y, z)
                     + (6 * Power(a, 6) * x + 4 * Power(a, 4) * x * Power(z, 2)
                        - 2 * Power(M, 2) * x * (Power(x, 2) - 3 * Power(y, 2) + Power(z, 2))
                              * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                        + a * M * y * Power(z, 2)
                              * (4 * Power(x, 2) + 3 * Power(y, 2) + 4 * Power(z, 2))
                        + Power(a, 3) * M * y
                              * (-2 * Power(x, 2) + 3 * Power(y, 2) + 12 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     + M
                           * (Power(a, 3) * (-2 * Power(x, 3) + 3 * x * Power(y, 2))
                              - a * x * Power(y, 2) * Power(z, 2)
                              - 2 * M * y * (-3 * Power(x, 2) + Power(y, 2))
                                    * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dx(x, y, z))
            + Power(a, 4) * Power(r(x, y, z), 7)
                  * (2 * a * M
                         * (2 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                            + a * Power(z, 2)
                                  * (Power(a, 2) + 5 * Power(x, 2) - 2 * Power(y, 2)
                                     + 6 * Power(z, 2)))
                     - M * z
                           * (a * x * y * (Power(a, 2) - 5 * Power(z, 2))
                              + 2 * M
                                    * (2 * Power(x, 4) - Power(y, 2) * (Power(y, 2) + Power(z, 2))
                                       + Power(x, 2) * (Power(y, 2) + 4 * Power(z, 2))))
                           * dr_dz(x, y, z)
                     + (M * x * (Power(a, 3) + 2 * M * x * y) * (Power(x, 2) + Power(y, 2))
                        + (-4 * Power(a, 3) * M * x + 12 * Power(a, 4) * y
                           + 10 * Power(M, 2) * Power(y, 3)
                           - 5 * a * M * x * (Power(x, 2) + Power(y, 2)))
                              * Power(z, 2)
                        + 2 * (-6 * a * M * x + 5 * Power(a, 2) * y + 6 * Power(M, 2) * y)
                              * Power(z, 4))
                           * dr_dy(x, y, z)
                     - 2 * Power(M, 2) * x
                           * (Power(y, 4) - 3 * Power(y, 2) * Power(z, 2)
                              + Power(x, 2) * (Power(y, 2) + 2 * Power(z, 2)))
                           * dr_dx(x, y, z))
            + Power(a, 2) * Power(r(x, y, z), 9)
                  * (2 * a * M
                         * (Power(a, 3) * (2 * Power(x, 2) - Power(y, 2))
                            + (4 * a * (Power(a, 2) + Power(x, 2)) - 2 * M * x * y
                               - a * Power(y, 2))
                                  * Power(z, 2)
                            + 4 * a * Power(z, 4))
                     + M * z
                           * (-4 * M * Power(x, 2) * (Power(x, 2) + Power(y, 2))
                              - 2 * M * (2 * Power(x, 2) + Power(y, 2)) * Power(z, 2)
                              + a * x * y * (Power(a, 2) + 3 * Power(z, 2)))
                           * dr_dz(x, y, z)
                     + 2
                           * (Power(a, 3) * M * x * (3 * Power(x, 2) + Power(y, 2))
                              + 6 * Power(a, 4) * y * Power(z, 2) + Power(a, 2) * y * Power(z, 4)
                              - a * M * x * Power(z, 2) * (Power(x, 2) - Power(y, 2) + Power(z, 2))
                              + Power(M, 2) * y
                                    * (-4 * Power(x, 4) - 2 * Power(x, 2) * Power(y, 2)
                                       + 2 * Power(y, 4) + 5 * Power(y, 2) * Power(z, 2)
                                       + 4 * Power(z, 4)))
                           * dr_dy(x, y, z)
                     + M
                           * (a * y * (-3 * Power(x, 2) + Power(y, 2)) * (a - z) * (a + z)
                              + M * x
                                    * (-4 * Power(x, 4) + 8 * Power(y, 4)
                                       + 6 * Power(y, 2) * Power(z, 2)
                                       + 4 * Power(x, 2) * (y - z) * (y + z)))
                           * dr_dx(x, y, z))
            + Power(r(x, y, z), 14)
                  * (2 * Power(a, 2) * (3 * Power(a, 2) + Power(z, 2) + 5 * a * x * dr_dy(x, y, z))
                     - M * Power(y, 2)
                           * (z * dr_dz(x, y, z) + y * dr_dy(x, y, z) + x * dr_dx(x, y, z)))
            + Power(a, 5) * z * Power(r(x, y, z), 6)
                  * (M
                         * (a * (-4 * Power(x, 2) + Power(y, 2)) * Power(z, 2)
                            - 2 * M * x * y * (Power(x, 2) + Power(y, 2) - Power(z, 2)))
                         * dr_dz(x, y, z)
                     + z
                           * (2 * a
                                  * (Power(a, 4) + 4 * a * M * x * y
                                     + 2 * Power(M, 2) * (x - y) * (x + y)
                                     + 3 * Power(a, 2) * Power(z, 2))
                              + (-4 * Power(a, 4) * x + 4 * Power(a, 3) * M * y
                                 - 6 * Power(a, 2) * x * Power(z, 2)
                                 + a * M * y
                                       * (10 * Power(x, 2) + 15 * Power(y, 2) + 36 * Power(z, 2))
                                 - 4 * Power(M, 2) * x
                                       * (Power(x, 2) + 3 * (Power(y, 2) + Power(z, 2))))
                                    * dr_dy(x, y, z)
                              + M
                                    * (-2 * a * Power(x, 3) - 6 * M * Power(x, 2) * y
                                       + 3 * a * x * Power(y, 2) + 2 * M * Power(y, 3))
                                    * dr_dx(x, y, z)))))
        / (Power(Power(a, 2) + Power(r(x, y, z), 2), 3)
           * (Power(a, 2) * Power(z, 2) + Power(r(x, y, z), 4))
           * Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
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
                        + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))));
  mixedK[1][2]
      = (M * r(x, y, z)
         * (3 * a * x * Power(r(x, y, z), 15) * dr_dz(x, y, z)
            + 2 * Power(a, 12) * Power(z, 5) * dr_dy(x, y, z)
            - 2 * Power(r(x, y, z), 16) * (y * dr_dz(x, y, z) + z * dr_dy(x, y, z))
            - 3 * Power(a, 10) * Power(z, 4) * r(x, y, z)
                  * (a * x * dr_dz(x, y, z) - 2 * M * z * dr_dy(x, y, z))
            + Power(a, 8) * Power(z, 3) * Power(r(x, y, z), 3)
                  * (-2 * a * (-2 * M * x + a * y) + (-7 * a * x + 8 * M * y) * z * dr_dz(x, y, z)
                     + 4 * M * (Power(x, 2) + Power(y, 2) + 6 * Power(z, 2)) * dr_dy(x, y, z))
            - 4 * Power(r(x, y, z), 14)
                  * (M * y * z + Power(a, 2) * (y * dr_dz(x, y, z) + 2 * z * dr_dy(x, y, z)))
            + 2 * Power(a, 8) * Power(z, 3) * Power(r(x, y, z), 2)
                  * (Power(a, 3) * x
                     + 2 * z
                           * (a * (-(M * x) + a * y) * dr_dz(x, y, z)
                              + (2 * Power(a, 2) + Power(M, 2)) * z * dr_dy(x, y, z)))
            + 2 * Power(a, 6) * Power(z, 3) * Power(r(x, y, z), 4)
                  * (Power(a, 2) * (3 * a * x - 4 * M * y)
                     + (-2 * a * M * x + 5 * Power(a, 2) * y + 2 * Power(M, 2) * y) * z
                           * dr_dz(x, y, z)
                     + (-(a * M * x * y) + 6 * Power(a, 2) * Power(z, 2)
                        + 2 * Power(M, 2) * (Power(x, 2) + Power(y, 2) + 4 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     + a * M * Power(y, 2) * dr_dx(x, y, z))
            + Power(a, 3) * z * Power(r(x, y, z), 7)
                  * (-2 * a
                         * (-2 * Power(a, 3) * M * x + Power(a, 4) * y - 2 * a * M * x * Power(z, 2)
                            + 3 * Power(a, 2) * y * Power(z, 2)
                            + 4 * Power(M, 2) * y * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                     + z
                           * (-2 * Power(a, 4) * x + 8 * Power(a, 3) * M * y
                              - Power(a, 2) * x * Power(z, 2)
                              + 8 * Power(M, 2) * x * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2))
                              + 4 * a * M * y * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                           * dr_dz(x, y, z)
                     - 4 * M
                           * ((Power(a, 3) + M * x * y) * (Power(x, 2) + Power(y, 2))
                              + (2 * Power(a, 3) - 3 * a * Power(x, 2) - 2 * M * x * y
                                 - 4 * a * Power(y, 2))
                                    * Power(z, 2)
                              - 6 * a * Power(z, 4))
                           * dr_dy(x, y, z)
                     + 4 * M
                           * (M * Power(y, 2) * (Power(x, 2) + Power(y, 2))
                              + x * (2 * M * x + a * y) * Power(z, 2))
                           * dr_dx(x, y, z))
            + 2 * Power(a, 2) * z * Power(r(x, y, z), 8)
                  * (a
                         * (3 * Power(a, 4) * x - 4 * Power(a, 3) * M * y
                            + Power(a, 2) * x * Power(z, 2) - 8 * a * M * y * Power(z, 2)
                            - 4 * Power(M, 2) * x * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                     + z
                           * (2 * Power(a, 3) * (-(M * x) + 2 * a * y)
                              + (2 * a * M * x + Power(a, 2) * y - 2 * Power(M, 2) * y)
                                    * Power(z, 2))
                           * dr_dz(x, y, z)
                     + (-Power(a, 6) - Power(a, 3) * M * x * y + 3 * a * M * x * y * Power(z, 2)
                        + Power(a, 2) * Power(z, 4)
                        + 2 * Power(M, 2)
                              * (-Power(x, 4) + 2 * Power(y, 4) + 5 * Power(y, 2) * Power(z, 2)
                                 + 4 * Power(z, 4) + Power(x, 2) * (Power(y, 2) + 3 * Power(z, 2))))
                           * dr_dy(x, y, z)
                     + M
                           * (y * (Power(a, 3) * y + 6 * M * x * (Power(x, 2) + Power(y, 2)))
                              + (2 * a * Power(x, 2) + 4 * M * x * y - a * Power(y, 2))
                                    * Power(z, 2))
                           * dr_dx(x, y, z))
            + 2 * z * Power(r(x, y, z), 10)
                  * (a
                         * (3 * Power(a, 4) * x - 10 * Power(a, 3) * M * y
                            - 2 * a * M * y * Power(z, 2)
                            - 4 * Power(M, 2) * x * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                     - 2 * z
                           * (-(Power(a, 3) * M * x) - Power(a, 4) * y
                              + Power(M, 2) * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dz(x, y, z)
                     + 2
                           * (-2 * Power(a, 6) + Power(a, 3) * M * x * y
                              + Power(M, 2) * (Power(x, 2) + Power(z, 2))
                                    * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dy(x, y, z)
                     - 2 * M * x
                           * (-(Power(a, 3) * x)
                              + M * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dx(x, y, z))
            + Power(a, 2) * Power(r(x, y, z), 11)
                  * (2 * a * (2 * M * x - 3 * a * y) * z
                     + (5 * Power(a, 3) * x + 2 * a * x * Power(z, 2)
                        + 4 * M * y * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                           * dr_dz(x, y, z)
                     + 4 * M * z
                           * (-((3 * Power(x, 2) + 2 * (Power(y, 2) + Power(z, 2)))
                                * dr_dy(x, y, z))
                              + x * y * dr_dx(x, y, z)))
            + Power(a, 5) * Power(z, 2) * Power(r(x, y, z), 5)
                  * ((-2 * Power(a, 4) * x - 5 * Power(a, 2) * x * Power(z, 2)
                      + 8 * Power(M, 2) * x * Power(z, 2)
                      + 2 * a * M * y * (Power(x, 2) + Power(y, 2) + 8 * Power(z, 2)))
                         * dr_dz(x, y, z)
                     + 2 * z
                           * (a * (4 * a * M * x - 3 * Power(a, 2) * y - 4 * Power(M, 2) * y)
                              + M
                                    * (-Power(a, 3) + 5 * a * Power(x, 2) - 2 * M * x * y
                                       + 8 * a * Power(y, 2) + 18 * a * Power(z, 2))
                                    * dr_dy(x, y, z)
                              + M * y * (3 * a * x + 2 * M * y) * dr_dx(x, y, z)))
            - 2 * a * Power(r(x, y, z), 12)
                  * ((Power(a, 3) * y - 2 * M * x * Power(z, 2)) * dr_dz(x, y, z)
                     + z
                           * (a * (-(a * x) + 8 * M * y)
                              + (6 * Power(a, 3) - 3 * M * x * y) * dr_dy(x, y, z)
                              + M * (-2 * Power(x, 2) + Power(y, 2)) * dr_dx(x, y, z)))
            + 2 * Power(a, 4) * z * Power(r(x, y, z), 6)
                  * (Power(a, 5) * x
                     + a * (3 * Power(a, 2) * x - 4 * Power(M, 2) * x - 10 * a * M * y)
                           * Power(z, 2)
                     + 2 * z
                           * ((-(Power(a, 3) * M * x) + Power(a, 4) * y + a * M * x * Power(z, 2)
                               + 2 * Power(a, 2) * y * Power(z, 2)
                               + Power(M, 2) * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                  * dr_dz(x, y, z)
                              + z
                                    * (a * M * x * y + 2 * Power(a, 2) * Power(z, 2)
                                       + Power(M, 2)
                                             * (2 * Power(x, 2) + 5 * Power(y, 2)
                                                + 6 * Power(z, 2)))
                                    * dr_dy(x, y, z)
                              + M * x * (a * x + 3 * M * y) * z * dr_dx(x, y, z)))
            + a * Power(r(x, y, z), 9)
                  * ((Power(a, 6) * x + 2 * Power(a, 4) * x * Power(z, 2)
                      + 2 * a * M * y * (Power(x, 2) + Power(y, 2)) * Power(z, 2)
                      + 8 * Power(M, 2) * x * Power(z, 2)
                            * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                      + 2 * Power(a, 3) * M * y * (Power(x, 2) + Power(y, 2) + 8 * Power(z, 2)))
                         * dr_dz(x, y, z)
                     + 2 * z
                           * (-(a
                                * (-4 * Power(a, 3) * M * x + 3 * Power(a, 4) * y
                                   + Power(a, 2) * y * Power(z, 2)
                                   + 4 * Power(M, 2) * y
                                         * (Power(x, 2) + Power(y, 2) + Power(z, 2))))
                              + M
                                    * (6 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                       + a * Power(z, 2)
                                             * (3 * Power(x, 2) + 2 * Power(y, 2) + 3 * Power(z, 2))
                                       - Power(a, 3)
                                             * (7 * Power(x, 2) + 4 * Power(y, 2)
                                                + 6 * Power(z, 2)))
                                    * dr_dy(x, y, z)
                              + M
                                    * (a * x * y * (3 * Power(a, 2) - Power(z, 2))
                                       + 2 * M * (2 * Power(x, 2) - Power(y, 2))
                                             * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                    * dr_dx(x, y, z)))
            + Power(r(x, y, z), 13)
                  * ((7 * Power(a, 3) * x + 2 * M * y * (Power(x, 2) + Power(y, 2)))
                         * dr_dz(x, y, z)
                     - 2 * z
                           * (M * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2)) * dr_dy(x, y, z)
                              + y * (Power(a, 2) + M * x * dr_dx(x, y, z))))))
        / (Power(Power(a, 2) + Power(r(x, y, z), 2), 2)
           * (Power(a, 2) * Power(z, 2) + Power(r(x, y, z), 4))
           * Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
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
                        + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))));
  mixedK[2][0]
      = (M * r(x, y, z)
         * (-3 * a * y * Power(r(x, y, z), 15) * dr_dz(x, y, z)
            + 2 * Power(a, 12) * Power(z, 5) * dr_dx(x, y, z)
            - 2 * Power(r(x, y, z), 16) * (x * dr_dz(x, y, z) + z * dr_dx(x, y, z))
            + Power(a, 10) * Power(z, 4) * r(x, y, z)
                  * (3 * a * y * dr_dz(x, y, z) + 2 * M * z * dr_dx(x, y, z))
            + Power(a, 8) * Power(z, 3) * Power(r(x, y, z), 3)
                  * (-2 * a * (a * x + 2 * M * y) + (8 * M * x + 7 * a * y) * z * dr_dz(x, y, z)
                     + 2 * M * x * y * dr_dy(x, y, z)
                     - 2 * M * (Power(x, 2) + 2 * Power(y, 2) - 4 * Power(z, 2)) * dr_dx(x, y, z))
            + Power(r(x, y, z), 13)
                  * (-((7 * Power(a, 3) * y
                        + 2 * M * x * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                       * dr_dz(x, y, z))
                     - 2 * x * z * (Power(a, 2) + M * y * dr_dy(x, y, z))
                     + 2 * M * z * (Power(y, 2) + Power(z, 2)) * dr_dx(x, y, z))
            - 2 * Power(a, 4) * z * Power(r(x, y, z), 6)
                  * (Power(a, 2)
                         * (a * y * (Power(a, 2) + 3 * Power(z, 2))
                            + 2 * M * x * (Power(x, 2) + Power(y, 2) + 6 * Power(z, 2)))
                     - a * z
                           * (2 * Power(a, 3) * x - Power(a, 2) * M * y + 4 * a * x * Power(z, 2)
                              + 5 * M * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dz(x, y, z)
                     - 2 * M
                           * (a * (x - y) * (x + y) * Power(z, 2)
                              + M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dy(x, y, z)
                     + 2
                           * (Power(M, 2) * Power(y, 2) * (Power(x, 2) + Power(y, 2))
                              + M * (2 * a * x * y + M * (Power(x, 2) + 2 * Power(y, 2)))
                                    * Power(z, 2)
                              - 2 * Power(a, 2) * Power(z, 4))
                           * dr_dx(x, y, z))
            + Power(a, 3) * z * Power(r(x, y, z), 7)
                  * (-2 * Power(a, 2)
                         * (Power(a, 3) * x + 4 * M * y * (Power(x, 2) + Power(y, 2))
                            + 3 * a * x * Power(z, 2) + 6 * M * y * Power(z, 2))
                     + z
                           * (2 * Power(a, 4) * y + 14 * a * M * x * (Power(x, 2) + Power(y, 2))
                              + 4 * Power(M, 2) * y * (Power(x, 2) + Power(y, 2))
                              + a * (16 * M * x + a * y) * Power(z, 2))
                           * dr_dz(x, y, z)
                     + 2 * M
                           * (a * x * y * (a - z) * (a + z)
                              + 2 * M * (x - y) * (x + y)
                                    * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     - 2 * M
                           * (a * Power(z, 2) * (Power(x, 2) - 4 * Power(z, 2))
                              + Power(a, 3) * (Power(x, 2) + 2 * Power(y, 2) - 4 * Power(z, 2))
                              + 4 * M * x * y * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                           * dr_dx(x, y, z))
            - 2 * a * Power(r(x, y, z), 12)
                  * ((Power(a, 3) * x + M * y * (2 * (Power(x, 2) + Power(y, 2)) + 3 * Power(z, 2)))
                         * dr_dz(x, y, z)
                     + z
                           * (a * (6 * M * x + a * y)
                              + M * (-Power(x, 2) + Power(y, 2)) * dr_dy(x, y, z)
                              + 2 * (3 * Power(a, 3) + M * x * y) * dr_dx(x, y, z)))
            - 4 * Power(r(x, y, z), 14)
                  * (M * x * z + Power(a, 2) * (x * dr_dz(x, y, z) + 2 * z * dr_dx(x, y, z)))
            + 2 * Power(a, 9) * Power(z, 3) * Power(r(x, y, z), 2)
                  * ((2 * a * x + 3 * M * y) * z * dr_dz(x, y, z)
                     + a * (-(a * y) + 4 * Power(z, 2) * dr_dx(x, y, z)))
            - Power(a, 2) * Power(r(x, y, z), 11)
                  * ((5 * Power(a, 3) * y + 2 * a * y * Power(z, 2)
                      + 2 * M * x * (Power(x, 2) + Power(y, 2) + 4 * Power(z, 2)))
                         * dr_dz(x, y, z)
                     + 2 * z
                           * (3 * Power(a, 2) * x + M * x * y * dr_dy(x, y, z)
                              + M * (Power(x, 2) - 4 * Power(z, 2)) * dr_dx(x, y, z)))
            + 2 * Power(a, 6) * Power(z, 2) * Power(r(x, y, z), 4)
                  * (a
                         * (3 * M * y * (Power(x, 2) + Power(y, 2))
                            + (5 * a * x + 7 * M * y) * Power(z, 2))
                         * dr_dz(x, y, z)
                     + z
                           * (-(Power(a, 2) * (4 * M * x + 3 * a * y))
                              + M * (2 * M * x * y + a * (x - y) * (x + y)) * dr_dy(x, y, z)
                              + (-2 * M * y * (a * x + M * y) + 6 * Power(a, 2) * Power(z, 2))
                                    * dr_dx(x, y, z)))
            - 2 * Power(r(x, y, z), 10)
                  * (-((2 * Power(a, 4) * x * Power(z, 2)
                        + 2 * Power(M, 2) * x * (Power(x, 2) + Power(y, 2))
                              * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                        - Power(a, 3) * M * y * (3 * (Power(x, 2) + Power(y, 2)) + 7 * Power(z, 2)))
                       * dr_dz(x, y, z))
                     + z
                           * (6 * Power(a, 4) * M * x + 3 * Power(a, 5) * y
                              + 2 * Power(a, 2) * M * x
                                    * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2))
                              + 2 * M
                                    * (Power(a, 3) * (-Power(x, 2) + Power(y, 2))
                                       + M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                    * dr_dy(x, y, z)
                              + 2
                                    * (2 * Power(a, 6) + 2 * Power(a, 3) * M * x * y
                                       + Power(M, 2) * Power(x, 2)
                                             * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                    * dr_dx(x, y, z)))
            - a * Power(r(x, y, z), 9)
                  * ((Power(a, 6) * y + 4 * Power(a, 3) * M * x * Power(z, 2)
                      + 2 * Power(a, 4) * y * Power(z, 2)
                      - 4 * Power(M, 2) * y * (Power(x, 2) + Power(y, 2))
                            * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                      - 2 * a * M * x * Power(z, 2)
                            * (3 * (Power(x, 2) + Power(y, 2)) + 2 * Power(z, 2)))
                         * dr_dz(x, y, z)
                     + 2 * z
                           * (Power(a, 2)
                                  * (3 * Power(a, 3) * x + a * x * Power(z, 2)
                                     + 2 * M * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                              + M
                                    * (a * x * y * (-Power(a, 2) + Power(z, 2))
                                       - 2 * M * (x - y) * (x + y)
                                             * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                    * dr_dy(x, y, z)
                              + M
                                    * (Power(a, 3)
                                           * (2 * Power(x, 2) + 3 * Power(y, 2) - 6 * Power(z, 2))
                                       - a * Power(z, 2) * (Power(y, 2) + Power(z, 2))
                                       + 4 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                    * dr_dx(x, y, z)))
            + Power(a, 5) * z * Power(r(x, y, z), 5)
                  * (-6 * Power(a, 3) * x * Power(z, 2)
                     - 4 * Power(a, 2) * M * y * (Power(x, 2) + Power(y, 2) + 3 * Power(z, 2))
                     + a * z
                           * (2 * Power(a, 3) * y + 8 * M * x * (Power(x, 2) + Power(y, 2))
                              + 5 * (4 * M * x + a * y) * Power(z, 2))
                           * dr_dz(x, y, z)
                     + 2 * M * Power(z, 2)
                           * ((a * x * y + 2 * M * (x - y) * (x + y)) * dr_dy(x, y, z)
                              + (Power(a, 3) - 4 * M * x * y
                                 + a * (-2 * Power(x, 2) - 3 * Power(y, 2) + 6 * Power(z, 2)))
                                    * dr_dx(x, y, z)))
            - 2 * Power(a, 2) * Power(r(x, y, z), 8)
                  * (-((-(Power(a, 3) * M * y * (Power(x, 2) + Power(y, 2)))
                        + (4 * Power(a, 4) * x - 5 * Power(a, 3) * M * y
                           + 2 * Power(M, 2) * x * (Power(x, 2) + Power(y, 2))
                           + 2 * a * M * y * (Power(x, 2) + Power(y, 2)))
                              * Power(z, 2)
                        + a * (a * x + M * y) * Power(z, 4))
                       * dr_dz(x, y, z))
                     + z
                           * (Power(a, 2)
                                  * (2 * Power(a, 2) * M * x + 3 * Power(a, 3) * y
                                     + a * y * Power(z, 2)
                                     + 4 * M * x * (Power(x, 2) + Power(y, 2) + 3 * Power(z, 2)))
                              + M
                                    * (Power(a, 3) * (-Power(x, 2) + Power(y, 2))
                                       + (2 * M * x * y + a * (-Power(x, 2) + Power(y, 2)))
                                             * Power(z, 2))
                                    * dr_dy(x, y, z)
                              + (Power(a, 6) + 2 * Power(a, 3) * M * x * y
                                 + 2 * a * M * x * y * Power(z, 2) - Power(a, 2) * Power(z, 4)
                                 + 2 * Power(M, 2)
                                       * (Power(Power(x, 2) + Power(y, 2), 2)
                                          + (2 * Power(x, 2) + Power(y, 2)) * Power(z, 2)))
                                    * dr_dx(x, y, z)))))
        / (Power(Power(a, 2) + Power(r(x, y, z), 2), 2)
           * (Power(a, 2) * Power(z, 2) + Power(r(x, y, z), 4))
           * Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
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
                        + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))));
  mixedK[2][1]
      = (M * r(x, y, z)
         * (3 * a * x * Power(r(x, y, z), 15) * dr_dz(x, y, z)
            + 2 * Power(a, 12) * Power(z, 5) * dr_dy(x, y, z)
            - 2 * Power(r(x, y, z), 16) * (y * dr_dz(x, y, z) + z * dr_dy(x, y, z))
            + Power(a, 10) * Power(z, 4) * r(x, y, z)
                  * (-3 * a * x * dr_dz(x, y, z) + 2 * M * z * dr_dy(x, y, z))
            - 4 * Power(r(x, y, z), 14)
                  * (M * y * z + Power(a, 2) * (y * dr_dz(x, y, z) + 2 * z * dr_dy(x, y, z)))
            + 2 * Power(a, 9) * Power(z, 3) * Power(r(x, y, z), 2)
                  * ((-3 * M * x + 2 * a * y) * z * dr_dz(x, y, z)
                     + a * (a * x + 4 * Power(z, 2) * dr_dy(x, y, z)))
            + Power(a, 8) * Power(z, 3) * Power(r(x, y, z), 3)
                  * (4 * a * M * x - 2 * Power(a, 2) * y
                     + (-7 * a * x + 8 * M * y) * z * dr_dz(x, y, z)
                     - 2 * M * (2 * Power(x, 2) + Power(y, 2) - 4 * Power(z, 2)) * dr_dy(x, y, z)
                     + 2 * M * x * y * dr_dx(x, y, z))
            + 2 * Power(a, 4) * z * Power(r(x, y, z), 6)
                  * (Power(a, 2)
                         * (Power(a, 3) * x + 3 * a * x * Power(z, 2)
                            - 2 * M * y * (Power(x, 2) + Power(y, 2) + 6 * Power(z, 2)))
                     + a * z
                           * (Power(a, 2) * M * x + 2 * Power(a, 3) * y + 4 * a * y * Power(z, 2)
                              - 5 * M * x * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dz(x, y, z)
                     - 2
                           * (Power(M, 2) * Power(x, 2) * (Power(x, 2) + Power(y, 2))
                              + M * (-2 * a * x * y + M * (2 * Power(x, 2) + Power(y, 2)))
                                    * Power(z, 2)
                              - 2 * Power(a, 2) * Power(z, 4))
                           * dr_dy(x, y, z)
                     + 2 * M
                           * (a * (x - y) * (x + y) * Power(z, 2)
                              + M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                           * dr_dx(x, y, z))
            + Power(a, 2) * Power(r(x, y, z), 11)
                  * ((5 * Power(a, 3) * x + 2 * a * x * Power(z, 2)
                      - 2 * M * y * (Power(x, 2) + Power(y, 2) + 4 * Power(z, 2)))
                         * dr_dz(x, y, z)
                     - 2 * z
                           * (3 * Power(a, 2) * y
                              + M * (Power(y, 2) - 4 * Power(z, 2)) * dr_dy(x, y, z)
                              + M * x * y * dr_dx(x, y, z)))
            + 2 * a * Power(r(x, y, z), 12)
                  * ((-(Power(a, 3) * y)
                      + M * x * (2 * (Power(x, 2) + Power(y, 2)) + 3 * Power(z, 2)))
                         * dr_dz(x, y, z)
                     + z
                           * (a * (a * x - 6 * M * y)
                              + (-6 * Power(a, 3) + 2 * M * x * y) * dr_dy(x, y, z)
                              + M * (x - y) * (x + y) * dr_dx(x, y, z)))
            + 2 * Power(a, 6) * Power(z, 2) * Power(r(x, y, z), 4)
                  * (a
                         * (-3 * M * x * (Power(x, 2) + Power(y, 2))
                            + (-7 * M * x + 5 * a * y) * Power(z, 2))
                         * dr_dz(x, y, z)
                     + z
                           * (Power(a, 2) * (3 * a * x - 4 * M * y)
                              + (2 * M * x * (-(M * x) + a * y) + 6 * Power(a, 2) * Power(z, 2))
                                    * dr_dy(x, y, z)
                              + M * (2 * M * x * y + a * (x - y) * (x + y)) * dr_dx(x, y, z)))
            + Power(a, 5) * z * Power(r(x, y, z), 5)
                  * (2 * Power(a, 2)
                         * (-3 * a * y * Power(z, 2)
                            + 2 * M * x * (Power(x, 2) + Power(y, 2) + 3 * Power(z, 2)))
                     + a * z
                           * (-2 * Power(a, 3) * x - 5 * a * x * Power(z, 2)
                              + 4 * M * y * (2 * (Power(x, 2) + Power(y, 2)) + 5 * Power(z, 2)))
                           * dr_dz(x, y, z)
                     + 2 * M * Power(z, 2)
                           * ((Power(a, 3) + 4 * M * x * y
                               + a * (-3 * Power(x, 2) - 2 * Power(y, 2) + 6 * Power(z, 2)))
                                  * dr_dy(x, y, z)
                              + (a * x * y + 2 * M * (x - y) * (x + y)) * dr_dx(x, y, z)))
            - 2 * Power(a, 2) * Power(r(x, y, z), 8)
                  * (-((Power(a, 3) * M * x * (Power(x, 2) + Power(y, 2))
                        + (5 * Power(a, 3) * M * x + 4 * Power(a, 4) * y
                           - 2 * a * M * x * (Power(x, 2) + Power(y, 2))
                           + 2 * Power(M, 2) * y * (Power(x, 2) + Power(y, 2)))
                              * Power(z, 2)
                        + a * (-(M * x) + a * y) * Power(z, 4))
                       * dr_dz(x, y, z))
                     + z
                           * (Power(a, 2)
                                  * (-3 * Power(a, 3) * x + 2 * Power(a, 2) * M * y
                                     - a * x * Power(z, 2)
                                     + 4 * M * y * (Power(x, 2) + Power(y, 2) + 3 * Power(z, 2)))
                              + (Power(a, 6) - 2 * Power(a, 3) * M * x * y
                                 - 2 * a * M * x * y * Power(z, 2) - Power(a, 2) * Power(z, 4)
                                 + 2 * Power(M, 2)
                                       * (Power(Power(x, 2) + Power(y, 2), 2)
                                          + (Power(x, 2) + 2 * Power(y, 2)) * Power(z, 2)))
                                    * dr_dy(x, y, z)
                              + M
                                    * (Power(a, 3) * (-Power(x, 2) + Power(y, 2))
                                       + (2 * M * x * y + a * (-Power(x, 2) + Power(y, 2)))
                                             * Power(z, 2))
                                    * dr_dx(x, y, z)))
            + 2 * Power(r(x, y, z), 10)
                  * ((2 * Power(a, 4) * y * Power(z, 2)
                      + 2 * Power(M, 2) * y * (Power(x, 2) + Power(y, 2))
                            * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                      + Power(a, 3) * M * x * (3 * (Power(x, 2) + Power(y, 2)) + 7 * Power(z, 2)))
                         * dr_dz(x, y, z)
                     + z
                           * (3 * Power(a, 5) * x - 6 * Power(a, 4) * M * y
                              - 2 * Power(a, 2) * M * y
                                    * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2))
                              - 2
                                    * (2 * Power(a, 6) - 2 * Power(a, 3) * M * x * y
                                       + Power(M, 2) * Power(y, 2)
                                             * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                    * dr_dy(x, y, z)
                              - 2 * M
                                    * (Power(a, 3) * (-Power(x, 2) + Power(y, 2))
                                       + M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                    * dr_dx(x, y, z)))
            + a * Power(r(x, y, z), 9)
                  * ((Power(a, 6) * x + 2 * Power(a, 4) * x * Power(z, 2)
                      - 4 * Power(a, 3) * M * y * Power(z, 2)
                      - 4 * Power(M, 2) * x * (Power(x, 2) + Power(y, 2))
                            * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                      + 2 * a * M * y * Power(z, 2)
                            * (3 * (Power(x, 2) + Power(y, 2)) + 2 * Power(z, 2)))
                         * dr_dz(x, y, z)
                     + 2 * z
                           * (Power(a, 2)
                                  * (-(a * y * (3 * Power(a, 2) + Power(z, 2)))
                                     + 2 * M * x * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                              + M
                                    * (a * Power(z, 2) * (Power(x, 2) + Power(z, 2))
                                       + 4 * M * x * y * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                       + Power(a, 3)
                                             * (-3 * Power(x, 2) - 2 * Power(y, 2)
                                                + 6 * Power(z, 2)))
                                    * dr_dy(x, y, z)
                              + M
                                    * (a * x * y * (a - z) * (a + z)
                                       + 2 * M * (x - y) * (x + y)
                                             * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                    * dr_dx(x, y, z)))
            + Power(a, 3) * z * Power(r(x, y, z), 7)
                  * (-(z
                       * (2 * Power(a, 4) * x + 4 * Power(M, 2) * x * (Power(x, 2) + Power(y, 2))
                          + Power(a, 2) * x * Power(z, 2)
                          - 2 * a * M * y * (7 * (Power(x, 2) + Power(y, 2)) + 8 * Power(z, 2)))
                       * dr_dz(x, y, z))
                     + 2
                           * (Power(a, 2)
                                  * (4 * M * x * (Power(x, 2) + Power(y, 2))
                                     + 6 * M * x * Power(z, 2)
                                     - a * y * (Power(a, 2) + 3 * Power(z, 2)))
                              + M
                                    * (-(Power(a, 3)
                                         * (2 * Power(x, 2) + Power(y, 2) - 4 * Power(z, 2)))
                                       + 4 * M * x * y
                                             * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2))
                                       + a * (-(Power(y, 2) * Power(z, 2)) + 4 * Power(z, 4)))
                                    * dr_dy(x, y, z)
                              + M
                                    * (a * x * y * (a - z) * (a + z)
                                       + 2 * M * (x - y) * (x + y)
                                             * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                                    * dr_dx(x, y, z)))
            + Power(r(x, y, z), 13)
                  * ((7 * Power(a, 3) * x
                      - 2 * M * y * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                         * dr_dz(x, y, z)
                     + 2 * z
                           * (M * (Power(x, 2) + Power(z, 2)) * dr_dy(x, y, z)
                              - y * (Power(a, 2) + M * x * dr_dx(x, y, z))))))
        / (Power(Power(a, 2) + Power(r(x, y, z), 2), 2)
           * (Power(a, 2) * Power(z, 2) + Power(r(x, y, z), 4))
           * Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
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
                        + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))));
  mixedK[2][2]
      = (2 * M * r(x, y, z)
         * (3 * Power(a, 2) * Power(r(x, y, z), 13) + Power(r(x, y, z), 15)
            + 2 * Power(a, 10) * Power(z, 5) * dr_dz(x, y, z)
            - 2 * z * Power(r(x, y, z), 14) * dr_dz(x, y, z)
            + 2 * Power(r(x, y, z), 12)
                  * (M * (Power(x, 2) + Power(y, 2)) - 3 * Power(a, 2) * z * dr_dz(x, y, z))
            - Power(a, 8) * Power(z, 4) * r(x, y, z) * (Power(a, 2) - 3 * M * z * dr_dz(x, y, z))
            + Power(a, 3) * Power(z, 2) * Power(r(x, y, z), 5)
                  * (-3 * Power(a, 3) * Power(z, 2)
                     + a * M * z
                           * (-Power(a, 2) + 5 * (Power(x, 2) + Power(y, 2)) + 9 * Power(z, 2))
                           * dr_dz(x, y, z)
                     + 2 * M
                           * (-(a * y * Power(z, 2))
                              + M * x * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                           * dr_dy(x, y, z)
                     - 2 * M
                           * (a * x * Power(z, 2)
                              + M * y * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                           * dr_dx(x, y, z))
            + Power(r(x, y, z), 11)
                  * (3 * Power(a, 4)
                     - M * Power(z, 2)
                           * (z * dr_dz(x, y, z) + y * dr_dy(x, y, z) + x * dr_dx(x, y, z)))
            + a * Power(r(x, y, z), 10)
                  * (4 * a * M * (Power(x, 2) + Power(y, 2)) - 6 * Power(a, 3) * z * dr_dz(x, y, z)
                     + M * Power(z, 2) * (x * dr_dy(x, y, z) - y * dr_dx(x, y, z)))
            + Power(a, 7) * Power(z, 4) * Power(r(x, y, z), 2)
                  * (6 * a * z * dr_dz(x, y, z)
                     + M * (-2 * a + x * dr_dy(x, y, z) - y * dr_dx(x, y, z)))
            - 2 * Power(a, 4) * Power(z, 4) * Power(r(x, y, z), 4)
                  * (-3 * Power(a, 2) * z * dr_dz(x, y, z)
                     + M
                           * (3 * Power(a, 2) + (-(a * x) + M * y) * dr_dy(x, y, z)
                              + (M * x + a * y) * dr_dx(x, y, z)))
            + Power(a, 5) * Power(z, 3) * Power(r(x, y, z), 3)
                  * (a * M * (Power(x, 2) + Power(y, 2) + 9 * Power(z, 2)) * dr_dz(x, y, z)
                     + z
                           * (-3 * Power(a, 3) + M * (2 * M * x - a * y) * dr_dy(x, y, z)
                              - M * (a * x + 2 * M * y) * dr_dx(x, y, z)))
            - a * z * Power(r(x, y, z), 7)
                  * (Power(a, 3) * Power(z, 3)
                     + M
                           * (a
                                  * (-4 * (Power(x, 2) + Power(y, 2)) * Power(z, 2)
                                     - 3 * Power(z, 4)
                                     + 3 * Power(a, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                  * dr_dz(x, y, z)
                              + z
                                    * (a * y * (Power(a, 2) + Power(z, 2))
                                       - 2 * M * x * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                    * dr_dy(x, y, z)
                              + z
                                    * (Power(a, 3) * x + a * x * Power(z, 2)
                                       + 2 * M * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                    * dr_dx(x, y, z)))
            + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 6)
                  * (2 * z * (Power(M, 2) * (Power(x, 2) + Power(y, 2)) + Power(a, 2) * Power(z, 2))
                         * dr_dz(x, y, z)
                     + M
                           * (-6 * Power(a, 2) * Power(z, 2)
                              + (Power(a, 3) * x + a * x * Power(z, 2)
                                 - 2 * M * y * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                                    * dr_dy(x, y, z)
                              - (a * y * (Power(a, 2) + Power(z, 2))
                                 + 2 * M * x * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2)))
                                    * dr_dx(x, y, z)))
            + Power(a, 2) * Power(r(x, y, z), 9)
                  * (Power(a, 4)
                     - M * z
                           * (3 * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * dr_dz(x, y, z)
                              + 2 * z * (y * dr_dy(x, y, z) + x * dr_dx(x, y, z))))
            - 2 * Power(r(x, y, z), 8)
                  * (Power(a, 2) * M * (-(Power(a, 2) * (Power(x, 2) + Power(y, 2))) + Power(z, 4))
                     + z
                           * ((Power(a, 6)
                               - Power(M, 2) * (Power(x, 2) + Power(y, 2))
                                     * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                  * dr_dz(x, y, z)
                              + M * z
                                    * ((-(Power(a, 3) * x)
                                        + M * y * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                           * dr_dy(x, y, z)
                                       + (Power(a, 3) * y
                                          + M * x * (Power(x, 2) + Power(y, 2) + Power(z, 2)))
                                             * dr_dx(x, y, z))))))
        / ((Power(a, 2) + Power(r(x, y, z), 2)) * (Power(a, 2) * Power(z, 2) + Power(r(x, y, z), 4))
           * Power(Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
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
                        + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))));

  return mixedK;
}