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

  auto r = r_KS(a, x, y, z);
  auto dr_dx = d_r_KS_dx(a, x, y, z);
  auto dr_dy = d_r_KS_dy(a, x, y, z);
  auto dr_dz = d_r_KS_dz(a, x, y, z);

  metric_server::spatial_matrix mixedK{};

  mixedK[0][0]
      = (2 * M * Power(r, 2)
         * ((M * Power(r, 3)
             * (Power(a, 2) * x * y + a * (Power(x, 2) - Power(y, 2)) * r - x * y * Power(r, 2))
             * (-4 * Power(r, 4) * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                    * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                + 3 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                - 2 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 3)
                      * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dy
                + r * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (a + x * dr_dy)
                - (2 * M * Power(r, 3) * Power(a * y + x * r, 2)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (3 * Power(a, 5) * y * Power(z, 2) * dr_dy
                      - Power(a, 3) * y * Power(z, 2) * Power(r, 2) * dr_dy
                      - Power(a, 3) * y * Power(r, 4) * dr_dy - 5 * a * y * Power(r, 6) * dr_dy
                      + Power(r, 7) * (2 * a - 3 * x * dr_dy)
                      + Power(a, 2) * Power(z, 2) * Power(r, 3) * (2 * a + x * dr_dy)
                      + Power(a, 2) * Power(r, 5) * (2 * a + x * dr_dy)
                      + Power(a, 4) * Power(z, 2) * r * (2 * a + 5 * x * dr_dy)))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                - 4 * Power(r, 4) * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                + 3 * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                + 2 * (a * x - y * r) * Power(Power(a, 2) + Power(r, 2), 3)
                      * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dx
                + r * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (-a + y * dr_dx)
                + (2 * M * Power(r, 3) * Power(a * x - y * r, 2)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (3 * Power(a, 5) * x * Power(z, 2) * dr_dx
                      - Power(a, 3) * x * Power(z, 2) * Power(r, 2) * dr_dx
                      - Power(a, 3) * x * Power(r, 4) * dr_dx - 5 * a * x * Power(r, 6) * dr_dx
                      + Power(a, 4) * Power(z, 2) * r * (2 * a - 5 * y * dr_dx)
                      + Power(a, 2) * Power(z, 2) * Power(r, 3) * (2 * a - y * dr_dx)
                      + Power(a, 2) * Power(r, 5) * (2 * a - y * dr_dx)
                      + Power(r, 7) * (2 * a + 3 * y * dr_dx)))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                + (2 * M * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                   * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                      + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                      + Power(a, 2) * Power(r, 4) + Power(r, 6))
                   * (2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy - x * dr_dx)
                      + 2 * Power(a, 7) * Power(z, 3) * (-(y * dr_dy) + x * dr_dx)
                      - 3 * Power(a, 6) * Power(z, 2) * r
                            * (x * y * dr_dz + x * z * dr_dy + y * z * dr_dx)
                      + 2 * Power(a, 5) * z * Power(r, 2)
                            * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                               - y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                      + 2 * Power(a, 4) * z * Power(r, 3)
                            * (a * Power(x, 2) - a * Power(y, 2) + 3 * x * y * z * dr_dz
                               - 2 * x * Power(z, 2) * dr_dy - 2 * y * Power(z, 2) * dr_dx)
                      - 2 * Power(a, 2) * z * Power(r, 6)
                            * (-3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                      + 4 * a * Power(r, 8)
                            * ((Power(x, 2) - Power(y, 2)) * dr_dz + z * (y * dr_dy - x * dr_dx))
                      + 2 * Power(a, 2) * Power(r, 7)
                            * (3 * x * y * dr_dz + 2 * z * (x * dr_dy + y * dr_dx))
                      + Power(r, 9) * (-3 * x * y * dr_dz + 3 * z * (x * dr_dy + y * dr_dx))
                      + Power(a, 2) * Power(r, 5)
                            * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                               + z
                                     * (2 * a * (Power(x, 2) - Power(y, 2))
                                        + x * (Power(a, 2) - Power(z, 2)) * dr_dy
                                        + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))))
                / Power(Power(a, 2) + Power(r, 2), 6)
            + ((Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
                + Power(r, 8))
               * (-4 * Power(r, 4) * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                  + 3 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                  - 2 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 3)
                        * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dx
                  + r * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (r + x * dr_dx)
                  - (M * Power(r, 3) * Power(a * y + x * r, 2)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (2 * Power(r, 8) + 3 * Power(a, 5) * y * Power(z, 2) * dr_dx
                        + 5 * Power(a, 4) * x * Power(z, 2) * r * dr_dx
                        + Power(a, 2) * x * Power(z, 2) * Power(r, 3) * dr_dx
                        + Power(a, 2) * x * Power(r, 5) * dr_dx - 3 * x * Power(r, 7) * dr_dx
                        + a * Power(r, 6) * (2 * a - 5 * y * dr_dx)
                        + Power(a, 3) * Power(z, 2) * Power(r, 2) * (2 * a - y * dr_dx)
                        + Power(a, 2) * Power(r, 4) * (2 * Power(z, 2) - a * y * dr_dx)))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  + (M * Power(r, 3) * (-(a * x) + y * r)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (2 * Power(a, 4) * y * Power(z, 2) * Power(r, 3)
                        + (4 * Power(a, 4) * y - 2 * Power(a, 2) * y * Power(z, 2)) * Power(r, 5)
                        - 2 * y * Power(r, 9)
                        + 3 * Power(a, 6) * y * Power(z, 2) * (y * dr_dy + 2 * x * dr_dx)
                        + 2 * a * Power(r, 7)
                              * (a * y - 4 * x * y * dr_dy
                                 - 4 * (Power(x, 2) - Power(y, 2)) * dr_dx)
                        + 4 * Power(a, 5) * Power(z, 2) * r
                              * (a * y + 2 * x * y * dr_dy
                                 + 2 * (Power(x, 2) - Power(y, 2)) * dr_dx)
                        + Power(a, 2) * Power(r, 4)
                              * (6 * a * x * Power(z, 2)
                                 + (-(Power(a, 2) * Power(y, 2)) + Power(x, 2) * Power(z, 2))
                                       * dr_dy
                                 - 2 * x * y * (Power(a, 2) + Power(z, 2)) * dr_dx)
                        + Power(a, 2) * Power(r, 6)
                              * ((Power(x, 2) - 5 * Power(y, 2)) * dr_dy
                                 + 6 * x * (a - 2 * y * dr_dx))
                        + Power(a, 4) * Power(z, 2) * Power(r, 2)
                              * ((5 * Power(x, 2) - Power(y, 2)) * dr_dy
                                 + 6 * x * (a - 2 * y * dr_dx))
                        + 3 * x * Power(r, 8) * (-(x * dr_dy) + 2 * (a + y * dr_dx))))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  - (M * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                     * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                        + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                        + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     * (2 * z * Power(r, 10) + 4 * Power(a, 7) * y * Power(z, 3) * dr_dx
                        + 2 * Power(a, 2) * z * Power(r, 6)
                              * (Power(a, 2) + Power(x, 2) + Power(z, 2) - 6 * a * y * dr_dx)
                        + 2 * Power(a, 4) * z * Power(r, 4)
                              * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2) - 2 * a * y * dr_dx)
                        + 3 * x * Power(r, 9) * (x * dr_dz - 2 * z * dr_dx)
                        + 3 * Power(a, 6) * Power(z, 2) * r
                              * (-(Power(y, 2) * dr_dz) + 2 * x * z * dr_dx)
                        - Power(a, 2) * Power(r, 7)
                              * ((Power(x, 2) - 5 * Power(y, 2)) * dr_dz + 8 * x * z * dr_dx)
                        + 2 * Power(a, 5) * z * Power(r, 2)
                              * (a * (Power(y, 2) + Power(z, 2)) - 4 * x * y * z * dr_dz
                                 + 2 * y * Power(z, 2) * dr_dx)
                        + 4 * a * Power(r, 8) * (2 * x * y * dr_dz + z * (a - 2 * y * dr_dx))
                        + Power(a, 4) * z * Power(r, 3)
                              * ((-5 * Power(x, 2) + Power(y, 2)) * z * dr_dz
                                 + 4 * x * (a * y + 2 * Power(z, 2) * dr_dx))
                        + Power(a, 2) * Power(r, 5)
                              * ((Power(a, 2) * Power(y, 2) - Power(x, 2) * Power(z, 2)) * dr_dz
                                 + 2 * x * z * (2 * a * y + (-Power(a, 2) + Power(z, 2)) * dr_dx))))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))))
                  / Power(Power(a, 2) + Power(r, 2), 6)
            - M * z * r * (a * y + x * r)
                  * ((M * Power(r, 4) * Power(a * y + x * r, 2)
                      * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                         + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                         + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                         + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                         + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                      * (-3 * Power(a, 5) * y * Power(z, 2) * dr_dz
                         - Power(a, 2) * x * Power(r, 5) * dr_dz + 5 * a * y * Power(r, 6) * dr_dz
                         + 3 * x * Power(r, 7) * dr_dz
                         + Power(a, 2) * Power(r, 4) * (2 * x * z + a * y * dr_dz)
                         + Power(a, 4) * z * r * (2 * a * y - 5 * x * z * dr_dz)
                         + Power(a, 2) * z * Power(r, 3) * (2 * a * y - x * z * dr_dz)
                         + Power(a, 3) * z * Power(r, 2) * (2 * a * x + y * z * dr_dz)))
                         / (Power(Power(a, 2) + Power(r, 2), 5)
                            * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                               + Power(a, 2) * Power(z, 2) * Power(r, 2)
                               + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                               + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     - 4 * z * Power(r, 4) * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                     + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                     - (M * Power(z, 3) * r * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                        * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                        * dr_dx)
                           / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                              + Power(a, 2) * Power(z, 2) * Power(r, 2)
                              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                              + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     + (M * Power(r, 3) * (-(a * x) + y * r)
                        * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                           + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                           + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                           + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                           + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                           + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                        * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                           - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                           + 2 * Power(a, 5) * z * Power(r, 2)
                                 * (-(a * x * y) + 2 * (Power(x, 2) - Power(y, 2)) * z * dr_dz
                                    + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                           + 2 * Power(a, 4) * z * Power(r, 3)
                                 * (-(a * Power(x, 2)) + a * Power(y, 2) + 2 * a * Power(z, 2)
                                    - 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                    - 2 * y * Power(z, 2) * dr_dx)
                           + 2 * Power(a, 2) * z * Power(r, 6)
                                 * (-3 * a * y * dr_dy + x * (y - 3 * a * dr_dx))
                           - 4 * a * Power(r, 8)
                                 * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                    + z * (y * dr_dy + x * dr_dx))
                           + Power(a, 6) * Power(z, 2) * r
                                 * (3 * x * y * dr_dz + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                           + 2 * Power(a, 2) * Power(r, 7)
                                 * (-3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                           + Power(r, 9)
                                 * (3 * x * y * dr_dz + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                           + Power(a, 2) * Power(r, 5)
                                 * (-(x * y * (Power(a, 2) + Power(z, 2)) * dr_dz)
                                    + z
                                          * (2 * a
                                                 * (Power(a, 2) - Power(x, 2) + Power(y, 2)
                                                    + Power(z, 2))
                                             + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                             + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                           / (Power(Power(a, 2) + Power(r, 2), 5)
                              * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                                 + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                 + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                                 + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     + r
                           * ((x * r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                  / (Power(a, 2) + Power(r, 2))
                              + (3 * (a * y + x * r)
                                 * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                    / (Power(a, 2) + Power(r, 2))
                              - (2 * (a * y + x * r)
                                 * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dz)
                                    / Power(Power(a, 2) + Power(r, 2), 2)
                              - (2 * r * (a * y + x * r) * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                                 * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz))
                                    / (Power(a, 2) + Power(r, 2))
                              + (M * Power(r, 3) * Power(a * y + x * r, 2)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (-3 * Power(a, 5) * y * Power(z, 2) * dr_dz
                                    - Power(a, 2) * x * Power(r, 5) * dr_dz
                                    + 5 * a * y * Power(r, 6) * dr_dz + 3 * x * Power(r, 7) * dr_dz
                                    + Power(a, 2) * Power(r, 4) * (2 * x * z + a * y * dr_dz)
                                    + Power(a, 4) * z * r * (2 * a * y - 5 * x * z * dr_dz)
                                    + Power(a, 2) * z * Power(r, 3) * (2 * a * y - x * z * dr_dz)
                                    + Power(a, 3) * z * Power(r, 2) * (2 * a * x + y * z * dr_dz)))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                              - (M * Power(z, 3) * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                                 * (Power(a, 4) * Power(z, 2)
                                    + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                    + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                                    + Power(a, 2) * Power(r, 4) + Power(r, 6))
                                 * dr_dx)
                                    / (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6))
                              + (M * Power(r, 2) * (-(a * x) + y * r)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                                    - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                                    + 2 * Power(a, 5) * z * Power(r, 2)
                                          * (-(a * x * y)
                                             + 2 * (Power(x, 2) - Power(y, 2)) * z * dr_dz
                                             + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                                    + 2 * Power(a, 4) * z * Power(r, 3)
                                          * (-(a * Power(x, 2)) + a * Power(y, 2)
                                             + 2 * a * Power(z, 2) - 3 * x * y * z * dr_dz
                                             + 2 * x * Power(z, 2) * dr_dy
                                             - 2 * y * Power(z, 2) * dr_dx)
                                    + 2 * Power(a, 2) * z * Power(r, 6)
                                          * (-3 * a * y * dr_dy + x * (y - 3 * a * dr_dx))
                                    - 4 * a * Power(r, 8)
                                          * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                             + z * (y * dr_dy + x * dr_dx))
                                    + Power(a, 6) * Power(z, 2) * r
                                          * (3 * x * y * dr_dz
                                             + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                                    + 2 * Power(a, 2) * Power(r, 7)
                                          * (-3 * x * y * dr_dz
                                             + 2 * z * (a - x * dr_dy + y * dr_dx))
                                    + Power(r, 9)
                                          * (3 * x * y * dr_dz
                                             + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                                    + Power(a, 2) * Power(r, 5)
                                          * (-(x * y * (Power(a, 2) + Power(z, 2)) * dr_dz)
                                             + z
                                                   * (2 * a
                                                          * (Power(a, 2) - Power(x, 2) + Power(y, 2)
                                                             + Power(z, 2))
                                                      + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                                      + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))))))
        / (Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 3)
           * Sqrt(
               (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                + Power(a, 2) * Power(z, 2) * Power(r, 2)
                + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                + Power(a, 2) * Power(r, 4) + Power(r, 6))
               * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                  + Power(a, 2) * Power(z, 2) * Power(r, 2)
                  + 2 * M * (-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                  + Power(a, 2) * Power(r, 4) - 2 * M * Power(r, 5) + Power(r, 6))));

  mixedK[0][1]
      = (M * Power(r, 2)
         * ((4 * M * Power(r, 3)
             * (Power(a, 2) * x * y + a * (Power(x, 2) - Power(y, 2)) * r - x * y * Power(r, 2))
             * (-4 * Power(r, 4) * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                    * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                + 3 * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                + 2 * (a * x - y * r) * Power(Power(a, 2) + Power(r, 2), 3)
                      * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dy
                + r * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (r + y * dr_dy)
                - (M * Power(r, 3) * Power(a * x - y * r, 2)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (2 * Power(r, 8) - 3 * Power(a, 5) * x * Power(z, 2) * dr_dy
                      + 5 * Power(a, 4) * y * Power(z, 2) * r * dr_dy
                      + Power(a, 2) * y * Power(z, 2) * Power(r, 3) * dr_dy
                      + Power(a, 2) * y * Power(r, 5) * dr_dy - 3 * y * Power(r, 7) * dr_dy
                      + Power(a, 3) * Power(z, 2) * Power(r, 2) * (2 * a + x * dr_dy)
                      + a * Power(r, 6) * (2 * a + 5 * x * dr_dy)
                      + Power(a, 2) * Power(r, 4) * (2 * Power(z, 2) + a * x * dr_dy)))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                - (M * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                   * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                      + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                      + Power(a, 2) * Power(r, 4) + Power(r, 6))
                   * (2 * z * Power(r, 10) - 4 * Power(a, 7) * x * Power(z, 3) * dr_dy
                      + 2 * Power(a, 4) * z * Power(r, 4)
                            * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2) + 2 * a * x * dr_dy)
                      + 2 * Power(a, 2) * z * Power(r, 6)
                            * (Power(a, 2) + Power(y, 2) + Power(z, 2) + 6 * a * x * dr_dy)
                      + 3 * y * Power(r, 9) * (y * dr_dz - 2 * z * dr_dy)
                      + Power(a, 2) * Power(r, 7)
                            * ((5 * Power(x, 2) - Power(y, 2)) * dr_dz - 8 * y * z * dr_dy)
                      + 3 * Power(a, 6) * Power(z, 2) * r
                            * (-(Power(x, 2) * dr_dz) + 2 * y * z * dr_dy)
                      + 2 * Power(a, 5) * z * Power(r, 2)
                            * (a * (Power(x, 2) + Power(z, 2)) + 4 * x * y * z * dr_dz
                               - 2 * x * Power(z, 2) * dr_dy)
                      + Power(a, 4) * z * Power(r, 3)
                            * (-4 * a * x * y + (Power(x, 2) - 5 * Power(y, 2)) * z * dr_dz
                               + 8 * y * Power(z, 2) * dr_dy)
                      + 4 * a * Power(r, 8) * (-2 * x * y * dr_dz + z * (a + 2 * x * dr_dy))
                      + Power(a, 2) * Power(r, 5)
                            * ((Power(a, 2) * Power(x, 2) - Power(y, 2) * Power(z, 2)) * dr_dz
                               + 2 * y * z * (-2 * a * x + (-Power(a, 2) + Power(z, 2)) * dr_dy))))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                - (M * Power(r, 3) * (a * y + x * r)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (-2 * Power(a, 4) * x * Power(z, 2) * Power(r, 3)
                      + 2 * Power(a, 2) * x * (-2 * Power(a, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * x * Power(r, 9)
                      - 3 * Power(a, 6) * x * Power(z, 2) * (2 * y * dr_dy + x * dr_dx)
                      + 3 * y * Power(r, 8) * (2 * a - 2 * x * dr_dy + y * dr_dx)
                      + Power(a, 4) * Power(z, 2) * Power(r, 2)
                            * (6 * a * y + 12 * x * y * dr_dy
                               + (Power(x, 2) - 5 * Power(y, 2)) * dr_dx)
                      + Power(a, 2) * Power(r, 6)
                            * (6 * a * y + 12 * x * y * dr_dy
                               + (5 * Power(x, 2) - Power(y, 2)) * dr_dx)
                      + Power(a, 2) * Power(r, 4)
                            * (6 * a * y * Power(z, 2)
                               + 2 * x * y * (Power(a, 2) + Power(z, 2)) * dr_dy
                               + (Power(a, 2) * Power(x, 2) - Power(y, 2) * Power(z, 2)) * dr_dx)
                      - 4 * Power(a, 5) * Power(z, 2) * r
                            * (2 * (Power(x, 2) - Power(y, 2)) * dr_dy + x * (a - 2 * y * dr_dx))
                      - 2 * a * Power(r, 7)
                            * (-4 * (Power(x, 2) - Power(y, 2)) * dr_dy + x * (a + 4 * y * dr_dx))))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))))
                / Power(Power(a, 2) + Power(r, 2), 6)
            + ((Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
                + Power(r, 8))
               * (-4 * Power(r, 4) * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                  + 3 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                  - 2 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 3)
                        * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dy
                  + r * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (a + x * dr_dy)
                  - (2 * M * Power(r, 3) * Power(a * y + x * r, 2)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (3 * Power(a, 5) * y * Power(z, 2) * dr_dy
                        - Power(a, 3) * y * Power(z, 2) * Power(r, 2) * dr_dy
                        - Power(a, 3) * y * Power(r, 4) * dr_dy - 5 * a * y * Power(r, 6) * dr_dy
                        + Power(r, 7) * (2 * a - 3 * x * dr_dy)
                        + Power(a, 2) * Power(z, 2) * Power(r, 3) * (2 * a + x * dr_dy)
                        + Power(a, 2) * Power(r, 5) * (2 * a + x * dr_dy)
                        + Power(a, 4) * Power(z, 2) * r * (2 * a + 5 * x * dr_dy)))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  - 4 * Power(r, 4) * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                        * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                  + 3 * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                  + 2 * (a * x - y * r) * Power(Power(a, 2) + Power(r, 2), 3)
                        * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dx
                  + r * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (-a + y * dr_dx)
                  + (2 * M * Power(r, 3) * Power(a * x - y * r, 2)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (3 * Power(a, 5) * x * Power(z, 2) * dr_dx
                        - Power(a, 3) * x * Power(z, 2) * Power(r, 2) * dr_dx
                        - Power(a, 3) * x * Power(r, 4) * dr_dx - 5 * a * x * Power(r, 6) * dr_dx
                        + Power(a, 4) * Power(z, 2) * r * (2 * a - 5 * y * dr_dx)
                        + Power(a, 2) * Power(z, 2) * Power(r, 3) * (2 * a - y * dr_dx)
                        + Power(a, 2) * Power(r, 5) * (2 * a - y * dr_dx)
                        + Power(r, 7) * (2 * a + 3 * y * dr_dx)))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  + (2 * M * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                     * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                        + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                        + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     * (2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy - x * dr_dx)
                        + 2 * Power(a, 7) * Power(z, 3) * (-(y * dr_dy) + x * dr_dx)
                        - 3 * Power(a, 6) * Power(z, 2) * r
                              * (x * y * dr_dz + x * z * dr_dy + y * z * dr_dx)
                        + 2 * Power(a, 5) * z * Power(r, 2)
                              * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                 - y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                        + 2 * Power(a, 4) * z * Power(r, 3)
                              * (a * Power(x, 2) - a * Power(y, 2) + 3 * x * y * z * dr_dz
                                 - 2 * x * Power(z, 2) * dr_dy - 2 * y * Power(z, 2) * dr_dx)
                        - 2 * Power(a, 2) * z * Power(r, 6)
                              * (-3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                        + 4 * a * Power(r, 8)
                              * ((Power(x, 2) - Power(y, 2)) * dr_dz + z * (y * dr_dy - x * dr_dx))
                        + 2 * Power(a, 2) * Power(r, 7)
                              * (3 * x * y * dr_dz + 2 * z * (x * dr_dy + y * dr_dx))
                        + Power(r, 9) * (-3 * x * y * dr_dz + 3 * z * (x * dr_dy + y * dr_dx))
                        + Power(a, 2) * Power(r, 5)
                              * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                 + z
                                       * (2 * a * (Power(x, 2) - Power(y, 2))
                                          + x * (Power(a, 2) - Power(z, 2)) * dr_dy
                                          + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))))
                  / Power(Power(a, 2) + Power(r, 2), 6)
            - 2 * M * z * r * (a * y + x * r)
                  * ((M * Power(r, 4) * Power(a * x - y * r, 2)
                      * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                         + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                         + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                         + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                         + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                      * (3 * Power(a, 5) * x * Power(z, 2) * dr_dz
                         - Power(a, 2) * y * Power(r, 5) * dr_dz - 5 * a * x * Power(r, 6) * dr_dz
                         + 3 * y * Power(r, 7) * dr_dz
                         + Power(a, 2) * Power(r, 4) * (2 * y * z - a * x * dr_dz)
                         + Power(a, 3) * z * Power(r, 2) * (2 * a * y - x * z * dr_dz)
                         - Power(a, 2) * z * Power(r, 3) * (2 * a * x + y * z * dr_dz)
                         - Power(a, 4) * z * r * (2 * a * x + 5 * y * z * dr_dz)))
                         / (Power(Power(a, 2) + Power(r, 2), 5)
                            * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                               + Power(a, 2) * Power(z, 2) * Power(r, 2)
                               + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                               + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     - 4 * z * Power(r, 4) * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                     + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                     - (M * Power(z, 3) * r * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                        * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                        * dr_dy)
                           / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                              + Power(a, 2) * Power(z, 2) * Power(r, 2)
                              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                              + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     - (M * Power(r, 3) * (a * y + x * r)
                        * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                           + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                           + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                           + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                           + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                           + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                        * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                           - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                           + 2 * Power(a, 5) * z * Power(r, 2)
                                 * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                    + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                           + 2 * Power(a, 4) * z * Power(r, 3)
                                 * (a * Power(x, 2) - a * Power(y, 2) + 2 * a * Power(z, 2)
                                    + 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                    - 2 * y * Power(z, 2) * dr_dx)
                           - 2 * Power(a, 2) * z * Power(r, 6)
                                 * (3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                           + 4 * a * Power(r, 8)
                                 * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                    - z * (y * dr_dy + x * dr_dx))
                           + Power(a, 6) * Power(z, 2) * r
                                 * (-3 * x * y * dr_dz
                                    + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                           + 2 * Power(a, 2) * Power(r, 7)
                                 * (3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                           + Power(r, 9)
                                 * (-3 * x * y * dr_dz
                                    + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                           + Power(a, 2) * Power(r, 5)
                                 * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                    + z
                                          * (2 * a
                                                 * (Power(a, 2) + Power(x, 2) - Power(y, 2)
                                                    + Power(z, 2))
                                             + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                             + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                           / (Power(Power(a, 2) + Power(r, 2), 5)
                              * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                                 + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                 + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                                 + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     + r
                           * ((y * r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                  / (Power(a, 2) + Power(r, 2))
                              + (3 * (-(a * x) + y * r)
                                 * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                    / (Power(a, 2) + Power(r, 2))
                              + (2 * (a * x - y * r)
                                 * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dz)
                                    / Power(Power(a, 2) + Power(r, 2), 2)
                              - (2 * r * (-(a * x) + y * r)
                                 * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                                 * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz))
                                    / (Power(a, 2) + Power(r, 2))
                              + (M * Power(r, 3) * Power(a * x - y * r, 2)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (3 * Power(a, 5) * x * Power(z, 2) * dr_dz
                                    - Power(a, 2) * y * Power(r, 5) * dr_dz
                                    - 5 * a * x * Power(r, 6) * dr_dz + 3 * y * Power(r, 7) * dr_dz
                                    + Power(a, 2) * Power(r, 4) * (2 * y * z - a * x * dr_dz)
                                    + Power(a, 3) * z * Power(r, 2) * (2 * a * y - x * z * dr_dz)
                                    - Power(a, 2) * z * Power(r, 3) * (2 * a * x + y * z * dr_dz)
                                    - Power(a, 4) * z * r * (2 * a * x + 5 * y * z * dr_dz)))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                              - (M * Power(z, 3) * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                                 * (Power(a, 4) * Power(z, 2)
                                    + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                    + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                                    + Power(a, 2) * Power(r, 4) + Power(r, 6))
                                 * dr_dy)
                                    / (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6))
                              - (M * Power(r, 2) * (a * y + x * r)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                                    - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                                    + 2 * Power(a, 5) * z * Power(r, 2)
                                          * (a * x * y
                                             + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                             + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                                    + 2 * Power(a, 4) * z * Power(r, 3)
                                          * (a * Power(x, 2) - a * Power(y, 2) + 2 * a * Power(z, 2)
                                             + 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                             - 2 * y * Power(z, 2) * dr_dx)
                                    - 2 * Power(a, 2) * z * Power(r, 6)
                                          * (3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                                    + 4 * a * Power(r, 8)
                                          * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                             - z * (y * dr_dy + x * dr_dx))
                                    + Power(a, 6) * Power(z, 2) * r
                                          * (-3 * x * y * dr_dz
                                             + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                                    + 2 * Power(a, 2) * Power(r, 7)
                                          * (3 * x * y * dr_dz
                                             + 2 * z * (a - x * dr_dy + y * dr_dx))
                                    + Power(r, 9)
                                          * (-3 * x * y * dr_dz
                                             + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                                    + Power(a, 2) * Power(r, 5)
                                          * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                             + z
                                                   * (2 * a
                                                          * (Power(a, 2) + Power(x, 2) - Power(y, 2)
                                                             + Power(z, 2))
                                                      + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                                      + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))))))
        / (Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 3)
           * Sqrt(
               (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                + Power(a, 2) * Power(z, 2) * Power(r, 2)
                + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                + Power(a, 2) * Power(r, 4) + Power(r, 6))
               * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                  + Power(a, 2) * Power(z, 2) * Power(r, 2)
                  + 2 * M * (-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                  + Power(a, 2) * Power(r, 4) - 2 * M * Power(r, 5) + Power(r, 6))));

  mixedK[0][2]
      = (M * r
         * (-4 * M * z * Power(r, 2) * (a * y + x * r)
                * (r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
                   + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz
                   - 2 * z * r * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                         * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz)
                   - (M * Power(z, 2) * r
                      * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                      * (2 * Power(r, 5) + Power(a, 2) * Power(z, 3) * dr_dz
                         - 3 * z * Power(r, 4) * dr_dz))
                         / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                            + Power(a, 2) * Power(z, 2) * Power(r, 2)
                            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                            + Power(a, 2) * Power(r, 4) + Power(r, 6))
                   + (M * Power(r, 2) * (-(a * x) + y * r)
                      * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                         + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                         + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                         + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                         + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                      * (2 * Power(a, 4) * y * Power(z, 2) * Power(r, 3) - 2 * y * Power(r, 9)
                         + 4 * Power(a, 5) * x * Power(z, 3) * r * dr_dz
                         + 2 * Power(a, 2) * z * Power(r, 5) * (y * z - 2 * a * x * dr_dz)
                         - 2 * a * Power(r, 7) * (a * y + 4 * x * z * dr_dz)
                         + Power(a, 6) * Power(z, 4) * dr_dy
                         + 2 * Power(a, 2) * Power(r, 6)
                               * (a * x + y * z * dr_dz - 3 * Power(z, 2) * dr_dy)
                         + Power(r, 8) * (2 * a * x + 6 * y * z * dr_dz - 3 * Power(z, 2) * dr_dy)
                         - 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                               * (a * x + 3 * y * z * dr_dz - Power(z, 2) * dr_dy)
                         + Power(a, 2) * Power(z, 2) * Power(r, 4)
                               * (-2 * a * x - 2 * y * z * dr_dz
                                  + (-3 * Power(a, 2) + Power(z, 2)) * dr_dy)))
                         / (Power(Power(a, 2) + Power(r, 2), 4)
                            * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                               + Power(a, 2) * Power(z, 2) * Power(r, 2)
                               + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                               + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                   + (M * Power(r, 2) * (a * y + x * r)
                      * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                         + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                         + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                         + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                         + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                      * (2 * Power(a, 4) * x * Power(z, 2) * Power(r, 3) - 2 * x * Power(r, 9)
                         - 4 * Power(a, 5) * y * Power(z, 3) * r * dr_dz
                         + 2 * Power(a, 2) * z * Power(r, 5) * (x * z + 2 * a * y * dr_dz)
                         - 2 * a * Power(r, 7) * (a * x - 4 * y * z * dr_dz)
                         + Power(a, 6) * Power(z, 4) * dr_dx
                         + Power(r, 8) * (-2 * a * y + 6 * x * z * dr_dz - 3 * Power(z, 2) * dr_dx)
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                               * (a * y - 3 * x * z * dr_dz + Power(z, 2) * dr_dx)
                         - 2 * Power(a, 2) * Power(r, 6)
                               * (a * y - x * z * dr_dz + 3 * Power(z, 2) * dr_dx)
                         + Power(a, 2) * Power(z, 2) * Power(r, 4)
                               * (2 * a * y - 2 * x * z * dr_dz
                                  + (-3 * Power(a, 2) + Power(z, 2)) * dr_dx)))
                         / (Power(Power(a, 2) + Power(r, 2), 4)
                            * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                               + Power(a, 2) * Power(z, 2) * Power(r, 2)
                               + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                               + Power(a, 2) * Power(r, 4) + Power(r, 6))))
            + (2 * M * Power(r, 3)
               * (Power(a, 2) * x * y + a * (Power(x, 2) - Power(y, 2)) * r - x * y * Power(r, 2))
               * ((M * Power(r, 4) * Power(a * x - y * r, 2)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (3 * Power(a, 5) * x * Power(z, 2) * dr_dz
                      - Power(a, 2) * y * Power(r, 5) * dr_dz - 5 * a * x * Power(r, 6) * dr_dz
                      + 3 * y * Power(r, 7) * dr_dz
                      + Power(a, 2) * Power(r, 4) * (2 * y * z - a * x * dr_dz)
                      + Power(a, 3) * z * Power(r, 2) * (2 * a * y - x * z * dr_dz)
                      - Power(a, 2) * z * Power(r, 3) * (2 * a * x + y * z * dr_dz)
                      - Power(a, 4) * z * r * (2 * a * x + 5 * y * z * dr_dz)))
                      / (Power(Power(a, 2) + Power(r, 2), 5)
                         * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                            + Power(a, 2) * Power(z, 2) * Power(r, 2)
                            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                            + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                  - 4 * z * Power(r, 4) * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                  + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                  - (M * Power(z, 3) * r * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                     * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                        + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                        + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     * dr_dy)
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  - (M * Power(r, 3) * (a * y + x * r)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                        - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                        + 2 * Power(a, 5) * z * Power(r, 2)
                              * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                 + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                        + 2 * Power(a, 4) * z * Power(r, 3)
                              * (a * Power(x, 2) - a * Power(y, 2) + 2 * a * Power(z, 2)
                                 + 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                 - 2 * y * Power(z, 2) * dr_dx)
                        - 2 * Power(a, 2) * z * Power(r, 6)
                              * (3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                        + 4 * a * Power(r, 8)
                              * ((Power(x, 2) - Power(y, 2)) * dr_dz - z * (y * dr_dy + x * dr_dx))
                        + Power(a, 6) * Power(z, 2) * r
                              * (-3 * x * y * dr_dz + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                        + 2 * Power(a, 2) * Power(r, 7)
                              * (3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                        + Power(r, 9)
                              * (-3 * x * y * dr_dz + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                        + Power(a, 2) * Power(r, 5)
                              * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                 + z
                                       * (2 * a
                                              * (Power(a, 2) + Power(x, 2) - Power(y, 2)
                                                 + Power(z, 2))
                                          + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                          + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                        / (Power(Power(a, 2) + Power(r, 2), 5)
                           * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                              + Power(a, 2) * Power(z, 2) * Power(r, 2)
                              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                              + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                  + r
                        * ((y * r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                               / (Power(a, 2) + Power(r, 2))
                           + (3 * (-(a * x) + y * r)
                              * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                 / (Power(a, 2) + Power(r, 2))
                           + (2 * (a * x - y * r)
                              * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dz)
                                 / Power(Power(a, 2) + Power(r, 2), 2)
                           - (2 * r * (-(a * x) + y * r) * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                              * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz))
                                 / (Power(a, 2) + Power(r, 2))
                           + (M * Power(r, 3) * Power(a * x - y * r, 2)
                              * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                                 + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                 + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2))
                                       * Power(r, 3)
                                 + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                                 + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                                 + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                              * (3 * Power(a, 5) * x * Power(z, 2) * dr_dz
                                 - Power(a, 2) * y * Power(r, 5) * dr_dz
                                 - 5 * a * x * Power(r, 6) * dr_dz + 3 * y * Power(r, 7) * dr_dz
                                 + Power(a, 2) * Power(r, 4) * (2 * y * z - a * x * dr_dz)
                                 + Power(a, 3) * z * Power(r, 2) * (2 * a * y - x * z * dr_dz)
                                 - Power(a, 2) * z * Power(r, 3) * (2 * a * x + y * z * dr_dz)
                                 - Power(a, 4) * z * r * (2 * a * x + 5 * y * z * dr_dz)))
                                 / (Power(Power(a, 2) + Power(r, 2), 5)
                                    * (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                           - (M * Power(z, 3) * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                              * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                 + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                                 + Power(a, 2) * Power(r, 4) + Power(r, 6))
                              * dr_dy)
                                 / (Power(a, 4) * Power(z, 2)
                                    + 2 * Power(a, 2) * M * Power(z, 2) * r
                                    + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                    + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                          * Power(r, 3)
                                    + Power(a, 2) * Power(r, 4) + Power(r, 6))
                           - (M * Power(r, 2) * (a * y + x * r)
                              * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                                 + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                 + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2))
                                       * Power(r, 3)
                                 + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                                 + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                                 + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                              * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                                 - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                                 + 2 * Power(a, 5) * z * Power(r, 2)
                                       * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                          + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                                 + 2 * Power(a, 4) * z * Power(r, 3)
                                       * (a * Power(x, 2) - a * Power(y, 2) + 2 * a * Power(z, 2)
                                          + 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                          - 2 * y * Power(z, 2) * dr_dx)
                                 - 2 * Power(a, 2) * z * Power(r, 6)
                                       * (3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                                 + 4 * a * Power(r, 8)
                                       * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                          - z * (y * dr_dy + x * dr_dx))
                                 + Power(a, 6) * Power(z, 2) * r
                                       * (-3 * x * y * dr_dz
                                          + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                                 + 2 * Power(a, 2) * Power(r, 7)
                                       * (3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                                 + Power(r, 9)
                                       * (-3 * x * y * dr_dz
                                          + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                                 + Power(a, 2) * Power(r, 5)
                                       * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                          + z
                                                * (2 * a
                                                       * (Power(a, 2) + Power(x, 2) - Power(y, 2)
                                                          + Power(z, 2))
                                                   + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                                   + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                                 / (Power(Power(a, 2) + Power(r, 2), 5)
                                    * (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6))))))
                  / (Power(a, 2) + Power(r, 2))
            + ((Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
                + Power(r, 8))
               * ((M * Power(r, 4) * Power(a * y + x * r, 2)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (-3 * Power(a, 5) * y * Power(z, 2) * dr_dz
                      - Power(a, 2) * x * Power(r, 5) * dr_dz + 5 * a * y * Power(r, 6) * dr_dz
                      + 3 * x * Power(r, 7) * dr_dz
                      + Power(a, 2) * Power(r, 4) * (2 * x * z + a * y * dr_dz)
                      + Power(a, 4) * z * r * (2 * a * y - 5 * x * z * dr_dz)
                      + Power(a, 2) * z * Power(r, 3) * (2 * a * y - x * z * dr_dz)
                      + Power(a, 3) * z * Power(r, 2) * (2 * a * x + y * z * dr_dz)))
                      / (Power(Power(a, 2) + Power(r, 2), 5)
                         * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                            + Power(a, 2) * Power(z, 2) * Power(r, 2)
                            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                            + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                  - 4 * z * Power(r, 4) * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                  + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                  - (M * Power(z, 3) * r * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                     * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                        + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                        + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     * dr_dx)
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  + (M * Power(r, 3) * (-(a * x) + y * r)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                        - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                        + 2 * Power(a, 5) * z * Power(r, 2)
                              * (-(a * x * y) + 2 * (Power(x, 2) - Power(y, 2)) * z * dr_dz
                                 + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                        + 2 * Power(a, 4) * z * Power(r, 3)
                              * (-(a * Power(x, 2)) + a * Power(y, 2) + 2 * a * Power(z, 2)
                                 - 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                 - 2 * y * Power(z, 2) * dr_dx)
                        + 2 * Power(a, 2) * z * Power(r, 6)
                              * (-3 * a * y * dr_dy + x * (y - 3 * a * dr_dx))
                        - 4 * a * Power(r, 8)
                              * ((Power(x, 2) - Power(y, 2)) * dr_dz + z * (y * dr_dy + x * dr_dx))
                        + Power(a, 6) * Power(z, 2) * r
                              * (3 * x * y * dr_dz + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                        + 2 * Power(a, 2) * Power(r, 7)
                              * (-3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                        + Power(r, 9)
                              * (3 * x * y * dr_dz + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                        + Power(a, 2) * Power(r, 5)
                              * (-(x * y * (Power(a, 2) + Power(z, 2)) * dr_dz)
                                 + z
                                       * (2 * a
                                              * (Power(a, 2) - Power(x, 2) + Power(y, 2)
                                                 + Power(z, 2))
                                          + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                          + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                        / (Power(Power(a, 2) + Power(r, 2), 5)
                           * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                              + Power(a, 2) * Power(z, 2) * Power(r, 2)
                              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                              + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                  + r
                        * ((x * r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                               / (Power(a, 2) + Power(r, 2))
                           + (3 * (a * y + x * r)
                              * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                 / (Power(a, 2) + Power(r, 2))
                           - (2 * (a * y + x * r)
                              * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dz)
                                 / Power(Power(a, 2) + Power(r, 2), 2)
                           - (2 * r * (a * y + x * r) * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                              * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz))
                                 / (Power(a, 2) + Power(r, 2))
                           + (M * Power(r, 3) * Power(a * y + x * r, 2)
                              * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                                 + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                 + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2))
                                       * Power(r, 3)
                                 + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                                 + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                                 + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                              * (-3 * Power(a, 5) * y * Power(z, 2) * dr_dz
                                 - Power(a, 2) * x * Power(r, 5) * dr_dz
                                 + 5 * a * y * Power(r, 6) * dr_dz + 3 * x * Power(r, 7) * dr_dz
                                 + Power(a, 2) * Power(r, 4) * (2 * x * z + a * y * dr_dz)
                                 + Power(a, 4) * z * r * (2 * a * y - 5 * x * z * dr_dz)
                                 + Power(a, 2) * z * Power(r, 3) * (2 * a * y - x * z * dr_dz)
                                 + Power(a, 3) * z * Power(r, 2) * (2 * a * x + y * z * dr_dz)))
                                 / (Power(Power(a, 2) + Power(r, 2), 5)
                                    * (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                           - (M * Power(z, 3) * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                              * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                 + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                                 + Power(a, 2) * Power(r, 4) + Power(r, 6))
                              * dr_dx)
                                 / (Power(a, 4) * Power(z, 2)
                                    + 2 * Power(a, 2) * M * Power(z, 2) * r
                                    + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                    + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                          * Power(r, 3)
                                    + Power(a, 2) * Power(r, 4) + Power(r, 6))
                           + (M * Power(r, 2) * (-(a * x) + y * r)
                              * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                                 + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                 + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2))
                                       * Power(r, 3)
                                 + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                                 + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                                 + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                              * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                                 - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                                 + 2 * Power(a, 5) * z * Power(r, 2)
                                       * (-(a * x * y) + 2 * (Power(x, 2) - Power(y, 2)) * z * dr_dz
                                          + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                                 + 2 * Power(a, 4) * z * Power(r, 3)
                                       * (-(a * Power(x, 2)) + a * Power(y, 2) + 2 * a * Power(z, 2)
                                          - 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                          - 2 * y * Power(z, 2) * dr_dx)
                                 + 2 * Power(a, 2) * z * Power(r, 6)
                                       * (-3 * a * y * dr_dy + x * (y - 3 * a * dr_dx))
                                 - 4 * a * Power(r, 8)
                                       * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                          + z * (y * dr_dy + x * dr_dx))
                                 + Power(a, 6) * Power(z, 2) * r
                                       * (3 * x * y * dr_dz
                                          + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                                 + 2 * Power(a, 2) * Power(r, 7)
                                       * (-3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                                 + Power(r, 9)
                                       * (3 * x * y * dr_dz
                                          + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                                 + Power(a, 2) * Power(r, 5)
                                       * (-(x * y * (Power(a, 2) + Power(z, 2)) * dr_dz)
                                          + z
                                                * (2 * a
                                                       * (Power(a, 2) - Power(x, 2) + Power(y, 2)
                                                          + Power(z, 2))
                                                   + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                                   + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                                 / (Power(Power(a, 2) + Power(r, 2), 5)
                                    * (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6))))))
                  / (Power(a, 2) + Power(r, 2))))
        / (Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 3)
           * Sqrt(
               (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                + Power(a, 2) * Power(z, 2) * Power(r, 2)
                + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                + Power(a, 2) * Power(r, 4) + Power(r, 6))
               * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                  + Power(a, 2) * Power(z, 2) * Power(r, 2)
                  + 2 * M * (-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                  + Power(a, 2) * Power(r, 4) - 2 * M * Power(r, 5) + Power(r, 6))));

  mixedK[1][0]
      = (M * Power(r, 2)
         * (((Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
              + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
              + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
              + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
              + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
              + Power(r, 8))
             * (-4 * Power(r, 4) * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                    * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                + 3 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                - 2 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 3)
                      * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dy
                + r * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (a + x * dr_dy)
                - (2 * M * Power(r, 3) * Power(a * y + x * r, 2)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (3 * Power(a, 5) * y * Power(z, 2) * dr_dy
                      - Power(a, 3) * y * Power(z, 2) * Power(r, 2) * dr_dy
                      - Power(a, 3) * y * Power(r, 4) * dr_dy - 5 * a * y * Power(r, 6) * dr_dy
                      + Power(r, 7) * (2 * a - 3 * x * dr_dy)
                      + Power(a, 2) * Power(z, 2) * Power(r, 3) * (2 * a + x * dr_dy)
                      + Power(a, 2) * Power(r, 5) * (2 * a + x * dr_dy)
                      + Power(a, 4) * Power(z, 2) * r * (2 * a + 5 * x * dr_dy)))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                - 4 * Power(r, 4) * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                + 3 * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                + 2 * (a * x - y * r) * Power(Power(a, 2) + Power(r, 2), 3)
                      * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dx
                + r * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (-a + y * dr_dx)
                + (2 * M * Power(r, 3) * Power(a * x - y * r, 2)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (3 * Power(a, 5) * x * Power(z, 2) * dr_dx
                      - Power(a, 3) * x * Power(z, 2) * Power(r, 2) * dr_dx
                      - Power(a, 3) * x * Power(r, 4) * dr_dx - 5 * a * x * Power(r, 6) * dr_dx
                      + Power(a, 4) * Power(z, 2) * r * (2 * a - 5 * y * dr_dx)
                      + Power(a, 2) * Power(z, 2) * Power(r, 3) * (2 * a - y * dr_dx)
                      + Power(a, 2) * Power(r, 5) * (2 * a - y * dr_dx)
                      + Power(r, 7) * (2 * a + 3 * y * dr_dx)))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                + (2 * M * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                   * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                      + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                      + Power(a, 2) * Power(r, 4) + Power(r, 6))
                   * (2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy - x * dr_dx)
                      + 2 * Power(a, 7) * Power(z, 3) * (-(y * dr_dy) + x * dr_dx)
                      - 3 * Power(a, 6) * Power(z, 2) * r
                            * (x * y * dr_dz + x * z * dr_dy + y * z * dr_dx)
                      + 2 * Power(a, 5) * z * Power(r, 2)
                            * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                               - y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                      + 2 * Power(a, 4) * z * Power(r, 3)
                            * (a * Power(x, 2) - a * Power(y, 2) + 3 * x * y * z * dr_dz
                               - 2 * x * Power(z, 2) * dr_dy - 2 * y * Power(z, 2) * dr_dx)
                      - 2 * Power(a, 2) * z * Power(r, 6)
                            * (-3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                      + 4 * a * Power(r, 8)
                            * ((Power(x, 2) - Power(y, 2)) * dr_dz + z * (y * dr_dy - x * dr_dx))
                      + 2 * Power(a, 2) * Power(r, 7)
                            * (3 * x * y * dr_dz + 2 * z * (x * dr_dy + y * dr_dx))
                      + Power(r, 9) * (-3 * x * y * dr_dz + 3 * z * (x * dr_dy + y * dr_dx))
                      + Power(a, 2) * Power(r, 5)
                            * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                               + z
                                     * (2 * a * (Power(x, 2) - Power(y, 2))
                                        + x * (Power(a, 2) - Power(z, 2)) * dr_dy
                                        + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))))
                / Power(Power(a, 2) + Power(r, 2), 6)
            + (4 * M * Power(r, 3)
               * (Power(a, 2) * x * y + a * (Power(x, 2) - Power(y, 2)) * r - x * y * Power(r, 2))
               * (-4 * Power(r, 4) * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                  + 3 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                  - 2 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 3)
                        * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dx
                  + r * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (r + x * dr_dx)
                  - (M * Power(r, 3) * Power(a * y + x * r, 2)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (2 * Power(r, 8) + 3 * Power(a, 5) * y * Power(z, 2) * dr_dx
                        + 5 * Power(a, 4) * x * Power(z, 2) * r * dr_dx
                        + Power(a, 2) * x * Power(z, 2) * Power(r, 3) * dr_dx
                        + Power(a, 2) * x * Power(r, 5) * dr_dx - 3 * x * Power(r, 7) * dr_dx
                        + a * Power(r, 6) * (2 * a - 5 * y * dr_dx)
                        + Power(a, 3) * Power(z, 2) * Power(r, 2) * (2 * a - y * dr_dx)
                        + Power(a, 2) * Power(r, 4) * (2 * Power(z, 2) - a * y * dr_dx)))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  + (M * Power(r, 3) * (-(a * x) + y * r)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (2 * Power(a, 4) * y * Power(z, 2) * Power(r, 3)
                        + (4 * Power(a, 4) * y - 2 * Power(a, 2) * y * Power(z, 2)) * Power(r, 5)
                        - 2 * y * Power(r, 9)
                        + 3 * Power(a, 6) * y * Power(z, 2) * (y * dr_dy + 2 * x * dr_dx)
                        + 2 * a * Power(r, 7)
                              * (a * y - 4 * x * y * dr_dy
                                 - 4 * (Power(x, 2) - Power(y, 2)) * dr_dx)
                        + 4 * Power(a, 5) * Power(z, 2) * r
                              * (a * y + 2 * x * y * dr_dy
                                 + 2 * (Power(x, 2) - Power(y, 2)) * dr_dx)
                        + Power(a, 2) * Power(r, 4)
                              * (6 * a * x * Power(z, 2)
                                 + (-(Power(a, 2) * Power(y, 2)) + Power(x, 2) * Power(z, 2))
                                       * dr_dy
                                 - 2 * x * y * (Power(a, 2) + Power(z, 2)) * dr_dx)
                        + Power(a, 2) * Power(r, 6)
                              * ((Power(x, 2) - 5 * Power(y, 2)) * dr_dy
                                 + 6 * x * (a - 2 * y * dr_dx))
                        + Power(a, 4) * Power(z, 2) * Power(r, 2)
                              * ((5 * Power(x, 2) - Power(y, 2)) * dr_dy
                                 + 6 * x * (a - 2 * y * dr_dx))
                        + 3 * x * Power(r, 8) * (-(x * dr_dy) + 2 * (a + y * dr_dx))))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  - (M * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                     * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                        + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                        + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     * (2 * z * Power(r, 10) + 4 * Power(a, 7) * y * Power(z, 3) * dr_dx
                        + 2 * Power(a, 2) * z * Power(r, 6)
                              * (Power(a, 2) + Power(x, 2) + Power(z, 2) - 6 * a * y * dr_dx)
                        + 2 * Power(a, 4) * z * Power(r, 4)
                              * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2) - 2 * a * y * dr_dx)
                        + 3 * x * Power(r, 9) * (x * dr_dz - 2 * z * dr_dx)
                        + 3 * Power(a, 6) * Power(z, 2) * r
                              * (-(Power(y, 2) * dr_dz) + 2 * x * z * dr_dx)
                        - Power(a, 2) * Power(r, 7)
                              * ((Power(x, 2) - 5 * Power(y, 2)) * dr_dz + 8 * x * z * dr_dx)
                        + 2 * Power(a, 5) * z * Power(r, 2)
                              * (a * (Power(y, 2) + Power(z, 2)) - 4 * x * y * z * dr_dz
                                 + 2 * y * Power(z, 2) * dr_dx)
                        + 4 * a * Power(r, 8) * (2 * x * y * dr_dz + z * (a - 2 * y * dr_dx))
                        + Power(a, 4) * z * Power(r, 3)
                              * ((-5 * Power(x, 2) + Power(y, 2)) * z * dr_dz
                                 + 4 * x * (a * y + 2 * Power(z, 2) * dr_dx))
                        + Power(a, 2) * Power(r, 5)
                              * ((Power(a, 2) * Power(y, 2) - Power(x, 2) * Power(z, 2)) * dr_dz
                                 + 2 * x * z * (2 * a * y + (-Power(a, 2) + Power(z, 2)) * dr_dx))))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))))
                  / Power(Power(a, 2) + Power(r, 2), 6)
            + 2 * M * z * r * (a * x - y * r)
                  * ((M * Power(r, 4) * Power(a * y + x * r, 2)
                      * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                         + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                         + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                         + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                         + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                      * (-3 * Power(a, 5) * y * Power(z, 2) * dr_dz
                         - Power(a, 2) * x * Power(r, 5) * dr_dz + 5 * a * y * Power(r, 6) * dr_dz
                         + 3 * x * Power(r, 7) * dr_dz
                         + Power(a, 2) * Power(r, 4) * (2 * x * z + a * y * dr_dz)
                         + Power(a, 4) * z * r * (2 * a * y - 5 * x * z * dr_dz)
                         + Power(a, 2) * z * Power(r, 3) * (2 * a * y - x * z * dr_dz)
                         + Power(a, 3) * z * Power(r, 2) * (2 * a * x + y * z * dr_dz)))
                         / (Power(Power(a, 2) + Power(r, 2), 5)
                            * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                               + Power(a, 2) * Power(z, 2) * Power(r, 2)
                               + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                               + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     - 4 * z * Power(r, 4) * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                     + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                     - (M * Power(z, 3) * r * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                        * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                        * dr_dx)
                           / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                              + Power(a, 2) * Power(z, 2) * Power(r, 2)
                              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                              + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     + (M * Power(r, 3) * (-(a * x) + y * r)
                        * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                           + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                           + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                           + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                           + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                           + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                        * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                           - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                           + 2 * Power(a, 5) * z * Power(r, 2)
                                 * (-(a * x * y) + 2 * (Power(x, 2) - Power(y, 2)) * z * dr_dz
                                    + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                           + 2 * Power(a, 4) * z * Power(r, 3)
                                 * (-(a * Power(x, 2)) + a * Power(y, 2) + 2 * a * Power(z, 2)
                                    - 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                    - 2 * y * Power(z, 2) * dr_dx)
                           + 2 * Power(a, 2) * z * Power(r, 6)
                                 * (-3 * a * y * dr_dy + x * (y - 3 * a * dr_dx))
                           - 4 * a * Power(r, 8)
                                 * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                    + z * (y * dr_dy + x * dr_dx))
                           + Power(a, 6) * Power(z, 2) * r
                                 * (3 * x * y * dr_dz + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                           + 2 * Power(a, 2) * Power(r, 7)
                                 * (-3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                           + Power(r, 9)
                                 * (3 * x * y * dr_dz + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                           + Power(a, 2) * Power(r, 5)
                                 * (-(x * y * (Power(a, 2) + Power(z, 2)) * dr_dz)
                                    + z
                                          * (2 * a
                                                 * (Power(a, 2) - Power(x, 2) + Power(y, 2)
                                                    + Power(z, 2))
                                             + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                             + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                           / (Power(Power(a, 2) + Power(r, 2), 5)
                              * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                                 + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                 + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                                 + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     + r
                           * ((x * r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                  / (Power(a, 2) + Power(r, 2))
                              + (3 * (a * y + x * r)
                                 * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                    / (Power(a, 2) + Power(r, 2))
                              - (2 * (a * y + x * r)
                                 * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dz)
                                    / Power(Power(a, 2) + Power(r, 2), 2)
                              - (2 * r * (a * y + x * r) * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                                 * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz))
                                    / (Power(a, 2) + Power(r, 2))
                              + (M * Power(r, 3) * Power(a * y + x * r, 2)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (-3 * Power(a, 5) * y * Power(z, 2) * dr_dz
                                    - Power(a, 2) * x * Power(r, 5) * dr_dz
                                    + 5 * a * y * Power(r, 6) * dr_dz + 3 * x * Power(r, 7) * dr_dz
                                    + Power(a, 2) * Power(r, 4) * (2 * x * z + a * y * dr_dz)
                                    + Power(a, 4) * z * r * (2 * a * y - 5 * x * z * dr_dz)
                                    + Power(a, 2) * z * Power(r, 3) * (2 * a * y - x * z * dr_dz)
                                    + Power(a, 3) * z * Power(r, 2) * (2 * a * x + y * z * dr_dz)))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                              - (M * Power(z, 3) * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                                 * (Power(a, 4) * Power(z, 2)
                                    + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                    + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                                    + Power(a, 2) * Power(r, 4) + Power(r, 6))
                                 * dr_dx)
                                    / (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6))
                              + (M * Power(r, 2) * (-(a * x) + y * r)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                                    - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                                    + 2 * Power(a, 5) * z * Power(r, 2)
                                          * (-(a * x * y)
                                             + 2 * (Power(x, 2) - Power(y, 2)) * z * dr_dz
                                             + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                                    + 2 * Power(a, 4) * z * Power(r, 3)
                                          * (-(a * Power(x, 2)) + a * Power(y, 2)
                                             + 2 * a * Power(z, 2) - 3 * x * y * z * dr_dz
                                             + 2 * x * Power(z, 2) * dr_dy
                                             - 2 * y * Power(z, 2) * dr_dx)
                                    + 2 * Power(a, 2) * z * Power(r, 6)
                                          * (-3 * a * y * dr_dy + x * (y - 3 * a * dr_dx))
                                    - 4 * a * Power(r, 8)
                                          * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                             + z * (y * dr_dy + x * dr_dx))
                                    + Power(a, 6) * Power(z, 2) * r
                                          * (3 * x * y * dr_dz
                                             + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                                    + 2 * Power(a, 2) * Power(r, 7)
                                          * (-3 * x * y * dr_dz
                                             + 2 * z * (a - x * dr_dy + y * dr_dx))
                                    + Power(r, 9)
                                          * (3 * x * y * dr_dz
                                             + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                                    + Power(a, 2) * Power(r, 5)
                                          * (-(x * y * (Power(a, 2) + Power(z, 2)) * dr_dz)
                                             + z
                                                   * (2 * a
                                                          * (Power(a, 2) - Power(x, 2) + Power(y, 2)
                                                             + Power(z, 2))
                                                      + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                                      + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))))))
        / (Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 3)
           * Sqrt(
               (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                + Power(a, 2) * Power(z, 2) * Power(r, 2)
                + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                + Power(a, 2) * Power(r, 4) + Power(r, 6))
               * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                  + Power(a, 2) * Power(z, 2) * Power(r, 2)
                  + 2 * M * (-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                  + Power(a, 2) * Power(r, 4) - 2 * M * Power(r, 5) + Power(r, 6))));

  mixedK[1][1]
      = (2 * M * Power(r, 2)
         * (((Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
              + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
              + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
              + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
              + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
              + Power(r, 8))
             * (-4 * Power(r, 4) * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                    * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                + 3 * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                + 2 * (a * x - y * r) * Power(Power(a, 2) + Power(r, 2), 3)
                      * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dy
                + r * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (r + y * dr_dy)
                - (M * Power(r, 3) * Power(a * x - y * r, 2)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (2 * Power(r, 8) - 3 * Power(a, 5) * x * Power(z, 2) * dr_dy
                      + 5 * Power(a, 4) * y * Power(z, 2) * r * dr_dy
                      + Power(a, 2) * y * Power(z, 2) * Power(r, 3) * dr_dy
                      + Power(a, 2) * y * Power(r, 5) * dr_dy - 3 * y * Power(r, 7) * dr_dy
                      + Power(a, 3) * Power(z, 2) * Power(r, 2) * (2 * a + x * dr_dy)
                      + a * Power(r, 6) * (2 * a + 5 * x * dr_dy)
                      + Power(a, 2) * Power(r, 4) * (2 * Power(z, 2) + a * x * dr_dy)))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                - (M * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                   * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                      + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                      + Power(a, 2) * Power(r, 4) + Power(r, 6))
                   * (2 * z * Power(r, 10) - 4 * Power(a, 7) * x * Power(z, 3) * dr_dy
                      + 2 * Power(a, 4) * z * Power(r, 4)
                            * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2) + 2 * a * x * dr_dy)
                      + 2 * Power(a, 2) * z * Power(r, 6)
                            * (Power(a, 2) + Power(y, 2) + Power(z, 2) + 6 * a * x * dr_dy)
                      + 3 * y * Power(r, 9) * (y * dr_dz - 2 * z * dr_dy)
                      + Power(a, 2) * Power(r, 7)
                            * ((5 * Power(x, 2) - Power(y, 2)) * dr_dz - 8 * y * z * dr_dy)
                      + 3 * Power(a, 6) * Power(z, 2) * r
                            * (-(Power(x, 2) * dr_dz) + 2 * y * z * dr_dy)
                      + 2 * Power(a, 5) * z * Power(r, 2)
                            * (a * (Power(x, 2) + Power(z, 2)) + 4 * x * y * z * dr_dz
                               - 2 * x * Power(z, 2) * dr_dy)
                      + Power(a, 4) * z * Power(r, 3)
                            * (-4 * a * x * y + (Power(x, 2) - 5 * Power(y, 2)) * z * dr_dz
                               + 8 * y * Power(z, 2) * dr_dy)
                      + 4 * a * Power(r, 8) * (-2 * x * y * dr_dz + z * (a + 2 * x * dr_dy))
                      + Power(a, 2) * Power(r, 5)
                            * ((Power(a, 2) * Power(x, 2) - Power(y, 2) * Power(z, 2)) * dr_dz
                               + 2 * y * z * (-2 * a * x + (-Power(a, 2) + Power(z, 2)) * dr_dy))))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                - (M * Power(r, 3) * (a * y + x * r)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (-2 * Power(a, 4) * x * Power(z, 2) * Power(r, 3)
                      + 2 * Power(a, 2) * x * (-2 * Power(a, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * x * Power(r, 9)
                      - 3 * Power(a, 6) * x * Power(z, 2) * (2 * y * dr_dy + x * dr_dx)
                      + 3 * y * Power(r, 8) * (2 * a - 2 * x * dr_dy + y * dr_dx)
                      + Power(a, 4) * Power(z, 2) * Power(r, 2)
                            * (6 * a * y + 12 * x * y * dr_dy
                               + (Power(x, 2) - 5 * Power(y, 2)) * dr_dx)
                      + Power(a, 2) * Power(r, 6)
                            * (6 * a * y + 12 * x * y * dr_dy
                               + (5 * Power(x, 2) - Power(y, 2)) * dr_dx)
                      + Power(a, 2) * Power(r, 4)
                            * (6 * a * y * Power(z, 2)
                               + 2 * x * y * (Power(a, 2) + Power(z, 2)) * dr_dy
                               + (Power(a, 2) * Power(x, 2) - Power(y, 2) * Power(z, 2)) * dr_dx)
                      - 4 * Power(a, 5) * Power(z, 2) * r
                            * (2 * (Power(x, 2) - Power(y, 2)) * dr_dy + x * (a - 2 * y * dr_dx))
                      - 2 * a * Power(r, 7)
                            * (-4 * (Power(x, 2) - Power(y, 2)) * dr_dy + x * (a + 4 * y * dr_dx))))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))))
                / Power(Power(a, 2) + Power(r, 2), 6)
            + (M * Power(r, 3)
               * (Power(a, 2) * x * y + a * (Power(x, 2) - Power(y, 2)) * r - x * y * Power(r, 2))
               * (-4 * Power(r, 4) * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                  + 3 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                  - 2 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 3)
                        * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dy
                  + r * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (a + x * dr_dy)
                  - (2 * M * Power(r, 3) * Power(a * y + x * r, 2)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (3 * Power(a, 5) * y * Power(z, 2) * dr_dy
                        - Power(a, 3) * y * Power(z, 2) * Power(r, 2) * dr_dy
                        - Power(a, 3) * y * Power(r, 4) * dr_dy - 5 * a * y * Power(r, 6) * dr_dy
                        + Power(r, 7) * (2 * a - 3 * x * dr_dy)
                        + Power(a, 2) * Power(z, 2) * Power(r, 3) * (2 * a + x * dr_dy)
                        + Power(a, 2) * Power(r, 5) * (2 * a + x * dr_dy)
                        + Power(a, 4) * Power(z, 2) * r * (2 * a + 5 * x * dr_dy)))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  - 4 * Power(r, 4) * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                        * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                  + 3 * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                  + 2 * (a * x - y * r) * Power(Power(a, 2) + Power(r, 2), 3)
                        * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dx
                  + r * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (-a + y * dr_dx)
                  + (2 * M * Power(r, 3) * Power(a * x - y * r, 2)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (3 * Power(a, 5) * x * Power(z, 2) * dr_dx
                        - Power(a, 3) * x * Power(z, 2) * Power(r, 2) * dr_dx
                        - Power(a, 3) * x * Power(r, 4) * dr_dx - 5 * a * x * Power(r, 6) * dr_dx
                        + Power(a, 4) * Power(z, 2) * r * (2 * a - 5 * y * dr_dx)
                        + Power(a, 2) * Power(z, 2) * Power(r, 3) * (2 * a - y * dr_dx)
                        + Power(a, 2) * Power(r, 5) * (2 * a - y * dr_dx)
                        + Power(r, 7) * (2 * a + 3 * y * dr_dx)))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  + (2 * M * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                     * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                        + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                        + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     * (2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy - x * dr_dx)
                        + 2 * Power(a, 7) * Power(z, 3) * (-(y * dr_dy) + x * dr_dx)
                        - 3 * Power(a, 6) * Power(z, 2) * r
                              * (x * y * dr_dz + x * z * dr_dy + y * z * dr_dx)
                        + 2 * Power(a, 5) * z * Power(r, 2)
                              * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                 - y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                        + 2 * Power(a, 4) * z * Power(r, 3)
                              * (a * Power(x, 2) - a * Power(y, 2) + 3 * x * y * z * dr_dz
                                 - 2 * x * Power(z, 2) * dr_dy - 2 * y * Power(z, 2) * dr_dx)
                        - 2 * Power(a, 2) * z * Power(r, 6)
                              * (-3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                        + 4 * a * Power(r, 8)
                              * ((Power(x, 2) - Power(y, 2)) * dr_dz + z * (y * dr_dy - x * dr_dx))
                        + 2 * Power(a, 2) * Power(r, 7)
                              * (3 * x * y * dr_dz + 2 * z * (x * dr_dy + y * dr_dx))
                        + Power(r, 9) * (-3 * x * y * dr_dz + 3 * z * (x * dr_dy + y * dr_dx))
                        + Power(a, 2) * Power(r, 5)
                              * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                 + z
                                       * (2 * a * (Power(x, 2) - Power(y, 2))
                                          + x * (Power(a, 2) - Power(z, 2)) * dr_dy
                                          + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))))
                  / Power(Power(a, 2) + Power(r, 2), 6)
            + M * z * r * (a * x - y * r)
                  * ((M * Power(r, 4) * Power(a * x - y * r, 2)
                      * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                         + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                         + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                         + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                         + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                      * (3 * Power(a, 5) * x * Power(z, 2) * dr_dz
                         - Power(a, 2) * y * Power(r, 5) * dr_dz - 5 * a * x * Power(r, 6) * dr_dz
                         + 3 * y * Power(r, 7) * dr_dz
                         + Power(a, 2) * Power(r, 4) * (2 * y * z - a * x * dr_dz)
                         + Power(a, 3) * z * Power(r, 2) * (2 * a * y - x * z * dr_dz)
                         - Power(a, 2) * z * Power(r, 3) * (2 * a * x + y * z * dr_dz)
                         - Power(a, 4) * z * r * (2 * a * x + 5 * y * z * dr_dz)))
                         / (Power(Power(a, 2) + Power(r, 2), 5)
                            * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                               + Power(a, 2) * Power(z, 2) * Power(r, 2)
                               + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                               + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     - 4 * z * Power(r, 4) * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                     + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                     - (M * Power(z, 3) * r * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                        * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                        * dr_dy)
                           / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                              + Power(a, 2) * Power(z, 2) * Power(r, 2)
                              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                              + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     - (M * Power(r, 3) * (a * y + x * r)
                        * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                           + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                           + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                           + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                           + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                           + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                        * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                           - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                           + 2 * Power(a, 5) * z * Power(r, 2)
                                 * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                    + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                           + 2 * Power(a, 4) * z * Power(r, 3)
                                 * (a * Power(x, 2) - a * Power(y, 2) + 2 * a * Power(z, 2)
                                    + 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                    - 2 * y * Power(z, 2) * dr_dx)
                           - 2 * Power(a, 2) * z * Power(r, 6)
                                 * (3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                           + 4 * a * Power(r, 8)
                                 * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                    - z * (y * dr_dy + x * dr_dx))
                           + Power(a, 6) * Power(z, 2) * r
                                 * (-3 * x * y * dr_dz
                                    + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                           + 2 * Power(a, 2) * Power(r, 7)
                                 * (3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                           + Power(r, 9)
                                 * (-3 * x * y * dr_dz
                                    + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                           + Power(a, 2) * Power(r, 5)
                                 * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                    + z
                                          * (2 * a
                                                 * (Power(a, 2) + Power(x, 2) - Power(y, 2)
                                                    + Power(z, 2))
                                             + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                             + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                           / (Power(Power(a, 2) + Power(r, 2), 5)
                              * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                                 + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                 + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                                 + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     + r
                           * ((y * r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                  / (Power(a, 2) + Power(r, 2))
                              + (3 * (-(a * x) + y * r)
                                 * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                    / (Power(a, 2) + Power(r, 2))
                              + (2 * (a * x - y * r)
                                 * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dz)
                                    / Power(Power(a, 2) + Power(r, 2), 2)
                              - (2 * r * (-(a * x) + y * r)
                                 * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                                 * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz))
                                    / (Power(a, 2) + Power(r, 2))
                              + (M * Power(r, 3) * Power(a * x - y * r, 2)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (3 * Power(a, 5) * x * Power(z, 2) * dr_dz
                                    - Power(a, 2) * y * Power(r, 5) * dr_dz
                                    - 5 * a * x * Power(r, 6) * dr_dz + 3 * y * Power(r, 7) * dr_dz
                                    + Power(a, 2) * Power(r, 4) * (2 * y * z - a * x * dr_dz)
                                    + Power(a, 3) * z * Power(r, 2) * (2 * a * y - x * z * dr_dz)
                                    - Power(a, 2) * z * Power(r, 3) * (2 * a * x + y * z * dr_dz)
                                    - Power(a, 4) * z * r * (2 * a * x + 5 * y * z * dr_dz)))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                              - (M * Power(z, 3) * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                                 * (Power(a, 4) * Power(z, 2)
                                    + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                    + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                                    + Power(a, 2) * Power(r, 4) + Power(r, 6))
                                 * dr_dy)
                                    / (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6))
                              - (M * Power(r, 2) * (a * y + x * r)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                                    - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                                    + 2 * Power(a, 5) * z * Power(r, 2)
                                          * (a * x * y
                                             + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                             + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                                    + 2 * Power(a, 4) * z * Power(r, 3)
                                          * (a * Power(x, 2) - a * Power(y, 2) + 2 * a * Power(z, 2)
                                             + 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                             - 2 * y * Power(z, 2) * dr_dx)
                                    - 2 * Power(a, 2) * z * Power(r, 6)
                                          * (3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                                    + 4 * a * Power(r, 8)
                                          * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                             - z * (y * dr_dy + x * dr_dx))
                                    + Power(a, 6) * Power(z, 2) * r
                                          * (-3 * x * y * dr_dz
                                             + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                                    + 2 * Power(a, 2) * Power(r, 7)
                                          * (3 * x * y * dr_dz
                                             + 2 * z * (a - x * dr_dy + y * dr_dx))
                                    + Power(r, 9)
                                          * (-3 * x * y * dr_dz
                                             + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                                    + Power(a, 2) * Power(r, 5)
                                          * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                             + z
                                                   * (2 * a
                                                          * (Power(a, 2) + Power(x, 2) - Power(y, 2)
                                                             + Power(z, 2))
                                                      + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                                      + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))))))
        / (Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 3)
           * Sqrt(
               (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                + Power(a, 2) * Power(z, 2) * Power(r, 2)
                + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                + Power(a, 2) * Power(r, 4) + Power(r, 6))
               * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                  + Power(a, 2) * Power(z, 2) * Power(r, 2)
                  + 2 * M * (-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                  + Power(a, 2) * Power(r, 4) - 2 * M * Power(r, 5) + Power(r, 6))));

  mixedK[1][2]
      = (M * r
         * (4 * M * z * Power(r, 2) * (a * x - y * r)
                * (r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
                   + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz
                   - 2 * z * r * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                         * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz)
                   - (M * Power(z, 2) * r
                      * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                      * (2 * Power(r, 5) + Power(a, 2) * Power(z, 3) * dr_dz
                         - 3 * z * Power(r, 4) * dr_dz))
                         / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                            + Power(a, 2) * Power(z, 2) * Power(r, 2)
                            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                            + Power(a, 2) * Power(r, 4) + Power(r, 6))
                   + (M * Power(r, 2) * (-(a * x) + y * r)
                      * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                         + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                         + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                         + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                         + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                      * (2 * Power(a, 4) * y * Power(z, 2) * Power(r, 3) - 2 * y * Power(r, 9)
                         + 4 * Power(a, 5) * x * Power(z, 3) * r * dr_dz
                         + 2 * Power(a, 2) * z * Power(r, 5) * (y * z - 2 * a * x * dr_dz)
                         - 2 * a * Power(r, 7) * (a * y + 4 * x * z * dr_dz)
                         + Power(a, 6) * Power(z, 4) * dr_dy
                         + 2 * Power(a, 2) * Power(r, 6)
                               * (a * x + y * z * dr_dz - 3 * Power(z, 2) * dr_dy)
                         + Power(r, 8) * (2 * a * x + 6 * y * z * dr_dz - 3 * Power(z, 2) * dr_dy)
                         - 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                               * (a * x + 3 * y * z * dr_dz - Power(z, 2) * dr_dy)
                         + Power(a, 2) * Power(z, 2) * Power(r, 4)
                               * (-2 * a * x - 2 * y * z * dr_dz
                                  + (-3 * Power(a, 2) + Power(z, 2)) * dr_dy)))
                         / (Power(Power(a, 2) + Power(r, 2), 4)
                            * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                               + Power(a, 2) * Power(z, 2) * Power(r, 2)
                               + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                               + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                   + (M * Power(r, 2) * (a * y + x * r)
                      * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                         + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                         + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                         + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                         + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                      * (2 * Power(a, 4) * x * Power(z, 2) * Power(r, 3) - 2 * x * Power(r, 9)
                         - 4 * Power(a, 5) * y * Power(z, 3) * r * dr_dz
                         + 2 * Power(a, 2) * z * Power(r, 5) * (x * z + 2 * a * y * dr_dz)
                         - 2 * a * Power(r, 7) * (a * x - 4 * y * z * dr_dz)
                         + Power(a, 6) * Power(z, 4) * dr_dx
                         + Power(r, 8) * (-2 * a * y + 6 * x * z * dr_dz - 3 * Power(z, 2) * dr_dx)
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                               * (a * y - 3 * x * z * dr_dz + Power(z, 2) * dr_dx)
                         - 2 * Power(a, 2) * Power(r, 6)
                               * (a * y - x * z * dr_dz + 3 * Power(z, 2) * dr_dx)
                         + Power(a, 2) * Power(z, 2) * Power(r, 4)
                               * (2 * a * y - 2 * x * z * dr_dz
                                  + (-3 * Power(a, 2) + Power(z, 2)) * dr_dx)))
                         / (Power(Power(a, 2) + Power(r, 2), 4)
                            * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                               + Power(a, 2) * Power(z, 2) * Power(r, 2)
                               + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                               + Power(a, 2) * Power(r, 4) + Power(r, 6))))
            + ((Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
                + Power(r, 8))
               * ((M * Power(r, 4) * Power(a * x - y * r, 2)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (3 * Power(a, 5) * x * Power(z, 2) * dr_dz
                      - Power(a, 2) * y * Power(r, 5) * dr_dz - 5 * a * x * Power(r, 6) * dr_dz
                      + 3 * y * Power(r, 7) * dr_dz
                      + Power(a, 2) * Power(r, 4) * (2 * y * z - a * x * dr_dz)
                      + Power(a, 3) * z * Power(r, 2) * (2 * a * y - x * z * dr_dz)
                      - Power(a, 2) * z * Power(r, 3) * (2 * a * x + y * z * dr_dz)
                      - Power(a, 4) * z * r * (2 * a * x + 5 * y * z * dr_dz)))
                      / (Power(Power(a, 2) + Power(r, 2), 5)
                         * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                            + Power(a, 2) * Power(z, 2) * Power(r, 2)
                            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                            + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                  - 4 * z * Power(r, 4) * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                  + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                  - (M * Power(z, 3) * r * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                     * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                        + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                        + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     * dr_dy)
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  - (M * Power(r, 3) * (a * y + x * r)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                        - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                        + 2 * Power(a, 5) * z * Power(r, 2)
                              * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                 + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                        + 2 * Power(a, 4) * z * Power(r, 3)
                              * (a * Power(x, 2) - a * Power(y, 2) + 2 * a * Power(z, 2)
                                 + 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                 - 2 * y * Power(z, 2) * dr_dx)
                        - 2 * Power(a, 2) * z * Power(r, 6)
                              * (3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                        + 4 * a * Power(r, 8)
                              * ((Power(x, 2) - Power(y, 2)) * dr_dz - z * (y * dr_dy + x * dr_dx))
                        + Power(a, 6) * Power(z, 2) * r
                              * (-3 * x * y * dr_dz + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                        + 2 * Power(a, 2) * Power(r, 7)
                              * (3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                        + Power(r, 9)
                              * (-3 * x * y * dr_dz + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                        + Power(a, 2) * Power(r, 5)
                              * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                 + z
                                       * (2 * a
                                              * (Power(a, 2) + Power(x, 2) - Power(y, 2)
                                                 + Power(z, 2))
                                          + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                          + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                        / (Power(Power(a, 2) + Power(r, 2), 5)
                           * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                              + Power(a, 2) * Power(z, 2) * Power(r, 2)
                              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                              + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                  + r
                        * ((y * r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                               / (Power(a, 2) + Power(r, 2))
                           + (3 * (-(a * x) + y * r)
                              * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                 / (Power(a, 2) + Power(r, 2))
                           + (2 * (a * x - y * r)
                              * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dz)
                                 / Power(Power(a, 2) + Power(r, 2), 2)
                           - (2 * r * (-(a * x) + y * r) * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                              * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz))
                                 / (Power(a, 2) + Power(r, 2))
                           + (M * Power(r, 3) * Power(a * x - y * r, 2)
                              * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                                 + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                 + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2))
                                       * Power(r, 3)
                                 + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                                 + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                                 + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                              * (3 * Power(a, 5) * x * Power(z, 2) * dr_dz
                                 - Power(a, 2) * y * Power(r, 5) * dr_dz
                                 - 5 * a * x * Power(r, 6) * dr_dz + 3 * y * Power(r, 7) * dr_dz
                                 + Power(a, 2) * Power(r, 4) * (2 * y * z - a * x * dr_dz)
                                 + Power(a, 3) * z * Power(r, 2) * (2 * a * y - x * z * dr_dz)
                                 - Power(a, 2) * z * Power(r, 3) * (2 * a * x + y * z * dr_dz)
                                 - Power(a, 4) * z * r * (2 * a * x + 5 * y * z * dr_dz)))
                                 / (Power(Power(a, 2) + Power(r, 2), 5)
                                    * (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                           - (M * Power(z, 3) * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                              * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                 + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                                 + Power(a, 2) * Power(r, 4) + Power(r, 6))
                              * dr_dy)
                                 / (Power(a, 4) * Power(z, 2)
                                    + 2 * Power(a, 2) * M * Power(z, 2) * r
                                    + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                    + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                          * Power(r, 3)
                                    + Power(a, 2) * Power(r, 4) + Power(r, 6))
                           - (M * Power(r, 2) * (a * y + x * r)
                              * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                                 + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                 + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2))
                                       * Power(r, 3)
                                 + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                                 + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                                 + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                              * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                                 - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                                 + 2 * Power(a, 5) * z * Power(r, 2)
                                       * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                          + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                                 + 2 * Power(a, 4) * z * Power(r, 3)
                                       * (a * Power(x, 2) - a * Power(y, 2) + 2 * a * Power(z, 2)
                                          + 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                          - 2 * y * Power(z, 2) * dr_dx)
                                 - 2 * Power(a, 2) * z * Power(r, 6)
                                       * (3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                                 + 4 * a * Power(r, 8)
                                       * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                          - z * (y * dr_dy + x * dr_dx))
                                 + Power(a, 6) * Power(z, 2) * r
                                       * (-3 * x * y * dr_dz
                                          + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                                 + 2 * Power(a, 2) * Power(r, 7)
                                       * (3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                                 + Power(r, 9)
                                       * (-3 * x * y * dr_dz
                                          + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                                 + Power(a, 2) * Power(r, 5)
                                       * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                          + z
                                                * (2 * a
                                                       * (Power(a, 2) + Power(x, 2) - Power(y, 2)
                                                          + Power(z, 2))
                                                   + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                                   + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                                 / (Power(Power(a, 2) + Power(r, 2), 5)
                                    * (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6))))))
                  / (Power(a, 2) + Power(r, 2))
            + (2 * M * Power(r, 3)
               * (Power(a, 2) * x * y + a * (Power(x, 2) - Power(y, 2)) * r - x * y * Power(r, 2))
               * ((M * Power(r, 4) * Power(a * y + x * r, 2)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (-3 * Power(a, 5) * y * Power(z, 2) * dr_dz
                      - Power(a, 2) * x * Power(r, 5) * dr_dz + 5 * a * y * Power(r, 6) * dr_dz
                      + 3 * x * Power(r, 7) * dr_dz
                      + Power(a, 2) * Power(r, 4) * (2 * x * z + a * y * dr_dz)
                      + Power(a, 4) * z * r * (2 * a * y - 5 * x * z * dr_dz)
                      + Power(a, 2) * z * Power(r, 3) * (2 * a * y - x * z * dr_dz)
                      + Power(a, 3) * z * Power(r, 2) * (2 * a * x + y * z * dr_dz)))
                      / (Power(Power(a, 2) + Power(r, 2), 5)
                         * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                            + Power(a, 2) * Power(z, 2) * Power(r, 2)
                            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                            + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                  - 4 * z * Power(r, 4) * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                  + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                  - (M * Power(z, 3) * r * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                     * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                        + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                        + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     * dr_dx)
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  + (M * Power(r, 3) * (-(a * x) + y * r)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                        - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                        + 2 * Power(a, 5) * z * Power(r, 2)
                              * (-(a * x * y) + 2 * (Power(x, 2) - Power(y, 2)) * z * dr_dz
                                 + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                        + 2 * Power(a, 4) * z * Power(r, 3)
                              * (-(a * Power(x, 2)) + a * Power(y, 2) + 2 * a * Power(z, 2)
                                 - 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                 - 2 * y * Power(z, 2) * dr_dx)
                        + 2 * Power(a, 2) * z * Power(r, 6)
                              * (-3 * a * y * dr_dy + x * (y - 3 * a * dr_dx))
                        - 4 * a * Power(r, 8)
                              * ((Power(x, 2) - Power(y, 2)) * dr_dz + z * (y * dr_dy + x * dr_dx))
                        + Power(a, 6) * Power(z, 2) * r
                              * (3 * x * y * dr_dz + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                        + 2 * Power(a, 2) * Power(r, 7)
                              * (-3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                        + Power(r, 9)
                              * (3 * x * y * dr_dz + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                        + Power(a, 2) * Power(r, 5)
                              * (-(x * y * (Power(a, 2) + Power(z, 2)) * dr_dz)
                                 + z
                                       * (2 * a
                                              * (Power(a, 2) - Power(x, 2) + Power(y, 2)
                                                 + Power(z, 2))
                                          + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                          + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                        / (Power(Power(a, 2) + Power(r, 2), 5)
                           * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                              + Power(a, 2) * Power(z, 2) * Power(r, 2)
                              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                              + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                  + r
                        * ((x * r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                               / (Power(a, 2) + Power(r, 2))
                           + (3 * (a * y + x * r)
                              * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                 / (Power(a, 2) + Power(r, 2))
                           - (2 * (a * y + x * r)
                              * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dz)
                                 / Power(Power(a, 2) + Power(r, 2), 2)
                           - (2 * r * (a * y + x * r) * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                              * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz))
                                 / (Power(a, 2) + Power(r, 2))
                           + (M * Power(r, 3) * Power(a * y + x * r, 2)
                              * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                                 + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                 + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2))
                                       * Power(r, 3)
                                 + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                                 + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                                 + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                              * (-3 * Power(a, 5) * y * Power(z, 2) * dr_dz
                                 - Power(a, 2) * x * Power(r, 5) * dr_dz
                                 + 5 * a * y * Power(r, 6) * dr_dz + 3 * x * Power(r, 7) * dr_dz
                                 + Power(a, 2) * Power(r, 4) * (2 * x * z + a * y * dr_dz)
                                 + Power(a, 4) * z * r * (2 * a * y - 5 * x * z * dr_dz)
                                 + Power(a, 2) * z * Power(r, 3) * (2 * a * y - x * z * dr_dz)
                                 + Power(a, 3) * z * Power(r, 2) * (2 * a * x + y * z * dr_dz)))
                                 / (Power(Power(a, 2) + Power(r, 2), 5)
                                    * (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                           - (M * Power(z, 3) * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                              * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                 + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                                 + Power(a, 2) * Power(r, 4) + Power(r, 6))
                              * dr_dx)
                                 / (Power(a, 4) * Power(z, 2)
                                    + 2 * Power(a, 2) * M * Power(z, 2) * r
                                    + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                    + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                          * Power(r, 3)
                                    + Power(a, 2) * Power(r, 4) + Power(r, 6))
                           + (M * Power(r, 2) * (-(a * x) + y * r)
                              * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                                 + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                 + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2))
                                       * Power(r, 3)
                                 + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                                 + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                                 + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                              * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                                 - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                                 + 2 * Power(a, 5) * z * Power(r, 2)
                                       * (-(a * x * y) + 2 * (Power(x, 2) - Power(y, 2)) * z * dr_dz
                                          + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                                 + 2 * Power(a, 4) * z * Power(r, 3)
                                       * (-(a * Power(x, 2)) + a * Power(y, 2) + 2 * a * Power(z, 2)
                                          - 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                          - 2 * y * Power(z, 2) * dr_dx)
                                 + 2 * Power(a, 2) * z * Power(r, 6)
                                       * (-3 * a * y * dr_dy + x * (y - 3 * a * dr_dx))
                                 - 4 * a * Power(r, 8)
                                       * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                          + z * (y * dr_dy + x * dr_dx))
                                 + Power(a, 6) * Power(z, 2) * r
                                       * (3 * x * y * dr_dz
                                          + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                                 + 2 * Power(a, 2) * Power(r, 7)
                                       * (-3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                                 + Power(r, 9)
                                       * (3 * x * y * dr_dz
                                          + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                                 + Power(a, 2) * Power(r, 5)
                                       * (-(x * y * (Power(a, 2) + Power(z, 2)) * dr_dz)
                                          + z
                                                * (2 * a
                                                       * (Power(a, 2) - Power(x, 2) + Power(y, 2)
                                                          + Power(z, 2))
                                                   + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                                   + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                                 / (Power(Power(a, 2) + Power(r, 2), 5)
                                    * (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6))))))
                  / (Power(a, 2) + Power(r, 2))))
        / (Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 3)
           * Sqrt(
               (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                + Power(a, 2) * Power(z, 2) * Power(r, 2)
                + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                + Power(a, 2) * Power(r, 4) + Power(r, 6))
               * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                  + Power(a, 2) * Power(z, 2) * Power(r, 2)
                  + 2 * M * (-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                  + Power(a, 2) * Power(r, 4) - 2 * M * Power(r, 5) + Power(r, 6))));

  mixedK[2][0]
      = (M * r
         * ((2 * M * z * Power(r, 3) * (a * x - y * r)
             * (-4 * Power(r, 4) * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                    * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                + 3 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                - 2 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 3)
                      * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dy
                + r * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (a + x * dr_dy)
                - (2 * M * Power(r, 3) * Power(a * y + x * r, 2)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (3 * Power(a, 5) * y * Power(z, 2) * dr_dy
                      - Power(a, 3) * y * Power(z, 2) * Power(r, 2) * dr_dy
                      - Power(a, 3) * y * Power(r, 4) * dr_dy - 5 * a * y * Power(r, 6) * dr_dy
                      + Power(r, 7) * (2 * a - 3 * x * dr_dy)
                      + Power(a, 2) * Power(z, 2) * Power(r, 3) * (2 * a + x * dr_dy)
                      + Power(a, 2) * Power(r, 5) * (2 * a + x * dr_dy)
                      + Power(a, 4) * Power(z, 2) * r * (2 * a + 5 * x * dr_dy)))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                - 4 * Power(r, 4) * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                + 3 * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                + 2 * (a * x - y * r) * Power(Power(a, 2) + Power(r, 2), 3)
                      * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dx
                + r * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (-a + y * dr_dx)
                + (2 * M * Power(r, 3) * Power(a * x - y * r, 2)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (3 * Power(a, 5) * x * Power(z, 2) * dr_dx
                      - Power(a, 3) * x * Power(z, 2) * Power(r, 2) * dr_dx
                      - Power(a, 3) * x * Power(r, 4) * dr_dx - 5 * a * x * Power(r, 6) * dr_dx
                      + Power(a, 4) * Power(z, 2) * r * (2 * a - 5 * y * dr_dx)
                      + Power(a, 2) * Power(z, 2) * Power(r, 3) * (2 * a - y * dr_dx)
                      + Power(a, 2) * Power(r, 5) * (2 * a - y * dr_dx)
                      + Power(r, 7) * (2 * a + 3 * y * dr_dx)))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                + (2 * M * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                   * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                      + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                      + Power(a, 2) * Power(r, 4) + Power(r, 6))
                   * (2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy - x * dr_dx)
                      + 2 * Power(a, 7) * Power(z, 3) * (-(y * dr_dy) + x * dr_dx)
                      - 3 * Power(a, 6) * Power(z, 2) * r
                            * (x * y * dr_dz + x * z * dr_dy + y * z * dr_dx)
                      + 2 * Power(a, 5) * z * Power(r, 2)
                            * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                               - y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                      + 2 * Power(a, 4) * z * Power(r, 3)
                            * (a * Power(x, 2) - a * Power(y, 2) + 3 * x * y * z * dr_dz
                               - 2 * x * Power(z, 2) * dr_dy - 2 * y * Power(z, 2) * dr_dx)
                      - 2 * Power(a, 2) * z * Power(r, 6)
                            * (-3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                      + 4 * a * Power(r, 8)
                            * ((Power(x, 2) - Power(y, 2)) * dr_dz + z * (y * dr_dy - x * dr_dx))
                      + 2 * Power(a, 2) * Power(r, 7)
                            * (3 * x * y * dr_dz + 2 * z * (x * dr_dy + y * dr_dx))
                      + Power(r, 9) * (-3 * x * y * dr_dz + 3 * z * (x * dr_dy + y * dr_dx))
                      + Power(a, 2) * Power(r, 5)
                            * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                               + z
                                     * (2 * a * (Power(x, 2) - Power(y, 2))
                                        + x * (Power(a, 2) - Power(z, 2)) * dr_dy
                                        + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))))
                / Power(Power(a, 2) + Power(r, 2), 5)
            - (4 * M * z * Power(r, 3) * (a * y + x * r)
               * (-4 * Power(r, 4) * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                  + 3 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                  - 2 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 3)
                        * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dx
                  + r * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (r + x * dr_dx)
                  - (M * Power(r, 3) * Power(a * y + x * r, 2)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (2 * Power(r, 8) + 3 * Power(a, 5) * y * Power(z, 2) * dr_dx
                        + 5 * Power(a, 4) * x * Power(z, 2) * r * dr_dx
                        + Power(a, 2) * x * Power(z, 2) * Power(r, 3) * dr_dx
                        + Power(a, 2) * x * Power(r, 5) * dr_dx - 3 * x * Power(r, 7) * dr_dx
                        + a * Power(r, 6) * (2 * a - 5 * y * dr_dx)
                        + Power(a, 3) * Power(z, 2) * Power(r, 2) * (2 * a - y * dr_dx)
                        + Power(a, 2) * Power(r, 4) * (2 * Power(z, 2) - a * y * dr_dx)))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  + (M * Power(r, 3) * (-(a * x) + y * r)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (2 * Power(a, 4) * y * Power(z, 2) * Power(r, 3)
                        + (4 * Power(a, 4) * y - 2 * Power(a, 2) * y * Power(z, 2)) * Power(r, 5)
                        - 2 * y * Power(r, 9)
                        + 3 * Power(a, 6) * y * Power(z, 2) * (y * dr_dy + 2 * x * dr_dx)
                        + 2 * a * Power(r, 7)
                              * (a * y - 4 * x * y * dr_dy
                                 - 4 * (Power(x, 2) - Power(y, 2)) * dr_dx)
                        + 4 * Power(a, 5) * Power(z, 2) * r
                              * (a * y + 2 * x * y * dr_dy
                                 + 2 * (Power(x, 2) - Power(y, 2)) * dr_dx)
                        + Power(a, 2) * Power(r, 4)
                              * (6 * a * x * Power(z, 2)
                                 + (-(Power(a, 2) * Power(y, 2)) + Power(x, 2) * Power(z, 2))
                                       * dr_dy
                                 - 2 * x * y * (Power(a, 2) + Power(z, 2)) * dr_dx)
                        + Power(a, 2) * Power(r, 6)
                              * ((Power(x, 2) - 5 * Power(y, 2)) * dr_dy
                                 + 6 * x * (a - 2 * y * dr_dx))
                        + Power(a, 4) * Power(z, 2) * Power(r, 2)
                              * ((5 * Power(x, 2) - Power(y, 2)) * dr_dy
                                 + 6 * x * (a - 2 * y * dr_dx))
                        + 3 * x * Power(r, 8) * (-(x * dr_dy) + 2 * (a + y * dr_dx))))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  - (M * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                     * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                        + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                        + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     * (2 * z * Power(r, 10) + 4 * Power(a, 7) * y * Power(z, 3) * dr_dx
                        + 2 * Power(a, 2) * z * Power(r, 6)
                              * (Power(a, 2) + Power(x, 2) + Power(z, 2) - 6 * a * y * dr_dx)
                        + 2 * Power(a, 4) * z * Power(r, 4)
                              * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2) - 2 * a * y * dr_dx)
                        + 3 * x * Power(r, 9) * (x * dr_dz - 2 * z * dr_dx)
                        + 3 * Power(a, 6) * Power(z, 2) * r
                              * (-(Power(y, 2) * dr_dz) + 2 * x * z * dr_dx)
                        - Power(a, 2) * Power(r, 7)
                              * ((Power(x, 2) - 5 * Power(y, 2)) * dr_dz + 8 * x * z * dr_dx)
                        + 2 * Power(a, 5) * z * Power(r, 2)
                              * (a * (Power(y, 2) + Power(z, 2)) - 4 * x * y * z * dr_dz
                                 + 2 * y * Power(z, 2) * dr_dx)
                        + 4 * a * Power(r, 8) * (2 * x * y * dr_dz + z * (a - 2 * y * dr_dx))
                        + Power(a, 4) * z * Power(r, 3)
                              * ((-5 * Power(x, 2) + Power(y, 2)) * z * dr_dz
                                 + 4 * x * (a * y + 2 * Power(z, 2) * dr_dx))
                        + Power(a, 2) * Power(r, 5)
                              * ((Power(a, 2) * Power(y, 2) - Power(x, 2) * Power(z, 2)) * dr_dz
                                 + 2 * x * z * (2 * a * y + (-Power(a, 2) + Power(z, 2)) * dr_dx))))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))))
                  / Power(Power(a, 2) + Power(r, 2), 5)
            + (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
               + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3) + Power(a, 2) * Power(r, 4)
               + Power(r, 6))
                  * ((M * Power(r, 4) * Power(a * y + x * r, 2)
                      * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                         + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                         + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                         + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                         + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                      * (-3 * Power(a, 5) * y * Power(z, 2) * dr_dz
                         - Power(a, 2) * x * Power(r, 5) * dr_dz + 5 * a * y * Power(r, 6) * dr_dz
                         + 3 * x * Power(r, 7) * dr_dz
                         + Power(a, 2) * Power(r, 4) * (2 * x * z + a * y * dr_dz)
                         + Power(a, 4) * z * r * (2 * a * y - 5 * x * z * dr_dz)
                         + Power(a, 2) * z * Power(r, 3) * (2 * a * y - x * z * dr_dz)
                         + Power(a, 3) * z * Power(r, 2) * (2 * a * x + y * z * dr_dz)))
                         / (Power(Power(a, 2) + Power(r, 2), 5)
                            * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                               + Power(a, 2) * Power(z, 2) * Power(r, 2)
                               + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                               + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     - 4 * z * Power(r, 4) * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                     + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                     - (M * Power(z, 3) * r * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                        * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                        * dr_dx)
                           / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                              + Power(a, 2) * Power(z, 2) * Power(r, 2)
                              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                              + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     + (M * Power(r, 3) * (-(a * x) + y * r)
                        * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                           + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                           + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                           + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                           + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                           + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                        * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                           - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                           + 2 * Power(a, 5) * z * Power(r, 2)
                                 * (-(a * x * y) + 2 * (Power(x, 2) - Power(y, 2)) * z * dr_dz
                                    + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                           + 2 * Power(a, 4) * z * Power(r, 3)
                                 * (-(a * Power(x, 2)) + a * Power(y, 2) + 2 * a * Power(z, 2)
                                    - 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                    - 2 * y * Power(z, 2) * dr_dx)
                           + 2 * Power(a, 2) * z * Power(r, 6)
                                 * (-3 * a * y * dr_dy + x * (y - 3 * a * dr_dx))
                           - 4 * a * Power(r, 8)
                                 * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                    + z * (y * dr_dy + x * dr_dx))
                           + Power(a, 6) * Power(z, 2) * r
                                 * (3 * x * y * dr_dz + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                           + 2 * Power(a, 2) * Power(r, 7)
                                 * (-3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                           + Power(r, 9)
                                 * (3 * x * y * dr_dz + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                           + Power(a, 2) * Power(r, 5)
                                 * (-(x * y * (Power(a, 2) + Power(z, 2)) * dr_dz)
                                    + z
                                          * (2 * a
                                                 * (Power(a, 2) - Power(x, 2) + Power(y, 2)
                                                    + Power(z, 2))
                                             + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                             + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                           / (Power(Power(a, 2) + Power(r, 2), 5)
                              * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                                 + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                 + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                                 + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     + r
                           * ((x * r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                  / (Power(a, 2) + Power(r, 2))
                              + (3 * (a * y + x * r)
                                 * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                    / (Power(a, 2) + Power(r, 2))
                              - (2 * (a * y + x * r)
                                 * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dz)
                                    / Power(Power(a, 2) + Power(r, 2), 2)
                              - (2 * r * (a * y + x * r) * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                                 * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz))
                                    / (Power(a, 2) + Power(r, 2))
                              + (M * Power(r, 3) * Power(a * y + x * r, 2)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (-3 * Power(a, 5) * y * Power(z, 2) * dr_dz
                                    - Power(a, 2) * x * Power(r, 5) * dr_dz
                                    + 5 * a * y * Power(r, 6) * dr_dz + 3 * x * Power(r, 7) * dr_dz
                                    + Power(a, 2) * Power(r, 4) * (2 * x * z + a * y * dr_dz)
                                    + Power(a, 4) * z * r * (2 * a * y - 5 * x * z * dr_dz)
                                    + Power(a, 2) * z * Power(r, 3) * (2 * a * y - x * z * dr_dz)
                                    + Power(a, 3) * z * Power(r, 2) * (2 * a * x + y * z * dr_dz)))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                              - (M * Power(z, 3) * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                                 * (Power(a, 4) * Power(z, 2)
                                    + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                    + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                                    + Power(a, 2) * Power(r, 4) + Power(r, 6))
                                 * dr_dx)
                                    / (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6))
                              + (M * Power(r, 2) * (-(a * x) + y * r)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                                    - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                                    + 2 * Power(a, 5) * z * Power(r, 2)
                                          * (-(a * x * y)
                                             + 2 * (Power(x, 2) - Power(y, 2)) * z * dr_dz
                                             + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                                    + 2 * Power(a, 4) * z * Power(r, 3)
                                          * (-(a * Power(x, 2)) + a * Power(y, 2)
                                             + 2 * a * Power(z, 2) - 3 * x * y * z * dr_dz
                                             + 2 * x * Power(z, 2) * dr_dy
                                             - 2 * y * Power(z, 2) * dr_dx)
                                    + 2 * Power(a, 2) * z * Power(r, 6)
                                          * (-3 * a * y * dr_dy + x * (y - 3 * a * dr_dx))
                                    - 4 * a * Power(r, 8)
                                          * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                             + z * (y * dr_dy + x * dr_dx))
                                    + Power(a, 6) * Power(z, 2) * r
                                          * (3 * x * y * dr_dz
                                             + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                                    + 2 * Power(a, 2) * Power(r, 7)
                                          * (-3 * x * y * dr_dz
                                             + 2 * z * (a - x * dr_dy + y * dr_dx))
                                    + Power(r, 9)
                                          * (3 * x * y * dr_dz
                                             + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                                    + Power(a, 2) * Power(r, 5)
                                          * (-(x * y * (Power(a, 2) + Power(z, 2)) * dr_dz)
                                             + z
                                                   * (2 * a
                                                          * (Power(a, 2) - Power(x, 2) + Power(y, 2)
                                                             + Power(z, 2))
                                                      + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                                      + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))))))
        / (Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 3)
           * Sqrt(
               (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                + Power(a, 2) * Power(z, 2) * Power(r, 2)
                + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                + Power(a, 2) * Power(r, 4) + Power(r, 6))
               * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                  + Power(a, 2) * Power(z, 2) * Power(r, 2)
                  + 2 * M * (-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                  + Power(a, 2) * Power(r, 4) - 2 * M * Power(r, 5) + Power(r, 6))));

  mixedK[2][1]
      = (M * r
         * ((4 * M * z * Power(r, 3) * (a * x - y * r)
             * (-4 * Power(r, 4) * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                    * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                + 3 * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                + 2 * (a * x - y * r) * Power(Power(a, 2) + Power(r, 2), 3)
                      * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dy
                + r * Power(Power(a, 2) + Power(r, 2), 4)
                      * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (r + y * dr_dy)
                - (M * Power(r, 3) * Power(a * x - y * r, 2)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (2 * Power(r, 8) - 3 * Power(a, 5) * x * Power(z, 2) * dr_dy
                      + 5 * Power(a, 4) * y * Power(z, 2) * r * dr_dy
                      + Power(a, 2) * y * Power(z, 2) * Power(r, 3) * dr_dy
                      + Power(a, 2) * y * Power(r, 5) * dr_dy - 3 * y * Power(r, 7) * dr_dy
                      + Power(a, 3) * Power(z, 2) * Power(r, 2) * (2 * a + x * dr_dy)
                      + a * Power(r, 6) * (2 * a + 5 * x * dr_dy)
                      + Power(a, 2) * Power(r, 4) * (2 * Power(z, 2) + a * x * dr_dy)))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                - (M * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                   * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                      + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                      + Power(a, 2) * Power(r, 4) + Power(r, 6))
                   * (2 * z * Power(r, 10) - 4 * Power(a, 7) * x * Power(z, 3) * dr_dy
                      + 2 * Power(a, 4) * z * Power(r, 4)
                            * (Power(x, 2) + Power(y, 2) + 2 * Power(z, 2) + 2 * a * x * dr_dy)
                      + 2 * Power(a, 2) * z * Power(r, 6)
                            * (Power(a, 2) + Power(y, 2) + Power(z, 2) + 6 * a * x * dr_dy)
                      + 3 * y * Power(r, 9) * (y * dr_dz - 2 * z * dr_dy)
                      + Power(a, 2) * Power(r, 7)
                            * ((5 * Power(x, 2) - Power(y, 2)) * dr_dz - 8 * y * z * dr_dy)
                      + 3 * Power(a, 6) * Power(z, 2) * r
                            * (-(Power(x, 2) * dr_dz) + 2 * y * z * dr_dy)
                      + 2 * Power(a, 5) * z * Power(r, 2)
                            * (a * (Power(x, 2) + Power(z, 2)) + 4 * x * y * z * dr_dz
                               - 2 * x * Power(z, 2) * dr_dy)
                      + Power(a, 4) * z * Power(r, 3)
                            * (-4 * a * x * y + (Power(x, 2) - 5 * Power(y, 2)) * z * dr_dz
                               + 8 * y * Power(z, 2) * dr_dy)
                      + 4 * a * Power(r, 8) * (-2 * x * y * dr_dz + z * (a + 2 * x * dr_dy))
                      + Power(a, 2) * Power(r, 5)
                            * ((Power(a, 2) * Power(x, 2) - Power(y, 2) * Power(z, 2)) * dr_dz
                               + 2 * y * z * (-2 * a * x + (-Power(a, 2) + Power(z, 2)) * dr_dy))))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                - (M * Power(r, 3) * (a * y + x * r)
                   * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                      + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                      + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                      + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                      + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                   * (-2 * Power(a, 4) * x * Power(z, 2) * Power(r, 3)
                      + 2 * Power(a, 2) * x * (-2 * Power(a, 2) + Power(z, 2)) * Power(r, 5)
                      + 2 * x * Power(r, 9)
                      - 3 * Power(a, 6) * x * Power(z, 2) * (2 * y * dr_dy + x * dr_dx)
                      + 3 * y * Power(r, 8) * (2 * a - 2 * x * dr_dy + y * dr_dx)
                      + Power(a, 4) * Power(z, 2) * Power(r, 2)
                            * (6 * a * y + 12 * x * y * dr_dy
                               + (Power(x, 2) - 5 * Power(y, 2)) * dr_dx)
                      + Power(a, 2) * Power(r, 6)
                            * (6 * a * y + 12 * x * y * dr_dy
                               + (5 * Power(x, 2) - Power(y, 2)) * dr_dx)
                      + Power(a, 2) * Power(r, 4)
                            * (6 * a * y * Power(z, 2)
                               + 2 * x * y * (Power(a, 2) + Power(z, 2)) * dr_dy
                               + (Power(a, 2) * Power(x, 2) - Power(y, 2) * Power(z, 2)) * dr_dx)
                      - 4 * Power(a, 5) * Power(z, 2) * r
                            * (2 * (Power(x, 2) - Power(y, 2)) * dr_dy + x * (a - 2 * y * dr_dx))
                      - 2 * a * Power(r, 7)
                            * (-4 * (Power(x, 2) - Power(y, 2)) * dr_dy + x * (a + 4 * y * dr_dx))))
                      / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                         + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))))
                / Power(Power(a, 2) + Power(r, 2), 5)
            - (2 * M * z * Power(r, 3) * (a * y + x * r)
               * (-4 * Power(r, 4) * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                      * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                  + 3 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                  - 2 * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 3)
                        * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dy
                  + r * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (a + x * dr_dy)
                  - (2 * M * Power(r, 3) * Power(a * y + x * r, 2)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (3 * Power(a, 5) * y * Power(z, 2) * dr_dy
                        - Power(a, 3) * y * Power(z, 2) * Power(r, 2) * dr_dy
                        - Power(a, 3) * y * Power(r, 4) * dr_dy - 5 * a * y * Power(r, 6) * dr_dy
                        + Power(r, 7) * (2 * a - 3 * x * dr_dy)
                        + Power(a, 2) * Power(z, 2) * Power(r, 3) * (2 * a + x * dr_dy)
                        + Power(a, 2) * Power(r, 5) * (2 * a + x * dr_dy)
                        + Power(a, 4) * Power(z, 2) * r * (2 * a + 5 * x * dr_dy)))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  - 4 * Power(r, 4) * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                        * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                  + 3 * (-(a * x) + y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                  + 2 * (a * x - y * r) * Power(Power(a, 2) + Power(r, 2), 3)
                        * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dx
                  + r * Power(Power(a, 2) + Power(r, 2), 4)
                        * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * (-a + y * dr_dx)
                  + (2 * M * Power(r, 3) * Power(a * x - y * r, 2)
                     * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                        + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                        + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                        + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                        + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                        + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                     * (3 * Power(a, 5) * x * Power(z, 2) * dr_dx
                        - Power(a, 3) * x * Power(z, 2) * Power(r, 2) * dr_dx
                        - Power(a, 3) * x * Power(r, 4) * dr_dx - 5 * a * x * Power(r, 6) * dr_dx
                        + Power(a, 4) * Power(z, 2) * r * (2 * a - 5 * y * dr_dx)
                        + Power(a, 2) * Power(z, 2) * Power(r, 3) * (2 * a - y * dr_dx)
                        + Power(a, 2) * Power(r, 5) * (2 * a - y * dr_dx)
                        + Power(r, 7) * (2 * a + 3 * y * dr_dx)))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                  + (2 * M * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                     * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                        + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                        + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     * (2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy - x * dr_dx)
                        + 2 * Power(a, 7) * Power(z, 3) * (-(y * dr_dy) + x * dr_dx)
                        - 3 * Power(a, 6) * Power(z, 2) * r
                              * (x * y * dr_dz + x * z * dr_dy + y * z * dr_dx)
                        + 2 * Power(a, 5) * z * Power(r, 2)
                              * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                 - y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                        + 2 * Power(a, 4) * z * Power(r, 3)
                              * (a * Power(x, 2) - a * Power(y, 2) + 3 * x * y * z * dr_dz
                                 - 2 * x * Power(z, 2) * dr_dy - 2 * y * Power(z, 2) * dr_dx)
                        - 2 * Power(a, 2) * z * Power(r, 6)
                              * (-3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                        + 4 * a * Power(r, 8)
                              * ((Power(x, 2) - Power(y, 2)) * dr_dz + z * (y * dr_dy - x * dr_dx))
                        + 2 * Power(a, 2) * Power(r, 7)
                              * (3 * x * y * dr_dz + 2 * z * (x * dr_dy + y * dr_dx))
                        + Power(r, 9) * (-3 * x * y * dr_dz + 3 * z * (x * dr_dy + y * dr_dx))
                        + Power(a, 2) * Power(r, 5)
                              * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                 + z
                                       * (2 * a * (Power(x, 2) - Power(y, 2))
                                          + x * (Power(a, 2) - Power(z, 2)) * dr_dy
                                          + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))))
                  / Power(Power(a, 2) + Power(r, 2), 5)
            + (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
               + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3) + Power(a, 2) * Power(r, 4)
               + Power(r, 6))
                  * ((M * Power(r, 4) * Power(a * x - y * r, 2)
                      * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                         + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                         + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                         + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                         + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                      * (3 * Power(a, 5) * x * Power(z, 2) * dr_dz
                         - Power(a, 2) * y * Power(r, 5) * dr_dz - 5 * a * x * Power(r, 6) * dr_dz
                         + 3 * y * Power(r, 7) * dr_dz
                         + Power(a, 2) * Power(r, 4) * (2 * y * z - a * x * dr_dz)
                         + Power(a, 3) * z * Power(r, 2) * (2 * a * y - x * z * dr_dz)
                         - Power(a, 2) * z * Power(r, 3) * (2 * a * x + y * z * dr_dz)
                         - Power(a, 4) * z * r * (2 * a * x + 5 * y * z * dr_dz)))
                         / (Power(Power(a, 2) + Power(r, 2), 5)
                            * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                               + Power(a, 2) * Power(z, 2) * Power(r, 2)
                               + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                               + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     - 4 * z * Power(r, 4) * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                     + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                     - (M * Power(z, 3) * r * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                        * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                        * dr_dy)
                           / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                              + Power(a, 2) * Power(z, 2) * Power(r, 2)
                              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                              + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     - (M * Power(r, 3) * (a * y + x * r)
                        * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                           + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                           + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                           + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                           + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                           + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                        * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                           - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                           + 2 * Power(a, 5) * z * Power(r, 2)
                                 * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                    + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                           + 2 * Power(a, 4) * z * Power(r, 3)
                                 * (a * Power(x, 2) - a * Power(y, 2) + 2 * a * Power(z, 2)
                                    + 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                    - 2 * y * Power(z, 2) * dr_dx)
                           - 2 * Power(a, 2) * z * Power(r, 6)
                                 * (3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                           + 4 * a * Power(r, 8)
                                 * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                    - z * (y * dr_dy + x * dr_dx))
                           + Power(a, 6) * Power(z, 2) * r
                                 * (-3 * x * y * dr_dz
                                    + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                           + 2 * Power(a, 2) * Power(r, 7)
                                 * (3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                           + Power(r, 9)
                                 * (-3 * x * y * dr_dz
                                    + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                           + Power(a, 2) * Power(r, 5)
                                 * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                    + z
                                          * (2 * a
                                                 * (Power(a, 2) + Power(x, 2) - Power(y, 2)
                                                    + Power(z, 2))
                                             + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                             + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                           / (Power(Power(a, 2) + Power(r, 2), 5)
                              * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                                 + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                 + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                                 + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     + r
                           * ((y * r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                  / (Power(a, 2) + Power(r, 2))
                              + (3 * (-(a * x) + y * r)
                                 * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                    / (Power(a, 2) + Power(r, 2))
                              + (2 * (a * x - y * r)
                                 * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dz)
                                    / Power(Power(a, 2) + Power(r, 2), 2)
                              - (2 * r * (-(a * x) + y * r)
                                 * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                                 * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz))
                                    / (Power(a, 2) + Power(r, 2))
                              + (M * Power(r, 3) * Power(a * x - y * r, 2)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (3 * Power(a, 5) * x * Power(z, 2) * dr_dz
                                    - Power(a, 2) * y * Power(r, 5) * dr_dz
                                    - 5 * a * x * Power(r, 6) * dr_dz + 3 * y * Power(r, 7) * dr_dz
                                    + Power(a, 2) * Power(r, 4) * (2 * y * z - a * x * dr_dz)
                                    + Power(a, 3) * z * Power(r, 2) * (2 * a * y - x * z * dr_dz)
                                    - Power(a, 2) * z * Power(r, 3) * (2 * a * x + y * z * dr_dz)
                                    - Power(a, 4) * z * r * (2 * a * x + 5 * y * z * dr_dz)))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                              - (M * Power(z, 3) * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                                 * (Power(a, 4) * Power(z, 2)
                                    + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                    + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                                    + Power(a, 2) * Power(r, 4) + Power(r, 6))
                                 * dr_dy)
                                    / (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6))
                              - (M * Power(r, 2) * (a * y + x * r)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                                    - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                                    + 2 * Power(a, 5) * z * Power(r, 2)
                                          * (a * x * y
                                             + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                             + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                                    + 2 * Power(a, 4) * z * Power(r, 3)
                                          * (a * Power(x, 2) - a * Power(y, 2) + 2 * a * Power(z, 2)
                                             + 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                             - 2 * y * Power(z, 2) * dr_dx)
                                    - 2 * Power(a, 2) * z * Power(r, 6)
                                          * (3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                                    + 4 * a * Power(r, 8)
                                          * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                             - z * (y * dr_dy + x * dr_dx))
                                    + Power(a, 6) * Power(z, 2) * r
                                          * (-3 * x * y * dr_dz
                                             + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                                    + 2 * Power(a, 2) * Power(r, 7)
                                          * (3 * x * y * dr_dz
                                             + 2 * z * (a - x * dr_dy + y * dr_dx))
                                    + Power(r, 9)
                                          * (-3 * x * y * dr_dz
                                             + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                                    + Power(a, 2) * Power(r, 5)
                                          * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                             + z
                                                   * (2 * a
                                                          * (Power(a, 2) + Power(x, 2) - Power(y, 2)
                                                             + Power(z, 2))
                                                      + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                                      + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))))))
        / (Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 3)
           * Sqrt(
               (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                + Power(a, 2) * Power(z, 2) * Power(r, 2)
                + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                + Power(a, 2) * Power(r, 4) + Power(r, 6))
               * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                  + Power(a, 2) * Power(z, 2) * Power(r, 2)
                  + 2 * M * (-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                  + Power(a, 2) * Power(r, 4) - 2 * M * Power(r, 5) + Power(r, 6))));

  mixedK[2][2]
      = (2 * M * r
         * ((Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
             + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3) + Power(a, 2) * Power(r, 4)
             + Power(r, 6))
                * (r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
                   + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz
                   - 2 * z * r * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                         * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz)
                   - (M * Power(z, 2) * r
                      * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                         + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                         + Power(a, 2) * Power(r, 4) + Power(r, 6))
                      * (2 * Power(r, 5) + Power(a, 2) * Power(z, 3) * dr_dz
                         - 3 * z * Power(r, 4) * dr_dz))
                         / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                            + Power(a, 2) * Power(z, 2) * Power(r, 2)
                            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                            + Power(a, 2) * Power(r, 4) + Power(r, 6))
                   + (M * Power(r, 2) * (-(a * x) + y * r)
                      * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                         + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                         + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                         + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                         + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                      * (2 * Power(a, 4) * y * Power(z, 2) * Power(r, 3) - 2 * y * Power(r, 9)
                         + 4 * Power(a, 5) * x * Power(z, 3) * r * dr_dz
                         + 2 * Power(a, 2) * z * Power(r, 5) * (y * z - 2 * a * x * dr_dz)
                         - 2 * a * Power(r, 7) * (a * y + 4 * x * z * dr_dz)
                         + Power(a, 6) * Power(z, 4) * dr_dy
                         + 2 * Power(a, 2) * Power(r, 6)
                               * (a * x + y * z * dr_dz - 3 * Power(z, 2) * dr_dy)
                         + Power(r, 8) * (2 * a * x + 6 * y * z * dr_dz - 3 * Power(z, 2) * dr_dy)
                         - 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                               * (a * x + 3 * y * z * dr_dz - Power(z, 2) * dr_dy)
                         + Power(a, 2) * Power(z, 2) * Power(r, 4)
                               * (-2 * a * x - 2 * y * z * dr_dz
                                  + (-3 * Power(a, 2) + Power(z, 2)) * dr_dy)))
                         / (Power(Power(a, 2) + Power(r, 2), 4)
                            * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                               + Power(a, 2) * Power(z, 2) * Power(r, 2)
                               + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                               + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                   + (M * Power(r, 2) * (a * y + x * r)
                      * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                         + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                         + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                         + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                         + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                      * (2 * Power(a, 4) * x * Power(z, 2) * Power(r, 3) - 2 * x * Power(r, 9)
                         - 4 * Power(a, 5) * y * Power(z, 3) * r * dr_dz
                         + 2 * Power(a, 2) * z * Power(r, 5) * (x * z + 2 * a * y * dr_dz)
                         - 2 * a * Power(r, 7) * (a * x - 4 * y * z * dr_dz)
                         + Power(a, 6) * Power(z, 4) * dr_dx
                         + Power(r, 8) * (-2 * a * y + 6 * x * z * dr_dz - 3 * Power(z, 2) * dr_dx)
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                               * (a * y - 3 * x * z * dr_dz + Power(z, 2) * dr_dx)
                         - 2 * Power(a, 2) * Power(r, 6)
                               * (a * y - x * z * dr_dz + 3 * Power(z, 2) * dr_dx)
                         + Power(a, 2) * Power(z, 2) * Power(r, 4)
                               * (2 * a * y - 2 * x * z * dr_dz
                                  + (-3 * Power(a, 2) + Power(z, 2)) * dr_dx)))
                         / (Power(Power(a, 2) + Power(r, 2), 4)
                            * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                               + Power(a, 2) * Power(z, 2) * Power(r, 2)
                               + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                               + Power(a, 2) * Power(r, 4) + Power(r, 6))))
            + M * z * Power(r, 2) * (a * x - y * r)
                  * ((M * Power(r, 4) * Power(a * x - y * r, 2)
                      * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                         + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                         + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                         + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                         + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                      * (3 * Power(a, 5) * x * Power(z, 2) * dr_dz
                         - Power(a, 2) * y * Power(r, 5) * dr_dz - 5 * a * x * Power(r, 6) * dr_dz
                         + 3 * y * Power(r, 7) * dr_dz
                         + Power(a, 2) * Power(r, 4) * (2 * y * z - a * x * dr_dz)
                         + Power(a, 3) * z * Power(r, 2) * (2 * a * y - x * z * dr_dz)
                         - Power(a, 2) * z * Power(r, 3) * (2 * a * x + y * z * dr_dz)
                         - Power(a, 4) * z * r * (2 * a * x + 5 * y * z * dr_dz)))
                         / (Power(Power(a, 2) + Power(r, 2), 5)
                            * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                               + Power(a, 2) * Power(z, 2) * Power(r, 2)
                               + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                               + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     - 4 * z * Power(r, 4) * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                     + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dy
                     - (M * Power(z, 3) * r * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                        * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                        * dr_dy)
                           / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                              + Power(a, 2) * Power(z, 2) * Power(r, 2)
                              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                              + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     - (M * Power(r, 3) * (a * y + x * r)
                        * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                           + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                           + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                           + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                           + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                           + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                        * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                           - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                           + 2 * Power(a, 5) * z * Power(r, 2)
                                 * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                    + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                           + 2 * Power(a, 4) * z * Power(r, 3)
                                 * (a * Power(x, 2) - a * Power(y, 2) + 2 * a * Power(z, 2)
                                    + 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                    - 2 * y * Power(z, 2) * dr_dx)
                           - 2 * Power(a, 2) * z * Power(r, 6)
                                 * (3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                           + 4 * a * Power(r, 8)
                                 * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                    - z * (y * dr_dy + x * dr_dx))
                           + Power(a, 6) * Power(z, 2) * r
                                 * (-3 * x * y * dr_dz
                                    + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                           + 2 * Power(a, 2) * Power(r, 7)
                                 * (3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                           + Power(r, 9)
                                 * (-3 * x * y * dr_dz
                                    + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                           + Power(a, 2) * Power(r, 5)
                                 * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                    + z
                                          * (2 * a
                                                 * (Power(a, 2) + Power(x, 2) - Power(y, 2)
                                                    + Power(z, 2))
                                             + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                             + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                           / (Power(Power(a, 2) + Power(r, 2), 5)
                              * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                                 + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                 + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                                 + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     + r
                           * ((y * r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                  / (Power(a, 2) + Power(r, 2))
                              + (3 * (-(a * x) + y * r)
                                 * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                    / (Power(a, 2) + Power(r, 2))
                              + (2 * (a * x - y * r)
                                 * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dz)
                                    / Power(Power(a, 2) + Power(r, 2), 2)
                              - (2 * r * (-(a * x) + y * r)
                                 * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                                 * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz))
                                    / (Power(a, 2) + Power(r, 2))
                              + (M * Power(r, 3) * Power(a * x - y * r, 2)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (3 * Power(a, 5) * x * Power(z, 2) * dr_dz
                                    - Power(a, 2) * y * Power(r, 5) * dr_dz
                                    - 5 * a * x * Power(r, 6) * dr_dz + 3 * y * Power(r, 7) * dr_dz
                                    + Power(a, 2) * Power(r, 4) * (2 * y * z - a * x * dr_dz)
                                    + Power(a, 3) * z * Power(r, 2) * (2 * a * y - x * z * dr_dz)
                                    - Power(a, 2) * z * Power(r, 3) * (2 * a * x + y * z * dr_dz)
                                    - Power(a, 4) * z * r * (2 * a * x + 5 * y * z * dr_dz)))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                              - (M * Power(z, 3) * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                                 * (Power(a, 4) * Power(z, 2)
                                    + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                    + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                                    + Power(a, 2) * Power(r, 4) + Power(r, 6))
                                 * dr_dy)
                                    / (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6))
                              - (M * Power(r, 2) * (a * y + x * r)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                                    - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                                    + 2 * Power(a, 5) * z * Power(r, 2)
                                          * (a * x * y
                                             + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                                             + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                                    + 2 * Power(a, 4) * z * Power(r, 3)
                                          * (a * Power(x, 2) - a * Power(y, 2) + 2 * a * Power(z, 2)
                                             + 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                             - 2 * y * Power(z, 2) * dr_dx)
                                    - 2 * Power(a, 2) * z * Power(r, 6)
                                          * (3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                                    + 4 * a * Power(r, 8)
                                          * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                             - z * (y * dr_dy + x * dr_dx))
                                    + Power(a, 6) * Power(z, 2) * r
                                          * (-3 * x * y * dr_dz
                                             + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                                    + 2 * Power(a, 2) * Power(r, 7)
                                          * (3 * x * y * dr_dz
                                             + 2 * z * (a - x * dr_dy + y * dr_dx))
                                    + Power(r, 9)
                                          * (-3 * x * y * dr_dz
                                             + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                                    + Power(a, 2) * Power(r, 5)
                                          * (x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                                             + z
                                                   * (2 * a
                                                          * (Power(a, 2) + Power(x, 2) - Power(y, 2)
                                                             + Power(z, 2))
                                                      + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                                      + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))))
            - M * z * Power(r, 2) * (a * y + x * r)
                  * ((M * Power(r, 4) * Power(a * y + x * r, 2)
                      * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                         + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                         + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                         + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                         + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                         + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                      * (-3 * Power(a, 5) * y * Power(z, 2) * dr_dz
                         - Power(a, 2) * x * Power(r, 5) * dr_dz + 5 * a * y * Power(r, 6) * dr_dz
                         + 3 * x * Power(r, 7) * dr_dz
                         + Power(a, 2) * Power(r, 4) * (2 * x * z + a * y * dr_dz)
                         + Power(a, 4) * z * r * (2 * a * y - 5 * x * z * dr_dz)
                         + Power(a, 2) * z * Power(r, 3) * (2 * a * y - x * z * dr_dz)
                         + Power(a, 3) * z * Power(r, 2) * (2 * a * x + y * z * dr_dz)))
                         / (Power(Power(a, 2) + Power(r, 2), 5)
                            * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                               + Power(a, 2) * Power(z, 2) * Power(r, 2)
                               + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                               + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     - 4 * z * Power(r, 4) * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                     + 2 * z * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dx
                     - (M * Power(z, 3) * r * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                        * (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6))
                        * dr_dx)
                           / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                              + Power(a, 2) * Power(z, 2) * Power(r, 2)
                              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                              + Power(a, 2) * Power(r, 4) + Power(r, 6))
                     + (M * Power(r, 3) * (-(a * x) + y * r)
                        * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                           + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                           + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                           + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                           + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                           + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                        * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                           - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                           + 2 * Power(a, 5) * z * Power(r, 2)
                                 * (-(a * x * y) + 2 * (Power(x, 2) - Power(y, 2)) * z * dr_dz
                                    + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                           + 2 * Power(a, 4) * z * Power(r, 3)
                                 * (-(a * Power(x, 2)) + a * Power(y, 2) + 2 * a * Power(z, 2)
                                    - 3 * x * y * z * dr_dz + 2 * x * Power(z, 2) * dr_dy
                                    - 2 * y * Power(z, 2) * dr_dx)
                           + 2 * Power(a, 2) * z * Power(r, 6)
                                 * (-3 * a * y * dr_dy + x * (y - 3 * a * dr_dx))
                           - 4 * a * Power(r, 8)
                                 * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                    + z * (y * dr_dy + x * dr_dx))
                           + Power(a, 6) * Power(z, 2) * r
                                 * (3 * x * y * dr_dz + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                           + 2 * Power(a, 2) * Power(r, 7)
                                 * (-3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                           + Power(r, 9)
                                 * (3 * x * y * dr_dz + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                           + Power(a, 2) * Power(r, 5)
                                 * (-(x * y * (Power(a, 2) + Power(z, 2)) * dr_dz)
                                    + z
                                          * (2 * a
                                                 * (Power(a, 2) - Power(x, 2) + Power(y, 2)
                                                    + Power(z, 2))
                                             + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                             + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                           / (Power(Power(a, 2) + Power(r, 2), 5)
                              * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                                 + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                 + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                                 + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                     + r
                           * ((x * r * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                  / (Power(a, 2) + Power(r, 2))
                              + (3 * (a * y + x * r)
                                 * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2) * dr_dz)
                                    / (Power(a, 2) + Power(r, 2))
                              - (2 * (a * y + x * r)
                                 * Power(Power(a, 2) * Power(z, 2) * r + Power(r, 5), 2) * dr_dz)
                                    / Power(Power(a, 2) + Power(r, 2), 2)
                              - (2 * r * (a * y + x * r) * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                                 * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz))
                                    / (Power(a, 2) + Power(r, 2))
                              + (M * Power(r, 3) * Power(a * y + x * r, 2)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (-3 * Power(a, 5) * y * Power(z, 2) * dr_dz
                                    - Power(a, 2) * x * Power(r, 5) * dr_dz
                                    + 5 * a * y * Power(r, 6) * dr_dz + 3 * x * Power(r, 7) * dr_dz
                                    + Power(a, 2) * Power(r, 4) * (2 * x * z + a * y * dr_dz)
                                    + Power(a, 4) * z * r * (2 * a * y - 5 * x * z * dr_dz)
                                    + Power(a, 2) * z * Power(r, 3) * (2 * a * y - x * z * dr_dz)
                                    + Power(a, 3) * z * Power(r, 2) * (2 * a * x + y * z * dr_dz)))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))
                              - (M * Power(z, 3) * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4))
                                 * (Power(a, 4) * Power(z, 2)
                                    + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                    + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3)
                                    + Power(a, 2) * Power(r, 4) + Power(r, 6))
                                 * dr_dx)
                                    / (Power(a, 4) * Power(z, 2)
                                       + 2 * Power(a, 2) * M * Power(z, 2) * r
                                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                             * Power(r, 3)
                                       + Power(a, 2) * Power(r, 4) + Power(r, 6))
                              + (M * Power(r, 2) * (-(a * x) + y * r)
                                 * (Power(a, 6) * Power(z, 2)
                                    + 2 * Power(a, 4) * M * Power(z, 2) * r
                                    + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                                    + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2))
                                          * Power(r, 3)
                                    + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2))
                                          * Power(r, 4)
                                    + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                                    + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                                 * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                                    - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                                    + 2 * Power(a, 5) * z * Power(r, 2)
                                          * (-(a * x * y)
                                             + 2 * (Power(x, 2) - Power(y, 2)) * z * dr_dz
                                             + y * Power(z, 2) * dr_dy + x * Power(z, 2) * dr_dx)
                                    + 2 * Power(a, 4) * z * Power(r, 3)
                                          * (-(a * Power(x, 2)) + a * Power(y, 2)
                                             + 2 * a * Power(z, 2) - 3 * x * y * z * dr_dz
                                             + 2 * x * Power(z, 2) * dr_dy
                                             - 2 * y * Power(z, 2) * dr_dx)
                                    + 2 * Power(a, 2) * z * Power(r, 6)
                                          * (-3 * a * y * dr_dy + x * (y - 3 * a * dr_dx))
                                    - 4 * a * Power(r, 8)
                                          * ((Power(x, 2) - Power(y, 2)) * dr_dz
                                             + z * (y * dr_dy + x * dr_dx))
                                    + Power(a, 6) * Power(z, 2) * r
                                          * (3 * x * y * dr_dz
                                             + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                                    + 2 * Power(a, 2) * Power(r, 7)
                                          * (-3 * x * y * dr_dz
                                             + 2 * z * (a - x * dr_dy + y * dr_dx))
                                    + Power(r, 9)
                                          * (3 * x * y * dr_dz
                                             + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                                    + Power(a, 2) * Power(r, 5)
                                          * (-(x * y * (Power(a, 2) + Power(z, 2)) * dr_dz)
                                             + z
                                                   * (2 * a
                                                          * (Power(a, 2) - Power(x, 2) + Power(y, 2)
                                                             + Power(z, 2))
                                                      + x * (-Power(a, 2) + Power(z, 2)) * dr_dy
                                                      + y * (Power(a, 2) - Power(z, 2)) * dr_dx))))
                                    / (Power(Power(a, 2) + Power(r, 2), 5)
                                       * (Power(a, 4) * Power(z, 2)
                                          + 2 * Power(a, 2) * M * Power(z, 2) * r
                                          + Power(a, 2) * Power(z, 2) * Power(r, 2)
                                          + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                                                * Power(r, 3)
                                          + Power(a, 2) * Power(r, 4) + Power(r, 6)))))))
        / (Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 3)
           * Sqrt(
               (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                + Power(a, 2) * Power(z, 2) * Power(r, 2)
                + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                + Power(a, 2) * Power(r, 4) + Power(r, 6))
               * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                  + Power(a, 2) * Power(z, 2) * Power(r, 2)
                  + 2 * M * (-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                  + Power(a, 2) * Power(r, 4) - 2 * M * Power(r, 5) + Power(r, 6))));

  return mixedK;
}