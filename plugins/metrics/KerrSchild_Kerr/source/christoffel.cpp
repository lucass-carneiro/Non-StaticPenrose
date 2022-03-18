#include "KerrSchild_Kerr.hpp"
#include "aux_functions.hpp"

using grlensing::KerrSchild_Kerr;
using grlensing::metric_server;
using ksk_aux::d_r_KS_dx;
using ksk_aux::d_r_KS_dy;
using ksk_aux::d_r_KS_dz;
using ksk_aux::Power;
using ksk_aux::r_KS;

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::spatial_christoffel(double, double x,
                                                                               double y, double z)
    -> metric_server::chirstofell_t {

  auto r = r_KS(a, x, y, z);
  auto dr_dx = d_r_KS_dx(a, x, y, z);
  auto dr_dy = d_r_KS_dy(a, x, y, z);
  auto dr_dz = d_r_KS_dz(a, x, y, z);

  metric_server::chirstofell_t Gamma{};

  Gamma[0][0][0] = -(
      (M * Power(r, 2) * (a * y + x * r)
       * ((Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
           + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
           + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
           + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
           + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
           + Power(r, 8))
              * (-2 * Power(r, 2) * (Power(a, 2) + Power(r, 2))
                     * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                 + (-3 * Power(a, 5) * y * Power(z, 2) - 5 * Power(a, 4) * x * Power(z, 2) * r
                    + Power(a, 3) * y * Power(z, 2) * Power(r, 2)
                    - Power(a, 2) * x * Power(z, 2) * Power(r, 3) + Power(a, 3) * y * Power(r, 4)
                    - Power(a, 2) * x * Power(r, 5) + 5 * a * y * Power(r, 6) + 3 * x * Power(r, 7))
                       * dr_dx)
          + 2 * M * z * r * (Power(a, 2) + Power(r, 2))
                * (4 * Power(r, 3) * Power(a * y + x * r, 2)
                       * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   - 2 * x * Power(r, 2) * (a * y + x * r) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   - 3 * r * Power(a * y + x * r, 2) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   + 2 * Power(r, 2) * Power(a * y + x * r, 2) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz)
                   - 8 * z * Power(r, 4) * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 2)
                         * dr_dx
                   - 4 * z * Power(r, 2) * (a * y + x * r) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                   + 4 * z * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 2)
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                   + 2 * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * (r + x * dr_dx))
          - 2 * M * Power(r, 3) * (a * x - y * r)
                * (-2 * Power(a, 4) * y * Power(z, 2) * Power(r, 3)
                   + 2 * Power(a, 2) * y * (-2 * Power(a, 2) + Power(z, 2)) * Power(r, 5)
                   + 2 * y * Power(r, 9)
                   - 3 * Power(a, 6) * y * Power(z, 2) * (y * dr_dy + 2 * x * dr_dx)
                   - 4 * Power(a, 5) * Power(z, 2) * r
                         * (a * y + 2 * x * y * dr_dy + 2 * (x - y) * (x + y) * dr_dx)
                   - 2 * a * Power(r, 7)
                         * (a * y - 4 * x * y * dr_dy + 4 * (-x + y) * (x + y) * dr_dx)
                   + Power(a, 2) * Power(r, 4)
                         * (-6 * a * x * Power(z, 2) + (a * y - x * z) * (a * y + x * z) * dr_dy
                            + 2 * x * y * (Power(a, 2) + Power(z, 2)) * dr_dx)
                   + Power(a, 2) * Power(r, 6)
                         * (-((Power(x, 2) - 5 * Power(y, 2)) * dr_dy)
                            - 6 * x * (a - 2 * y * dr_dx))
                   + Power(a, 4) * Power(z, 2) * Power(r, 2)
                         * ((-5 * Power(x, 2) + Power(y, 2)) * dr_dy - 6 * x * (a - 2 * y * dr_dx))
                   + 3 * x * Power(r, 8) * (x * dr_dy - 2 * (a + y * dr_dx)))))
      / (Power(Power(a, 2) + Power(r, 2), 4) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
         * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
            + Power(a, 2) * Power(z, 2) * Power(r, 2)
            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
            + Power(a, 2) * Power(r, 4) + Power(r, 6))));

  Gamma[0][0][1]
      = (M * Power(r, 2) * (a * y + x * r)
         * (-((Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
               + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
               + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
               + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
               + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
               + Power(r, 8))
              * (-2 * a * r * (Power(a, 2) + Power(r, 2))
                     * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                 + (-3 * Power(a, 5) * y * Power(z, 2) - 5 * Power(a, 4) * x * Power(z, 2) * r
                    + Power(a, 3) * y * Power(z, 2) * Power(r, 2)
                    - Power(a, 2) * x * Power(z, 2) * Power(r, 3) + Power(a, 3) * y * Power(r, 4)
                    - Power(a, 2) * x * Power(r, 5) + 5 * a * y * Power(r, 6) + 3 * x * Power(r, 7))
                       * dr_dy))
            + 2 * M * Power(r, 3) * Power(a * x - y * r, 2)
                  * (2 * a * r * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                     + (3 * Power(a, 5) * x * Power(z, 2) - 5 * Power(a, 4) * y * Power(z, 2) * r
                        - Power(a, 3) * x * Power(z, 2) * Power(r, 2)
                        - Power(a, 2) * y * Power(z, 2) * Power(r, 3)
                        - Power(a, 3) * x * Power(r, 4) - Power(a, 2) * y * Power(r, 5)
                        - 5 * a * x * Power(r, 6) + 3 * y * Power(r, 7))
                           * dr_dx)
            + 2 * M * z * r * (Power(a, 2) + Power(r, 2))
                  * (2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy - x * dr_dx)
                     + 2 * Power(a, 7) * Power(z, 3) * (-(y * dr_dy) + x * dr_dx)
                     + 3 * Power(r, 9) * (-(x * y * dr_dz) + x * z * dr_dy + y * z * dr_dx)
                     - 3 * Power(a, 6) * Power(z, 2) * r
                           * (x * y * dr_dz + x * z * dr_dy + y * z * dr_dx)
                     - 2 * Power(a, 2) * z * Power(r, 6)
                           * (-3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                     + 4 * a * Power(r, 8)
                           * ((x - y) * (x + y) * dr_dz + z * (y * dr_dy - x * dr_dx))
                     + 2 * Power(a, 2) * Power(r, 7)
                           * (3 * x * y * dr_dz + 2 * z * (x * dr_dy + y * dr_dx))
                     + 2 * Power(a, 4) * z * Power(r, 3)
                           * (a * (x - y) * (x + y) + 3 * x * y * z * dr_dz
                              - 2 * Power(z, 2) * (x * dr_dy + y * dr_dx))
                     + Power(a, 2) * Power(r, 5)
                           * (2 * a * (x - y) * (x + y) * z
                              + x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                              + (a - z) * z * (a + z) * (x * dr_dy + y * dr_dx))
                     + 2 * Power(a, 5) * z * Power(r, 2)
                           * (a * x * y
                              + z
                                    * (2 * (-Power(x, 2) + Power(y, 2)) * dr_dz - y * z * dr_dy
                                       + x * z * dr_dx)))))
        / (Power(Power(a, 2) + Power(r, 2), 4) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
           * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
              + Power(a, 2) * Power(z, 2) * Power(r, 2)
              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
              + Power(a, 2) * Power(r, 4) + Power(r, 6)));

  Gamma[0][0][2] = -(
      (M * Power(r, 2) * (a * y + x * r)
       * ((Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
           + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
           + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
           + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
           + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
           + Power(r, 8))
              * (2 * Power(a, 2) * z * r * (a * y + x * r) * (Power(a, 2) + Power(r, 2))
                 + (-3 * Power(a, 5) * y * Power(z, 2) - 5 * Power(a, 4) * x * Power(z, 2) * r
                    + Power(a, 3) * y * Power(z, 2) * Power(r, 2)
                    - Power(a, 2) * x * Power(z, 2) * Power(r, 3) + Power(a, 3) * y * Power(r, 4)
                    - Power(a, 2) * x * Power(r, 5) + 5 * a * y * Power(r, 6) + 3 * x * Power(r, 7))
                       * dr_dz)
          + 2 * M * Power(z, 3) * Power(Power(a, 2) + Power(r, 2), 4)
                * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4)) * dr_dx
          - 2 * M * Power(r, 2) * (a * x - y * r)
                * (-2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                   + 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                   + 2 * Power(a, 2) * z * Power(r, 6)
                         * (-(x * y) + 3 * a * y * dr_dy + 3 * a * x * dr_dx)
                   + 4 * a * Power(r, 8) * ((x - y) * (x + y) * dr_dz + z * (y * dr_dy + x * dr_dx))
                   - Power(a, 6) * Power(z, 2) * r
                         * (3 * x * y * dr_dz + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                   + Power(a, 2) * Power(r, 5)
                         * (-2 * a * z * (Power(a, 2) - Power(x, 2) + Power(y, 2) + Power(z, 2))
                            + x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                            + (a - z) * z * (a + z) * (x * dr_dy - y * dr_dx))
                   + 2 * Power(a, 2) * Power(r, 7)
                         * (3 * x * y * dr_dz - 2 * z * (a - x * dr_dy + y * dr_dx))
                   - Power(r, 9) * (3 * x * y * dr_dz + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                   + 2 * Power(a, 4) * z * Power(r, 3)
                         * (a * (Power(x, 2) - Power(y, 2) - 2 * Power(z, 2))
                            + z * (3 * x * y * dr_dz - 2 * x * z * dr_dy + 2 * y * z * dr_dx))
                   + 2 * Power(a, 5) * z * Power(r, 2)
                         * (a * x * y
                            - z * (2 * (x - y) * (x + y) * dr_dz + z * (y * dr_dy + x * dr_dx))))))
      / (Power(Power(a, 2) + Power(r, 2), 4) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
         * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
            + Power(a, 2) * Power(z, 2) * Power(r, 2)
            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
            + Power(a, 2) * Power(r, 4) + Power(r, 6))));

  Gamma[0][1][1] = -(
      (M * Power(r, 2)
       * (2 * M * Power(r, 3) * (a * y + x * r) * Power(a * x - y * r, 2)
              * (2 * Power(r, 2) * (Power(a, 2) + Power(r, 2))
                     * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                 + (-3 * Power(a, 5) * x * Power(z, 2) + 5 * Power(a, 4) * y * Power(z, 2) * r
                    + Power(a, 3) * x * Power(z, 2) * Power(r, 2)
                    + Power(a, 2) * y * Power(z, 2) * Power(r, 3) + Power(a, 3) * x * Power(r, 4)
                    + Power(a, 2) * y * Power(r, 5) + 5 * a * x * Power(r, 6) - 3 * y * Power(r, 7))
                       * dr_dy)
          + 2 * M * z * r * (a * y + x * r) * (Power(a, 2) + Power(r, 2))
                * (4 * Power(r, 3) * Power(a * x - y * r, 2)
                       * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   + 2 * y * Power(r, 2) * (a * x - y * r) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   - 3 * r * Power(a * x - y * r, 2) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   + 2 * Power(r, 2) * Power(a * x - y * r, 2) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz)
                   + 8 * z * Power(r, 4) * (a * x - y * r) * Power(Power(a, 2) + Power(r, 2), 2)
                         * dr_dy
                   + 4 * z * Power(r, 2) * (a * x - y * r) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                   - 4 * z * (a * x - y * r) * Power(Power(a, 2) + Power(r, 2), 2)
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                   + 2 * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * (r + y * dr_dy))
          - (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
             + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
             + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
             + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
             + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
             + Power(r, 8))
                * (-2 * Power(a, 4) * x * Power(z, 2) * Power(r, 3)
                   + 2 * Power(a, 2) * x * (-2 * Power(a, 2) + Power(z, 2)) * Power(r, 5)
                   + 2 * x * Power(r, 9)
                   - 3 * Power(a, 6) * x * Power(z, 2) * (2 * y * dr_dy + x * dr_dx)
                   + 3 * y * Power(r, 8) * (2 * a - 2 * x * dr_dy + y * dr_dx)
                   + Power(a, 4) * Power(z, 2) * Power(r, 2)
                         * (6 * y * (a + 2 * x * dr_dy) + (Power(x, 2) - 5 * Power(y, 2)) * dr_dx)
                   + Power(a, 2) * Power(r, 6)
                         * (6 * y * (a + 2 * x * dr_dy) + (5 * Power(x, 2) - Power(y, 2)) * dr_dx)
                   + Power(a, 2) * Power(r, 4)
                         * (6 * a * y * Power(z, 2)
                            + 2 * x * y * (Power(a, 2) + Power(z, 2)) * dr_dy
                            + (a * x - y * z) * (a * x + y * z) * dr_dx)
                   - 4 * Power(a, 5) * Power(z, 2) * r
                         * (2 * (x - y) * (x + y) * dr_dy + x * (a - 2 * y * dr_dx))
                   - 2 * a * Power(r, 7)
                         * (4 * (-Power(x, 2) + Power(y, 2)) * dr_dy + x * (a + 4 * y * dr_dx)))))
      / (Power(Power(a, 2) + Power(r, 2), 4) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
         * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
            + Power(a, 2) * Power(z, 2) * Power(r, 2)
            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
            + Power(a, 2) * Power(r, 4) + Power(r, 6))));

  Gamma[0][1][2] = -(
      (M * r
       * (2 * M * Power(r, 4) * (a * y + x * r) * Power(a * x - y * r, 2)
              * (2 * Power(a, 2) * z * r * (a * x - y * r) * (Power(a, 2) + Power(r, 2))
                 + (-3 * Power(a, 5) * x * Power(z, 2) + 5 * Power(a, 4) * y * Power(z, 2) * r
                    + Power(a, 3) * x * Power(z, 2) * Power(r, 2)
                    + Power(a, 2) * y * Power(z, 2) * Power(r, 3) + Power(a, 3) * x * Power(r, 4)
                    + Power(a, 2) * y * Power(r, 5) + 5 * a * x * Power(r, 6) - 3 * y * Power(r, 7))
                       * dr_dz)
          + 2 * M * Power(z, 3) * r * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 4)
                * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4)) * dr_dy
          - (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
             + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
             + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
             + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
             + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
             + Power(r, 8))
                * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                   - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                   + Power(r, 9) * (2 * a * z - 3 * x * (y * dr_dz + z * dr_dy) + 3 * y * z * dr_dx)
                   - 2 * Power(a, 2) * z * Power(r, 6)
                         * (3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                   + 4 * a * Power(r, 8) * ((x - y) * (x + y) * dr_dz - z * (y * dr_dy + x * dr_dx))
                   + 2 * Power(a, 5) * z * Power(r, 2)
                         * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                            + Power(z, 2) * (y * dr_dy + x * dr_dx))
                   + Power(a, 6) * Power(z, 2) * r
                         * (-3 * x * y * dr_dz + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                   + Power(a, 2) * Power(r, 5)
                         * (2 * a * z * (Power(a, 2) + Power(x, 2) - Power(y, 2) + Power(z, 2))
                            + x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                            + z * (-Power(a, 2) + Power(z, 2)) * (x * dr_dy - y * dr_dx))
                   + 2 * Power(a, 2) * Power(r, 7)
                         * (3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                   + 2 * Power(a, 4) * z * Power(r, 3)
                         * (a * (Power(x, 2) - Power(y, 2) + 2 * Power(z, 2))
                            + z * (3 * x * y * dr_dz + 2 * x * z * dr_dy - 2 * y * z * dr_dx)))))
      / (Power(Power(a, 2) + Power(r, 2), 4) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
         * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
            + Power(a, 2) * Power(z, 2) * Power(r, 2)
            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
            + Power(a, 2) * Power(r, 4) + Power(r, 6))));

  Gamma[0][2][2] = -(
      (M
       * (2 * M * Power(z, 2) * Power(r, 2) * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 3)
              * (2 * Power(r, 5) + z * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4)) * dr_dz)
          + 2 * M * Power(r, 3) * (a * y + x * r) * (a * x - y * r)
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
                            + (-3 * Power(a, 2) + Power(z, 2)) * dr_dy))
          - (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
             + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
             + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
             + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
             + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
             + Power(r, 8))
                * (2 * r
                       * (r * (a * y + x * r) * (Power(a, 2) + Power(r, 2))
                              * (-(Power(a, 2) * Power(z, 2)) + Power(r, 4))
                          + z
                                * (2 * Power(a, 5) * y * Power(z, 2)
                                   + 3 * Power(a, 4) * x * Power(z, 2) * r
                                   + Power(a, 2) * x * Power(z, 2) * Power(r, 3)
                                   - 2 * Power(a, 3) * y * Power(r, 4)
                                   - Power(a, 2) * x * Power(r, 5) - 4 * a * y * Power(r, 6)
                                   - 3 * x * Power(r, 7))
                                * dr_dz)
                   + Power(z, 2) * Power(Power(a, 2) + Power(r, 2), 2)
                         * (-(Power(a, 2) * Power(z, 2)) + 3 * Power(r, 4)) * dr_dx)))
      / (Power(Power(a, 2) + Power(r, 2), 3) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
         * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
            + Power(a, 2) * Power(z, 2) * Power(r, 2)
            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
            + Power(a, 2) * Power(r, 4) + Power(r, 6))));

  Gamma[1][0][0] = -(
      (M * Power(r, 2)
       * (2 * M * Power(r, 3) * Power(a * y + x * r, 2) * (-(a * x) + y * r)
              * (2 * Power(r, 2) * (Power(a, 2) + Power(r, 2))
                     * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                 + (3 * Power(a, 5) * y * Power(z, 2) + 5 * Power(a, 4) * x * Power(z, 2) * r
                    - Power(a, 3) * y * Power(z, 2) * Power(r, 2)
                    + Power(a, 2) * x * Power(z, 2) * Power(r, 3) - Power(a, 3) * y * Power(r, 4)
                    + Power(a, 2) * x * Power(r, 5) - 5 * a * y * Power(r, 6) - 3 * x * Power(r, 7))
                       * dr_dx)
          - 2 * M * z * r * (a * x - y * r) * (Power(a, 2) + Power(r, 2))
                * (4 * Power(r, 3) * Power(a * y + x * r, 2)
                       * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   - 2 * x * Power(r, 2) * (a * y + x * r) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   - 3 * r * Power(a * y + x * r, 2) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   + 2 * Power(r, 2) * Power(a * y + x * r, 2) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz)
                   - 8 * z * Power(r, 4) * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 2)
                         * dr_dx
                   - 4 * z * Power(r, 2) * (a * y + x * r) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                   + 4 * z * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 2)
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                   + 2 * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * (r + x * dr_dx))
          - (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
             + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
             + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
             + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
             + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
             + Power(r, 8))
                * (-2 * Power(a, 4) * y * Power(z, 2) * Power(r, 3)
                   + 2 * Power(a, 2) * y * (-2 * Power(a, 2) + Power(z, 2)) * Power(r, 5)
                   + 2 * y * Power(r, 9)
                   - 3 * Power(a, 6) * y * Power(z, 2) * (y * dr_dy + 2 * x * dr_dx)
                   - 4 * Power(a, 5) * Power(z, 2) * r
                         * (a * y + 2 * x * y * dr_dy + 2 * (x - y) * (x + y) * dr_dx)
                   - 2 * a * Power(r, 7)
                         * (a * y - 4 * x * y * dr_dy + 4 * (-x + y) * (x + y) * dr_dx)
                   + Power(a, 2) * Power(r, 4)
                         * (-6 * a * x * Power(z, 2) + (a * y - x * z) * (a * y + x * z) * dr_dy
                            + 2 * x * y * (Power(a, 2) + Power(z, 2)) * dr_dx)
                   + Power(a, 2) * Power(r, 6)
                         * (-((Power(x, 2) - 5 * Power(y, 2)) * dr_dy)
                            - 6 * x * (a - 2 * y * dr_dx))
                   + Power(a, 4) * Power(z, 2) * Power(r, 2)
                         * ((-5 * Power(x, 2) + Power(y, 2)) * dr_dy - 6 * x * (a - 2 * y * dr_dx))
                   + 3 * x * Power(r, 8) * (x * dr_dy - 2 * (a + y * dr_dx)))))
      / (Power(Power(a, 2) + Power(r, 2), 4) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
         * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
            + Power(a, 2) * Power(z, 2) * Power(r, 2)
            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
            + Power(a, 2) * Power(r, 4) + Power(r, 6))));

  Gamma[1][0][1] = -(
      (M * Power(r, 2)
       * (2 * M * Power(r, 3) * Power(a * y + x * r, 2) * (-(a * x) + y * r)
              * (2 * a * r * (Power(a, 2) + Power(r, 2)) * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                 + (3 * Power(a, 5) * y * Power(z, 2) + 5 * Power(a, 4) * x * Power(z, 2) * r
                    - Power(a, 3) * y * Power(z, 2) * Power(r, 2)
                    + Power(a, 2) * x * Power(z, 2) * Power(r, 3) - Power(a, 3) * y * Power(r, 4)
                    + Power(a, 2) * x * Power(r, 5) - 5 * a * y * Power(r, 6) - 3 * x * Power(r, 7))
                       * dr_dy)
          + (a * x - y * r)
                * (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                   + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                   + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                   + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                   + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                   + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                * (-2 * a * r * (Power(a, 2) + Power(r, 2))
                       * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                   + (-3 * Power(a, 5) * x * Power(z, 2) + 5 * Power(a, 4) * y * Power(z, 2) * r
                      + Power(a, 3) * x * Power(z, 2) * Power(r, 2)
                      + Power(a, 2) * y * Power(z, 2) * Power(r, 3) + Power(a, 3) * x * Power(r, 4)
                      + Power(a, 2) * y * Power(r, 5) + 5 * a * x * Power(r, 6)
                      - 3 * y * Power(r, 7))
                         * dr_dx)
          + 2 * M * z * r * (a * x - y * r) * (Power(a, 2) + Power(r, 2))
                * (2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy - x * dr_dx)
                   + 2 * Power(a, 7) * Power(z, 3) * (-(y * dr_dy) + x * dr_dx)
                   + 3 * Power(r, 9) * (-(x * y * dr_dz) + x * z * dr_dy + y * z * dr_dx)
                   - 3 * Power(a, 6) * Power(z, 2) * r
                         * (x * y * dr_dz + x * z * dr_dy + y * z * dr_dx)
                   - 2 * Power(a, 2) * z * Power(r, 6)
                         * (-3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                   + 4 * a * Power(r, 8) * ((x - y) * (x + y) * dr_dz + z * (y * dr_dy - x * dr_dx))
                   + 2 * Power(a, 2) * Power(r, 7)
                         * (3 * x * y * dr_dz + 2 * z * (x * dr_dy + y * dr_dx))
                   + 2 * Power(a, 4) * z * Power(r, 3)
                         * (a * (x - y) * (x + y) + 3 * x * y * z * dr_dz
                            - 2 * Power(z, 2) * (x * dr_dy + y * dr_dx))
                   + Power(a, 2) * Power(r, 5)
                         * (2 * a * (x - y) * (x + y) * z
                            + x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                            + (a - z) * z * (a + z) * (x * dr_dy + y * dr_dx))
                   + 2 * Power(a, 5) * z * Power(r, 2)
                         * (a * x * y
                            + z
                                  * (2 * (-Power(x, 2) + Power(y, 2)) * dr_dz - y * z * dr_dy
                                     + x * z * dr_dx)))))
      / (Power(Power(a, 2) + Power(r, 2), 4) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
         * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
            + Power(a, 2) * Power(z, 2) * Power(r, 2)
            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
            + Power(a, 2) * Power(r, 4) + Power(r, 6))));

  Gamma[1][0][2]
      = (M * r
         * (2 * M * Power(r, 4) * Power(a * y + x * r, 2) * (-(a * x) + y * r)
                * (2 * Power(a, 2) * z * r * (a * y + x * r) * (Power(a, 2) + Power(r, 2))
                   + (-3 * Power(a, 5) * y * Power(z, 2) - 5 * Power(a, 4) * x * Power(z, 2) * r
                      + Power(a, 3) * y * Power(z, 2) * Power(r, 2)
                      - Power(a, 2) * x * Power(z, 2) * Power(r, 3) + Power(a, 3) * y * Power(r, 4)
                      - Power(a, 2) * x * Power(r, 5) + 5 * a * y * Power(r, 6)
                      + 3 * x * Power(r, 7))
                         * dr_dz)
            + 2 * M * Power(z, 3) * r * (a * x - y * r) * Power(Power(a, 2) + Power(r, 2), 4)
                  * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4)) * dr_dx
            + (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
               + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
               + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
               + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
               + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
               + Power(r, 8))
                  * (-2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                     + 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                     + 2 * Power(a, 2) * z * Power(r, 6)
                           * (-(x * y) + 3 * a * y * dr_dy + 3 * a * x * dr_dx)
                     + 4 * a * Power(r, 8)
                           * ((x - y) * (x + y) * dr_dz + z * (y * dr_dy + x * dr_dx))
                     - Power(a, 6) * Power(z, 2) * r
                           * (3 * x * y * dr_dz + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                     + Power(a, 2) * Power(r, 5)
                           * (-2 * a * z * (Power(a, 2) - Power(x, 2) + Power(y, 2) + Power(z, 2))
                              + x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                              + (a - z) * z * (a + z) * (x * dr_dy - y * dr_dx))
                     + 2 * Power(a, 2) * Power(r, 7)
                           * (3 * x * y * dr_dz - 2 * z * (a - x * dr_dy + y * dr_dx))
                     - Power(r, 9)
                           * (3 * x * y * dr_dz + z * (2 * a - 3 * x * dr_dy + 3 * y * dr_dx))
                     + 2 * Power(a, 4) * z * Power(r, 3)
                           * (a * (Power(x, 2) - Power(y, 2) - 2 * Power(z, 2))
                              + z * (3 * x * y * dr_dz - 2 * x * z * dr_dy + 2 * y * z * dr_dx))
                     + 2 * Power(a, 5) * z * Power(r, 2)
                           * (a * x * y
                              - z
                                    * (2 * (x - y) * (x + y) * dr_dz
                                       + z * (y * dr_dy + x * dr_dx))))))
        / (Power(Power(a, 2) + Power(r, 2), 4) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
           * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
              + Power(a, 2) * Power(z, 2) * Power(r, 2)
              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
              + Power(a, 2) * Power(r, 4) + Power(r, 6)));

  Gamma[1][1][1]
      = (M * Power(r, 2) * (a * x - y * r)
         * (2 * M * z * r * (Power(a, 2) + Power(r, 2))
                * (2 * z * Power(r, 2) * Power(Power(a, 2) + Power(r, 2), 2)
                       * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                   + 4 * Power(r, 3) * Power(a * x - y * r, 2)
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   + 2 * y * Power(r, 2) * (a * x - y * r) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   - 3 * r * Power(a * x - y * r, 2) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   + 2 * Power(r, 2) * Power(a * x - y * r, 2) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz)
                   + 2 * z * (Power(a, 2) + Power(r, 2))
                         * (-2 * Power(a, 5) * x * Power(z, 2)
                            + 3 * Power(a, 4) * y * Power(z, 2) * r
                            + Power(a, 2) * y * Power(z, 2) * Power(r, 3)
                            + 2 * Power(a, 3) * x * Power(r, 4) - Power(a, 2) * y * Power(r, 5)
                            + 4 * a * x * Power(r, 6) - 3 * y * Power(r, 7))
                         * dr_dy)
            - (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
               + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
               + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
               + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
               + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
               + Power(r, 8))
                  * (2 * Power(r, 2) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                     + (-3 * Power(a, 5) * x * Power(z, 2) + 5 * Power(a, 4) * y * Power(z, 2) * r
                        + Power(a, 3) * x * Power(z, 2) * Power(r, 2)
                        + Power(a, 2) * y * Power(z, 2) * Power(r, 3)
                        + Power(a, 3) * x * Power(r, 4) + Power(a, 2) * y * Power(r, 5)
                        + 5 * a * x * Power(r, 6) - 3 * y * Power(r, 7))
                           * dr_dy)
            + 2 * M * Power(r, 3) * (a * y + x * r)
                  * (-2 * Power(a, 4) * x * Power(z, 2) * Power(r, 3)
                     + 2 * Power(a, 2) * x * (-2 * Power(a, 2) + Power(z, 2)) * Power(r, 5)
                     + 2 * x * Power(r, 9)
                     - 3 * Power(a, 6) * x * Power(z, 2) * (2 * y * dr_dy + x * dr_dx)
                     + 3 * y * Power(r, 8) * (2 * a - 2 * x * dr_dy + y * dr_dx)
                     + Power(a, 4) * Power(z, 2) * Power(r, 2)
                           * (6 * y * (a + 2 * x * dr_dy) + (Power(x, 2) - 5 * Power(y, 2)) * dr_dx)
                     + Power(a, 2) * Power(r, 6)
                           * (6 * y * (a + 2 * x * dr_dy) + (5 * Power(x, 2) - Power(y, 2)) * dr_dx)
                     + Power(a, 2) * Power(r, 4)
                           * (6 * a * y * Power(z, 2)
                              + 2 * x * y * (Power(a, 2) + Power(z, 2)) * dr_dy
                              + (a * x - y * z) * (a * x + y * z) * dr_dx)
                     - 4 * Power(a, 5) * Power(z, 2) * r
                           * (2 * (x - y) * (x + y) * dr_dy + x * (a - 2 * y * dr_dx))
                     - 2 * a * Power(r, 7)
                           * (4 * (-Power(x, 2) + Power(y, 2)) * dr_dy + x * (a + 4 * y * dr_dx)))))
        / (Power(Power(a, 2) + Power(r, 2), 4) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
           * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
              + Power(a, 2) * Power(z, 2) * Power(r, 2)
              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
              + Power(a, 2) * Power(r, 4) + Power(r, 6)));

  Gamma[1][1][2]
      = (M * Power(r, 2) * (a * x - y * r)
         * (-((Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
               + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
               + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
               + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
               + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5) + 2 * Power(a, 2) * Power(r, 6)
               + Power(r, 8))
              * (2 * Power(a, 2) * z * r * (a * x - y * r) * (Power(a, 2) + Power(r, 2))
                 + (-3 * Power(a, 5) * x * Power(z, 2) + 5 * Power(a, 4) * y * Power(z, 2) * r
                    + Power(a, 3) * x * Power(z, 2) * Power(r, 2)
                    + Power(a, 2) * y * Power(z, 2) * Power(r, 3) + Power(a, 3) * x * Power(r, 4)
                    + Power(a, 2) * y * Power(r, 5) + 5 * a * x * Power(r, 6) - 3 * y * Power(r, 7))
                       * dr_dz))
            + 2 * M * Power(z, 3) * Power(Power(a, 2) + Power(r, 2), 4)
                  * (Power(a, 2) * Power(z, 2) - 3 * Power(r, 4)) * dr_dy
            + 2 * M * Power(r, 2) * (a * y + x * r)
                  * (2 * Power(a, 7) * Power(z, 3) * (y * dr_dy + x * dr_dx)
                     - 2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy + x * dr_dx)
                     + Power(r, 9)
                           * (2 * a * z - 3 * x * (y * dr_dz + z * dr_dy) + 3 * y * z * dr_dx)
                     - 2 * Power(a, 2) * z * Power(r, 6)
                           * (3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                     + 4 * a * Power(r, 8)
                           * ((x - y) * (x + y) * dr_dz - z * (y * dr_dy + x * dr_dx))
                     + 2 * Power(a, 5) * z * Power(r, 2)
                           * (a * x * y + 2 * (-Power(x, 2) + Power(y, 2)) * z * dr_dz
                              + Power(z, 2) * (y * dr_dy + x * dr_dx))
                     + Power(a, 6) * Power(z, 2) * r
                           * (-3 * x * y * dr_dz + z * (2 * a + 3 * x * dr_dy - 3 * y * dr_dx))
                     + Power(a, 2) * Power(r, 5)
                           * (2 * a * z * (Power(a, 2) + Power(x, 2) - Power(y, 2) + Power(z, 2))
                              + x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                              + z * (-Power(a, 2) + Power(z, 2)) * (x * dr_dy - y * dr_dx))
                     + 2 * Power(a, 2) * Power(r, 7)
                           * (3 * x * y * dr_dz + 2 * z * (a - x * dr_dy + y * dr_dx))
                     + 2 * Power(a, 4) * z * Power(r, 3)
                           * (a * (Power(x, 2) - Power(y, 2) + 2 * Power(z, 2))
                              + z * (3 * x * y * dr_dz + 2 * x * z * dr_dy - 2 * y * z * dr_dx)))))
        / (Power(Power(a, 2) + Power(r, 2), 4) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
           * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
              + Power(a, 2) * Power(z, 2) * Power(r, 2)
              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
              + Power(a, 2) * Power(r, 4) + Power(r, 6)));

  Gamma[1][2][2]
      = (M
         * (2 * y * Power(r, 15) + 4 * a * Power(r, 13) * (a * y + 2 * x * z * dr_dz)
            - Power(a, 10) * Power(z, 6) * dr_dy
            - 2 * Power(a, 8) * Power(z, 5) * r * (2 * a * x * dr_dz + M * z * dr_dy)
            + Power(r, 14) * (-2 * a * x - 6 * y * z * dr_dz + 3 * Power(z, 2) * dr_dy)
            + Power(a, 2) * Power(r, 12)
                  * (-4 * a * x - 8 * y * z * dr_dz + 9 * Power(z, 2) * dr_dy)
            + 2 * Power(a, 2) * z * Power(r, 9)
                  * (2
                         * (Power(a, 3) * x + 2 * a * x * Power(z, 2)
                            + M * y * (Power(x, 2) + Power(y, 2) - Power(z, 2)))
                         * dr_dz
                     + 3 * M * z * (Power(x, 2) + Power(y, 2) + 3 * Power(z, 2)) * dr_dy)
            + 2 * Power(a, 4) * Power(z, 3) * Power(r, 5)
                  * (-2 * a * (-2 * M * x + a * y) * z
                     + 2 * M * y * (Power(x, 2) + Power(y, 2) + 3 * Power(z, 2)) * dr_dz
                     - M * z * (-3 * Power(a, 2) + Power(x, 2) + Power(y, 2) + 3 * Power(z, 2))
                           * dr_dy)
            + Power(a, 7) * Power(z, 4) * Power(r, 2)
                  * (6 * (-(M * x) + a * y) * z * dr_dz + a * (2 * a * x - 3 * Power(z, 2) * dr_dy))
            - 2 * Power(a, 6) * Power(z, 4) * Power(r, 3)
                  * (a * (-2 * M * x + a * y) + (2 * a * x - 5 * M * y) * z * dr_dz
                     + M * (Power(y, 2) + 3 * Power(z, 2)) * dr_dy + M * x * y * dr_dx)
            + Power(a, 5) * Power(z, 4) * Power(r, 4)
                  * (4 * a * (a * x - M * y) + (-4 * M * x * z + 8 * a * y * z) * dr_dz
                     + (2 * Power(a, 3) - 4 * M * x * y - 3 * a * Power(z, 2)) * dr_dy
                     - 2 * M * (x - y) * (x + y) * dr_dx)
            + 2 * Power(a, 2) * Power(z, 2) * Power(r, 7)
                  * (a * (2 * M * x - a * y) * Power(z, 2)
                     + z * (4 * Power(a, 3) * x + Power(a, 2) * M * y + M * y * Power(z, 2)) * dr_dz
                     + M
                           * (-(Power(z, 2) * (Power(x, 2) + Power(z, 2)))
                              + 3 * Power(a, 2) * (Power(y, 2) + 3 * Power(z, 2)))
                           * dr_dy
                     + M * x * y * (3 * Power(a, 2) + Power(z, 2)) * dr_dx)
            + 2 * Power(r, 11)
                  * (Power(a, 4) * y + (6 * Power(a, 3) * x * z - 3 * M * y * Power(z, 3)) * dr_dz
                     + 3 * M * Power(z, 2) * ((Power(x, 2) + Power(z, 2)) * dr_dy - x * y * dr_dx))
            + Power(a, 3) * Power(z, 3) * Power(r, 6)
                  * (2
                         * (Power(a, 2) * M * x + 2 * Power(a, 3) * y + a * y * Power(z, 2)
                            + M * x * (2 * (Power(x, 2) + Power(y, 2)) + Power(z, 2)))
                         * dr_dz
                     + z
                           * (2 * a * (a * x - 4 * M * y)
                              + (6 * Power(a, 3) - 4 * M * x * y - a * Power(z, 2)) * dr_dy
                              - 2 * M * (x - y) * (x + y) * dr_dx))
            + Power(a, 2) * Power(z, 2) * Power(r, 8)
                  * (-4 * M * y * Power(z, 2)
                     + 3 * a
                           * (4 * M * x * z * dr_dz
                              + (Power(a, 3) + 4 * M * x * y + 2 * a * Power(z, 2)) * dr_dy
                              + 2 * M * (x - y) * (x + y) * dr_dx))
            + a * Power(r, 10)
                  * (-2 * Power(a, 4) * x
                     + z
                           * (2
                                  * (-(Power(a, 3) * y) + 2 * M * x * (Power(x, 2) + Power(y, 2))
                                     + (5 * M * x - 2 * a * y) * Power(z, 2))
                                  * dr_dz
                              + z * (9 * Power(a, 3) + 12 * M * x * y + 2 * a * Power(z, 2)) * dr_dy
                              + 6 * M * (x - y) * (x + y) * z * dr_dx))))
        / (Power(Power(a, 2) + Power(r, 2), 2) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
           * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
              + Power(a, 2) * Power(z, 2) * Power(r, 2)
              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
              + Power(a, 2) * Power(r, 4) + Power(r, 6)));

  Gamma[2][0][0]
      = (M * r
         * (2 * M * z * Power(r, 3) * Power(a * y + x * r, 2)
                * (-2 * Power(r, 2) * (Power(a, 2) + Power(r, 2))
                       * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                   + (-3 * Power(a, 5) * y * Power(z, 2) - 5 * Power(a, 4) * x * Power(z, 2) * r
                      + Power(a, 3) * y * Power(z, 2) * Power(r, 2)
                      - Power(a, 2) * x * Power(z, 2) * Power(r, 3) + Power(a, 3) * y * Power(r, 4)
                      - Power(a, 2) * x * Power(r, 5) + 5 * a * y * Power(r, 6)
                      + 3 * x * Power(r, 7))
                         * dr_dx)
            + (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
               + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3) + Power(a, 2) * Power(r, 4)
               + Power(r, 6))
                  * (4 * Power(r, 3) * Power(a * y + x * r, 2)
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                     - 2 * x * Power(r, 2) * (a * y + x * r) * (Power(a, 2) + Power(r, 2))
                           * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                     - 3 * r * Power(a * y + x * r, 2) * (Power(a, 2) + Power(r, 2))
                           * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                     + 2 * Power(r, 2) * Power(a * y + x * r, 2) * (Power(a, 2) + Power(r, 2))
                           * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz)
                     - 8 * z * Power(r, 4) * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 2)
                           * dr_dx
                     - 4 * z * Power(r, 2) * (a * y + x * r) * (Power(a, 2) + Power(r, 2))
                           * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                     + 4 * z * (a * y + x * r) * Power(Power(a, 2) + Power(r, 2), 2)
                           * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx
                     + 2 * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                           * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * (r + x * dr_dx))
            + 2 * M * z * Power(r, 3) * (a * x - y * r)
                  * (-2 * Power(a, 4) * y * Power(z, 2) * Power(r, 3)
                     + 2 * Power(a, 2) * y * (-2 * Power(a, 2) + Power(z, 2)) * Power(r, 5)
                     + 2 * y * Power(r, 9)
                     - 3 * Power(a, 6) * y * Power(z, 2) * (y * dr_dy + 2 * x * dr_dx)
                     - 4 * Power(a, 5) * Power(z, 2) * r
                           * (a * y + 2 * x * y * dr_dy + 2 * (x - y) * (x + y) * dr_dx)
                     - 2 * a * Power(r, 7)
                           * (a * y - 4 * x * y * dr_dy + 4 * (-x + y) * (x + y) * dr_dx)
                     + Power(a, 2) * Power(r, 4)
                           * (-6 * a * x * Power(z, 2) + (a * y - x * z) * (a * y + x * z) * dr_dy
                              + 2 * x * y * (Power(a, 2) + Power(z, 2)) * dr_dx)
                     + Power(a, 2) * Power(r, 6)
                           * (-((Power(x, 2) - 5 * Power(y, 2)) * dr_dy)
                              - 6 * x * (a - 2 * y * dr_dx))
                     + Power(a, 4) * Power(z, 2) * Power(r, 2)
                           * ((-5 * Power(x, 2) + Power(y, 2)) * dr_dy
                              - 6 * x * (a - 2 * y * dr_dx))
                     + 3 * x * Power(r, 8) * (x * dr_dy - 2 * (a + y * dr_dx)))))
        / (Power(Power(a, 2) + Power(r, 2), 3) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
           * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
              + Power(a, 2) * Power(z, 2) * Power(r, 2)
              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
              + Power(a, 2) * Power(r, 4) + Power(r, 6)));

  Gamma[2][0][1]
      = (M * r
         * (2 * M * z * Power(r, 3) * Power(a * y + x * r, 2)
                * (-2 * a * r * (Power(a, 2) + Power(r, 2))
                       * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                   + (-3 * Power(a, 5) * y * Power(z, 2) - 5 * Power(a, 4) * x * Power(z, 2) * r
                      + Power(a, 3) * y * Power(z, 2) * Power(r, 2)
                      - Power(a, 2) * x * Power(z, 2) * Power(r, 3) + Power(a, 3) * y * Power(r, 4)
                      - Power(a, 2) * x * Power(r, 5) + 5 * a * y * Power(r, 6)
                      + 3 * x * Power(r, 7))
                         * dr_dy)
            + 2 * M * z * Power(r, 3) * Power(a * x - y * r, 2)
                  * (2 * a * r * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                     - (-3 * Power(a, 5) * x * Power(z, 2) + 5 * Power(a, 4) * y * Power(z, 2) * r
                        + Power(a, 3) * x * Power(z, 2) * Power(r, 2)
                        + Power(a, 2) * y * Power(z, 2) * Power(r, 3)
                        + Power(a, 3) * x * Power(r, 4) + Power(a, 2) * y * Power(r, 5)
                        + 5 * a * x * Power(r, 6) - 3 * y * Power(r, 7))
                           * dr_dx)
            - (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
               + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3) + Power(a, 2) * Power(r, 4)
               + Power(r, 6))
                  * (2 * Power(a, 5) * z * Power(r, 4) * (y * dr_dy - x * dr_dx)
                     + 2 * Power(a, 7) * Power(z, 3) * (-(y * dr_dy) + x * dr_dx)
                     + 3 * Power(r, 9) * (-(x * y * dr_dz) + x * z * dr_dy + y * z * dr_dx)
                     - 3 * Power(a, 6) * Power(z, 2) * r
                           * (x * y * dr_dz + x * z * dr_dy + y * z * dr_dx)
                     - 2 * Power(a, 2) * z * Power(r, 6)
                           * (-3 * a * y * dr_dy + x * (y + 3 * a * dr_dx))
                     + 4 * a * Power(r, 8)
                           * ((x - y) * (x + y) * dr_dz + z * (y * dr_dy - x * dr_dx))
                     + 2 * Power(a, 2) * Power(r, 7)
                           * (3 * x * y * dr_dz + 2 * z * (x * dr_dy + y * dr_dx))
                     + 2 * Power(a, 4) * z * Power(r, 3)
                           * (a * (x - y) * (x + y) + 3 * x * y * z * dr_dz
                              - 2 * Power(z, 2) * (x * dr_dy + y * dr_dx))
                     + Power(a, 2) * Power(r, 5)
                           * (2 * a * (x - y) * (x + y) * z
                              + x * y * (Power(a, 2) + Power(z, 2)) * dr_dz
                              + (a - z) * z * (a + z) * (x * dr_dy + y * dr_dx))
                     + 2 * Power(a, 5) * z * Power(r, 2)
                           * (a * x * y
                              + z
                                    * (2 * (-Power(x, 2) + Power(y, 2)) * dr_dz - y * z * dr_dy
                                       + x * z * dr_dx)))))
        / (Power(Power(a, 2) + Power(r, 2), 3) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
           * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
              + Power(a, 2) * Power(z, 2) * Power(r, 2)
              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
              + Power(a, 2) * Power(r, 4) + Power(r, 6)));

  Gamma[2][0][2]
      = (M * z
         * (Power(a, 10) * Power(z, 5) * dr_dx + 3 * Power(a, 8) * Power(z, 5) * Power(r, 2) * dr_dx
            - 9 * Power(a, 2) * z * Power(r, 12) * dr_dx - 3 * z * Power(r, 14) * dr_dx
            - 2 * Power(a, 6) * M * Power(z, 3) * Power(r, 3)
                  * (2 * x * y * dr_dy + (x - y) * (x + y) * dr_dx)
            + 2 * Power(a, 2) * M * z * Power(r, 9)
                  * (2 * a * y + 3 * x * y * dr_dy - (2 * Power(x, 2) + 5 * Power(y, 2)) * dr_dx)
            + Power(a, 2) * Power(r, 8)
                  * (2 * M
                         * (2 * x * (-Power(a, 2) + Power(x, 2) + Power(y, 2)) * z
                            + a * y * (Power(x, 2) + Power(y, 2)) * dr_dz
                            + a * (Power(x, 2) - 2 * Power(y, 2)) * z * dr_dy)
                     - 3 * a * z * (Power(a, 3) + 2 * M * x * y + 2 * a * Power(z, 2)) * dr_dx)
            + 2 * Power(a, 2) * M * z * Power(r, 7)
                  * (2 * a * y * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                     - x * (Power(x, 2) + Power(y, 2)) * z * dr_dz
                     + x * y * (2 * Power(a, 2) + Power(z, 2)) * dr_dy
                     + (-(Power(a, 2) * (Power(x, 2) + 3 * Power(y, 2)))
                        + Power(x, 2) * Power(z, 2))
                           * dr_dx)
            + 2 * M * Power(r, 11)
                  * (2 * a * y * z + 3 * x * (Power(x, 2) + Power(y, 2)) * dr_dz
                     - 3 * x * z * (y * dr_dy + x * dr_dx))
            + 2 * Power(a, 4) * M * z * Power(r, 5)
                  * (2 * a * y * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                     + z
                           * (-4 * x * (Power(x, 2) + Power(y, 2)) * dr_dz + 3 * x * y * z * dr_dy
                              + (2 * Power(x, 2) - Power(y, 2)) * z * dr_dx))
            + Power(a, 3) * z * Power(r, 6)
                  * (4 * a * M * x * (Power(x, 2) + Power(y, 2) - Power(z, 2))
                     + Power(z, 2)
                           * (-2 * M * Power(x, 2) * dr_dy
                              + (-6 * Power(a, 3) + 2 * M * x * y + a * Power(z, 2)) * dr_dx))
            + a * Power(r, 10)
                  * (8 * M * y * (Power(x, 2) + Power(y, 2)) * dr_dz
                     + z
                           * (-4 * a * M * x + M * (6 * Power(x, 2) - 8 * Power(y, 2)) * dr_dy
                              - (9 * Power(a, 3) + 14 * M * x * y + 2 * a * Power(z, 2)) * dr_dx))
            + Power(a, 5) * Power(z, 2) * Power(r, 4)
                  * (-6 * M * y * (Power(x, 2) + Power(y, 2)) * dr_dz
                     + z
                           * (-4 * a * M * x + M * (-6 * Power(x, 2) + 4 * Power(y, 2)) * dr_dy
                              + (-2 * Power(a, 3) + 10 * M * x * y + 3 * a * Power(z, 2))
                                    * dr_dx))))
        / (Power(Power(a, 2) + Power(r, 2), 2) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
           * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
              + Power(a, 2) * Power(z, 2) * Power(r, 2)
              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
              + Power(a, 2) * Power(r, 4) + Power(r, 6)));

  Gamma[2][1][1] = -(
      (M * r
       * (2 * M * z * Power(r, 3) * Power(a * x - y * r, 2)
              * (2 * Power(r, 2) * (Power(a, 2) + Power(r, 2))
                     * (Power(a, 2) * Power(z, 2) + Power(r, 4))
                 + (-3 * Power(a, 5) * x * Power(z, 2) + 5 * Power(a, 4) * y * Power(z, 2) * r
                    + Power(a, 3) * x * Power(z, 2) * Power(r, 2)
                    + Power(a, 2) * y * Power(z, 2) * Power(r, 3) + Power(a, 3) * x * Power(r, 4)
                    + Power(a, 2) * y * Power(r, 5) + 5 * a * x * Power(r, 6) - 3 * y * Power(r, 7))
                       * dr_dy)
          - (Power(a, 4) * Power(z, 2) + Power(a, 2) * Power(z, 2) * Power(r, 2)
             + 2 * M * (Power(x, 2) + Power(y, 2)) * Power(r, 3) + Power(a, 2) * Power(r, 4)
             + Power(r, 6))
                * (4 * Power(r, 3) * Power(a * x - y * r, 2)
                       * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   + 2 * y * Power(r, 2) * (a * x - y * r) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   - 3 * r * Power(a * x - y * r, 2) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dz
                   + 2 * Power(r, 2) * Power(a * x - y * r, 2) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * z + 2 * Power(r, 3) * dr_dz)
                   + 8 * z * Power(r, 4) * (a * x - y * r) * Power(Power(a, 2) + Power(r, 2), 2)
                         * dr_dy
                   + 4 * z * Power(r, 2) * (a * x - y * r) * (Power(a, 2) + Power(r, 2))
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                   - 4 * z * (a * x - y * r) * Power(Power(a, 2) + Power(r, 2), 2)
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy
                   + 2 * z * r * Power(Power(a, 2) + Power(r, 2), 2)
                         * (Power(a, 2) * Power(z, 2) + Power(r, 4)) * (r + y * dr_dy))
          + 2 * M * z * Power(r, 3) * (a * y + x * r)
                * (-2 * Power(a, 4) * x * Power(z, 2) * Power(r, 3)
                   + 2 * Power(a, 2) * x * (-2 * Power(a, 2) + Power(z, 2)) * Power(r, 5)
                   + 2 * x * Power(r, 9)
                   - 3 * Power(a, 6) * x * Power(z, 2) * (2 * y * dr_dy + x * dr_dx)
                   + 3 * y * Power(r, 8) * (2 * a - 2 * x * dr_dy + y * dr_dx)
                   + Power(a, 4) * Power(z, 2) * Power(r, 2)
                         * (6 * y * (a + 2 * x * dr_dy) + (Power(x, 2) - 5 * Power(y, 2)) * dr_dx)
                   + Power(a, 2) * Power(r, 6)
                         * (6 * y * (a + 2 * x * dr_dy) + (5 * Power(x, 2) - Power(y, 2)) * dr_dx)
                   + Power(a, 2) * Power(r, 4)
                         * (6 * a * y * Power(z, 2)
                            + 2 * x * y * (Power(a, 2) + Power(z, 2)) * dr_dy
                            + (a * x - y * z) * (a * x + y * z) * dr_dx)
                   - 4 * Power(a, 5) * Power(z, 2) * r
                         * (2 * (x - y) * (x + y) * dr_dy + x * (a - 2 * y * dr_dx))
                   - 2 * a * Power(r, 7)
                         * (4 * (-Power(x, 2) + Power(y, 2)) * dr_dy + x * (a + 4 * y * dr_dx)))))
      / (Power(Power(a, 2) + Power(r, 2), 3) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
         * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
            + Power(a, 2) * Power(z, 2) * Power(r, 2)
            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
            + Power(a, 2) * Power(r, 4) + Power(r, 6))));

  Gamma[2][1][2] = -(
      (M * z
       * (-(Power(a, 10) * Power(z, 5) * dr_dy)
          - 3 * Power(a, 8) * Power(z, 5) * Power(r, 2) * dr_dy
          + 9 * Power(a, 2) * z * Power(r, 12) * dr_dy + 3 * z * Power(r, 14) * dr_dy
          + 2 * Power(a, 2) * M * z * Power(r, 9)
                * (2 * a * x + (5 * Power(x, 2) + 2 * Power(y, 2)) * dr_dy - 3 * x * y * dr_dx)
          + 2 * Power(a, 6) * M * Power(z, 3) * Power(r, 3)
                * ((-Power(x, 2) + Power(y, 2)) * dr_dy + 2 * x * y * dr_dx)
          + 2 * Power(a, 2) * M * z * Power(r, 7)
                * (2 * a * x * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                   + y * (Power(x, 2) + Power(y, 2)) * z * dr_dz
                   + (Power(a, 2) * (3 * Power(x, 2) + Power(y, 2)) - Power(y, 2) * Power(z, 2))
                         * dr_dy
                   - x * y * (2 * Power(a, 2) + Power(z, 2)) * dr_dx)
          + 2 * M * Power(r, 11)
                * (2 * a * x * z - 3 * y * (Power(x, 2) + Power(y, 2)) * dr_dz
                   + 3 * y * z * (y * dr_dy + x * dr_dx))
          + Power(a, 5) * Power(z, 2) * Power(r, 4)
                * (-6 * M * x * (Power(x, 2) + Power(y, 2)) * dr_dz
                   + z
                         * (4 * a * M * y
                            + (2 * Power(a, 3) + 10 * M * x * y - 3 * a * Power(z, 2)) * dr_dy
                            + M * (4 * Power(x, 2) - 6 * Power(y, 2)) * dr_dx))
          + Power(a, 2) * Power(r, 8)
                * (2 * a * M * x * (Power(x, 2) + Power(y, 2)) * dr_dz
                   + z
                         * (-4 * M * y * (-Power(a, 2) + Power(x, 2) + Power(y, 2))
                            + 3 * a * (Power(a, 3) - 2 * M * x * y + 2 * a * Power(z, 2)) * dr_dy
                            + 2 * a * M * (-2 * Power(x, 2) + Power(y, 2)) * dr_dx))
          + a * Power(r, 10)
                * (8 * M * x * (Power(x, 2) + Power(y, 2)) * dr_dz
                   + z
                         * (4 * a * M * y
                            + (9 * Power(a, 3) - 14 * M * x * y + 2 * a * Power(z, 2)) * dr_dy
                            + M * (-8 * Power(x, 2) + 6 * Power(y, 2)) * dr_dx))
          + 2 * Power(a, 4) * M * z * Power(r, 5)
                * (2 * a * x * (Power(x, 2) + Power(y, 2) + Power(z, 2))
                   + z
                         * (4 * y * (Power(x, 2) + Power(y, 2)) * dr_dz
                            + (Power(x, 2) - 2 * Power(y, 2)) * z * dr_dy - 3 * x * y * z * dr_dx))
          + Power(a, 3) * z * Power(r, 6)
                * (Power(z, 2) * (6 * Power(a, 3) + 2 * M * x * y - a * Power(z, 2)) * dr_dy
                   - 2 * M * y
                         * (2 * a * (Power(x, 2) + Power(y, 2) - Power(z, 2))
                            + y * Power(z, 2) * dr_dx))))
      / (Power(Power(a, 2) + Power(r, 2), 2) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
         * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
            + Power(a, 2) * Power(z, 2) * Power(r, 2)
            + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
            + Power(a, 2) * Power(r, 4) + Power(r, 6))));

  Gamma[2][2][2]
      = (M * z
         * (4 * Power(a, 2) * Power(r, 11) + 2 * Power(r, 13) + Power(a, 8) * Power(z, 5) * dr_dz
            - 6 * Power(a, 2) * z * Power(r, 10) * dr_dz - 3 * z * Power(r, 12) * dr_dz
            + 2 * Power(a, 5) * Power(z, 4) * Power(r, 2)
                  * (a * z * dr_dz - M * x * dr_dy + M * y * dr_dx)
            + 2 * Power(a, 4) * M * Power(z, 3) * Power(r, 3)
                  * (-3 * (Power(x, 2) + Power(y, 2)) * dr_dz + z * (y * dr_dy + x * dr_dx))
            + 2 * Power(a, 2) * z * Power(r, 7)
                  * (2 * Power(a, 2) * z + M * (Power(x, 2) + Power(y, 2)) * dr_dz
                     - 3 * M * z * (y * dr_dy + x * dr_dx))
            - a * z * Power(r, 8)
                  * ((3 * Power(a, 3) + 2 * a * Power(z, 2)) * dr_dz
                     + 6 * M * z * (-(x * dr_dy) + y * dr_dx))
            + Power(a, 3) * Power(z, 2) * Power(r, 4)
                  * (4 * a * M * (Power(x, 2) + Power(y, 2))
                     + a * z * (-2 * Power(a, 2) + Power(z, 2)) * dr_dz
                     + 2 * M * Power(z, 2) * (-(x * dr_dy) + y * dr_dx))
            + 2 * Power(a, 2) * Power(z, 2) * Power(r, 6)
                  * (2 * M * (Power(x, 2) + Power(y, 2))
                     + a * (-2 * a * z * dr_dz + 3 * M * x * dr_dy - 3 * M * y * dr_dx))
            + 2 * Power(r, 9)
                  * (Power(a, 2) * (Power(a, 2) + Power(z, 2))
                     + 3 * M * z
                           * ((Power(x, 2) + Power(y, 2)) * dr_dz - z * (y * dr_dy + x * dr_dx)))
            + 2 * Power(a, 2) * Power(z, 2) * Power(r, 5)
                  * (Power(a, 4)
                     + M * z
                           * (-((Power(x, 2) + Power(y, 2)) * dr_dz)
                              + z * (y * dr_dy + x * dr_dx)))))
        / ((Power(a, 2) + Power(r, 2)) * Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2)
           * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
              + Power(a, 2) * Power(z, 2) * Power(r, 2)
              + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
              + Power(a, 2) * Power(r, 4) + Power(r, 6)));

  Gamma[0][1][0] = Gamma[0][0][1];
  Gamma[0][2][0] = Gamma[0][0][2];
  Gamma[0][2][1] = Gamma[0][1][2];
  Gamma[1][1][0] = Gamma[1][0][1];
  Gamma[1][2][0] = Gamma[1][0][2];
  Gamma[1][2][1] = Gamma[1][1][2];
  Gamma[2][1][0] = Gamma[2][0][1];
  Gamma[2][2][0] = Gamma[2][0][2];
  Gamma[2][2][1] = Gamma[2][1][2];

  return Gamma;
}