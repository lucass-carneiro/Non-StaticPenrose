#include "KerrSchild_Kerr.hpp"
#include "aux_functions.hpp"

using grlensing::KerrSchild_Kerr;
using grlensing::metric_server;
using ksk_aux::Power;
using ksk_aux::r_KS;

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::ll_smetric(double, double x, double y,
                                                                      double z)
    -> metric_server::spatial_matrix {

  auto r = r_KS(a, x, y, z);

  metric_server::spatial_matrix metric{};

  metric[0][0]
      = 1
        + (2 * M * Power(r, 3) * Power(a * y + x * r, 2))
              / (Power(Power(a, 2) + Power(r, 2), 2) * (Power(a, 2) * Power(z, 2) + Power(r, 4)));

  metric[0][1]
      = (2 * M * Power(r, 3) * (a * y + x * r) * (-(a * x) + y * r))
        / (Power(Power(a, 2) + Power(r, 2), 2) * (Power(a, 2) * Power(z, 2) + Power(r, 4)));

  metric[0][2] = (2 * M * z * Power(r, 2) * (a * y + x * r))
                 / ((Power(a, 2) + Power(r, 2)) * (Power(a, 2) * Power(z, 2) + Power(r, 4)));

  metric[1][1]
      = 1
        + (2 * M * Power(r, 3) * Power(a * x - y * r, 2))
              / (Power(Power(a, 2) + Power(r, 2), 2) * (Power(a, 2) * Power(z, 2) + Power(r, 4)));

  metric[1][2] = (2 * M * z * Power(r, 2) * (-(a * x) + y * r))
                 / ((Power(a, 2) + Power(r, 2)) * (Power(a, 2) * Power(z, 2) + Power(r, 4)));

  metric[2][2] = 1 + (2 * M * Power(z, 2) * r) / (Power(a, 2) * Power(z, 2) + Power(r, 4));

  metric[1][0] = metric[0][1];
  metric[2][0] = metric[0][2];
  metric[2][1] = metric[1][2];

  return metric;
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::uu_smetric(double, double x, double y,
                                                                      double z)
    -> metric_server::spatial_matrix {

  auto r = r_KS(a, x, y, z);

  metric_server::spatial_matrix imetric{};

  imetric[0][0] = (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                   + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                   + 2 * Power(a, 2) * M * (Power(x, 2) + 2 * Power(z, 2)) * Power(r, 3)
                   + a * (Power(a, 3) - 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                   + 2 * M * (Power(y, 2) + Power(z, 2)) * Power(r, 5)
                   + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                  / ((Power(a, 2) + Power(r, 2))
                     * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                        + Power(a, 2) * Power(z, 2) * Power(r, 2)
                        + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                        + Power(a, 2) * Power(r, 4) + Power(r, 6)));

  imetric[0][1] = (2 * M * Power(r, 3) * (a * y + x * r) * (a * x - y * r))
                  / ((Power(a, 2) + Power(r, 2))
                     * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                        + Power(a, 2) * Power(z, 2) * Power(r, 2)
                        + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                        + Power(a, 2) * Power(r, 4) + Power(r, 6)));

  imetric[0][2] = (-2 * M * z * Power(r, 2) * (a * y + x * r))
                  / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                     + Power(a, 2) * Power(z, 2) * Power(r, 2)
                     + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                     + Power(a, 2) * Power(r, 4) + Power(r, 6));

  imetric[1][1] = (Power(a, 6) * Power(z, 2) + 2 * Power(a, 4) * M * Power(z, 2) * r
                   + 2 * Power(a, 4) * Power(z, 2) * Power(r, 2)
                   + 2 * Power(a, 2) * M * (Power(y, 2) + 2 * Power(z, 2)) * Power(r, 3)
                   + a * (Power(a, 3) + 4 * M * x * y + a * Power(z, 2)) * Power(r, 4)
                   + 2 * M * (Power(x, 2) + Power(z, 2)) * Power(r, 5)
                   + 2 * Power(a, 2) * Power(r, 6) + Power(r, 8))
                  / ((Power(a, 2) + Power(r, 2))
                     * (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                        + Power(a, 2) * Power(z, 2) * Power(r, 2)
                        + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                        + Power(a, 2) * Power(r, 4) + Power(r, 6)));

  imetric[1][2] = (2 * M * z * Power(r, 2) * (a * x - y * r))
                  / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                     + Power(a, 2) * Power(z, 2) * Power(r, 2)
                     + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                     + Power(a, 2) * Power(r, 4) + Power(r, 6));

  imetric[2][2] = 1
                  - (2 * M * Power(z, 2) * r * (Power(a, 2) + Power(r, 2)))
                        / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                           + Power(a, 2) * Power(z, 2) * Power(r, 2)
                           + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                           + Power(a, 2) * Power(r, 4) + Power(r, 6));

  imetric[1][0] = imetric[0][1];
  imetric[2][0] = imetric[0][2];
  imetric[2][1] = imetric[1][2];

  return imetric;
}