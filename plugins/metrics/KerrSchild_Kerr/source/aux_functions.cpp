#include "aux_functions.hpp"

auto ksk_aux::r_KS(double a, double x, double y, double z) -> double {
  return Sqrt(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)
              + Sqrt(4 * Power(a, 2) * Power(z, 2)
                     + Power(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2), 2)))
         / Sqrt(2);
}

auto ksk_aux::d_r_KS_dx(double a, double x, double y, double z) -> double {
  return (x
          * Sqrt((-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)
                  + Sqrt(4 * Power(a, 2) * Power(z, 2)
                         + Power(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2), 2)))
                 / (4 * Power(a, 2) * Power(z, 2)
                    + Power(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2), 2))))
         / Sqrt(2);
}

auto ksk_aux::d_r_KS_dy(double a, double x, double y, double z) -> double {
  return (y
          * Sqrt((-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)
                  + Sqrt(4 * Power(a, 2) * Power(z, 2)
                         + Power(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2), 2)))
                 / (4 * Power(a, 2) * Power(z, 2)
                    + Power(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2), 2))))
         / Sqrt(2);
}

auto ksk_aux::d_r_KS_dz(double a, double x, double y, double z) -> double {
  return (z
          * (1
             + (Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2))
                   / Sqrt(4 * Power(a, 2) * Power(z, 2)
                          + Power(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2), 2))))
         / (Sqrt(2)
            * Sqrt(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)
                   + Sqrt(4 * Power(a, 2) * Power(z, 2)
                          + Power(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2), 2))));
}