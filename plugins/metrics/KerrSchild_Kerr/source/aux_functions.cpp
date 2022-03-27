#include "aux_functions.hpp"

auto ksk_aux::r_KS(double a, double x, double y, double z) -> double {
  double part1 = Power(x, 2) + Power(y, 2) + Power(z, 2) - Power(a, 2);
  return Sqrt(part1 + Sqrt(4 * Power(a, 2) * Power(z, 2) + Power(part1, 2))) / Sqrt(2);
}

auto ksk_aux::d_r_KS_dx(double a, double x, double y, double z) -> double {
  return (x
          * Sqrt(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)
                 + Sqrt(4 * Power(a, 2) * Power(z, 2)
                        + Power(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2), 2))))
         / (Sqrt(2)
            * Sqrt(4 * Power(a, 2) * Power(z, 2)
                   + Power(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2), 2)));
}

auto ksk_aux::d_r_KS_dy(double a, double x, double y, double z) -> double {
  return (y
          * Sqrt(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2)
                 + Sqrt(4 * Power(a, 2) * Power(z, 2)
                        + Power(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2), 2))))
         / (Sqrt(2)
            * Sqrt(4 * Power(a, 2) * Power(z, 2)
                   + Power(-Power(a, 2) + Power(x, 2) + Power(y, 2) + Power(z, 2), 2)));
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

auto ksk_aux::H_KS(double M, double a, double r, double z) -> double {
  return (M * r) / (Power(r, 2) + Power(a, 2) * Power(z / r, 2));
}

auto ksk_aux::d_H_KS_dx(double M, double a, double dr_dx, double r, double z) -> double {
  return -((M * Power(r, 2) * (-3 * Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dx)
           / Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2));
}

auto ksk_aux::d_H_KS_dy(double M, double a, double dr_dy, double r, double z) -> double {
  return -((M * Power(r, 2) * (-3 * Power(a, 2) * Power(z, 2) + Power(r, 4)) * dr_dy)
           / Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2));
}

auto ksk_aux::d_H_KS_dz(double M, double a, double dr_dz, double r, double z) -> double {
  return (M * Power(r, 2)
          * (-2 * Power(a, 2) * z * r + (3 * Power(a, 2) * Power(z, 2) - Power(r, 4)) * dr_dz))
         / Power(Power(a, 2) * Power(z, 2) + Power(r, 4), 2);
}

auto ksk_aux::l1_KS(double a, double r, double x, double y) -> double {
  return (r * x + a * y) / (Power(r, 2) + Power(a, 2));
}

auto ksk_aux::d_l1_KS_dx(double a, double dr_dx, double r, double x, double y) -> double {
  return (r * (Power(a, 2) + Power(r, 2)) + (Power(a, 2) * x - r * (2 * a * y + x * r)) * dr_dx)
         / Power(Power(a, 2) + Power(r, 2), 2);
}
auto ksk_aux::d_l1_KS_dy(double a, double dr_dy, double r, double x, double y) -> double {
  return (Power(a, 3) + a * Power(r, 2) + (Power(a, 2) * x - r * (2 * a * y + x * r)) * dr_dy)
         / Power(Power(a, 2) + Power(r, 2), 2);
}
auto ksk_aux::d_l1_KS_dz(double a, double dr_dz, double r, double x, double y) -> double {
  return ((Power(a, 2) * x - r * (2 * a * y + x * r)) * dr_dz)
         / Power(Power(a, 2) + Power(r, 2), 2);
}

auto ksk_aux::l2_KS(double a, double r, double x, double y) -> double {
  return (r * y - a * x) / (Power(r, 2) + Power(a, 2));
}

auto ksk_aux::d_l2_KS_dx(double a, double dr_dx, double r, double x, double y) -> double {
  return (-(a * (Power(a, 2) + Power(r, 2)))
          + (Power(a, 2) * y + 2 * a * x * r - y * Power(r, 2)) * dr_dx)
         / Power(Power(a, 2) + Power(r, 2), 2);
}

auto ksk_aux::d_l2_KS_dy(double a, double dr_dy, double r, double x, double y) -> double {
  return (r * (Power(a, 2) + Power(r, 2))
          + (Power(a, 2) * y + 2 * a * x * r - y * Power(r, 2)) * dr_dy)
         / Power(Power(a, 2) + Power(r, 2), 2);
}

auto ksk_aux::d_l2_KS_dz(double a, double dr_dz, double r, double x, double y) -> double {
  return ((Power(a, 2) * y + 2 * a * x * r - y * Power(r, 2)) * dr_dz)
         / Power(Power(a, 2) + Power(r, 2), 2);
}

auto ksk_aux::l3_KS(double r, double z) -> double { return z / r; }

auto ksk_aux::d_l3_KS_dx(double r, double dr_dx, double z) -> double {
  return -((z * dr_dx) / Power(r, 2));
}

auto ksk_aux::d_l3_KS_dy(double r, double dr_dy, double z) -> double {
  return -((z * dr_dy) / Power(r, 2));
}

auto ksk_aux::d_l3_KS_dz(double r, double dr_dz, double z) -> double {
  return (r - z * dr_dz) / Power(r, 2);
}