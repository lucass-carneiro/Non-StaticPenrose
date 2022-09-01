#ifndef GRLENSING_SKS_METRIC_HPP
#define GRLENSING_SKS_METRIC_HPP

#include "api_macros.hpp"

namespace grlensing {

class SKS final : public metric_server::adm_metric {
public:
  GRLENSING_SKS_METRIC_API void load_parameters(const YAML::Node &) final;

  GRLENSING_SKS_METRIC_API auto lapse(double, double, double, double) -> double final;

  GRLENSING_SKS_METRIC_API auto l_shift(double, double, double, double)
      -> metric_server::spatial_vector final;

  GRLENSING_SKS_METRIC_API auto u_shift(double, double, double, double)
      -> metric_server::spatial_vector final;

  GRLENSING_SKS_METRIC_API auto ll_smetric(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_SKS_METRIC_API auto uu_smetric(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_SKS_METRIC_API auto ll_extrinsic(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_SKS_METRIC_API auto ul_extrinsic(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_SKS_METRIC_API auto spatial_christoffel(double, double, double, double)
      -> metric_server::chirstofell_t final;

  GRLENSING_SKS_METRIC_API auto grad_lapse(double, double, double, double)
      -> metric_server::spatial_vector final;

  GRLENSING_SKS_METRIC_API auto grad_ushift(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_SKS_METRIC_API auto name() -> std::string_view final;

private:
  double M1 = 0.5;
  double M2 = 0.5;
  double a1 = 0.0;
  double a2 = 0.0;
  double b = 10.0;

  [[nodiscard]] auto sx(unsigned bhIdx, double t) const noexcept -> double;
  [[nodiscard]] auto sy(unsigned bhIdx, double t) const noexcept -> double;

  [[nodiscard]] auto dsx_dt(unsigned bhIdx, double t) const noexcept -> double;
  [[nodiscard]] auto dsy_dt(unsigned bhIdx, double t) const noexcept -> double;

  [[nodiscard]] auto d2sx_dt2(unsigned bhIdx, double t) const noexcept -> double;
  [[nodiscard]] auto d2sy_dt2(unsigned bhIdx, double t) const noexcept -> double;

  [[nodiscard]] auto T(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto X(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto Y(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto Z(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;

  [[nodiscard]] auto dT_dt(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto dT_dx(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto dT_dy(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto dT_dz(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;

  [[nodiscard]] auto dX_dt(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto dX_dx(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto dX_dy(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto dX_dz(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;

  [[nodiscard]] auto dY_dt(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto dY_dx(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto dY_dy(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto dY_dz(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;

  [[nodiscard]] auto dZ_dt(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto dZ_dx(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto dZ_dy(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;
  [[nodiscard]] auto dZ_dz(unsigned bhIdx, double t, double x, double y, double z) const noexcept
      -> double;

  [[nodiscard]] auto llgSKS_00(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto llgSKS_01(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto llgSKS_02(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto llgSKS_03(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto llgSKS_11(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto llgSKS_12(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto llgSKS_13(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto llgSKS_22(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto llgSKS_23(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto llgSKS_33(double t, double x, double y, double z) const noexcept -> double;

  [[nodiscard]] auto dllgSKS_00_dt(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_00_dx(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_00_dy(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_00_dz(double t, double x, double y, double z) const noexcept -> double;

  [[nodiscard]] auto dllgSKS_01_dt(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_01_dx(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_01_dy(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_01_dz(double t, double x, double y, double z) const noexcept -> double;

  [[nodiscard]] auto dllgSKS_02_dt(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_02_dx(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_02_dy(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_02_dz(double t, double x, double y, double z) const noexcept -> double;

  [[nodiscard]] auto dllgSKS_03_dt(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_03_dx(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_03_dy(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_03_dz(double t, double x, double y, double z) const noexcept -> double;

  [[nodiscard]] auto dllgSKS_11_dt(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_11_dx(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_11_dy(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_11_dz(double t, double x, double y, double z) const noexcept -> double;

  [[nodiscard]] auto dllgSKS_12_dt(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_12_dx(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_12_dy(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_12_dz(double t, double x, double y, double z) const noexcept -> double;

  [[nodiscard]] auto dllgSKS_13_dt(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_13_dx(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_13_dy(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_13_dz(double t, double x, double y, double z) const noexcept -> double;

  [[nodiscard]] auto dllgSKS_22_dt(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_22_dx(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_22_dy(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_22_dz(double t, double x, double y, double z) const noexcept -> double;

  [[nodiscard]] auto dllgSKS_23_dt(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_23_dx(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_23_dy(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_23_dz(double t, double x, double y, double z) const noexcept -> double;

  [[nodiscard]] auto dllgSKS_33_dt(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_33_dx(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_33_dy(double t, double x, double y, double z) const noexcept -> double;
  [[nodiscard]] auto dllgSKS_33_dz(double t, double x, double y, double z) const noexcept -> double;
};

extern "C" GRLENSING_SKS_METRIC_API auto get_engine_version() -> unsigned;

extern "C" GRLENSING_SKS_METRIC_API void register_plugin(kernel &kernel);

} // namespace grlensing

#endif // GRLENSING_SKS_METRIC_HPP ";