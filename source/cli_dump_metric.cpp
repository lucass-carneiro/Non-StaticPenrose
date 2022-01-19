#include "cli.hpp"
#include "log.hpp"
#include "mpi_index_map_3D.hpp"

#include <concepts>
#include <memory>
#include <string>
#include <string_view>

using writer_ptr = std::unique_ptr<grlensing::storage_server::trajectory_writer>;
using metric_ptr = std::unique_ptr<grlensing::metric_server::adm_metric>;
using std_path = std::filesystem::path;

template <std::unsigned_integral index_t, grlensing::ordering_type ordering, bool bound_checked>
static void
dump_kernel(grlensing::kernel &kernel, const std::string &metric_name, const std::string &basename,
            const std::string &extension, double rf, double points,
            const grlensing::mpi_index_map_3D<index_t, ordering, bound_checked> &coords_map,
            void (*operation)(double, double, double, const metric_ptr &, const writer_ptr &)) {

  using grlensing::log;
  using grlensing::LogEvent;

  const auto &writer = kernel.get_storage_server().get_writer(extension);
  const auto &metric = kernel.get_metric_server().get_metric(metric_name);

  const std::filesystem::path data_file(metric_name + std::string("_") + basename + std::string("_")
                                        + std::to_string(MPI::COMM_WORLD.Get_rank()) + extension);
  const std::filesystem::path metadata_file(
      metric_name + std::string("_") + basename + std::string("_metadata_")
      + std::to_string(MPI::COMM_WORLD.Get_rank()) + extension);

  writer->open_file(data_file, metadata_file);

  if (points <= 1) {
    log<LogEvent::error>("The number of points in a metric dump must be at least 2");
    throw std::runtime_error("Invalid parameter error");
  }

  auto r0 = -rf;
  const auto dr = 2 * rf / points;

  const auto glor = coords_map.get_glor();
  const auto start_ijk = coords_map.global_linear_to_matrix(glor.first);
  const auto end_ijk = coords_map.global_linear_to_matrix(glor.second);

  for (auto i = std::get<0>(start_ijk); i <= std::get<0>(end_ijk); i++) {
    for (auto j = std::get<1>(start_ijk); j <= std::get<1>(end_ijk); j++) {
      for (auto k = std::get<2>(start_ijk); k <= std::get<2>(end_ijk); k++) {
        const auto x = r0 + i * dr;
        const auto y = r0 + j * dr;
        const auto z = r0 + k * dr;

        writer->push_real(x);
        writer->push_real(y);
        writer->push_real(z);

        operation(x, y, z, metric, writer);
      }
    }
  }

  writer->push_metadata("glor", std::size_t(glor.first));

  writer->push_metadata("start_i", std::size_t(std::get<0>(start_ijk)));
  writer->push_metadata("start_j", std::size_t(std::get<1>(start_ijk)));
  writer->push_metadata("start_k", std::size_t(std::get<2>(start_ijk)));

  writer->push_metadata("end_i", std::size_t(std::get<0>(end_ijk)));
  writer->push_metadata("end_j", std::size_t(std::get<1>(end_ijk)));
  writer->push_metadata("end_k", std::size_t(std::get<2>(end_ijk)));

  writer->close_file();
}

void grlensing::dump_metric(grlensing::kernel &kernel, const std::string &metric_name,
                            const YAML::Node &config_file) {
  using grlensing::log;
  using grlensing::LogEvent;

  log<LogEvent::message>("Dumping metric {:s}", metric_name);

  const auto extension = config_file["dump_metric_settings"]["extension"].as<std::string>();
  auto rf = config_file["dump_metric_settings"]["radius"].as<double>();
  auto points = config_file["dump_metric_settings"]["points"].as<unsigned>();

  mpi_index_map_3D coords_map(points + 1, points + 1, points + 1);

  // Lapse dumps
  dump_kernel(kernel, metric_name, "lapse", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                writer->push_final_real(metric->lapse(0.0, x, y, z));
              });

  // Lower shift dumps
  dump_kernel(kernel, metric_name, "l_shift", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                writer->push_real(metric->l_shift(0.0, x, y, z)[0]);
                writer->push_real(metric->l_shift(0.0, x, y, z)[1]);
                writer->push_final_real(metric->l_shift(0.0, x, y, z)[2]);
              });

  // Upper shift dumps
  dump_kernel(kernel, metric_name, "u_shift", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                writer->push_real(metric->u_shift(0.0, x, y, z)[0]);
                writer->push_real(metric->u_shift(0.0, x, y, z)[1]);
                writer->push_final_real(metric->u_shift(0.0, x, y, z)[2]);
              });

  // Lower metric dumps
  dump_kernel(kernel, metric_name, "ll_smetric", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                writer->push_real(metric->ll_smetric(0.0, x, y, z)[0][0]);
                writer->push_real(metric->ll_smetric(0.0, x, y, z)[0][1]);
                writer->push_real(metric->ll_smetric(0.0, x, y, z)[0][2]);

                writer->push_real(metric->ll_smetric(0.0, x, y, z)[1][0]);
                writer->push_real(metric->ll_smetric(0.0, x, y, z)[1][1]);
                writer->push_real(metric->ll_smetric(0.0, x, y, z)[1][2]);

                writer->push_real(metric->ll_smetric(0.0, x, y, z)[2][0]);
                writer->push_real(metric->ll_smetric(0.0, x, y, z)[2][1]);
                writer->push_final_real(metric->ll_smetric(0.0, x, y, z)[2][2]);
              });

  // Upper metric dumps
  dump_kernel(kernel, metric_name, "uu_smetric", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                writer->push_real(metric->uu_smetric(0.0, x, y, z)[0][0]);
                writer->push_real(metric->uu_smetric(0.0, x, y, z)[0][1]);
                writer->push_real(metric->uu_smetric(0.0, x, y, z)[0][2]);

                writer->push_real(metric->uu_smetric(0.0, x, y, z)[1][0]);
                writer->push_real(metric->uu_smetric(0.0, x, y, z)[1][1]);
                writer->push_real(metric->uu_smetric(0.0, x, y, z)[1][2]);

                writer->push_real(metric->uu_smetric(0.0, x, y, z)[2][0]);
                writer->push_real(metric->uu_smetric(0.0, x, y, z)[2][1]);
                writer->push_final_real(metric->uu_smetric(0.0, x, y, z)[2][2]);
              });

  // Lower extrinsic dumps
  dump_kernel(kernel, metric_name, "ll_extrinsic", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                writer->push_real(metric->ll_extrinsic(0.0, x, y, z)[0][0]);
                writer->push_real(metric->ll_extrinsic(0.0, x, y, z)[0][1]);
                writer->push_real(metric->ll_extrinsic(0.0, x, y, z)[0][2]);

                writer->push_real(metric->ll_extrinsic(0.0, x, y, z)[1][0]);
                writer->push_real(metric->ll_extrinsic(0.0, x, y, z)[1][1]);
                writer->push_real(metric->ll_extrinsic(0.0, x, y, z)[1][2]);

                writer->push_real(metric->ll_extrinsic(0.0, x, y, z)[2][0]);
                writer->push_real(metric->ll_extrinsic(0.0, x, y, z)[2][1]);
                writer->push_final_real(metric->ll_extrinsic(0.0, x, y, z)[2][2]);
              });

  // Mixed extrinsic dumps
  dump_kernel(kernel, metric_name, "ul_extrinsic", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                writer->push_real(metric->ul_extrinsic(0.0, x, y, z)[0][0]);
                writer->push_real(metric->ul_extrinsic(0.0, x, y, z)[0][1]);
                writer->push_real(metric->ul_extrinsic(0.0, x, y, z)[0][2]);

                writer->push_real(metric->ul_extrinsic(0.0, x, y, z)[1][0]);
                writer->push_real(metric->ul_extrinsic(0.0, x, y, z)[1][1]);
                writer->push_real(metric->ul_extrinsic(0.0, x, y, z)[1][2]);

                writer->push_real(metric->ul_extrinsic(0.0, x, y, z)[2][0]);
                writer->push_real(metric->ul_extrinsic(0.0, x, y, z)[2][1]);
                writer->push_final_real(metric->ul_extrinsic(0.0, x, y, z)[2][2]);
              });

  // Chirstoffel dump
  dump_kernel(kernel, metric_name, "spatial_christoffel", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[0][0][0]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[0][0][1]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[0][0][2]);

                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[0][1][0]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[0][1][1]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[0][1][2]);

                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[0][2][0]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[0][2][1]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[0][2][2]);

                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[1][0][0]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[1][0][1]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[1][0][2]);

                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[1][1][0]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[1][1][1]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[1][1][2]);

                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[1][2][0]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[1][2][1]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[1][2][2]);

                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[2][0][0]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[2][0][1]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[2][0][2]);

                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[2][1][0]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[2][1][1]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[2][1][2]);

                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[2][2][0]);
                writer->push_real(metric->spatial_christoffel(0.0, x, y, z)[2][2][1]);
                writer->push_final_real(metric->spatial_christoffel(0.0, x, y, z)[2][2][2]);
              });

  // Lapse gradient
  dump_kernel(kernel, metric_name, "grad_lapse", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                writer->push_real(metric->grad_lapse(0.0, x, y, z)[0]);
                writer->push_real(metric->grad_lapse(0.0, x, y, z)[1]);
                writer->push_final_real(metric->grad_lapse(0.0, x, y, z)[2]);
              });

  // Shift gradient dumps
  dump_kernel(kernel, metric_name, "grad_ushift", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                writer->push_real(metric->grad_ushift(0.0, x, y, z)[0][0]);
                writer->push_real(metric->grad_ushift(0.0, x, y, z)[0][1]);
                writer->push_real(metric->grad_ushift(0.0, x, y, z)[0][2]);

                writer->push_real(metric->grad_ushift(0.0, x, y, z)[1][0]);
                writer->push_real(metric->grad_ushift(0.0, x, y, z)[1][1]);
                writer->push_real(metric->grad_ushift(0.0, x, y, z)[1][2]);

                writer->push_real(metric->grad_ushift(0.0, x, y, z)[2][0]);
                writer->push_real(metric->grad_ushift(0.0, x, y, z)[2][1]);
                writer->push_final_real(metric->grad_ushift(0.0, x, y, z)[2][2]);
              });
}