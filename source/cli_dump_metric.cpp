#include "cli.hpp"
#include "log.hpp"
#include "mpi_index_map_3D.hpp"

#include <cmath>
#include <concepts>
#include <cstddef>
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

  const auto r0 = -rf;
  const auto dr = (rf - r0) / (points - 1);

  const auto [I_0, I_f] = coords_map.get_glor();
  for (auto I = I_0; I <= I_f; I++) {
    const auto [i, j, k] = coords_map.global_linear_to_global_matrix(I);

    const auto x = r0 + i * dr;
    const auto y = r0 + j * dr;
    const auto z = r0 + k * dr;

    writer->push_real(x);
    writer->push_real(y);
    writer->push_real(z);

    operation(x, y, z, metric, writer);
  }

  writer->push_metadata("glor start", std::size_t(I_0));
  writer->push_metadata("glor end", std::size_t(I_f));

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

  mpi_index_map_3D coords_map(points, points, points);

  // Lapse dumps
  dump_kernel(kernel, metric_name, "lapse", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                const auto value = metric->lapse(10.0, x, y, z);
                writer->push_final_real(std::isnan(value) ? 0.0 : value);
              });

  // Lower shift dumps
  dump_kernel(kernel, metric_name, "l_shift", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                auto value = metric->l_shift(0.0, x, y, z)[0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->l_shift(0.0, x, y, z)[1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->l_shift(0.0, x, y, z)[2];
                writer->push_final_real(std::isnan(value) ? 0.0 : value);
              });

  // Upper shift dumps
  dump_kernel(kernel, metric_name, "u_shift", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                auto value = metric->u_shift(0.0, x, y, z)[0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->u_shift(0.0, x, y, z)[1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->u_shift(0.0, x, y, z)[2];
                writer->push_final_real(std::isnan(value) ? 0.0 : value);
              });

  // Lower metric dumps
  dump_kernel(kernel, metric_name, "ll_smetric", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                auto value = metric->ll_smetric(0.0, x, y, z)[0][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_smetric(0.0, x, y, z)[0][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_smetric(0.0, x, y, z)[0][2];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_smetric(0.0, x, y, z)[1][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_smetric(0.0, x, y, z)[1][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_smetric(0.0, x, y, z)[1][2];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_smetric(0.0, x, y, z)[2][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_smetric(0.0, x, y, z)[2][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_smetric(0.0, x, y, z)[2][2];
                writer->push_final_real(std::isnan(value) ? 0.0 : value);
              });

  // Upper metric dumps
  dump_kernel(kernel, metric_name, "uu_smetric", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                auto value = metric->uu_smetric(0.0, x, y, z)[0][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->uu_smetric(0.0, x, y, z)[0][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->uu_smetric(0.0, x, y, z)[0][2];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->uu_smetric(0.0, x, y, z)[1][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->uu_smetric(0.0, x, y, z)[1][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->uu_smetric(0.0, x, y, z)[1][2];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->uu_smetric(0.0, x, y, z)[2][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->uu_smetric(0.0, x, y, z)[2][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->uu_smetric(0.0, x, y, z)[2][2];
                writer->push_final_real(std::isnan(value) ? 0.0 : value);
              });

  // Lower extrinsic dumps
  dump_kernel(kernel, metric_name, "ll_extrinsic", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                auto value = metric->ll_extrinsic(0.0, x, y, z)[0][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_extrinsic(0.0, x, y, z)[0][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_extrinsic(0.0, x, y, z)[0][2];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_extrinsic(0.0, x, y, z)[1][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_extrinsic(0.0, x, y, z)[1][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_extrinsic(0.0, x, y, z)[1][2];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_extrinsic(0.0, x, y, z)[2][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_extrinsic(0.0, x, y, z)[2][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ll_extrinsic(0.0, x, y, z)[2][2];
                writer->push_final_real(std::isnan(value) ? 0.0 : value);
              });

  // Mixed extrinsic dumps
  dump_kernel(kernel, metric_name, "ul_extrinsic", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                auto value = metric->ul_extrinsic(0.0, x, y, z)[0][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ul_extrinsic(0.0, x, y, z)[0][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ul_extrinsic(0.0, x, y, z)[0][2];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ul_extrinsic(0.0, x, y, z)[1][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ul_extrinsic(0.0, x, y, z)[1][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ul_extrinsic(0.0, x, y, z)[1][2];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ul_extrinsic(0.0, x, y, z)[2][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ul_extrinsic(0.0, x, y, z)[2][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->ul_extrinsic(0.0, x, y, z)[2][2];
                writer->push_final_real(std::isnan(value) ? 0.0 : value);
              });

  // Chirstoffel dump
  dump_kernel(kernel, metric_name, "spatial_christoffel", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                for (std::size_t i = 0; i < 3; i++) {
                  for (std::size_t j = 0; j < 3; j++) {
                    for (size_t k = 0; k < 3; k++) {
                      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
                      const auto value = metric->spatial_christoffel(0.0, x, y, z)[i][j][k];
                      if (i == 2 && j == 2 && k == 2) {
                        writer->push_final_real(std::isnan(value) ? 0.0 : value);
                      } else {
                        writer->push_real(std::isnan(value) ? 0.0 : value);
                      }
                    }
                  }
                }
              });

  // Lapse gradient
  dump_kernel(kernel, metric_name, "grad_lapse", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                auto value = metric->grad_lapse(0.0, x, y, z)[0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->grad_lapse(0.0, x, y, z)[1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->grad_lapse(0.0, x, y, z)[2];
                writer->push_final_real(std::isnan(value) ? 0.0 : value);
              });

  // Shift gradient dumps
  dump_kernel(kernel, metric_name, "grad_ushift", extension, rf, points, coords_map,
              [](double x, double y, double z, const metric_ptr &metric, // NOLINT
                 const writer_ptr &writer) -> void {
                auto value = metric->grad_ushift(0.0, x, y, z)[0][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->grad_ushift(0.0, x, y, z)[0][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->grad_ushift(0.0, x, y, z)[0][2];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->grad_ushift(0.0, x, y, z)[1][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->grad_ushift(0.0, x, y, z)[1][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->grad_ushift(0.0, x, y, z)[1][2];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->grad_ushift(0.0, x, y, z)[2][0];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->grad_ushift(0.0, x, y, z)[2][1];
                writer->push_real(std::isnan(value) ? 0.0 : value);

                value = metric->grad_ushift(0.0, x, y, z)[2][2];
                writer->push_final_real(std::isnan(value) ? 0.0 : value);
              });

  // Ricci scalar dump
  dump_kernel(
      kernel, metric_name, "ricci", extension, rf, points, coords_map,
      [](double x, double y, double z, const metric_ptr &metric, const writer_ptr &writer) -> void {
        /* To compute the Ricci scaler, we will use the Hamiltonian Constraint:
         * K_{ij} K^{ij} - K^2 = \gamma^{ik} \gamma^{jl} K_{ik} K_{kl} - (K^{i}_{i})^2
         */
        const auto uu_gamma = metric->uu_smetric(0.0, x, y, z);
        const auto ll_K = metric->ll_smetric(0.0, x, y, z);
        const auto ul_K = metric->ul_extrinsic(0.0, x, y, z);

        double part1 = 0.0;
        const auto part2 = ul_K(0, 0) + ul_K(1, 1) + ul_K(2, 2);

        for (std::size_t i = 0; i < 3; i++) {
          for (std::size_t j = 0; j < 3; j++) {
            for (std::size_t k = 0; k < 3; k++) {
              for (std::size_t l = 0; l < 3; l++) {
                part1 += uu_gamma(i, k) * uu_gamma(j, l) * ll_K(i, k) * ll_K(k, l);
              }
            }
          }
        }

        const auto value = part1 - part2 * part2;
        writer->push_final_real(std::isnan(value) ? 0.0 : value);
      });
}