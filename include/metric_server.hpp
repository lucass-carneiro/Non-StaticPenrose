#ifndef GRLENSING_METRIC_SERVER_HPP
#define GRLENSING_METRIC_SERVER_HPP

#include "api_macros.hpp"

#include <array>
#include <list>
#include <memory>
#include <string>
#include <string_view>
#include <yaml-cpp/yaml.h>

namespace grlensing {

class metric_server {
public:
  /**
   * Alias for a 3-vector.
   */
  using spatial_vector = std::array<double, 3>;

  /**
   * Alias for a 3-matrix
   */
  using spatial_matrix = std::array<std::array<double, 3>, 3>;

  /**
   * Alias used for storing Chirstoffel symbol components
   */
  using chirstofell_t = std::array<std::array<std::array<double, 3>, 3>, 3>;

  /**
   * Define a spacetime metric using it's ADM components and
   * their derivatives.
   */
  class adm_metric { // TODO: add remaining constructors
  public:
    virtual void load_parameters(const YAML::Node &) = 0;

    /**
     * The ADM lapse function ($\alpha$)
     *
     * @param t The time coordinate value.
     * @param x The first spatial coordinate value.
     * @param y The first spatial coordinate value.
     * @param z The first spatial coordinate value.
     * @return The value of the lapse at the given spacetime point.
     */
    virtual auto lapse(double, double, double, double) -> double = 0;

    /**
     * The ADM lower shift vector ($\beta_i$)
     *
     * @param t The time coordinate value.
     * @param x The first spatial coordinate value.
     * @param y The first spatial coordinate value.
     * @param z The first spatial coordinate value.
     * @return A vector containing the shift values.
     */
    virtual auto l_shift(double, double, double, double) -> spatial_vector = 0;

    /**
     * The ADM upper shift vector ($\beta^i$)
     *
     * @param t The time coordinate value.
     * @param x The first spatial coordinate value.
     * @param y The first spatial coordinate value.
     * @param z The first spatial coordinate value.
     * @return A vector containing the shift values.
     */
    virtual auto u_shift(double, double, double, double) -> spatial_vector = 0;

    /**
     * The ADM lower spatial metric ($\gamma_{ij}$)
     *
     * @param t The time coordinate value.
     * @param x The first spatial coordinate value.
     * @param y The first spatial coordinate value.
     * @param z The first spatial coordinate value.
     * @return A matrix containing the metric values.
     */
    virtual auto ll_smetric(double, double, double, double) -> spatial_matrix = 0;

    /**
     * The ADM upper spatial metric ($\gamma^{ij}$)
     *
     * @param t The time coordinate value.
     * @param x The first spatial coordinate value.
     * @param y The first spatial coordinate value.
     * @param z The first spatial coordinate value.
     * @return A matrix containing the metric values.
     */
    virtual auto uu_smetric(double, double, double, double) -> spatial_matrix = 0;

    /**
     * The ADM lower extrinsic curvature ($K_{ij}$)
     *
     * @param t The time coordinate value.
     * @param x The first spatial coordinate value.
     * @param y The first spatial coordinate value.
     * @param z The first spatial coordinate value.
     * @return A matrix containing the extrinsic curvature values.
     */
    virtual auto ll_extrinsic(double, double, double, double) -> spatial_matrix = 0;

    /**
     * The ADM upper/lower extrinsic curvature ($K^i_j$)
     *
     * @param t The time coordinate value.
     * @param x The first spatial coordinate value.
     * @param y The first spatial coordinate value.
     * @param z The first spatial coordinate value.
     * @return A matrix containing the extrinsic curvature values.
     */
    virtual auto ul_extrinsic(double, double, double, double) -> spatial_matrix = 0;

    /**
     * The ADM 3-Chrostoffel symbol (the symble with respect to the induced metric)
     *
     * @param t The time coordinate value.
     * @param x The first spatial coordinate value.
     * @param y The first spatial coordinate value.
     * @param z The first spatial coordinate value.
     * @return A 3D array with the values of the symbos.
     */
    virtual auto spatial_christoffel(double, double, double, double) -> chirstofell_t = 0;

    /**
     * The gradient of the lapse
     *
     * @param t The time coordinate value.
     * @param x The first spatial coordinate value.
     * @param y The first spatial coordinate value.
     * @param z The first spatial coordinate value.
     * @return A 3 dimentional vector where each element represents a derivative.
     */
    virtual auto grad_lapse(double, double, double, double) -> spatial_vector = 0;

    /**
     * The gradient of the shift (upper index) vector.
     *
     * @param t The time coordinate value.
     * @param x The first spatial coordinate value.
     * @param y The first spatial coordinate value.
     * @param z The first spatial coordinate value.
     * @return A 3 dimentional matrix where each row represents a derivative and column a $\beta^i$
     * component.
     */
    virtual auto grad_ushift(double, double, double, double) -> spatial_matrix = 0;

    /**
     * The display name of the archive reader
     *
     * @return A string with the name of the metric.
     */
    virtual auto name() -> std::string_view = 0;

    virtual ~adm_metric() = default;
  };

  /**
   * Alias for a unique pointer to a adm_metric object
   */
  using metric_ptr = std::unique_ptr<adm_metric>;

  /**
   * Allows plugins to add new metrics
   */
  GRLENSING_API void add_metric(metric_ptr metric) { metrics.push_back(std::move(metric)); }

  /**
   * Print the name of the stored metrics, as provided by the plugin
   * implementer.
   */
  GRLENSING_API void print_names() const;

  /**
   * Gets a metric from it's name
   *
   * @param name The name of the metric.
   * @return A constant reference to a unique_ptr wraping the metric.
   */
  GRLENSING_API auto get_metric(std::string_view name) const noexcept(false) -> const metric_ptr &;

private:
  /**
   * A list of metrics
   */
  using metric_list = std::list<metric_ptr>;

  /**
   * All available metrics
   */
  metric_list metrics;
};

/**
 * The module of the shift four vector.
 *
 * This quantity is simply the contraction $\beta_i \beta^i$.
 * computed directly from metric components
 *
 * @param metric The spacetime metric to use.
 * @param t The time coordinate value.
 * @param x The first spatial coordinate value.
 * @param y The first spatial coordinate value.
 * @param z The first spatial coordinate value.
 * @return The reconstructed $\beta_i \beta^i$ contraction.
 * component.
 */
auto shift_module(const metric_server::metric_ptr &, double, double, double, double) -> double;

/**
 * The (covariant) diagonal time component of the spacetime metric.
 *
 * This quantity is computed from a ADM metric by the relation
 * $g_{00} = \gamma_{ij} \beta^i \beta^j - \alpha^2$
 *
 * @param metric The spacetime metric to use.
 * @param t The time coordinate value.
 * @param x The first spatial coordinate value.
 * @param y The first spatial coordinate value.
 * @param z The first spatial coordinate value.
 * @return The reconstructed $g_{00}$ metric component.
 * component.
 */
auto ll_g_00(const metric_server::metric_ptr &metric, double t, double x, double y, double z)
    -> double;

/**
 * The (covariant) full pacetime metric reconstructed from the ADM components
 *
 * @param metric The spacetime ADM metric to use.
 * @param t The time coordinate value.
 * @param x The first spatial coordinate value.
 * @param y The first spatial coordinate value.
 * @param z The first spatial coordinate value.
 * @return The reconstructed 4D spacetime metric component.
 * component.
 */
auto ll_g(const metric_server::metric_ptr &metric, double t, double x, double y, double z)
    -> std::array<std::array<double, 4>, 4>;

/**
 * Computes the global energy from local quantities.
 *
 * The global energy is the energy as measured by an observer that uses the propper time to
 * parametrize it's trajectory, or equivalently the energy as measured by a static observer at
 * infinity
 *
 * @param metric The background metric where the trajectory is being integrated.
 * @param ti The current coordinate time being integrated.
 * @param Vi The current local velocities of the particle.
 * @param Xi The current local positions of the particle.
 * @param En The current local energy of the particle.
 * @return The global energy.
 */
auto compute_global_energy(const metric_server::metric_ptr &metric, double ti,
                           const metric_server::spatial_vector &Vi,
                           const metric_server::spatial_vector &Xi, double En) -> double;

/**
 * Reconstructs the upper index four momentum of the particle from local quantities.
 *
 * The four momentum is given by $p^\mu = E_l * (1/N, V^i - \beta^i/N)$
 *
 * @param metric The background metric where the trajectory is being integrated.
 * @param ti The current coordinate time being integrated.
 * @param Vi The current local velocities of the particle.
 * @param Xi The current local positions of the particle.
 * @param En The current local energy of the particle.
 * @return The global energy.
 */
auto reconstruct_u_p(const metric_server::metric_ptr &metric, double ti,
                     const metric_server::spatial_vector &Vi,
                     const metric_server::spatial_vector &Xi, double En) -> std::array<double, 4>;

auto decompose_u_p(const metric_server::metric_ptr &metric, const std::array<double, 4> &u_p,
                   double ti, const metric_server::spatial_vector &Xi)
    -> std::tuple<metric_server::spatial_vector, double>;

/**
 * TODO: Doc
 */
auto compute_mass(const metric_server::metric_ptr &metric, double ti,
                  const metric_server::spatial_vector &Vi, const metric_server::spatial_vector &Xi,
                  double En) -> double;

} // namespace grlensing

#endif // GRLENSING_METRIC_SERVER_HPP