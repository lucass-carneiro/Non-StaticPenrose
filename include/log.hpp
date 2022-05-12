/**
 * log.hpp - Handles messages and event informations to the user.
 * All log data is sent to stdout. Color can be controlled with
 * the CMake option WITH_COLORS.
 */
#ifndef GRLENSING_LOG_HPP
#define GRLENSING_LOG_HPP

#include "options.hpp"

#include <array>
#include <ctime>
#include <fmt/color.h>
#include <fmt/core.h>
#include <mpi.h>
#include <string>

namespace grlensing {

static constexpr const std::size_t max_datetimestring_size = 20;

/**
 * Returns a formatted string containing the current (local) date and time
 *
 * This function uses ctime instead of the chrono machinary, but it is
 * thread safe. It uses a local stack allocated buffer to store the
 * date string, which can be adjusted at compile time.
 *
 * @param n The size of the internal buffer used to store the date and time.
 * @return A C++ string containing the current local date and time.
 */
template <size_t n> auto get_datetime_string() -> std::string {
  std::time_t t = std::time(nullptr);
  std::array<char, n> mbstr{};
  std::memset(mbstr.data(), 0, n);
  auto bytes_written = std::strftime(mbstr.data(), n, "%F %T", std::localtime(&t));

  if (bytes_written != 0)
    return std::string(mbstr.data());
  else
    return std::string{"YYYY-MM-DD HH:MM:SS"};
}

/**
 * The type of event to be logged.
 *
 * The type of event determines the string that is output to the screen as
 * well as the color of the log entry. These types serve as tags that are
 * passed to the general log function. Based on type information, another
 * function, specific to the type of the event being logged is called and
 * executed. This dispatch happens at compile time.
 */
enum class LogEvent { warning, error, message, missing_feature };

/**
 * Print a log entry to the screen - non variadic version used
 * internally.
 *
 * This functions is intended to be used internally. It is called
 * from the variadic version log and prints a log entry to the screen using
 * user provided formats.
 * This function reads from the USE_LOG_COLORS macro, which can be set during
 * CMake configuration process by passing -D LOG_COLORS=ON to the CMake command
 * line. This variable is set to ON by default. This should be disabled if using
 * terminals that do not support ANSI colors.
 *
 * @param e The type of event to log.
 * @param str A string containing fmt format information.
 * @param args The arguments of the format string.
 */
template <LogEvent e> void vlog(fmt::string_view str, fmt::format_args args);

/**
 * Print a regular log message to the screen - non variadic version
 * used internally.
 */
template <> inline void vlog<LogEvent::message>(fmt::string_view str, fmt::format_args args) {
#ifdef GRLENSING_USE_LOG_COLOR
  fmt::print(fg(fmt::color::steel_blue), "{:s} (process {:d}) - GRLensing message: ",
             get_datetime_string<max_datetimestring_size>(), MPI::COMM_WORLD.Get_rank());
#else
  fmt::print("{:s} (process {:d}) - GRLensing message: ",
             get_datetime_string<max_datetimestring_size>(), MPI::COMM_WORLD.Get_rank());
#endif
  fmt::vprint(str, args);
  fmt::print("\n");
}

/**
 * Print a regular log warning to the screen - non variadic version
 * used internally.
 */
template <> inline void vlog<LogEvent::warning>(fmt::string_view str, fmt::format_args args) {
#ifdef GRLENSING_USE_LOG_COLOR
  fmt::print(fg(fmt::color::yellow), "{:s} (process {:d}) - GRLensing warning: ",
             get_datetime_string<max_datetimestring_size>(), MPI::COMM_WORLD.Get_rank());
#else
  fmt::print("{:s} (process {:d}) - GRLensing warning: ",
             get_datetime_string<max_datetimestring_size>(), MPI::COMM_WORLD.Get_rank());
#endif
  fmt::vprint(str, args);
  fmt::print("\n");
}

/**
 * Print a missing feature to the screen - non variadic version
 * used internally.
 */
template <>
inline void vlog<LogEvent::missing_feature>(fmt::string_view str, fmt::format_args args) {
#ifdef GRLENSING_USE_LOG_COLOR
  fmt::print(fmt::emphasis::bold | fg(fmt::color::yellow),
             "{:s} (process {:d}) - GRLensing missing feature: ",
             get_datetime_string<max_datetimestring_size>(), MPI::COMM_WORLD.Get_rank());
#else
  fmt::print("{:s} (process {:d}) - GRLensing missing feature: ",
             get_datetime_string<max_datetimestring_size>(), MPI::COMM_WORLD.Get_rank());
#endif
  fmt::vprint(str, args);
  fmt::print("\n");
}

/**
 * Print a regular log error to the screen - non variadic version
 * used internally.
 */
template <> inline void vlog<LogEvent::error>(fmt::string_view str, fmt::format_args args) {
#ifdef GRLENSING_USE_LOG_COLOR
  fmt::print(fmt::emphasis::bold | fg(fmt::color::red), "{:s} (process {:d}) - GRLensing error: ",
             get_datetime_string<max_datetimestring_size>(), MPI::COMM_WORLD.Get_rank());
#else
  fmt::print("{:s} (process {:d}) - GRLensing error: ",
             get_datetime_string<max_datetimestring_size>(), MPI::COMM_WORLD.Get_rank());
#endif
  fmt::vprint(str, args);
  fmt::print("\n");
}

/**
 * Log an event.
 *
 * This function logs an event using fmt print syntax. @see vlog()
 *
 * @param e The type of event to log.
 * @param Args The types of the format arguments. This is automatically deduced by the compiler.
 * and should not be specified.
 * @param str A fmt-like format string.
 * @param args The arguments of the fmt-like format string.
 */
#ifdef GRLENSING_USE_LOG
template <LogEvent e, typename... Args> void log(fmt::string_view str, Args &&...args) {
  vlog<e>(str, fmt::make_format_args(args...));
}
#else
template <LogEvent e, typename... Args> void log(fmt::string_view, Args &&...) { return; }
#endif

} // namespace grlensing

#endif // GRLENSING_LOG_HPP
