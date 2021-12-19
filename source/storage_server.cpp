#include "storage_server.hpp"

GRLENSING_API auto grlensing::storage_server::open_archive(const std::filesystem::path &p)
    -> std::unique_ptr<grlensing::archive> {
  for (const auto &reader : archive_readers) {
    if (reader->can_open_archive(p)) {
      return reader->open_archive(p);
    }
  }
  log<LogEvent::error>("No plugin available to open file extension {:s}", p.extension().c_str());
  throw std::runtime_error("IO error");
}

GRLENSING_API void grlensing::storage_server::print_names() const {
  for (const auto &reader : archive_readers) {
    log<LogEvent::message>("Archive reader: {:s}", reader->name());
  }

  for (const auto &writer : writers) {
    log<LogEvent::message>("Trajectory file writer: {:s}", writer->name());
  }
}

GRLENSING_API auto grlensing::storage_server::get_writer(std::string_view extension) const
    noexcept(false) -> const std::unique_ptr<grlensing::storage_server::trajectory_writer> & {
  for (const auto &writer : writers) {
    if (writer->can_write_extension(extension))
      return writer;
  }
  log<LogEvent::error>("No plugin available to write file extension {:s}", extension);
  throw std::runtime_error("IO error");
}