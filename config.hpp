#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <filesystem>
#include <optional>
#include <string>

// Run settings for ek_solver, loaded from a TOML file or from positional args.
struct Config {
    int num_pores = 0;
    int num_reservoirs = 0;
    int num_horizontal_pores = 0;
    double dt = 0.;
    double max_time = 0.;
    int write_every = 0;
    bool restart = false;
    std::filesystem::path tables_dir;
    std::filesystem::path network_file;
    std::optional<std::filesystem::path> source_file;  // TOML file it came from, if any

    // Network file name as Network::setup expects it (relative to input_files/).
    std::string network_name_for_setup() const;
    // tables_dir as a string with a trailing '/', as the table loader concatenates names.
    std::string tables_prefix() const;
};

// Throw std::runtime_error with a clear message on any problem.
Config load_config_toml(const std::filesystem::path& file);
Config load_config_args(int argc, char* argv[]);
void validate_config(const Config& cfg);

#endif
