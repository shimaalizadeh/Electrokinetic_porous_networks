#include "config.hpp"

#include <stdexcept>
#include <cstdlib>

#define TOML_EXCEPTIONS 1
#include <toml.hpp>

namespace fs = std::filesystem;

namespace {

const fs::path kInputDir = "input_files";

template <typename T>
T required(const toml::table& tbl, const char* section, const char* key, const fs::path& file) {
    auto node = tbl[section][key];
    if (!node) {
        throw std::runtime_error(file.string() + ": missing required key [" + section + "]." + key);
    }
    auto value = node.template value<T>();
    if (!value) {
        throw std::runtime_error(file.string() + ": key [" + section + "]." + key + " has the wrong type");
    }
    return *value;
}

template <typename T>
T optional_or(const toml::table& tbl, const char* section, const char* key, T fallback, const fs::path& file) {
    auto node = tbl[section][key];
    if (!node) return fallback;
    auto value = node.template value<T>();
    if (!value) {
        throw std::runtime_error(file.string() + ": key [" + section + "]." + key + " has the wrong type");
    }
    return *value;
}

}  // namespace

std::string Config::network_name_for_setup() const {
    // Network::setup prepends "input_files/"; hand it the path relative to that dir.
    fs::path rel = network_file.lexically_normal().lexically_relative(kInputDir);
    if (rel.empty() || *rel.begin() == "..") {
        throw std::runtime_error("network_file must be inside input_files/ (got " + network_file.string() + ")");
    }
    return rel.string();
}

std::string Config::tables_prefix() const {
    std::string s = tables_dir.string();
    if (!s.empty() && s.back() != '/') s += '/';
    return s;
}

Config load_config_toml(const fs::path& file) {
    if (!fs::is_regular_file(file)) {
        throw std::runtime_error("config file not found: " + file.string());
    }
    toml::table tbl;
    try {
        tbl = toml::parse_file(file.string());
    } catch (const toml::parse_error& err) {
        throw std::runtime_error(file.string() + ": TOML parse error: " + std::string(err.description()));
    }

    Config cfg;
    cfg.source_file = file;
    cfg.num_pores = static_cast<int>(required<int64_t>(tbl, "network", "num_pores", file));
    cfg.num_reservoirs = static_cast<int>(required<int64_t>(tbl, "network", "num_reservoirs", file));
    cfg.num_horizontal_pores = static_cast<int>(required<int64_t>(tbl, "network", "num_horizontal_pores", file));
    cfg.dt = required<double>(tbl, "solver", "dt", file);
    cfg.max_time = required<double>(tbl, "solver", "max_time", file);
    cfg.write_every = static_cast<int>(required<int64_t>(tbl, "solver", "write_every", file));
    cfg.restart = optional_or<bool>(tbl, "solver", "restart", false, file);
    cfg.tables_dir = optional_or<std::string>(tbl, "io", "tables_dir", "tables/", file);
    cfg.network_file = required<std::string>(tbl, "io", "network_file", file);
    return cfg;
}

Config load_config_args(int argc, char* argv[]) {
    if (argc != 10) throw std::runtime_error("expected 9 positional arguments");
    Config cfg;
    cfg.num_pores = std::atoi(argv[1]);
    cfg.num_reservoirs = std::atoi(argv[2]);
    cfg.num_horizontal_pores = std::atoi(argv[3]);
    cfg.dt = std::atof(argv[4]);
    cfg.max_time = std::atof(argv[5]);
    cfg.write_every = std::atoi(argv[6]);
    cfg.restart = std::atoi(argv[7]) != 0;
    cfg.tables_dir = argv[8];
    // Legacy behavior: the network file name is relative to input_files/.
    cfg.network_file = kInputDir / argv[9];
    return cfg;
}

void validate_config(const Config& cfg) {
    if (cfg.num_pores <= 0) throw std::runtime_error("num_pores must be > 0");
    if (cfg.num_reservoirs <= 0) throw std::runtime_error("num_reservoirs must be > 0");
    if (cfg.num_horizontal_pores < 0) throw std::runtime_error("num_horizontal_pores must be >= 0");
    if (cfg.dt <= 0.) throw std::runtime_error("dt must be > 0");
    if (cfg.max_time <= 0.) throw std::runtime_error("max_time must be > 0");
    if (cfg.write_every <= 0) throw std::runtime_error("write_every must be > 0");
    if (!fs::is_directory(cfg.tables_dir)) {
        throw std::runtime_error("tables_dir not found: " + cfg.tables_dir.string());
    }
    for (const char* name : {"slit_sigma_lambda0.in", "slit_sigma_lambdastar.in",
                             "circle_sigma_lambda0.in", "circle_sigma_lambdastar.in"}) {
        if (!fs::is_regular_file(cfg.tables_dir / name)) {
            throw std::runtime_error("missing table file: " + (cfg.tables_dir / name).string());
        }
    }
    if (cfg.restart) {
        for (const char* name : {"network_final.in", "channels_final.bin"}) {
            if (!fs::is_regular_file(kInputDir / name)) {
                throw std::runtime_error("restart = true but restart file is missing: " + (kInputDir / name).string());
            }
        }
    } else if (!fs::is_regular_file(cfg.network_file)) {
        throw std::runtime_error("network_file not found: " + cfg.network_file.string());
    }
    (void)cfg.network_name_for_setup();  // checks it lives under input_files/
}
