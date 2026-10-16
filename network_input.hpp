#ifndef NETWORK_INPUT_HPP
#define NETWORK_INPUT_HPP

// Network description, read either from a legacy Network.in or from a network.toml.
// Both readers fill the same NetworkInput; Network::apply_input copies it into the
// solver arrays, so the physics never depends on which file format was used.

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

struct CellProfile {                  // optional per-cell values (old Channels.in capability)
    std::vector<double> sigma_star;
    std::vector<double> lambda_D;
    std::vector<double> Sp;           // nondimensional area / Sp
    std::vector<double> C0;
};

struct PoreInput {
    int type = 0;                     // 0 = slit, 1 = circular
    int dead_end = 0;                 // legacy column 2 (Channel_blocked)
    int ends[2] = {0, 0};             // [inlet reservoir, outlet reservoir]
    double length = 0.;               // total length; dx = length / n_cells
    int n_cells = 0;
    double sigma_star = 0.;
    double lambda_D = 0.;
    double area = 0.;                 // -> nondimensional_Sp
    double c0 = 0.;
    std::optional<CellProfile> profile;
};

struct ReservoirInput {
    int type = 0;                     // 0 = slit, 1 = circular
    bool pressure_flag = false;       // raw legacy flag: 0 = KNOWN (fixed), 1 = solved for
    bool potential_flag = false;      // raw legacy flag: 0 = KNOWN (fixed), 1 = solved for
    double pressure = 0.;             // read but overwritten to 0 by the solver (legacy behavior)
    double potential = 0.;
    double c0 = 0.;
    double sigma_star = 0.;
    double lambda_D = 0.;
    double volume = 0.;
};

struct ConnectionInput { int from = 0, to = 0, pore = 0; };

struct NetworkInput {
    std::vector<PoreInput> pores;
    std::vector<ReservoirInput> reservoirs;
    std::vector<ConnectionInput> connections;   // legacy: read from file; TOML: derived from pore ends
};

// Both throw std::runtime_error with a clear message on any problem.
NetworkInput read_network_legacy(const std::filesystem::path& file, int num_pores, int num_reservoirs);
NetworkInput read_network_toml(const std::filesystem::path& file, int num_pores, int num_reservoirs);
NetworkInput read_network_file(const std::filesystem::path& file, int num_pores, int num_reservoirs);

#endif
