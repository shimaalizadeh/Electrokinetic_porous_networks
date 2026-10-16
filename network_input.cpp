#include "network_input.hpp"

#include <fstream>
#include <initializer_list>
#include <set>
#include <sstream>
#include <stdexcept>

#define TOML_EXCEPTIONS 1
#include <toml.hpp>

namespace fs = std::filesystem;
using std::runtime_error;
using std::string;
using std::to_string;

namespace {

// Built-in defaults: precedence is entry field > [defaults] table > these values.
// Only neutral values get a built-in; physics-defining fields stay required.
constexpr bool   kDefaultDeadEnd = false;  // pore: not a dead end
constexpr double kDefaultArea    = 1.0;    // pore: unit (reference) cross-section
constexpr double kDefaultC0      = 1.0;    // pore + reservoir: reference concentration
constexpr double kDefaultVolume  = 0.0;    // reservoir: end reservoir (no volume)
constexpr double kDefaultPressure = 0.0;   // reservoir: ignored anyway (solver forces 0)
// type defaults to "slit" in get_type().

[[noreturn]] void fail(const fs::path& file, const string& msg) {
    throw runtime_error(file.string() + ": " + msg);
}

// Checks shared by both formats.
void validate(const NetworkInput& net, const fs::path& file) {
    const int np = static_cast<int>(net.pores.size());
    const int nr = static_cast<int>(net.reservoirs.size());
    std::vector<int> degree(nr, 0);
    for (int i = 0; i < np; ++i) {
        const auto& p = net.pores[i];
        const string tag = "pore " + to_string(i);
        if (p.type != 0 && p.type != 1) fail(file, tag + ": type must be slit(0) or circular(1)");
        if (!(p.length > 0.)) fail(file, tag + ": length must be > 0");
        if (p.n_cells <= 0) fail(file, tag + ": n_cells must be > 0");
        for (int e : p.ends)
            if (e < 0 || e >= nr) fail(file, tag + ": end reservoir " + to_string(e) + " out of range 0.." + to_string(nr - 1));
        if (p.ends[0] == p.ends[1]) fail(file, tag + ": both ends are reservoir " + to_string(p.ends[0]));
        ++degree[p.ends[0]];
        ++degree[p.ends[1]];
        if (p.profile && static_cast<int>(p.profile->sigma_star.size()) != p.n_cells)
            fail(file, tag + ": profile has " + to_string(p.profile->sigma_star.size()) +
                           " rows but n_cells = " + to_string(p.n_cells));
    }
    for (int r = 0; r < nr; ++r) {
        if (net.reservoirs[r].type != 0 && net.reservoirs[r].type != 1)
            fail(file, "reservoir " + to_string(r) + ": type must be slit(0) or circular(1)");
        if (degree[r] == 0) fail(file, "reservoir " + to_string(r) + " is not connected to any pore");
    }
    for (const auto& c : net.connections) {
        if (c.from < 0 || c.from >= nr || c.to < 0 || c.to >= nr || c.pore < 0 || c.pore >= np)
            fail(file, "connectivity entry (" + to_string(c.from) + " " + to_string(c.to) + " " +
                           to_string(c.pore) + ") out of range");
    }
}

std::vector<ConnectionInput> connections_from_ends(const std::vector<PoreInput>& pores) {
    std::vector<ConnectionInput> c;
    for (int i = 0; i < static_cast<int>(pores.size()); ++i) c.push_back({pores[i].ends[0], pores[i].ends[1], i});
    return c;
}

// ---- per-cell CSV profile: header sigma_star,lambda_D,Sp,C0 ----
CellProfile read_profile_csv(const fs::path& file) {
    std::ifstream in(file);
    if (!in) throw runtime_error("profile file not found: " + file.string());
    auto trim = [](string s) {
        const char* ws = " \t\r";
        s.erase(0, s.find_first_not_of(ws));
        s.erase(s.find_last_not_of(ws) + 1);
        return s;
    };
    string line;
    int lineno = 0;
    std::vector<string> header;
    CellProfile prof;
    while (std::getline(in, line)) {
        ++lineno;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        std::vector<string> fields;
        std::stringstream ss(line);
        string f;
        while (std::getline(ss, f, ',')) fields.push_back(trim(f));
        if (header.empty()) {
            header = fields;
            if (header != std::vector<string>{"sigma_star", "lambda_D", "Sp", "C0"})
                throw runtime_error(file.string() + ": header must be exactly sigma_star,lambda_D,Sp,C0");
            continue;
        }
        if (fields.size() != 4)
            throw runtime_error(file.string() + ":" + to_string(lineno) + ": expected 4 values, got " + to_string(fields.size()));
        double v[4];
        for (int k = 0; k < 4; ++k) {
            size_t used = 0;
            try { v[k] = std::stod(fields[k], &used); } catch (...) { used = 0; }
            if (used == 0 || used != fields[k].size())
                throw runtime_error(file.string() + ":" + to_string(lineno) + ": bad number '" + fields[k] + "'");
        }
        prof.sigma_star.push_back(v[0]);
        prof.lambda_D.push_back(v[1]);
        prof.Sp.push_back(v[2]);
        prof.C0.push_back(v[3]);
    }
    if (prof.sigma_star.empty()) throw runtime_error(file.string() + ": profile has no data rows");
    return prof;
}

// ---- TOML helpers ----
const toml::node* lookup(const toml::table& item, const toml::table* defaults, const char* key) {
    if (auto n = item.get(key)) return n;
    if (defaults) return defaults->get(key);
    return nullptr;
}

double get_double(const toml::table& item, const toml::table* defaults, const char* key,
                  const fs::path& file, const string& tag, std::optional<double> builtin = std::nullopt) {
    const toml::node* n = lookup(item, defaults, key);
    if (!n && builtin) return *builtin;
    if (!n) fail(file, tag + ": missing '" + key + "'" + (defaults ? " (and no value in [defaults])" : ""));
    if (auto v = n->value_exact<double>()) return *v;
    if (auto v = n->value_exact<int64_t>()) return static_cast<double>(*v);
    fail(file, tag + ": '" + key + "' must be a number");
}

int64_t get_int(const toml::table& item, const toml::table* defaults, const char* key,
                const fs::path& file, const string& tag, std::optional<int64_t> fallback = std::nullopt) {
    const toml::node* n = lookup(item, defaults, key);
    if (!n) {
        if (fallback) return *fallback;
        fail(file, tag + ": missing '" + key + "'");
    }
    if (auto v = n->value_exact<int64_t>()) return *v;
    fail(file, tag + ": '" + key + "' must be an integer");
}

bool get_bool(const toml::table& item, const toml::table* defaults, const char* key,
              const fs::path& file, const string& tag, std::optional<bool> fallback = std::nullopt) {
    const toml::node* n = lookup(item, defaults, key);
    if (!n) {
        if (fallback) return *fallback;
        fail(file, tag + ": missing '" + key + "'");
    }
    if (auto v = n->value_exact<bool>()) return *v;
    fail(file, tag + ": '" + key + "' must be true or false");
}

int get_type(const toml::table& item, const toml::table* defaults, const fs::path& file, const string& tag) {
    const toml::node* n = lookup(item, defaults, "type");
    if (!n) return 0;
    if (auto s = n->value_exact<string>()) {
        if (*s == "slit") return 0;
        if (*s == "circular") return 1;
        fail(file, tag + ": type must be \"slit\" or \"circular\" (got \"" + *s + "\")");
    }
    fail(file, tag + ": type must be \"slit\" or \"circular\"");
}

void check_keys(const toml::table& t, std::initializer_list<const char*> allowed,
                const fs::path& file, const string& tag) {
    for (const auto& [k, v] : t) {
        (void)v;
        bool ok = false;
        for (const char* a : allowed) ok = ok || k.str() == a;
        if (!ok) fail(file, tag + ": unknown key '" + string(k.str()) + "' (typo?)");
    }
}

const toml::array& get_array_of_tables(const toml::table& tbl, const char* key, const fs::path& file) {
    const toml::array* arr = tbl[key].as_array();
    if (!arr || !arr->is_array_of_tables()) fail(file, string("missing [[") + key + "]] entries");
    return *arr;
}

}  // namespace

NetworkInput read_network_legacy(const fs::path& file, int num_pores, int num_reservoirs) {
    std::ifstream in(file);
    if (!in) fail(file, "cannot open network file");
    NetworkInput net;
    net.pores.resize(num_pores);
    net.reservoirs.resize(num_reservoirs);
    for (int i = 0; i < num_pores; ++i) {
        auto& p = net.pores[i];
        in >> p.type >> p.dead_end >> p.ends[0] >> p.ends[1] >> p.length >> p.n_cells >> p.sigma_star >> p.lambda_D >> p.area >> p.c0;
        if (!in) fail(file, "could not read pore row " + to_string(i) + " (expected " + to_string(num_pores) + " pore rows of 10 values)");
    }
    for (int i = 0; i < num_reservoirs; ++i) {
        auto& r = net.reservoirs[i];
        in >> r.type >> r.pressure_flag >> r.potential_flag >> r.pressure >> r.potential >> r.c0 >> r.sigma_star >> r.lambda_D >> r.volume;
        if (!in) fail(file, "could not read reservoir row " + to_string(i) + " (expected " + to_string(num_reservoirs) + " reservoir rows of 9 values)");
    }
    for (int i = 0; i < num_pores; ++i) {
        ConnectionInput c;
        in >> c.from >> c.to >> c.pore;
        if (!in) fail(file, "could not read connectivity row " + to_string(i) + " (expected " + to_string(num_pores) + " rows of 'i j pore')");
        net.connections.push_back(c);
    }
    validate(net, file);
    return net;
}

NetworkInput read_network_toml(const fs::path& file, int num_pores, int num_reservoirs) {
    if (!fs::is_regular_file(file)) fail(file, "cannot open network file");
    toml::table tbl;
    try {
        tbl = toml::parse_file(file.string());
    } catch (const toml::parse_error& err) {
        fail(file, "TOML parse error at line " + to_string(err.source().begin.line) + ": " + string(err.description()));
    }
    check_keys(tbl, {"defaults", "pores", "reservoirs"}, file, "top level");
    const toml::table* defaults = tbl["defaults"].as_table();
    if (defaults) check_keys(*defaults, {"type", "dead_end", "length", "n_cells", "sigma_star", "lambda_D", "area", "c0"}, file, "[defaults]");
    const auto& pores = get_array_of_tables(tbl, "pores", file);
    const auto& reservoirs = get_array_of_tables(tbl, "reservoirs", file);

    if (static_cast<int>(pores.size()) != num_pores)
        fail(file, "has " + to_string(pores.size()) + " pores but the run config says num_pores = " + to_string(num_pores));
    if (static_cast<int>(reservoirs.size()) != num_reservoirs)
        fail(file, "has " + to_string(reservoirs.size()) + " reservoirs but the run config says num_reservoirs = " + to_string(num_reservoirs));

    NetworkInput net;
    net.pores.resize(num_pores);
    net.reservoirs.resize(num_reservoirs);

    std::set<int64_t> seen;
    for (size_t k = 0; k < pores.size(); ++k) {
        const toml::table& t = *pores[k].as_table();
        string tag = "[[pores]] entry " + to_string(k);
        int64_t id = get_int(t, nullptr, "id", file, tag);
        if (id < 0 || id >= num_pores) fail(file, tag + ": id " + to_string(id) + " out of range 0.." + to_string(num_pores - 1));
        if (!seen.insert(id).second) fail(file, tag + ": duplicate pore id " + to_string(id));
        tag = "pore " + to_string(id);
        check_keys(t, {"id", "ends", "profile", "type", "dead_end", "length", "n_cells", "sigma_star", "lambda_D", "area", "c0"}, file, tag);
        PoreInput& p = net.pores[id];
        p.type = get_type(t, defaults, file, tag);
        p.dead_end = get_bool(t, defaults, "dead_end", file, tag, kDefaultDeadEnd) ? 1 : 0;
        const toml::array* ends = t["ends"].as_array();
        if (!ends || ends->size() != 2 || !(*ends)[0].is_integer() || !(*ends)[1].is_integer())
            fail(file, tag + ": 'ends' must be [inlet_reservoir, outlet_reservoir]");
        p.ends[0] = static_cast<int>(*(*ends)[0].value<int64_t>());
        p.ends[1] = static_cast<int>(*(*ends)[1].value<int64_t>());
        p.length = get_double(t, defaults, "length", file, tag);
        p.n_cells = static_cast<int>(get_int(t, defaults, "n_cells", file, tag));
        p.sigma_star = get_double(t, defaults, "sigma_star", file, tag);
        p.lambda_D = get_double(t, defaults, "lambda_D", file, tag);
        p.area = get_double(t, defaults, "area", file, tag, kDefaultArea);
        p.c0 = get_double(t, defaults, "c0", file, tag, kDefaultC0);
        if (auto prof = t["profile"].value<string>()) {
            fs::path pf = fs::path(*prof).is_absolute() ? fs::path(*prof) : file.parent_path() / *prof;
            p.profile = read_profile_csv(pf);
        }
    }

    seen.clear();
    for (size_t k = 0; k < reservoirs.size(); ++k) {
        const toml::table& t = *reservoirs[k].as_table();
        string tag = "[[reservoirs]] entry " + to_string(k);
        int64_t id = get_int(t, nullptr, "id", file, tag);
        if (id < 0 || id >= num_reservoirs) fail(file, tag + ": id " + to_string(id) + " out of range 0.." + to_string(num_reservoirs - 1));
        if (!seen.insert(id).second) fail(file, tag + ": duplicate reservoir id " + to_string(id));
        tag = "reservoir " + to_string(id);
        check_keys(t, {"id", "type", "pressure_known", "potential_known", "pressure", "potential",
                       "c0", "sigma_star", "lambda_D", "volume"}, file, tag);
        ReservoirInput& r = net.reservoirs[id];
        r.type = get_type(t, nullptr, file, tag);
        // TOML uses the physical meaning; the solver's raw flag is the opposite (0 = known).
        r.pressure_flag = !get_bool(t, nullptr, "pressure_known", file, tag);
        r.potential_flag = !get_bool(t, nullptr, "potential_known", file, tag);
        r.pressure = get_double(t, nullptr, "pressure", file, tag, kDefaultPressure);
        r.potential = get_double(t, nullptr, "potential", file, tag);
        r.c0 = get_double(t, nullptr, "c0", file, tag, kDefaultC0);
        r.sigma_star = get_double(t, defaults, "sigma_star", file, tag);
        r.lambda_D = get_double(t, defaults, "lambda_D", file, tag);
        r.volume = get_double(t, nullptr, "volume", file, tag, kDefaultVolume);
    }

    net.connections = connections_from_ends(net.pores);
    validate(net, file);
    return net;
}

NetworkInput read_network_file(const fs::path& file, int num_pores, int num_reservoirs) {
    if (file.extension() == ".toml") return read_network_toml(file, num_pores, num_reservoirs);
    return read_network_legacy(file, num_pores, num_reservoirs);
}
