#include "SteeringParams.h"

#include "nlohmann/json.hpp"

#include <fstream>
#include <regex>

namespace mkfit {

// Begin AUTO code, some members commented out.

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(mkfit::IterationLayerConfig,
  /* float */   m_select_min_dphi,
  /* float */   m_select_max_dphi,
  /* float */   m_select_min_dq,
  /* float */   m_select_max_dq,
  // /* function<void(const Track&,const float,const float,float&,float&)> */   m_dynamic_windows,
  /* float */   m_qf_treg,
  /* float */   m_phif_treg,
  /* float */   m_phif_lpt_brl,
  /* float */   m_phif_lpt_treg,
  /* float */   m_phif_lpt_ec
)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(mkfit::IterationParams,
  /* int */   nlayers_per_seed,
  /* int */   maxCandsPerSeed,
  /* int */   maxHolesPerCand,
  /* int */   maxConsecHoles,
  /* float */   chi2Cut,
  /* float */   chi2CutOverlap,
  /* float */   pTCutOverlap,
  /* float */   c_ptthr_hpt,
  /* float */   c_drmax_bh,
  /* float */   c_dzmax_bh,
  /* float */   c_drmax_eh,
  /* float */   c_dzmax_eh,
  /* float */   c_drmax_bl,
  /* float */   c_dzmax_bl,
  /* float */   c_drmax_el,
  /* float */   c_dzmax_el
)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(mkfit::IterationConfig,
  // /* int */   m_iteration_index,
  /* int */   m_track_algorithm,
  /* mkfit::IterationParams */   m_params,
  // /* int */   m_n_regions,
  // /* vector<int> */   m_region_order,
  // /* vector<mkfit::SteeringParams> */   m_steering_params,
  /* vector<mkfit::IterationLayerConfig> */   m_layer_configs
  // /* function<void(const TrackerInfo&,const TrackVec&,const EventOfHits&,IterationSeedPartition&)> */   m_partition_seeds
)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(mkfit::IterationsInfo,
  /* vector<mkfit::IterationConfig> */   m_iterations
)

// End AUTO code.


// ============================================================================
// ConfigJsonPatcher
// ============================================================================

ConfigJsonPatcher::ConfigJsonPatcher(bool verbose) : m_verbose(verbose) {}

ConfigJsonPatcher::~ConfigJsonPatcher() { release_json(); }

void ConfigJsonPatcher::release_json()
{
    if (m_owner) delete m_json;
    m_json  = nullptr;
    m_owner = false;
}

std::string ConfigJsonPatcher::get_abs_path() const
{
    std::string s;
    s.reserve(64);
    for (auto &p : m_path_stack) s += p;
    return s;
}

std::string ConfigJsonPatcher::exc_hdr(const char* func) const
{
    std::string s;
    s.reserve(128);
    s = "ConfigJsonPatcher";
    if (func) { s += "::"; s += func; }
    s += " '";
    s += get_abs_path();
    s += "' ";
    return s;
}

template<class T> void ConfigJsonPatcher::Load(const T& o)
{
    release_json();
    m_json  = new nlohmann::json;
    *m_json = o;
    m_owner = true;
    cd_top();
}
template void ConfigJsonPatcher::Load<IterationsInfo> (const IterationsInfo  &o);
template void ConfigJsonPatcher::Load<IterationConfig>(const IterationConfig &o);

template<class T> void ConfigJsonPatcher::Save(T& o)
{
    from_json(*m_json, o);
}
template void ConfigJsonPatcher::Save<IterationConfig>(IterationConfig &o);

// Must not bork the IterationConfig elements of IterationsInfo ... default
// deserializator apparently reinitializes the vectors with defaults c-tors.
template <> void ConfigJsonPatcher::Save<IterationsInfo> (IterationsInfo &o)
{
    auto &itc_arr = m_json->at("m_iterations");
    for (int i = 0; i < o.size(); ++i)
    {
        from_json(itc_arr[i], o[i]);
    }
}

void ConfigJsonPatcher::cd(const std::string& path)
{
    nlohmann::json::json_pointer jp(path);
    m_json_stack.push_back(m_current);
    m_path_stack.push_back(path);
    m_current = & m_current->at(jp);
}

void ConfigJsonPatcher::cd_up(const std::string& path)
{
    if (m_json_stack.empty()) throw std::runtime_error("JSON stack empty on cd_up");

    m_current = m_json_stack.back();
    m_json_stack.pop_back();
    m_path_stack.pop_back();
    if ( ! path.empty()) cd(path);
}

void ConfigJsonPatcher::cd_top(const std::string& path)
{
    m_current = m_json;
    m_json_stack.clear();
    m_path_stack.clear();
    if ( ! path.empty()) cd(path);
}

template<typename T>
void ConfigJsonPatcher::replace(const std::string& path, T val)
{
    nlohmann::json::json_pointer jp(path);
    m_current->at(jp) = val;
}
template void ConfigJsonPatcher::replace<int>   (const std::string& path, int    val);
template void ConfigJsonPatcher::replace<float> (const std::string& path, float  val);
template void ConfigJsonPatcher::replace<double>(const std::string& path, double val);

template<typename T>
void ConfigJsonPatcher::replace(int first, int last, const std::string& path, T val)
{
    nlohmann::json::json_pointer jp(path);
    for (int i = first; i <= last; ++i)
    {
        m_current->at(i).at(jp) = val;
    }
}
template void ConfigJsonPatcher::replace<int>   (int first, int last, const std::string& path, int    val);
template void ConfigJsonPatcher::replace<float> (int first, int last, const std::string& path, float  val);
template void ConfigJsonPatcher::replace<double>(int first, int last, const std::string& path, double val);

nlohmann::json& ConfigJsonPatcher::get(const std::string& path)
{
    nlohmann::json::json_pointer jp(path);
    return m_current->at(jp);
}

int ConfigJsonPatcher::replace(const nlohmann::json &j)
{
    if (j.is_null()) throw std::runtime_error(exc_hdr(__func__) + "null not expected");

    if (j.is_boolean() || j.is_number() || j.is_string())
    {
        throw std::runtime_error(exc_hdr(__func__) + "value not expected on this parsing level");
    }

    int n_replaced = 0;

    if (j.is_object())
    {
        std::regex index_range_re("^\\[(\\d+)..(\\d+)\\]$");

        for (auto& [key, value] : j.items())
        {
            std::smatch m;
            std::regex_search(key, m, index_range_re);

            if (m.size() == 3)
            {
                if ( ! m_current->is_array())  throw std::runtime_error(exc_hdr(__func__) + "array range encountered when current json is not an array");
                int first = std::stoi(m.str(1));
                int last  = std::stoi(m.str(2));
                for (int i = first; i <= last; ++i)
                {
                    std::string s("/");
                    s += std::to_string(i);
                    cd(s);
                    n_replaced += replace(value);
                    cd_up();
                }
            }
            else if (value.is_array() || value.is_object())
            {
                std::string s("/");
                s += key;
                cd(s);
                n_replaced += replace(value);
                cd_up();
            }
            else if (value.is_number() || value.is_boolean() || value.is_string())
            {
                std::string s("/");
                s += key;
                nlohmann::json::json_pointer jp(s);
                if (m_verbose)
                {
                    std::cout << "  " << get_abs_path() << s << ": " << m_current->at(jp) << " -> " << value << "\n";
                }
                m_current->at(jp) = value;
                ++n_replaced;
            }
            else
            {
                throw std::runtime_error(exc_hdr(__func__) + "unexpected value type");
            }
        }
    }
    else if (j.is_array())
    {
        for (auto& element : j)
        {
            if ( ! element.is_object()) throw std::runtime_error(exc_hdr(__func__) + "array elements expected to be objects");
            n_replaced += replace(element);
        }
    }
    else
    {
        throw std::runtime_error(exc_hdr(__func__) + "unexpected json type");
    }

    return n_replaced;
}

std::string ConfigJsonPatcher::dump(int indent)
{
    return m_json->dump(indent);
}


// ============================================================================
// ConfigJson_Patch_File steering function
// ============================================================================
/*
    See example JSON patcher input: "mkFit/config-parse/test.json"

    The file can contain several valid JSON dumps in sequence.

    '/' character can be used to descend more than one level at a time.

    A number can be used to specify an array index. This can be combined with
    the '/' syntax.

    "[first,last]" key (as string) can be used to denote a range of array
    elements. Such a key must not be combined with a '/' syntax.
*/

namespace
{
    // Open file, throw exception on failure.
    void open_ifstream(std::ifstream &ifs, const std::string &fname, const char *pfx=0)
    {
        ifs.open(fname);
        if ( ! ifs)
        {
            char m[2048];
            snprintf(m, 2048, "%s%sError opening %s: %m", pfx ? pfx : "", pfx ? " " : "", fname.c_str());
            throw std::runtime_error(m);
        }
    }

    // Skip white-space, return true if more characters are available, false if eof.
    bool skipws_ifstream(std::ifstream &ifs)
    {
        while (std::isspace(ifs.peek())) ifs.get();
        return ! ifs.eof();
    }
}

void ConfigJson_Patch_File(IterationsInfo &its_info, const std::string &fname)
{
    using nlohmann::json;

    std::ifstream  ifs;
    open_ifstream(ifs, fname, __func__);

    ConfigJsonPatcher cjp(Config::json_patch_verbose);
    cjp.Load(its_info);

    if (Config::json_patch_dump_before)
    {
        std::cout << cjp.dump(3) << "\n";
    }

    printf("%s begin reading from file %s.\n", __func__, fname.c_str());

    int n_read = 0, n_tot_replaced = 0;
    while (skipws_ifstream(ifs))
    {
        json j;
        ifs >> j;
        ++n_read;

        std::cout << " Read JSON entity " << n_read << " -- applying patch:\n";
        // std::cout << j.dump(3) << "\n";

        int n_replaced = cjp.replace(j);
        std::cout << " Replaced " << n_replaced << " entries.\n";

        cjp.cd_top();
        n_tot_replaced += n_replaced;
    }
    printf("%s read %d JSON entities from file %s, replaced %d parameters.\n",
           __func__, n_read, fname.c_str(), n_tot_replaced);

     if (Config::json_patch_dump_after)
    {
        std::cout << cjp.dump(3) << "\n";
    }

    if (n_read > 0)
    {
        cjp.Save(its_info);
    }

    ifs.close();
}


// ============================================================================
// Tests for ConfigJson stuff
// ============================================================================

void ConfigJson_Test_Direct(IterationConfig &it_cfg)
{
    using nlohmann::json;

    std::string lojz("/m_select_max_dphi");

    json j = it_cfg;
    std::cout << j << "\n";
    std::cout << j.dump(3) << "\n";

    std::cout << "Layer 43, m_select_max_dphi = " << j["/m_layer_configs/43/m_select_max_dphi"_json_pointer] << "\n";
    std::cout << "Patching it to pi ...\n";
    json p = R"([
        { "op": "replace", "path": "/m_layer_configs/43/m_select_max_dphi", "value": 3.141 }
    ])"_json;
    j = j.patch(p);
    std::cout << "Layer 43, m_select_max_dphi = " << j["/m_layer_configs/43/m_select_max_dphi"_json_pointer] << "\n";

    auto &jx = j["/m_layer_configs/60"_json_pointer];
    // jx["m_select_max_dphi"] = 99.876;
    json::json_pointer jp(lojz);
    jx[jp] = 99.876;

    // try loading it back, see what happens to vector m_layer_configs.

    from_json(j, it_cfg);
    printf("Layer 43 : m_select_max_dphi = %f, size_of_layer_vec=%d, m_n_regions=%d, size_of_steering_params=%d\n",
           it_cfg.m_layer_configs[43].m_select_max_dphi, (int) it_cfg.m_layer_configs.size(),
           it_cfg.m_n_regions, (int) it_cfg.m_steering_params.size());

    printf("Layer 60 : m_select_max_dphi = %f, size_of_layer_vec=%d, m_n_regions=%d, size_of_steering_params=%d\n",
           it_cfg.m_layer_configs[60].m_select_max_dphi, (int) it_cfg.m_layer_configs.size(),
           it_cfg.m_n_regions, (int) it_cfg.m_steering_params.size());

    // try accessing something that does not exist

    // std::cout << "Non-existent path " << j["/m_layer_configs/143/m_select_max_dphi"_json_pointer] << "\n";

    auto &x = j["/m_layer_configs"_json_pointer];
    std::cout << "Typename /m_layer_configs " << x.type_name() << "\n";
    auto &y = j["/m_layer_configs/143"_json_pointer];
    std::cout << "Typename /m_layer_configs/143 " << y.type_name() << ", is_null=" << y.is_null() << "\n";
}

void ConfigJson_Test_Patcher(IterationConfig &it_cfg)
{
    ConfigJsonPatcher cjp;
    cjp.Load(it_cfg);

    std::cout << cjp.dump(3) << "\n";

    {
        cjp.cd("/m_layer_configs/43/m_select_max_dphi");
        std::cout << "Layer 43, m_select_max_dphi = " << cjp.get("") << "\n";
        std::cout << "Setting it to pi ...\n";
        cjp.replace("", 3.141);
        cjp.cd_top();
        std::cout << "Layer 43, m_select_max_dphi = " << cjp.get("/m_layer_configs/43/m_select_max_dphi") << "\n";
    }
    {
        std::cout << "Replacing layer 60 m_select_max_dphi with full path\n";
        cjp.replace("/m_layer_configs/60/m_select_max_dphi", 99.876);
    }
    try
    {
        std::cout << "Trying to replace an non-existent array entry\n";
        cjp.replace("/m_layer_configs/1460/m_select_max_dphi", 666.666);
    }
    catch (std::exception &exc)
    {
        std::cout << "Caugth exception: " << exc.what() << "\n";
    }
    try
    {
        std::cout << "Trying to replace an non-existent object entry\n";
        cjp.replace("/m_layer_configs/1/moo_select_max_dphi", 666.666);
    }
    catch (std::exception &exc)
    {
        std::cout << "Caugth exception: " << exc.what() << "\n";
    }
    {
        std::cout << "Replacing m_select_max_dphi on layers 1 to 3 to 7.7\n";
        cjp.cd("/m_layer_configs");
        cjp.replace(1, 3, "/m_select_max_dphi", 7.7);
        cjp.cd_top();
    }

    // try getting it back into c++, see what happens to vector m_layer_configs.

    cjp.Save(it_cfg);

    printf("Layer 43: m_select_max_dphi = %f, size_of_layer_vec=%d, m_n_regions=%d, size_of_steering_params=%d\n",
           it_cfg.m_layer_configs[43].m_select_max_dphi, (int) it_cfg.m_layer_configs.size(),
           it_cfg.m_n_regions, (int) it_cfg.m_steering_params.size());

    printf("Layer 60: m_select_max_dphi = %f\n", it_cfg.m_layer_configs[60].m_select_max_dphi);
    for (int i = 0; i < 5; ++i)
        printf("Layer %2d: m_select_max_dphi = %f\n", i, it_cfg.m_layer_configs[i].m_select_max_dphi);

    // try accessing something that does not exist

    // std::cout << "Non-existent path " << j["/m_layer_configs/143/m_select_max_dphi"_json_pointer] << "\n";

    auto &j = cjp.get("");

    auto &x = j["/m_layer_configs"_json_pointer];
    std::cout << "Typename /m_layer_configs " << x.type_name() << "\n";
    auto &y = j["/m_layer_configs/143"_json_pointer];
    std::cout << "Typename /m_layer_configs/143 " << y.type_name() << ", is_null=" << y.is_null() << "\n";
}

}
