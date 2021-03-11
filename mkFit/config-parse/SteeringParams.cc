#include "SteeringParams.h"

#include "json.hpp"

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


ConfigJsonPatcher::ConfigJsonPatcher() {}

ConfigJsonPatcher::~ConfigJsonPatcher() { release_json(); }

void ConfigJsonPatcher::release_json()
{
    if (m_owner) delete m_json;
    m_json  = nullptr;
    m_owner = false;
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
template void ConfigJsonPatcher::Save<IterationsInfo> (IterationsInfo  &o);
template void ConfigJsonPatcher::Save<IterationConfig>(IterationConfig &o);

void ConfigJsonPatcher::cd(const std::string& path)
{
    nlohmann::json::json_pointer jp(path);
    m_current = & m_current->at(jp);
}

void ConfigJsonPatcher::cd_top(const std::string& path)
{
    m_current = m_json;
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
void ConfigJsonPatcher::replace(int beg, int end, const std::string& path, T val)
{
    nlohmann::json::json_pointer jp(path);
    for (int i = beg; i < end; ++i)
    {
        m_current->at(i).at(jp) = val;
    }
}
template void ConfigJsonPatcher::replace<int>   (int beg, int end, const std::string& path, int    val);
template void ConfigJsonPatcher::replace<float> (int beg, int end, const std::string& path, float  val);
template void ConfigJsonPatcher::replace<double>(int beg, int end, const std::string& path, double val);

nlohmann::json& ConfigJsonPatcher::get(const std::string& path)
{
    nlohmann::json::json_pointer jp(path);
    return m_current->at(jp);
}

std::string ConfigJsonPatcher::dump(int indent)
{
    return m_json->dump(indent);
}

// ============================================================================

void TestJson_Dump_Direct(IterationConfig &it_cfg)
{
    using nlohmann::json;
    // using nlohmann::json::json_pointer;

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

void TestJson_Dump_ConfigPatcher(IterationConfig &it_cfg)
{
    using nlohmann::json;
    // using nlohmann::json::json_pointer;

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
        cjp.replace(1, 4, "/m_select_max_dphi", 7.7);
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
