#include "parameters_sim.h"
#include <spdlog/spdlog.h>


void SimParams::Reset()
{
    SimulationStartUnixTime = -1;

    LatMin = 63.2540559;
    LonMin = 19.7794673;
    LatMax = 63.881440;
    LonMax = 21.877357;

    nPtsInitial = 0;
    GlenA = 1e-24;
    windDragCoeff_airDensity = 0.0025 * 1.2;
    waterDrag_waterDensity = 0.000005 * 1025;

    InitialTimeStep = 3.e-5;
    YoungsModulus = 5.e8;
    GridXTotal = 128;
    GridY = 55;
    ParticleViewSize = 2.5f;

    SimulationEndTime = 12;
    AnimationFramesRequested = 5000;

    PoissonsRatio = 0.3;
    Gravity = 9.81;
    Density = 916*0.2;  // surface density

    SimulationStep = 0;
    SimulationTime = 0;

    IceCompressiveStrength = 100e6;
    IceTensileStrength = 10e6;
    IceShearStrength = 1e6;
    IceTensileStrength2 = 10e6;

    DP_tan_phi = std::tan(30*pi/180.);
    DP_threshold_p = 0;

    tpb_P2G = 256;
    tpb_Upd = 512;
    tpb_G2P = 128;

    SetupType = 0;
    GrainVariability = 0.50;

    ComputeLame();
    ComputeCamClayParams2();
    ComputeHelperVariables();
    spdlog::info("SimParams reset");
}


std::map<std::string,std::string> SimParams::ParseFile(std::string fileName)
{
    spdlog::info("SimParams ParseFile {}",fileName);
    if(!std::filesystem::exists(fileName)) throw std::runtime_error("configuration file is not found");
    std::ifstream fileStream(fileName);
    std::string strConfigFile;
    strConfigFile.resize(std::filesystem::file_size(fileName));
    fileStream.read(strConfigFile.data(), strConfigFile.length());
    fileStream.close();


    rapidjson::Document doc;
    doc.Parse(strConfigFile.data());
    if(!doc.IsObject()) throw std::runtime_error("configuration file is not JSON");

    if(doc.HasMember("SimulationStartUnixTime")) SimulationStartUnixTime = doc["SimulationStartUnixTime"].GetInt64();

    if(doc.HasMember("LatMin")) LatMin = doc["LatMin"].GetDouble();
    if(doc.HasMember("LonMin")) LonMin = doc["LonMin"].GetDouble();
    if(doc.HasMember("LatMax")) LatMax = doc["LatMax"].GetDouble();
    if(doc.HasMember("LonMax")) LonMax = doc["LonMax"].GetDouble();

    if(doc.HasMember("SetupType")) SetupType = doc["SetupType"].GetInt();
    if(doc.HasMember("InitialTimeStep")) InitialTimeStep = doc["InitialTimeStep"].GetDouble();
    if(doc.HasMember("AnimationFramesRequested")) AnimationFramesRequested = doc["AnimationFramesRequested"].GetInt();

    if(doc.HasMember("YoungsModulus")) YoungsModulus = doc["YoungsModulus"].GetDouble();
    if(doc.HasMember("ParticleViewSize")) ParticleViewSize = doc["ParticleViewSize"].GetDouble();
    if(doc.HasMember("SimulationEndTime")) SimulationEndTime = doc["SimulationEndTime"].GetDouble();
    if(doc.HasMember("PoissonsRatio")) PoissonsRatio = doc["PoissonsRatio"].GetDouble();
    if(doc.HasMember("Gravity")) Gravity = doc["Gravity"].GetDouble();
    if(doc.HasMember("Density")) Density = doc["Density"].GetDouble();

    if(doc.HasMember("IceCompressiveStrength")) IceCompressiveStrength = doc["IceCompressiveStrength"].GetDouble();
    if(doc.HasMember("IceTensileStrength")) IceTensileStrength = doc["IceTensileStrength"].GetDouble();
    if(doc.HasMember("IceTensileStrength2")) IceTensileStrength2 = doc["IceTensileStrength2"].GetDouble();
    if(doc.HasMember("IceShearStrength")) IceShearStrength = doc["IceShearStrength"].GetDouble();

    if(doc.HasMember("DP_phi")) DP_tan_phi = std::tan(doc["DP_phi"].GetDouble()*pi/180);
    if(doc.HasMember("DP_threshold_p")) DP_threshold_p = doc["DP_threshold_p"].GetDouble();
    if(doc.HasMember("GrainVariability")) GrainVariability = doc["GrainVariability"].GetDouble();

    if(doc.HasMember("tpb_P2G")) tpb_P2G = doc["tpb_P2G"].GetInt();
    if(doc.HasMember("tpb_Upd")) tpb_Upd = doc["tpb_Upd"].GetInt();
    if(doc.HasMember("tpb_G2P")) tpb_G2P = doc["tpb_G2P"].GetInt();

//    ComputeCamClayParams2();
//    ComputeHelperVariables();


    std::map<std::string,std::string> result;

    if(!doc.HasMember("InputRawPoints"))
    {
        spdlog::critical("InputRawPoints entry is missing in JSON config file");
        throw std::runtime_error("config parameter missing");
    }

    result["InputRawPoints"] = doc["InputRawPoints"].GetString();
    spdlog::info("ParseFile; raw point data {}", result["InputRawPoints"]);

    if(doc.HasMember("InputWindData"))
    {
        result["InputWindData"] = doc["InputWindData"].GetString();
        spdlog::info("ParseFile; InputWindData file {}", result["InputWindData"]);
    }

    spdlog::info("SimParams::ParseFile done");

    return result;
}

void SimParams::ComputeLame()
{
    lambda = YoungsModulus*PoissonsRatio/((1+PoissonsRatio)*(1-2*PoissonsRatio));
    mu = YoungsModulus/(2*(1+PoissonsRatio));
    kappa = mu*2./3. + lambda;
}

void SimParams::ComputeHelperVariables()
{
    ParticleMass = ParticleVolume * Density;

    UpdateEveryNthStep = (int)(SimulationEndTime/(AnimationFramesRequested*InitialTimeStep));
    cellsize_inv = 1./cellsize; // cellsize itself is set when loading .h5 file
    Dp_inv = 4./(cellsize*cellsize);
    dt_vol_Dpinv = InitialTimeStep*ParticleVolume*Dp_inv;
    dt_Gravity = InitialTimeStep*Gravity;
    vmax = 0.5*cellsize/InitialTimeStep;
    vmax_squared = vmax*vmax;
}

void SimParams::ComputeCamClayParams2()
{
    ComputeLame();
    NACC_beta = IceTensileStrength/IceCompressiveStrength;
    const double &beta = NACC_beta;
    const double &q = IceShearStrength;
    const double &p0 = IceCompressiveStrength;
    NACC_M = (2*q*sqrt(1+2*beta))/(p0*(1+beta));
    NACC_Msq = NACC_M*NACC_M;
    spdlog::info("ComputeCamClayParams2() done");
}

void SimParams::ComputeIntegerBlockCoords()
{
    nxmin = floor(xmin/cellsize);
    nxmax = ceil(xmax/cellsize);
    nymin = floor(ymin/cellsize);
    nymax = ceil(ymax/cellsize);
}
