#include "parameters_sim.h"
#include <spdlog/spdlog.h>

#include <string>
#include <filesystem>

void SimParams::Reset()
{
    UseWindData = false;
    UseCurrentData = false;
    SimulationStartUnixTime = 0;
    DimensionHorizontal = 0;
    SaveSnapshots = false;

    nPtsInitial = 0;
    windDragCoeff_airDensity = 0.0025 * 1.2;
//    currentDragCoeff_waterDensity = 0.0025 * 1025;
    currentDragCoeff_waterDensity = 0.15 * 1025;

    InitialTimeStep = 3.e-5;
    YoungsModulus = 5.e8;
    ParticleViewSize = 2.5f;

    SimulationEndTime = 12;
    AnimationFramesRequested = 5000;

    PoissonsRatio = 0.3;
    SurfaceDensity = 916*0.2;  // surface density

    SimulationStep = 0;
    SimulationTime = 0;

    IceCompressiveStrength = 100e6;
    IceTensileStrength = 10e6;
    IceShearStrength = 1e6;
    IceTensileStrength2 = 10e6;

    DP_phi = 62;
    DP_threshold_p = 1e4;

    tpb_P2G = 256;
    tpb_Upd = 512;
    tpb_G2P = 128;

    FluentDataScale = 3;
    FluentDataOffsetX = 0;
    FluentDataOffsetY = 0;
    FrameTimeInterval = 1.0;

    ComputeLame();
    ComputeHelperVariables();
    spdlog::info("SimParams reset");
}



std::map<std::string,std::string> SimParams::ParseFile(std::string fileName)
{
    LOGR("SimParams ParseFile {}",fileName);
    if(!std::filesystem::exists(fileName)) throw std::runtime_error("configuration file is not found");
    std::ifstream fileStream(fileName);
    std::string strConfigFile;
    strConfigFile.resize(std::filesystem::file_size(fileName));
    fileStream.read(strConfigFile.data(), strConfigFile.length());
    fileStream.close();


    rapidjson::Document doc;
    doc.Parse(strConfigFile.data());
    if(!doc.IsObject()) throw std::runtime_error("configuration file is not JSON");

    // parse strings and save into "result"
    std::map<std::string,std::string> result;
    result["InputPNG"] = doc["InputPNG"].GetString();
    result["InputMap"] = doc["InputMap"].GetString();
    if(doc.HasMember("InputFlowVelocity")) result["InputFlowVelocity"] = doc["InputFlowVelocity"].GetString();
    if(doc.HasMember("SimulationTitle")) result["SimulationTitle"] = doc["SimulationTitle"].GetString();
    else result["SimulationTitle"] = "sim1";

    if(doc.HasMember("DimensionHorizontal")) DimensionHorizontal = doc["DimensionHorizontal"].GetDouble();
    if(doc.HasMember("SurfaceDensity")) SurfaceDensity = doc["SurfaceDensity"].GetDouble();
    if(doc.HasMember("SaveSnapshots")) SaveSnapshots = doc["SaveSnapshots"].GetBool();
    if(doc.HasMember("UseCurrentData")) UseCurrentData = doc["UseCurrentData"].GetBool();

    if(doc.HasMember("SimulationStartUnixTime")) SimulationStartUnixTime = doc["SimulationStartUnixTime"].GetInt64();
    if(doc.HasMember("SimulationEndTime")) SimulationEndTime = doc["SimulationEndTime"].GetDouble();

    if(doc.HasMember("InitialTimeStep")) InitialTimeStep = doc["InitialTimeStep"].GetDouble();
    if(doc.HasMember("ParticleViewSize")) ParticleViewSize = doc["ParticleViewSize"].GetDouble();
    if(doc.HasMember("AnimationFramesRequested")) AnimationFramesRequested = doc["AnimationFramesRequested"].GetInt();

    if(doc.HasMember("YoungsModulus")) YoungsModulus = doc["YoungsModulus"].GetDouble();
    if(doc.HasMember("PoissonsRatio")) PoissonsRatio = doc["PoissonsRatio"].GetDouble();

    if(doc.HasMember("IceCompressiveStrength")) IceCompressiveStrength = doc["IceCompressiveStrength"].GetDouble();
    if(doc.HasMember("IceTensileStrength")) IceTensileStrength = doc["IceTensileStrength"].GetDouble();
    if(doc.HasMember("IceTensileStrength2")) IceTensileStrength2 = doc["IceTensileStrength2"].GetDouble();
    if(doc.HasMember("IceShearStrength")) IceShearStrength = doc["IceShearStrength"].GetDouble();

    if(doc.HasMember("DP_phi")) DP_phi = doc["DP_phi"].GetDouble();
    if(doc.HasMember("DP_threshold_p")) DP_threshold_p = doc["DP_threshold_p"].GetDouble();

    if(doc.HasMember("tpb_P2G")) tpb_P2G = doc["tpb_P2G"].GetInt();
    if(doc.HasMember("tpb_Upd")) tpb_Upd = doc["tpb_Upd"].GetInt();
    if(doc.HasMember("tpb_G2P")) tpb_G2P = doc["tpb_G2P"].GetInt();

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
    ParticleMass = ParticleVolume * SurfaceDensity;

    UpdateEveryNthStep = (int)(SimulationEndTime/(AnimationFramesRequested*InitialTimeStep));
    cellsize_inv = 1./cellsize; // cellsize itself is set when loading .h5 file
    Dp_inv = 4./(cellsize*cellsize);
    dt_vol_Dpinv = InitialTimeStep*ParticleVolume*Dp_inv;
    vmax = 0.25*cellsize/InitialTimeStep;

    ComputeLame();
}




void SimParams::SaveParametersAsHDF5Attributes(H5::DataSet &dataset)
{
    // Create a scalar dataspace for the attributes
    H5::DataSpace att_dspace(H5S_SCALAR);

    // Save each member as an attribute

    // integration, time, initial
    dataset.createAttribute("nPtsInitial", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &nPtsInitial);
    dataset.createAttribute("SimulationStartUnixTime", H5::PredType::NATIVE_INT64, att_dspace).write(H5::PredType::NATIVE_INT64, &SimulationStartUnixTime);
    dataset.createAttribute("InitialTimeStep", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &InitialTimeStep);
    dataset.createAttribute("SimulationEndTime", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &SimulationEndTime);
    dataset.createAttribute("AnimationFramesRequested", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &AnimationFramesRequested);
    dataset.createAttribute("ParticleViewSize", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &ParticleViewSize);

    dataset.createAttribute("UseWindData", H5::PredType::NATIVE_HBOOL, att_dspace).write(H5::PredType::NATIVE_HBOOL, &UseWindData);
    dataset.createAttribute("UseCurrentData", H5::PredType::NATIVE_HBOOL, att_dspace).write(H5::PredType::NATIVE_HBOOL, &UseCurrentData);
    dataset.createAttribute("SaveSnapshots", H5::PredType::NATIVE_HBOOL, att_dspace).write(H5::PredType::NATIVE_HBOOL, &SaveSnapshots);

    // parameters
    dataset.createAttribute("windDragCoeff_airDensity", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &windDragCoeff_airDensity);
    dataset.createAttribute("currentDragCoeff_waterDensity", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &currentDragCoeff_waterDensity);
    dataset.createAttribute("SurfaceDensity", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &SurfaceDensity);
    dataset.createAttribute("PoissonsRatio", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &PoissonsRatio);
    dataset.createAttribute("YoungsModulus", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &YoungsModulus);
    dataset.createAttribute("IceCompressiveStrength", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &IceCompressiveStrength);
    dataset.createAttribute("IceTensileStrength", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &IceTensileStrength);
    dataset.createAttribute("IceShearStrength", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &IceShearStrength);
    dataset.createAttribute("IceTensileStrength2", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &IceTensileStrength2);
    dataset.createAttribute("DP_phi", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &DP_phi);
    dataset.createAttribute("DP_threshold_p", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &DP_threshold_p);
    dataset.createAttribute("ParticleVolume", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &ParticleVolume);

    // grid
    dataset.createAttribute("cellsize", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &cellsize);

    dataset.createAttribute("GridXTotal", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &GridXTotal);
    dataset.createAttribute("GridYTotal", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &GridYTotal);
    dataset.createAttribute("ModeledRegionOffsetX", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &ModeledRegionOffsetX);
    dataset.createAttribute("ModeledRegionOffsetY", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &ModeledRegionOffsetY);
    dataset.createAttribute("InitializationImageSizeX", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &InitializationImageSizeX);
    dataset.createAttribute("InitializationImageSizeY", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &InitializationImageSizeY);

    // flow render / interpolator
    dataset.createAttribute("FrameTimeInterval", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &FrameTimeInterval);
}

void SimParams::ReadParametersFromHDF5Attributes(H5::DataSet &dataset)
{
    // Read each attribute and assign it to the corresponding member variable

    // integration, time, initial
    dataset.openAttribute("nPtsInitial").read(H5::PredType::NATIVE_INT, &nPtsInitial);
    dataset.openAttribute("SimulationStartUnixTime").read(H5::PredType::NATIVE_INT64, &SimulationStartUnixTime);
    dataset.openAttribute("InitialTimeStep").read(H5::PredType::NATIVE_DOUBLE, &InitialTimeStep);
    dataset.openAttribute("SimulationEndTime").read(H5::PredType::NATIVE_DOUBLE, &SimulationEndTime);
    dataset.openAttribute("AnimationFramesRequested").read(H5::PredType::NATIVE_INT, &AnimationFramesRequested);
    dataset.openAttribute("ParticleViewSize").read(H5::PredType::NATIVE_DOUBLE, &ParticleViewSize);

    dataset.openAttribute("UseWindData").read(H5::PredType::NATIVE_HBOOL, &UseWindData);
    dataset.openAttribute("UseCurrentData").read(H5::PredType::NATIVE_HBOOL, &UseCurrentData);
    dataset.openAttribute("SaveSnapshots").read(H5::PredType::NATIVE_HBOOL, &SaveSnapshots);

    // parameters
    dataset.openAttribute("windDragCoeff_airDensity").read(H5::PredType::NATIVE_DOUBLE, &windDragCoeff_airDensity);
    dataset.openAttribute("currentDragCoeff_waterDensity").read(H5::PredType::NATIVE_DOUBLE, &currentDragCoeff_waterDensity);

    dataset.openAttribute("SurfaceDensity").read(H5::PredType::NATIVE_DOUBLE, &SurfaceDensity);
    dataset.openAttribute("PoissonsRatio").read(H5::PredType::NATIVE_DOUBLE, &PoissonsRatio);
    dataset.openAttribute("YoungsModulus").read(H5::PredType::NATIVE_DOUBLE, &YoungsModulus);
    dataset.openAttribute("IceCompressiveStrength").read(H5::PredType::NATIVE_DOUBLE, &IceCompressiveStrength);
    dataset.openAttribute("IceTensileStrength").read(H5::PredType::NATIVE_DOUBLE, &IceTensileStrength);
    dataset.openAttribute("IceShearStrength").read(H5::PredType::NATIVE_DOUBLE, &IceShearStrength);
    dataset.openAttribute("IceTensileStrength2").read(H5::PredType::NATIVE_DOUBLE, &IceTensileStrength2);
    dataset.openAttribute("DP_phi").read(H5::PredType::NATIVE_DOUBLE, &DP_phi);
    dataset.openAttribute("DP_threshold_p").read(H5::PredType::NATIVE_DOUBLE, &DP_threshold_p);
    dataset.openAttribute("ParticleVolume").read(H5::PredType::NATIVE_DOUBLE, &ParticleVolume);

    // grid
    dataset.openAttribute("cellsize").read(H5::PredType::NATIVE_DOUBLE, &cellsize);
    dataset.openAttribute("GridXTotal").read(H5::PredType::NATIVE_INT, &GridXTotal);
    dataset.openAttribute("GridYTotal").read(H5::PredType::NATIVE_INT, &GridYTotal);

    dataset.openAttribute("ModeledRegionOffsetX").read(H5::PredType::NATIVE_INT, &ModeledRegionOffsetX);
    dataset.openAttribute("ModeledRegionOffsetY").read(H5::PredType::NATIVE_INT, &ModeledRegionOffsetY);
    dataset.openAttribute("InitializationImageSizeX").read(H5::PredType::NATIVE_INT, &InitializationImageSizeX);
    dataset.openAttribute("InitializationImageSizeY").read(H5::PredType::NATIVE_INT, &InitializationImageSizeY);

    // flow render / interpolator
    dataset.openAttribute("FrameTimeInterval").read(H5::PredType::NATIVE_DOUBLE, &FrameTimeInterval);

    tpb_P2G = 256;
    tpb_Upd = 512;
    tpb_G2P = 128;
    ComputeLame();
    ComputeHelperVariables();
}


void SimParams::Printout()
{
    spdlog::info(fmt::format(fmt::runtime("Simulation Parameters:")));
    spdlog::info(fmt::format(fmt::runtime("SimulationStartUnixTime: {}"), SimulationStartUnixTime));
    spdlog::info(fmt::format(fmt::runtime("InitialTimeStep: {}, SimulationEndTime: {}"), InitialTimeStep, SimulationEndTime));
    spdlog::info(fmt::format(fmt::runtime("AnimationFramesRequested: {}"), AnimationFramesRequested));
    spdlog::info(fmt::format(fmt::runtime("SimulationStep: {}"), SimulationStep));
    spdlog::info(fmt::format(fmt::runtime("SimulationTime: {}"), SimulationTime));
    spdlog::info(fmt::format(fmt::runtime("UpdateEveryNthStep: {}"), UpdateEveryNthStep));
    spdlog::info(fmt::format(fmt::runtime("UseWindData: {}"), UseWindData));
    spdlog::info(fmt::format(fmt::runtime("UseCurrentData: {}"), UseCurrentData));

    // parameters
    spdlog::info("");
    spdlog::info(fmt::format(fmt::runtime("Parameters:")));
    spdlog::info(fmt::format(fmt::runtime("dt_vol_Dpinv: {}, vmax: {}"), dt_vol_Dpinv, vmax));
    spdlog::info(fmt::format(fmt::runtime("windDragCoeff_airDensity: {}"), windDragCoeff_airDensity));
    spdlog::info(fmt::format(fmt::runtime("lambda: {}, mu: {}, kappa: {}"), lambda, mu, kappa));
    spdlog::info(fmt::format(fmt::runtime("ParticleVolume: {}, ParticleViewSize: {}"), ParticleVolume, ParticleViewSize));
    spdlog::info(fmt::format(fmt::runtime("ParticleMass: {}"), ParticleMass));
    spdlog::info(fmt::format(fmt::runtime("DP_phi: {}, DP_threshold_p: {}"), DP_phi, DP_threshold_p));
    spdlog::info(fmt::format(fmt::runtime("SurfaceDensity: {}, PoissonsRatio: {}, YoungsModulus: {}"),
                             SurfaceDensity, PoissonsRatio, YoungsModulus));
    spdlog::info(fmt::format(fmt::runtime("IceCompressiveStrength: {}, IceTensileStrength: {}, IceShearStrength: {}, IceTensileStrength2: {}"),
                             IceCompressiveStrength, IceTensileStrength, IceShearStrength, IceTensileStrength2));

    // points
    spdlog::info("");
    spdlog::info(fmt::format(fmt::runtime("Points:")));
    spdlog::info(fmt::format(fmt::runtime("nPtsInitial: {}"), nPtsInitial));

    // grid
    spdlog::info("");
    spdlog::info(fmt::format(fmt::runtime("Grid:")));
    spdlog::info(fmt::format(fmt::runtime("cellsize: {}; cellsize*InitializationImageSizeX: {}"),
                             cellsize, cellsize * InitializationImageSizeX));
    spdlog::info(fmt::format(fmt::runtime("Sim grid: {} x {}"), GridXTotal, GridYTotal));
    spdlog::info(fmt::format(fmt::runtime("Original image: {} x {}"), InitializationImageSizeX, InitializationImageSizeY));
    spdlog::info(fmt::format(fmt::runtime("Offset of the modelled region: [{}, {}]"),
                             ModeledRegionOffsetX, ModeledRegionOffsetY));
}
