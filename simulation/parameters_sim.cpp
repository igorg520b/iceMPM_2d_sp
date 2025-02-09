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
    spdlog::info("SimParams reset; nPtsArrays {}", nPtsArrays);
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

    if(doc.HasMember("InitialTimeStep")) InitialTimeStep = doc["InitialTimeStep"].GetDouble();
    if(doc.HasMember("AnimationFramesRequested")) AnimationFramesRequested = doc["AnimationFramesRequested"].GetInt();

    if(doc.HasMember("YoungsModulus")) YoungsModulus = doc["YoungsModulus"].GetDouble();
    if(doc.HasMember("ParticleViewSize")) ParticleViewSize = doc["ParticleViewSize"].GetDouble();
    if(doc.HasMember("SimulationEndTime")) SimulationEndTime = doc["SimulationEndTime"].GetDouble();
    if(doc.HasMember("PoissonsRatio")) PoissonsRatio = doc["PoissonsRatio"].GetDouble();
    if(doc.HasMember("SurfaceDensity")) SurfaceDensity = doc["SurfaceDensity"].GetDouble();

    if(doc.HasMember("IceCompressiveStrength")) IceCompressiveStrength = doc["IceCompressiveStrength"].GetDouble();
    if(doc.HasMember("IceTensileStrength")) IceTensileStrength = doc["IceTensileStrength"].GetDouble();
    if(doc.HasMember("IceTensileStrength2")) IceTensileStrength2 = doc["IceTensileStrength2"].GetDouble();
    if(doc.HasMember("IceShearStrength")) IceShearStrength = doc["IceShearStrength"].GetDouble();

    if(doc.HasMember("DP_phi")) DP_phi = doc["DP_phi"].GetDouble();
    if(doc.HasMember("DP_threshold_p")) DP_threshold_p = doc["DP_threshold_p"].GetDouble();

    if(doc.HasMember("tpb_P2G")) tpb_P2G = doc["tpb_P2G"].GetInt();
    if(doc.HasMember("tpb_Upd")) tpb_Upd = doc["tpb_Upd"].GetInt();
    if(doc.HasMember("tpb_G2P")) tpb_G2P = doc["tpb_G2P"].GetInt();

    if(doc.HasMember("FluentDataScale")) FluentDataScale = doc["FluentDataScale"].GetDouble();
    if(doc.HasMember("FluentDataOffsetX")) FluentDataOffsetX = doc["FluentDataOffsetX"].GetDouble();
    if(doc.HasMember("FluentDataOffsetY")) FluentDataOffsetY = doc["FluentDataOffsetY"].GetDouble();
    if(doc.HasMember("FrameTimeInterval")) FrameTimeInterval = doc["FrameTimeInterval"].GetDouble();

    if(doc.HasMember("currentDragCoeff_waterDensity")) currentDragCoeff_waterDensity = doc["currentDragCoeff_waterDensity"].GetDouble();



    if(doc.HasMember("SaveSnapshots")) SaveSnapshots = doc["SaveSnapshots"].GetBool();


    std::map<std::string,std::string> result;

    if(!doc.HasMember("InputPNG") || !doc.HasMember("ModeledRegion"))
    {
        spdlog::critical("InputPNG and/or ModeledRegion entry is missing in JSON config file");
        throw std::runtime_error("config parameter missing");
    }
    result["InputPNG"] = doc["InputPNG"].GetString();
    result["ModeledRegion"] = doc["ModeledRegion"].GetString();
    spdlog::info("ParseFile; png map data {}", result["InputPNG"]);
    spdlog::info("ModeledRegion png {}", result["ModeledRegion"]);

    UseWindData = doc.HasMember("InputWindData");
    if(UseWindData)
    {
        std::string strFileName = doc["InputWindData"].GetString();
        result["InputWindData"] = strFileName;

        if (!std::filesystem::exists(strFileName))
        {
            spdlog::critical("file does not exist: {}", strFileName);
            throw std::runtime_error("file does not exist");
        }
    }

    UseCurrentData = doc.HasMember("InputCurrentData");
    if(UseCurrentData)
    {
        std::string strFileName = doc["InputCurrentData"].GetString();
        result["InputCurrentData"] = strFileName;
        if (!std::filesystem::exists(strFileName))
        {
            spdlog::critical("file does not exist: {}", strFileName);
            throw std::runtime_error("file does not exist");
        }
    }

    if(doc.HasMember("DimensionHorizontal")) DimensionHorizontal = doc["DimensionHorizontal"].GetDouble();
    else DimensionHorizontal = 0;

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
    spdlog::info("Simulation Parameters:");
    spdlog::info("SimulationStartUnixTime: {}", SimulationStartUnixTime);
    spdlog::info("InitialTimeStep: {}, SimulationEndTime: {}", InitialTimeStep, SimulationEndTime);
    spdlog::info("AnimationFramesRequested: {}", AnimationFramesRequested);
    spdlog::info("SimulationStep: {}", SimulationStep);
    spdlog::info("SimulationTime: {}", SimulationTime);
    spdlog::info("UpdateEveryNthStep: {}", UpdateEveryNthStep);
    spdlog::info("UseWindData: {}", UseWindData);
    spdlog::info("UseCurrentData: {}", UseCurrentData);

    // parameters
    spdlog::info("\n");
    spdlog::info("Parameters:");
    spdlog::info("dt_vol_Dpinv: {}, vmax: {}", dt_vol_Dpinv, vmax);
    spdlog::info("windDragCoeff_airDensity: {}", windDragCoeff_airDensity);
    spdlog::info("lambda: {}, mu: {}, kappa: {}", lambda, mu, kappa);
    spdlog::info("ParticleVolume: {}, ParticleViewSize: {}", ParticleVolume, ParticleViewSize);
    spdlog::info("ParticleMass: {}", ParticleMass);
    spdlog::info("DP_phi: {}, DP_threshold_p: {}", DP_phi, DP_threshold_p);
    spdlog::info("SurfaceDensity: {}, PoissonsRatio: {}, YoungsModulus: {}", SurfaceDensity, PoissonsRatio, YoungsModulus);
    spdlog::info("IceCompressiveStrength: {}, IceTensileStrength: {}, IceShearStrength: {}, IceTensileStrength2: {}",
                 IceCompressiveStrength, IceTensileStrength, IceShearStrength, IceTensileStrength2);

    // points
    spdlog::info("\n");
    spdlog::info("Points:");
    spdlog::info("nPtsInitial: {}", nPtsInitial);

    // grid
    spdlog::info("\n");
    spdlog::info("Grid:");
    spdlog::info("cellsize: {}; cellsize*InitializationImageSizeX", cellsize, cellsize*InitializationImageSizeX);
    spdlog::info("Sim grid: {} x {}", GridXTotal, GridYTotal);
    spdlog::info("Original image: {} x {}", InitializationImageSizeX, InitializationImageSizeY);
    spdlog::info("Offset of the modelled region: [{}, {}]", ModeledRegionOffsetX, ModeledRegionOffsetY);
}
