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
    windDragCoeff_airDensity = 0.0025 * 1.2;

    InitialTimeStep = 3.e-5;
    YoungsModulus = 5.e8;
    GridXTotal = 128;
    GridY = 55;
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
    DP_threshold_p = 0;

    tpb_P2G = 256;
    tpb_Upd = 512;
    tpb_G2P = 128;

    GrainVariability = 0.50;

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

    if(doc.HasMember("LatMin")) LatMin = doc["LatMin"].GetDouble();
    if(doc.HasMember("LonMin")) LonMin = doc["LonMin"].GetDouble();
    if(doc.HasMember("LatMax")) LatMax = doc["LatMax"].GetDouble();
    if(doc.HasMember("LonMax")) LonMax = doc["LonMax"].GetDouble();

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
    if(doc.HasMember("GrainVariability")) GrainVariability = doc["GrainVariability"].GetDouble();

    if(doc.HasMember("tpb_P2G")) tpb_P2G = doc["tpb_P2G"].GetInt();
    if(doc.HasMember("tpb_Upd")) tpb_Upd = doc["tpb_Upd"].GetInt();
    if(doc.HasMember("tpb_G2P")) tpb_G2P = doc["tpb_G2P"].GetInt();

    std::map<std::string,std::string> result;

    if(!doc.HasMember("InputPNG"))
    {
        spdlog::critical("InputPNG entry is missing in JSON config file");
        throw std::runtime_error("config parameter missing");
    }

    result["InputPNG"] = doc["InputPNG"].GetString();
    spdlog::info("ParseFile; png map data {}", result["InputPNG"]);

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
    ParticleMass = ParticleVolume * SurfaceDensity;

    UpdateEveryNthStep = (int)(SimulationEndTime/(AnimationFramesRequested*InitialTimeStep));
    cellsize_inv = 1./cellsize; // cellsize itself is set when loading .h5 file
    Dp_inv = 4./(cellsize*cellsize);
    dt_vol_Dpinv = InitialTimeStep*ParticleVolume*Dp_inv;
    vmax = 0.5*cellsize/InitialTimeStep;
    vmax_squared = vmax*vmax;

    ComputeLame();
}




void SimParams::SaveParametersAsHDF5Attributes(H5::DataSet &dataset) {
    // Create a scalar dataspace for the attributes
    H5::DataSpace att_dspace(H5S_SCALAR);

    // Save each member as an attribute
    dataset.createAttribute("nPtsInitial", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &nPtsInitial);
    dataset.createAttribute("SimulationStartUnixTime", H5::PredType::NATIVE_INT64, att_dspace).write(H5::PredType::NATIVE_INT64, &SimulationStartUnixTime);
    dataset.createAttribute("GridXTotal", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &GridXTotal);
    dataset.createAttribute("GridY", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &GridY);

    dataset.createAttribute("LatMin", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &LatMin);
    dataset.createAttribute("LatMax", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &LatMax);
    dataset.createAttribute("LonMin", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &LonMin);
    dataset.createAttribute("LonMax", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &LonMax);

    dataset.createAttribute("gridLatMin", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &gridLatMin);
    dataset.createAttribute("gridLonMin", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &gridLonMin);
    dataset.createAttribute("DimensionHorizontal", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &DimensionHorizontal);
    dataset.createAttribute("DimensionVertical", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &DimensionVertical);
    dataset.createAttribute("InitialTimeStep", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &InitialTimeStep);
    dataset.createAttribute("SimulationEndTime", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &SimulationEndTime);

    dataset.createAttribute("AnimationFramesRequested", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &AnimationFramesRequested);
    dataset.createAttribute("windDragCoeff_airDensity", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &windDragCoeff_airDensity);
    dataset.createAttribute("SurfaceDensity", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &SurfaceDensity);

    dataset.createAttribute("PoissonsRatio", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &PoissonsRatio);
    dataset.createAttribute("YoungsModulus", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &YoungsModulus);

    dataset.createAttribute("IceCompressiveStrength", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &IceCompressiveStrength);

    dataset.createAttribute("IceTensileStrength", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &IceTensileStrength);
    dataset.createAttribute("IceShearStrength", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &IceShearStrength);
    dataset.createAttribute("IceTensileStrength2", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &IceTensileStrength2);
    dataset.createAttribute("DP_phi", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &DP_phi);
    dataset.createAttribute("DP_threshold_p", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &DP_threshold_p);
    dataset.createAttribute("cellsize", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &cellsize);

    dataset.createAttribute("ParticleVolume", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &ParticleVolume);
    dataset.createAttribute("ParticleViewSize", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &ParticleViewSize);
    dataset.createAttribute("GrainVariability", H5::PredType::NATIVE_DOUBLE, att_dspace).write(H5::PredType::NATIVE_DOUBLE, &GrainVariability);
}

void SimParams::ReadParametersFromHDF5Attributes(H5::DataSet &dataset) {
    // Read each attribute and assign it to the corresponding member variable
    dataset.openAttribute("nPtsInitial").read(H5::PredType::NATIVE_INT, &nPtsInitial);
    dataset.openAttribute("SimulationStartUnixTime").read(H5::PredType::NATIVE_INT64, &SimulationStartUnixTime);
    dataset.openAttribute("GridXTotal").read(H5::PredType::NATIVE_INT, &GridXTotal);
    dataset.openAttribute("GridY").read(H5::PredType::NATIVE_INT, &GridY);

    dataset.openAttribute("LatMin").read(H5::PredType::NATIVE_DOUBLE, &LatMin);
    dataset.openAttribute("LatMax").read(H5::PredType::NATIVE_DOUBLE, &LatMax);
    dataset.openAttribute("LonMin").read(H5::PredType::NATIVE_DOUBLE, &LonMin);
    dataset.openAttribute("LonMax").read(H5::PredType::NATIVE_DOUBLE, &LonMax);

    dataset.openAttribute("gridLatMin").read(H5::PredType::NATIVE_DOUBLE, &gridLatMin);
    dataset.openAttribute("gridLonMin").read(H5::PredType::NATIVE_DOUBLE, &gridLonMin);

    dataset.openAttribute("DimensionHorizontal").read(H5::PredType::NATIVE_DOUBLE, &DimensionHorizontal);
    dataset.openAttribute("DimensionVertical").read(H5::PredType::NATIVE_DOUBLE, &DimensionVertical);
    dataset.openAttribute("InitialTimeStep").read(H5::PredType::NATIVE_DOUBLE, &InitialTimeStep);

    dataset.openAttribute("SimulationEndTime").read(H5::PredType::NATIVE_DOUBLE, &SimulationEndTime);
    dataset.openAttribute("AnimationFramesRequested").read(H5::PredType::NATIVE_INT, &AnimationFramesRequested);
    dataset.openAttribute("windDragCoeff_airDensity").read(H5::PredType::NATIVE_DOUBLE, &windDragCoeff_airDensity);
    dataset.openAttribute("SurfaceDensity").read(H5::PredType::NATIVE_DOUBLE, &SurfaceDensity);
    dataset.openAttribute("PoissonsRatio").read(H5::PredType::NATIVE_DOUBLE, &PoissonsRatio);
    dataset.openAttribute("YoungsModulus").read(H5::PredType::NATIVE_DOUBLE, &YoungsModulus);
    dataset.openAttribute("IceCompressiveStrength").read(H5::PredType::NATIVE_DOUBLE, &IceCompressiveStrength);
    dataset.openAttribute("IceTensileStrength").read(H5::PredType::NATIVE_DOUBLE, &IceTensileStrength);
    dataset.openAttribute("IceShearStrength").read(H5::PredType::NATIVE_DOUBLE, &IceShearStrength);
    dataset.openAttribute("IceTensileStrength2").read(H5::PredType::NATIVE_DOUBLE, &IceTensileStrength2);
    dataset.openAttribute("DP_phi").read(H5::PredType::NATIVE_DOUBLE, &DP_phi);
    dataset.openAttribute("DP_threshold_p").read(H5::PredType::NATIVE_DOUBLE, &DP_threshold_p);
    dataset.openAttribute("cellsize").read(H5::PredType::NATIVE_DOUBLE, &cellsize);
    dataset.openAttribute("ParticleVolume").read(H5::PredType::NATIVE_DOUBLE, &ParticleVolume);
    dataset.openAttribute("ParticleViewSize").read(H5::PredType::NATIVE_DOUBLE, &ParticleViewSize);
    dataset.openAttribute("GrainVariability").read(H5::PredType::NATIVE_DOUBLE, &GrainVariability);

    tpb_P2G = 256;
    tpb_Upd = 512;
    tpb_G2P = 128;
    ComputeLame();
    ComputeHelperVariables();

}


