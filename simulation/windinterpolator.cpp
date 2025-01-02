#include "windinterpolator.h"
#include <H5Cpp.h>
#include <spdlog/spdlog.h>
#include <iostream>

WindInterpolator::WindInterpolator()
{
    spdlog::info("WindInterpolator compile-time data storage {} bytes",allocatedLonExtent*allocatedLatExtent*4*sizeof(float));
}

void WindInterpolator::LoadWindData(std::string netCDF4FileName,
                                    double latMin, double latMax, double lonMin, double lonMax, int64_t sim_start_date)
{
    std::vector<double> latitudes, longitudes;
    std::vector<int64_t> valid_time;
    int idxLatMin, idxLatMax, idxLonMin, idxLonMax;

    simulation_start_date = sim_start_date;
    H5::H5File file(netCDF4FileName, H5F_ACC_RDONLY);

    // READ INDICES FROM HDF5
    {
        // Read latitude dataset
        H5::DataSet latDataset = file.openDataSet("latitude");
        H5::DataSpace latSpace = latDataset.getSpace();
        hsize_t latDims[1];
        latSpace.getSimpleExtentDims(latDims, nullptr);
        latitudes.resize(latDims[0]);
        latDataset.read(latitudes.data(), H5::PredType::NATIVE_DOUBLE);

        // Read longitude dataset
        H5::DataSet lonDataset = file.openDataSet("longitude");
        H5::DataSpace lonSpace = lonDataset.getSpace();
        hsize_t lonDims[1];
        lonSpace.getSimpleExtentDims(lonDims, nullptr);
        longitudes.resize(lonDims[0]);
        lonDataset.read(longitudes.data(), H5::PredType::NATIVE_DOUBLE);

        // read valid_time
        H5::DataSet vtds = file.openDataSet("valid_time");
        H5::DataSpace vtSpace = vtds.getSpace();
        hsize_t vtDims[1];
        vtSpace.getSimpleExtentDims(vtDims, nullptr);
        nTimeIntervals = vtDims[0];
        valid_time.resize(vtDims[0]);
        vtds.read(valid_time.data(), H5::PredType::NATIVE_INT64);
        valid_time_front = valid_time.front();
    }


    // VERIFY DELTAS
    {
        // Verify latitude deltas
        for (size_t i = 1; i < latitudes.size(); ++i) {
            double delta = latitudes[i - 1] - latitudes[i];
            if (std::abs(delta - gridCellSize) > 1e-9) {
                spdlog::error("Latitude delta verification failed at index {}: delta = {}", i, delta);
                throw std::runtime_error("Latitude delta verification failed");
            }
        }
        spdlog::info("Latitude deltas verified successfully.");

        // Verify longitude deltas
        for (size_t i = 1; i < longitudes.size(); ++i) {
            double delta = longitudes[i] - longitudes[i - 1];
            if (std::abs(delta - gridCellSize) > 1e-9) {
                spdlog::error("Longitude delta verification failed at index {}: delta = {}", i, delta);
                throw std::runtime_error("Longitude delta verification failed");
            }
        }
        spdlog::info("Longitude deltas verified successfully.");

        // Verify valid_time deltas
        for (size_t i = 1; i < valid_time.size(); ++i) {
            int64_t delta = valid_time[i] - valid_time[i - 1];
            if (delta != timeResolution) {
                spdlog::error("Valid_time delta verification failed at index {}: delta = {}", i, delta);
                throw std::runtime_error("Valid_time delta verification failed");
            }
        }
        spdlog::info("Valid_time deltas verified successfully.");
    }


    LatMin = latMin;
    LatMax = latMax;
    LonMin = lonMin;
    LonMax = lonMax;

    // FIND INDICES OF OVERLAPPING REGION
    {
        // Find indices for latitude (descending order)
        auto itLatMax = std::lower_bound(latitudes.begin(), latitudes.end(), LatMax, std::greater<double>());
        idxLatMax = std::distance(latitudes.begin(), itLatMax)-1;
        if (itLatMax == latitudes.end()) {
            spdlog::error("No elements in the latitude vector are greater than or equal to LatMax ({})", LatMax);
            throw std::runtime_error("No elements in the latitude vector are greater than or equal to LatMax");
        }
        else if(itLatMax == latitudes.begin())
        {
            throw std::runtime_error("LatMax: range is not covered");
        }

        auto itLatMin = std::upper_bound(latitudes.begin(), latitudes.end(), LatMin, std::greater<double>());
        idxLatMin = std::distance(latitudes.begin(), itLatMin);

        if (itLatMin == latitudes.end()) {
            spdlog::error("No elements in the latitude vector are less than or equal to LatMin ({})", LatMin);
            throw std::runtime_error("No elements in the latitude vector are less than or equal to LatMin");
        }

        // Find indices for longitude (ascending order)
        auto itLonMin = std::lower_bound(longitudes.begin(), longitudes.end(), LonMin);
        if (itLonMin == longitudes.end()) {
            spdlog::error("No elements in the longitude vector are greater than or equal to LonMin ({})", LonMin);
            throw std::runtime_error("No elements in the longitude vector are greater than or equal to LonMin");
        } else if (itLonMin == longitudes.begin()) {
            spdlog::error("No elements in the longitude vector are lower than LonMin ({})", LonMin);
            throw std::runtime_error("Region cannot be covered");
        }
        idxLonMin = std::distance(longitudes.begin(), itLonMin) - 1;

        auto itLonMax = std::upper_bound(longitudes.begin(), longitudes.end(), LonMax);
        if (itLonMax == longitudes.end()) {
            spdlog::error("No elements in the longitude vector are greater than LonMax ({})", LonMax);
            throw std::runtime_error("No elements in the longitude vector are greater than LonMax");
        }
        idxLonMax = std::distance(longitudes.begin(), itLonMax);

        // Ensure the identified region fully overlaps the defined region
        if (!(latitudes[idxLatMin] < LatMin && latitudes[idxLatMax] > LatMax &&
              longitudes[idxLonMin] < LonMin && longitudes[idxLonMax] > LonMax)) {
            spdlog::error("The identified region does not fully overlap the defined region.");
            throw std::runtime_error("The identified region does not fully overlap the defined region.");
        }

        extentLat = idxLatMin - idxLatMax + 1;  // latitudes are in descending order
        extentLon = idxLonMax - idxLonMin + 1;
        gridLatMin = latitudes[idxLatMin];
        gridLonMin = longitudes[idxLonMin];

        spdlog::info("are of inerest origin; lat {:05.4f}; lon {:05.4f}", LatMin, LonMin);
        spdlog::info("grid origin lat {:05.2f}; lon {:05.2f}", gridLatMin, gridLonMin);
        spdlog::info("grid size (lat x lon) {} x {}", extentLat, extentLon);
        spdlog::info("idxLatMin-Max [{}, {}]; idxLonMin-Max [{}, {}]", idxLatMin, idxLatMax, idxLonMin, idxLonMax);

        if(extentLat > allocatedLatExtent || extentLon > allocatedLonExtent)
        {
            spdlog::error("wind grid extent exceeds the compile-time allocated space");
            throw std::runtime_error("wind grid extent exceeds the compile-time allocated space");
        }
    }

    // READ U AND V FROM HDF5
    {
        // Open dataset and dataspace
        H5::DataSet uDataset = file.openDataSet("u");
        H5::DataSpace uSpace = uDataset.getSpace();

        hsize_t uDims[4];
        uSpace.getSimpleExtentDims(uDims, nullptr);

        if (uDims[1] != 1) {
            spdlog::error("The mysterious dimension (probably pressure values) must have size 1, but found {}.", uDims[1]);
            throw std::runtime_error("Invalid mysterious dimension size");
        }

        hsize_t offset[4] = {0, 0, (hsize_t)idxLatMax, (hsize_t)idxLonMin};
        hsize_t count[4] = {(hsize_t)nTimeIntervals, 1, (hsize_t)extentLat, (hsize_t)extentLon};
        spdlog::info("count[4]: {}, {}, {}, {}", count[0], count[1], count[2], count[3]);

        H5::DataSpace memSpace(4, count);
        uSpace.selectHyperslab(H5S_SELECT_SET, count, offset);

        u_data.resize(count[0] * count[2] * count[3]);
        spdlog::info("u_data size is now {}",u_data.size());
        uDataset.read(u_data.data(), H5::PredType::NATIVE_FLOAT, memSpace, uSpace);


        // similarly, laod the V-dataset
        H5::DataSet vDataset = file.openDataSet("v");
        H5::DataSpace vSpace = vDataset.getSpace();
        vSpace.selectHyperslab(H5S_SELECT_SET, count, offset);

        v_data.resize(count[0] * count[2] * count[3]);
        vDataset.read(v_data.data(), H5::PredType::NATIVE_FLOAT, memSpace, vSpace);
    }

    isInitialized = true;
    spdlog::info("WindInterpolator::LoadWindData done");
}



size_t WindInterpolator::getIndexInUV(size_t timeIndex, size_t idxLat, size_t idxLon)
{
    return timeIndex * extentLat * extentLon + (extentLat - idxLat - 1) * extentLon + idxLon;
}

bool WindInterpolator::setTime(const double simulation_time)
{
    if(!isInitialized) return false;
    uint64_t sim_time_seconds = static_cast<int64_t>(simulation_time) + simulation_start_date;
    int interval = (int)((sim_time_seconds - valid_time_front) / timeResolution);

    if(currentInterval != interval)
    {
        if(interval > nTimeIntervals-2 || interval < 0)
        {
            spdlog::error("WindInterpolator::setTime selected time is out of range");
            throw std::runtime_error("WindInteroplator time is out of range");
        }
        currentInterval = interval;

        // fill the grid[] array
        for(int idxLat = 0; idxLat < extentLat; idxLat++)
            for(int idxLon = 0; idxLon < extentLon; idxLon++)
            {
                grid[idxLat][idxLon][0] = u_data[getIndexInUV(interval, idxLat, idxLon)];
                grid[idxLat][idxLon][1] = v_data[getIndexInUV(interval, idxLat, idxLon)];
                grid[idxLat][idxLon][2] = u_data[getIndexInUV(interval+1, idxLat, idxLon)];
                grid[idxLat][idxLon][3] = v_data[getIndexInUV(interval+1, idxLat, idxLon)];

                gridv[idxLat][idxLon][0].x() = u_data[getIndexInUV(interval, idxLat, idxLon)];
                gridv[idxLat][idxLon][0].y() = v_data[getIndexInUV(interval, idxLat, idxLon)];
                gridv[idxLat][idxLon][1].x() = u_data[getIndexInUV(interval+1, idxLat, idxLon)];
                gridv[idxLat][idxLon][1].y() = v_data[getIndexInUV(interval+1, idxLat, idxLon)];
            }
        spdlog::info("UPDATE REQUIRED");
        return true;
    }
    return false;
}

float WindInterpolator::interpolationCoeffFromTime(const double simulation_time)
{
    if(!isInitialized) return 0;
    uint64_t sim_time_seconds = static_cast<int64_t>(simulation_time) + simulation_start_date;
    int interval = (int)((sim_time_seconds - valid_time_front) / timeResolution);
    float interpolationParam = (float)(simulation_time + (simulation_start_date - valid_time_front - interval*timeResolution))/(float)timeResolution;
    if(currentInterval != interval) throw std::runtime_error("WindInterpolator::interpolationCoeffFromTime");
    return interpolationParam;
}


Eigen::Vector2f WindInterpolator::interpolationResult(float lat, float lon, float tb)
{
    if(!isInitialized) return Eigen::Vector2f::Zero();

    // space
    int lat_cell = (int)((lat-gridLatMin)/gridCellSize);
    int lon_cell = (int)((lon-gridLonMin)/gridCellSize);

    // Compute local coordinates within the cell
    float localLon = lon - (gridLonMin + lon_cell * gridCellSize);
    float localLat = lat - (gridLatMin + lat_cell * gridCellSize);

    // Compute barycentric coordinates
    float ub = localLon / gridCellSize;
    float vb = localLat / gridCellSize;

//    spdlog::info("input values: lat {}; lon {}", lat, lon);
//    spdlog::info("cell index {}, {}", lat_cell, lon_cell);
//    spdlog::info("local coords v, u: {}, {}; tb {}", vb, ub, tb);

    Eigen::Vector2f cell_values0[2][2], cell_values1[2][2];
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
        {
            cell_values0[i][j] = Eigen::Vector2f(grid[lat_cell+i][lon_cell+j][0],grid[lat_cell+i][lon_cell+j][1]);
            cell_values1[i][j] = Eigen::Vector2f(grid[lat_cell+i][lon_cell+j][2],grid[lat_cell+i][lon_cell+j][3]);
        }
    Eigen::Vector2f ipVal[2];

    ipVal[0] =
        (1 - ub) * (1 - vb) * cell_values0[0][0] +
        ub * (1 - vb) * cell_values0[0][1] +
        (1 - ub) * vb * cell_values0[1][0] +
        ub * vb * cell_values0[1][1];

    ipVal[1] =
        (1 - ub) * (1 - vb) * cell_values1[0][0] +
        ub * (1 - vb) * cell_values1[0][1] +
        (1 - ub) * vb * cell_values1[1][0] +
        ub * vb * cell_values1[1][1];
/*

    for(int k=0;k<2;k++)
    {
        ipVal[k] = (1 - ub) * (1 - vb) * gridv[lat_cell][lon_cell][k] +
                   ub * (1 - vb) * gridv[lat_cell][lon_cell+1][k] +
                   (1 - ub) * vb * gridv[lat_cell+1][lon_cell][k] +
                   ub * vb * gridv[lat_cell+1][lon_cell+1][k];
    }

*/

    Eigen::Vector2f final_result = (1-tb)*ipVal[0] + tb*ipVal[1];
    return final_result;
}


void WindInterpolator::SaveToOwnHDF5(H5::H5File &file)
{
    hsize_t dims_wind[3] = {(hsize_t)nTimeIntervals,(hsize_t)extentLat, (hsize_t)extentLon};
    H5::DataSpace dataspace_wind(3, dims_wind);

    H5::DataSet u_ds = file.createDataSet("u_data", H5::PredType::NATIVE_FLOAT, dataspace_wind);
    u_ds.write(u_data.data(), H5::PredType::NATIVE_FLOAT);
    H5::DataSet v_dataset = file.createDataSet("v_data", H5::PredType::NATIVE_FLOAT, dataspace_wind);
    v_dataset.write(v_data.data(), H5::PredType::NATIVE_FLOAT);

    H5::DataSpace att_sp(H5S_SCALAR);
    u_ds.createAttribute("extentLat", H5::PredType::NATIVE_INT, att_sp).write(H5::PredType::NATIVE_INT, &extentLat);
    u_ds.createAttribute("extentLon", H5::PredType::NATIVE_INT, att_sp).write(H5::PredType::NATIVE_INT, &extentLon);
    u_ds.createAttribute("nTimeIntervals", H5::PredType::NATIVE_INT, att_sp).write(H5::PredType::NATIVE_INT, &nTimeIntervals);

    u_ds.createAttribute("gridLatMin", H5::PredType::NATIVE_DOUBLE, att_sp).write(H5::PredType::NATIVE_DOUBLE, &gridLatMin);
    u_ds.createAttribute("gridLonMin", H5::PredType::NATIVE_DOUBLE, att_sp).write(H5::PredType::NATIVE_DOUBLE, &gridLonMin);

    u_ds.createAttribute("valid_time_front", H5::PredType::NATIVE_INT64, att_sp).write(H5::PredType::NATIVE_INT64, &valid_time_front);
    u_ds.createAttribute("simulation_start_date", H5::PredType::NATIVE_INT64, att_sp).write(H5::PredType::NATIVE_INT64, &simulation_start_date);

    u_ds.createAttribute("LatMin", H5::PredType::NATIVE_DOUBLE, att_sp).write(H5::PredType::NATIVE_DOUBLE, &LatMin);
    u_ds.createAttribute("LatMax", H5::PredType::NATIVE_DOUBLE, att_sp).write(H5::PredType::NATIVE_DOUBLE, &LatMax);
    u_ds.createAttribute("LonMin", H5::PredType::NATIVE_DOUBLE, att_sp).write(H5::PredType::NATIVE_DOUBLE, &LonMin);
    u_ds.createAttribute("LonMax", H5::PredType::NATIVE_DOUBLE, att_sp).write(H5::PredType::NATIVE_DOUBLE, &LonMax);
}

void WindInterpolator::ReadFromOwnHDF5(H5::H5File &file)
{
    spdlog::info("WindInterpolator::ReadFromOwnHDF5");
    H5::DataSet uDataset = file.openDataSet("u_data");

    uDataset.openAttribute("extentLat").read(H5::PredType::NATIVE_INT, &extentLat);
    uDataset.openAttribute("extentLon").read(H5::PredType::NATIVE_INT, &extentLon);
    uDataset.openAttribute("nTimeIntervals").read(H5::PredType::NATIVE_INT, &nTimeIntervals);

    H5::DataSpace dataspace = uDataset.getSpace();
    // Get the dimensions of the dataset
    hsize_t dims[3];
    dataspace.getSimpleExtentDims(dims, nullptr);
    if(dims[0] != nTimeIntervals || dims[1] != extentLat || dims[2] != extentLon)
        throw std::runtime_error("WindInterpolator::Read dataset size mismatch");


    uDataset.openAttribute("gridLatMin").read(H5::PredType::NATIVE_DOUBLE, &gridLatMin);
    uDataset.openAttribute("gridLonMin").read(H5::PredType::NATIVE_DOUBLE, &gridLonMin);

    uDataset.openAttribute("valid_time_front").read(H5::PredType::NATIVE_INT64, &valid_time_front);
    uDataset.openAttribute("simulation_start_date").read(H5::PredType::NATIVE_INT64, &simulation_start_date);

    uDataset.openAttribute("LatMin").read(H5::PredType::NATIVE_DOUBLE, &LatMin);
    uDataset.openAttribute("LatMax").read(H5::PredType::NATIVE_DOUBLE, &LatMax);
    uDataset.openAttribute("LonMin").read(H5::PredType::NATIVE_DOUBLE, &LonMin);
    uDataset.openAttribute("LonMax").read(H5::PredType::NATIVE_DOUBLE, &LonMax);

    u_data.resize(nTimeIntervals*extentLat*extentLon);
    v_data.resize(nTimeIntervals*extentLat*extentLon);
    uDataset.read(u_data.data(),H5::PredType::NATIVE_FLOAT);
    file.openDataSet("v_data").read(v_data.data(),H5::PredType::NATIVE_FLOAT);

    isInitialized = true;
    currentInterval = -1;

    spdlog::info("extentLat x extentLon: {} x {}", extentLat, extentLon);
    spdlog::info("nTimeIntervals {}", nTimeIntervals);
    spdlog::info("LatMin - LatMax: [{}; {}]", LatMin, LatMax);
    spdlog::info("LonMin - LonMax: [{}; {}]", LonMin, LonMax);
    spdlog::info("valid_time_front {}", valid_time_front);
    spdlog::info("simulation_start_date {}", simulation_start_date);
    spdlog::info("gridLatMin; gridLonMin [{},{}]",gridLatMin, gridLonMin);
    spdlog::info("WindInterpolator::ReadFromOwnHDF5 done");
}


