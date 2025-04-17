 #ifndef WINDINTERPOLATOR_H
#define WINDINTERPOLATOR_H

#include <string_view>
#include <string>
#include <utility>
#include <vector>
#include <H5Cpp.h>
#include <Eigen/Core>

class WindInterpolator
{
public:
    WindInterpolator();









    void LoadWindData(std::string netCDF4FileName, double latMin, double latMax, double lonMin, double lonMax, int64_t sim_start_date);

    constexpr static int allocatedLatExtent = 10;
    constexpr static int allocatedLonExtent = 20;    // grid max dims must be known at compile time, because it goes into __constant__
    constexpr static double gridCellSize = 0.25; // in degrees
    constexpr static int timeResolution = 3600; // data updated every hour
    constexpr static size_t gridArraySize = allocatedLatExtent*allocatedLonExtent*4*sizeof(float);

    float grid[allocatedLatExtent][allocatedLonExtent][4];  // copy to __constant__ on the device

    int extentLat, extentLon, nTimeIntervals;   // actual extent
    double gridLatMin, gridLonMin;   // this is where the grid starts (copy to SimParams)
    int64_t valid_time_front;
    double LatMin, LatMax, LonMin, LonMax;  // region of interest

     // returns true if the contents of "grid" array changed;
    // interpolationParam is a value in [0,1) range between the two selected time frames
    bool setTime(const double simulation_time); // return true if transfer of "grid" to GPU is required
    float interpolationCoeffFromTime(const double simulation_time);
    Eigen::Vector2f interpolationResult(float lat, float lon, float tb);

    int currentInterval = -1;   // selected time interval in the "grid" array
    bool isInitialized = false;

    void SaveToOwnHDF5(H5::H5File &file);
    void ReadFromOwnHDF5(H5::H5File &file);

    // additional functionality for visualization
    void prepareVisualizationData(const double simulation_time);
    Eigen::Vector2f vis_interp(float lat, float lon);

    Eigen::Vector2f vis_coarse_grid[allocatedLatExtent][allocatedLonExtent]; // stores wind velocity grid

private:
    int64_t simulation_start_date;  // set once upon data load


    std::vector<float> u_data, v_data;  // raw data loaded from netCDF4

    size_t getIndexInUV(size_t timeIndex, size_t idxLat, size_t idxLon);



};

#endif // WINDINTERPOLATOR_H
