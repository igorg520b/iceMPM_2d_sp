#ifndef FLUENTINTERPOLATOR_H
#define FLUENTINTERPOLATOR_H

#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkLookupTable.h>
#include <vtkCellData.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkProbeFilter.h>
#include <vtkImageData.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkTransformFilter.h>
#include <vtkTransform.h>

#include <string>
#include <utility>
#include <string>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <atomic>
#include <optional>

#include <Eigen/Core>

#include "vtkFLUENTCFFCustomReader.h"
#include "parameters_sim.h"

class FluentInterpolator
{
public:
    FluentInterpolator();
    ~FluentInterpolator();
    SimParams *prms;
    bool is_initialized = false;
    double position;

    void PrepareFlowDataCache(std::string geometryFile);

    void TestLoad(double scale, double ox, double oy);

    bool SetTime(double t);

    Eigen::Vector2f getInterpolation(int i, int j) const;

    float* getFramePtr(int frame, int component);   // frame is either 0 (from) or 1 (to)

    // for setting up the scale
    vtkNew<vtkActor> actor_original;
/*
    vtkNew<vtkActor> actor;

    vtkNew<vtkDataSetMapper> mapper;
    vtkNew<vtkUnstructuredGrid> grid;
    vtkNew<vtkLookupTable> lut;
    vtkNew<vtkDataSetMapper> probeMapper;
*/
private:
    constexpr static std::string_view flow_cache_path = "_data/flow_cache";
    constexpr static int preloadedFrames = 3;   // 2 current and one extra for fast switching (circular buffer size)

    int file_count, interval_size;   // files from the scanned directory
    double currentTime;
    std::string geometryFilePrefix;
    std::string cachedFileName;
    std::vector<float> _data;

    int currentFrame = -1;
    int circularBufferIdx = -1; // between 0 and preloadedFrames-1; points to currentFrame

/*
    vtkNew<vtkFLUENTCFFCustomReader> fluentReader;
    vtkNew<vtkTransform> transform;
    vtkNew<vtkTransformFilter> transformFilter;
    vtkNew<vtkCellDataToPointData> filter_cd2pd;
    vtkNew<vtkImageData> imageData;
    vtkNew<vtkProbeFilter> probeFilter;
*/

    void ParseDataFrame(int frame, float *U, float *V);
    void LoadFrame(int frame, int slot);    // from HDF5 into a given slot of _data
    std::string CachedFileName();

    // async loading of next frame
    std::mutex loadMutex;
    std::condition_variable loadCV;
    std::atomic<bool> loadInProgress{false};
    std::optional<std::thread> loadThread;

    // for setting up the scale
    vtkNew<vtkDataSetMapper> mapper;

};

#endif // FLUENTINTERPOLATOR_H
