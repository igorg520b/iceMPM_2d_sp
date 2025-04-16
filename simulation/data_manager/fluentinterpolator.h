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

#include "parameters_sim.h"
#include "gpu_implementation5.h"

class FluentInterpolator
{
public:
    FluentInterpolator();
    ~FluentInterpolator();
    SimParams *prms;
    GPU_Implementation5 *gpu;
    bool is_initialized = false;
    double position;

    void PrepareFlowDataCache(std::string geometryFile);
    void PrepareFromCachedFile(std::string cachedFile);

    void TestLoad(double scale, double ox, double oy);

    bool SetTime(double t);

    Eigen::Vector2f getInterpolation(int i, int j) const;

    float* getFramePtr(int frame, int component);   // frame is either 0 (from) or 1 (to)


    std::string cachedFileName;

    // for setting up the scale
    vtkNew<vtkActor> actor_original;
private:
    constexpr static std::string_view flow_cache_path = "_data/flow_cache";
    constexpr static int preloadedFrames = 3;   // 2 current and one extra for fast switching (circular buffer size)

    int file_count, interval_size;   // files from the scanned directory
    double currentTime;
    std::string geometryFilePrefix;
    std::vector<float> _data;

    int currentFrame = -1;
    int circularBufferIdx = -1; // between 0 and preloadedFrames-1; points to currentFrame

    void ParseDataFrame(int frame, float *U, float *V);
    void LoadFrame(int frame, int slot);    // from HDF5 into a given slot of _data
    std::string CachedFileName();
    void applyDiffusion(float* V, int gx, int gy, float D, float dt, int steps);

    // async loading of next frame
    std::mutex loadMutex;
    std::condition_variable loadCV;
    std::atomic<bool> loadInProgress{false};
    std::optional<std::thread> loadThread;

    // for setting up the scale
    vtkNew<vtkDataSetMapper> mapper;

};

#endif // FLUENTINTERPOLATOR_H
