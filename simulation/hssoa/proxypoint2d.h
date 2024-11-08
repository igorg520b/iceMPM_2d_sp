#ifndef PROXYPOINT2D_H
#define PROXYPOINT2D_H

#include <Eigen/Core>
#include <spdlog/spdlog.h>
#include "parameters_sim.h"

struct ProxyPoint
{
    constexpr static unsigned nArrays = icy::SimParams::nPtsArrays;  // count of arrays in SOA
    bool isReference = false;
    unsigned pos, pitch;    // element # and capacity of each array in SOA
    float *soa;            // reference to SOA (assume contiguous space of size nArrays*pitch)
    float data[nArrays];    // local copy of the data when isReference==true

    ProxyPoint() { isReference = false; }

    ProxyPoint(const ProxyPoint &other);
    ProxyPoint& operator=(const ProxyPoint &other);

    // access data
    float getValue(size_t valueIdx) const;   // valueIdx < nArrays
    void setValue(size_t valueIdx, float value);
    Eigen::Vector2f getPos() const;
    Eigen::Vector2f getVelocity() const;
    bool getCrushedStatus();
    bool getDisabledStatus();
    uint16_t getGrain();
    int getCellIndex(float hinv, unsigned GridY);  // index of the grid cell at the point's location
    int getXIndex(float hinv) const;                     // x-index of the grid cell
};

#endif
