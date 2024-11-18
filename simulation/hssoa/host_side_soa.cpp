#include "host_side_soa.h"


std::pair<PointVector2r, PointVector2r> HostSideSOA::getBlockDimensions()
{
    PointVector2r result[2];
    for(int k=0;k<SimParams::dim;k++)
    {
        std::pair<SOAIterator, SOAIterator> it_res = std::minmax_element(begin(), end(),
                                                                         [k](ProxyPoint &p1, ProxyPoint &p2)
                                                                         {return p1.getValue(SimParams::posx+k)<p2.getValue(SimParams::posx+k);});
        result[0][k] = (*it_res.first).getValue(SimParams::posx+k);
        result[1][k] = (*it_res.second).getValue(SimParams::posx+k);
    }
    return {result[0], result[1]};
}

void HostSideSOA::offsetBlock(PointVector2r offset)
{
    for(SOAIterator it = begin(); it!=end(); ++it)
    {
        ProxyPoint &p = *it;
        PointVector2r pos = p.getPos() + offset;
        p.setValue(SimParams::posx, pos.x());
        p.setValue(SimParams::posx+1, pos.y());
    }
}

void HostSideSOA::convertToIntegerCellFormat(t_PointReal h)
{
    for(SOAIterator it = begin(); it!=end(); ++it)
    {
        ProxyPoint &p = *it;
        p.ConvertToIntegerCellFormat(h);
    }
}



void HostSideSOA::RemoveDisabledAndSort(int GridY)
{
    spdlog::info("RemoveDisabledAndSort; nPtsArrays {}", SimParams::nPtsArrays);
    unsigned size_before = size;
    SOAIterator it_result = std::remove_if(begin(), end(), [](ProxyPoint &p){return p.getDisabledStatus();});
    size = it_result.m_point.pos;
    spdlog::info("RemoveDisabledAndSort: {} removed; new size {}", size_before-size, size);
    std::sort(begin(), end(),
              [&](ProxyPoint &p1, ProxyPoint &p2)
              {return p1.getCellIndex(GridY)<p2.getCellIndex(GridY);});
    spdlog::info("RemoveDisabledAndSort done");
}


void HostSideSOA::Allocate(int pts_capacity, int gridTotal)
{
    spdlog::info("HostSideSOA::Allocate; pts {}; grid {}", pts_capacity, gridTotal);
    cudaFreeHost(host_buffer);
    this->capacity = pts_capacity;
    size_t allocation_size = sizeof(t_PointReal)*capacity*SimParams::nPtsArrays;
    cudaError_t err = cudaMallocHost(&host_buffer, allocation_size);
    if(err != cudaSuccess)
    {
        const char *description = cudaGetErrorString(err);
        spdlog::critical("allocating host buffer of size {}: {}",allocation_size,description);
        throw std::runtime_error("allocating host buffer for points");
    }
    size = 0;
    memset(host_buffer, 0, allocation_size);

    // grid buffer
    grid_status_buffer.resize(gridTotal);
    //grid_status_buffer
    spdlog::info("HSSOA allocate capacity {} pt; toal {} Gb", capacity, (double)allocation_size/(1024.*1024.*1024.));
}


void HostSideSOA::InitializeBlock()
{
    for(SOAIterator it = begin(); it!=end(); ++it)
    {
        ProxyPoint &p = *it;
        p.setValue(SimParams::idx_Jp_inv, 1.f);
        for(int i=0; i<SimParams::dim; i++)
                p.setValue(SimParams::Fe00+i*2+i, 1.f);
    }
    spdlog::info("HostSideSOA::InitializeBlock() done");
}




// ==================================================== SOAIterator

SOAIterator::SOAIterator(unsigned pos, t_PointReal *soa_data, unsigned pitch)
{
    m_point.isReference = true;
    m_point.pos = pos;
    m_point.soa = soa_data;
    m_point.pitch = pitch;
}

SOAIterator::SOAIterator(const SOAIterator& other)
{
    m_point.isReference = other.m_point.isReference;
    m_point.pos = other.m_point.pos;
    m_point.soa = other.m_point.soa;
    m_point.pitch = other.m_point.pitch;
}

SOAIterator& SOAIterator::operator=(const SOAIterator& other)
{
    m_point.isReference = other.m_point.isReference;
    m_point.pos = other.m_point.pos;
    m_point.soa = other.m_point.soa;
    m_point.pitch = other.m_point.pitch;
    return *this;
}

