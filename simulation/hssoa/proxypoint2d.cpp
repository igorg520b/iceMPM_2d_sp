#include "proxypoint2d.h"

// ====================================================== ProxyPoint
ProxyPoint::ProxyPoint(const ProxyPoint &other)
{
    isReference = false;
    *this = other;
}

ProxyPoint& ProxyPoint::operator=(const ProxyPoint &other)
{
    if(isReference)
    {
        // distribute into soa
        if(other.isReference)
        {
            for(int i=0;i<nArrays;i++) soa[pos + i*pitch] = other.soa[other.pos + i*other.pitch];
        }
        else
        {
            for(int i=0;i<nArrays;i++) soa[pos + i*pitch] = other.data[i];
        }
    }
    else
    {
        // local copy
        if(other.isReference)
        {
            for(int i=0;i<nArrays;i++) data[i] = other.soa[other.pos + i*other.pitch];
        }
        else
        {
            for(int i=0;i<nArrays;i++) data[i] = other.data[i];
        }
    }
    return *this;
}

Eigen::Vector2f ProxyPoint::getPos() const
{
    Eigen::Vector2f result;
    for(int i=0; i<icy::SimParams::dim;i++)
        result[i] = getValue(icy::SimParams::posx+i);
    return result;
}

Eigen::Vector2f ProxyPoint::getVelocity() const
{
    Eigen::Vector2f result;
    for(int i=0; i<icy::SimParams::dim;i++)
        result[i] = getValue(icy::SimParams::velx+i);
    return result;
}


float ProxyPoint::getValue(size_t valueIdx) const
{
    if(isReference)
        return soa[pos + pitch*valueIdx];
    else
        return data[valueIdx];
}


void ProxyPoint::setValue(size_t valueIdx, float value)
{
    if(isReference)
        soa[pos + pitch*valueIdx] = value;
    else
        data[valueIdx] = value;
}


bool ProxyPoint::getCrushedStatus()
{
    float dval = getValue(icy::SimParams::idx_utility_data);
    uint32_t val = *reinterpret_cast<uint32_t*>(&dval);
    return (val & 0x10000);
}

bool ProxyPoint::getDisabledStatus()
{
    float dval = getValue(icy::SimParams::idx_utility_data);
    uint32_t val = *reinterpret_cast<uint32_t*>(&dval);
    return (val & 0x20000) == 0x20000;
}


uint16_t ProxyPoint::getGrain()
{
    float dval = getValue(icy::SimParams::idx_utility_data);
    uint32_t val = *reinterpret_cast<uint32_t*>(&dval);
    return (val & 0xffff);
}

int ProxyPoint::getCellIndex(float hinv, unsigned GridY)
{
    Eigen::Vector2f v = getPos();
    Eigen::Vector2i idx = (v*hinv + Eigen::Vector2f::Constant(0.5)).cast<int>();
    return idx[0]*GridY + idx[1];
}

int ProxyPoint::getXIndex(float hinv) const
{
    float x = getValue(icy::SimParams::posx);
    int x_idx = (int)(x*hinv + 0.5f);
    return x_idx;
}
