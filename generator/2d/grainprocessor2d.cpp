#include "grainprocessor2d.h"
#include "poisson_disk_sampling.h"

#include <vector>
#include <unordered_map>
#include <filesystem>
#include <algorithm>
#include <iostream>
#include <algorithm>

#include <H5Cpp.h>
#include <gmsh.h>
#include <spdlog/spdlog.h>

#include <Eigen/Dense>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"



void GrainProcessor2D::generate_block_and_write()
{
    load_png();
//    LoadMSH();
    print_out_parameters();
    GenerateBlock();
//    IdentifyGrains();
    Write_HDF5();
    stbi_image_free(png_data);
}




void GrainProcessor2D::print_out_parameters()
{
    spdlog::info("requestedPointsPerCell {}", requestedPointsPerCell);
    spdlog::info("grid dimensions {} x {}", gridx, gridy);
    spdlog::info("proportion of landmass {}",(float)nLandNodes/(gridx*gridy));
    spdlog::info("expected points {}", (gridx*gridy-nLandNodes)*requestedPointsPerCell);
    spdlog::info("grain scale {}", scale);
    spdlog::info("output file: {}", outputFileName);
}


void GrainProcessor2D::load_png()
{
    spdlog::info("load_png");
    const char* filename = landPNGFileName.c_str();
    png_data = stbi_load(filename, &imgx, &imgy, &channels, 3);  // Request 1 channel (grayscale)

    if (!png_data)
    {
        std::cerr << "Failed to load image: " << filename << std::endl;
        throw std::runtime_error("png not loaded");
    }

    scale_img = (float)imgx / gridx;

    gridy = (int)(imgy/scale_img);
    grid_buffer.resize(gridx*gridy);
    spdlog::info("image ({} x {}) loaded; {} channels", imgx, imgy, channels);
    spdlog::info("grid size {} x {}", gridx, gridy);
    spdlog::info("grid_buffer size {}",grid_buffer.size());

    nLandNodes = 0;
    nWaterNodes = 0;
    for (int y = 0; y < gridy; ++y)
    {
        for (int x = 0; x < gridx; ++x)
        {
            int px = (int)(x * scale_img);
            int py = (int)(y * scale_img);
            px = std::min(px, imgx-1);
            py = std::min(py, imgy-1);

            int index_png = py * imgx + px;
            unsigned char &r = png_data[index_png*3 + 0];   // red channel
            unsigned char &g = png_data[index_png*3 + 1];   // red channel
            unsigned char &b = png_data[index_png*3 + 2];   // red channel
            bool is_land = (r > 128 && g < 128 & b < 128);
            if(!is_land) nWaterNodes++;
            else nLandNodes++;

            grid_buffer[(gridy-y-1) + x*gridy] = is_land ? 1 : 0;
        }
    }
}




void GrainProcessor2D::GenerateBlock()
{
    constexpr float magic_constant = 0.65;

    float h_img = 1.f/(float(imgx-1));

    float h = 1.f/(float)(gridx-1);
    float dx = 1.f;
    float dy = 1.f*(gridy-1)/(float)(gridx-1);
    float nPts = requestedPointsPerCell*gridx*gridy;
    const float kRadius = sqrt(magic_constant*dx*dy/nPts);

    const std::array<float, 2>kXMin {2*h, 2*h};
    const std::array<float, 2>kXMax {dx-2*h, dy-2*h};

    int pt_counters[5]{};
    {
    buffer = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax);
    spdlog::info("generating rectangle [{},{}] - [{},{}]",kXMin[0],kXMin[1],kXMax[0],kXMax[1]);


    // copy to buffer_categorized and categorize
    buffer_categorized.resize(buffer.size());

    for(int idx=0;idx<buffer.size();idx++)
    {
        float x = buffer[idx][0];
        float y = buffer[idx][1];

        int i = (int)(x/h + 0.5f);
        int j = (int)(y/h + 0.5f);
        if(i<0 || j< 0 || i>=gridx || j>=gridy)
        {
            spdlog::critical("error pt ({},{}); grid cell ({},{})",x,y,i,j);
            throw std::runtime_error("particle is out of grid bounds");
        }

        int category = -2;

        // if on land cell of the grid
        unsigned char pixel_land = grid_buffer[j + i*gridy];
        if(pixel_land) category = -1;
        else
        {
            // pick from png image
            int px = (int)(x/h_img + 0.5f);
            int py = (int)(y/h_img + 0.5f);
            int index_png = (imgy - py -1) * imgx + px;
            float r = (float)png_data[index_png*3 + 0]/255.;
            float g = (float)png_data[index_png*3 + 1]/255.;
            float b = (float)png_data[index_png*3 + 2]/255.;
            const float color[3] {r,g,b};
            category = categorizeColor(color);
        }
        pt_counters[category+1]++;

        buffer_categorized[idx] = {x,y,category};
    }
    }

    size_t nPtsCountInitial = buffer_categorized.size();
    spdlog::info("initial point cloud: {} pts", nPtsCountInitial);
    spdlog::info("pt_counters; land {} water {} crushed {} weakened {} solid {}",
                     pt_counters[0],pt_counters[1],pt_counters[2],pt_counters[3],pt_counters[4]);

    auto result = std::remove_if(buffer_categorized.begin(),buffer_categorized.end(),
                                 [&](std::tuple<float,float,int> &pt){
//        float x = std::get<0>(pt);
//        float y = std::get<1>(pt);
        int cat = std::get<2>(pt);
        return (cat == 0 || cat == -1); // remove points that are open water or land
    });
    buffer_categorized.erase(result,buffer_categorized.end());

    const size_t nPtsCountFinal = buffer_categorized.size();
    spdlog::info("after erasing land points: {} pts", nPtsCountFinal);


    volume = (kXMax[0]-kXMin[0])*(kXMax[1]-kXMin[1])*block_length*block_length;
    volume *= (float)nPtsCountFinal/(float)nPtsCountInitial;


    auto getidx = [&] (std::tuple<float,float,int> &pt)
    {
        float x = std::get<0>(pt);
        float y = std::get<1>(pt);

        int i = (int)(x/h + 0.5f);
        int j = (int)(y/h + 0.5f);
        int cellidx = j + gridy*i;
        return cellidx;
    };

    // sort, to facilitate subsequent loading
    std::sort(buffer_categorized.begin(), buffer_categorized.end(),
              [&](std::tuple<float,float,int> &p1, std::tuple<float,float,int> &p2)
              {return getidx(p1)<getidx(p2);});
    spdlog::info("points sortred by cell");

    coordinates[0].resize(buffer_categorized.size());
    coordinates[1].resize(buffer_categorized.size());
    llGrainID.resize(buffer_categorized.size());

    for(int i=0;i<buffer_categorized.size();i++)
    {
        std::tuple<float,float,int> &pt = buffer_categorized[i];
        float x = std::get<0>(pt);
        float y = std::get<1>(pt);

        coordinates[0][i] = x*block_length;
        coordinates[1][i] = y*block_length;

        int category = std::get<2>(pt);
        constexpr uint32_t status_crushed = 0x10000;
        constexpr uint32_t status_weakened = 0x40000;
        if(category == 3) llGrainID[i] = 0;
//        else if(category == 2) llGrainID[i] = status_weakened;
        else if(category == 2) llGrainID[i] = status_crushed;
        else if(category == 1) llGrainID[i] = status_crushed;

    }
    cellsize = block_length/(float)(gridx-1);
}


void GrainProcessor2D::Write_HDF5()
{
    spdlog::info("Write_HDF5()");
    H5::H5File file(outputFileName, H5F_ACC_TRUNC);

    hsize_t dims_grains[1] = {llGrainID.size()};
    H5::DataSpace dataspace_points_grains(1, dims_grains);
    hsize_t g_chunk_dims[1] = {1024*256};
    if(g_chunk_dims[0] > llGrainID.size()) g_chunk_dims[0] = std::max(llGrainID.size()/3,(size_t)1);
    spdlog::info("llGrain chunk {}; dims_grains {}", g_chunk_dims[0], dims_grains[0]);
    H5::DSetCreatPropList proplist1;
    proplist1.setChunk(1, g_chunk_dims);
    proplist1.setDeflate(7);
    H5::DataSet dataset_grainids = file.createDataSet("llGrainIDs", H5::PredType::NATIVE_UINT32, dataspace_points_grains, proplist1);
    dataset_grainids.write(llGrainID.data(), H5::PredType::NATIVE_UINT32);

    hsize_t att_dim = 1;
    H5::DataSpace att_dspace(1, &att_dim);
    H5::Attribute att_volume = dataset_grainids.createAttribute("volume", H5::PredType::NATIVE_FLOAT, att_dspace);
    att_volume.write(H5::PredType::NATIVE_FLOAT, &volume);
    spdlog::info("volume written {}",volume);


    H5::DataSet dataset_x = file.createDataSet("x", H5::PredType::NATIVE_FLOAT, dataspace_points_grains, proplist1);
    dataset_x.write(coordinates[0].data(), H5::PredType::NATIVE_FLOAT);

    H5::DataSet dataset_y = file.createDataSet("y", H5::PredType::NATIVE_FLOAT, dataspace_points_grains, proplist1);
    dataset_y.write(coordinates[1].data(), H5::PredType::NATIVE_FLOAT);



    // write grid info
    hsize_t dims_grid[1] = {grid_buffer.size()};

    H5::DataSpace dataspace_grid(1, dims_grid);
    hsize_t g_chunk_grid[1] = {1024*128};
    if(g_chunk_grid[0] > grid_buffer.size()) g_chunk_grid[0] = std::max(grid_buffer.size()/3,(size_t)1);
    H5::DSetCreatPropList proplist2;

    proplist2.setChunk(1, g_chunk_grid);
    proplist2.setDeflate(7);
    H5::DataSet dataset_grid = file.createDataSet("GridLand", H5::PredType::NATIVE_UINT8, dataspace_grid, proplist2);
    dataset_grid.write(grid_buffer.data(), H5::PredType::NATIVE_UINT8);

    dataset_grid.createAttribute("GridX", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &gridx);
    dataset_grid.createAttribute("GridY", H5::PredType::NATIVE_INT, att_dspace).write(H5::PredType::NATIVE_INT, &gridy);
    dataset_grid.createAttribute("cellsize", H5::PredType::NATIVE_FLOAT, att_dspace).write(H5::PredType::NATIVE_FLOAT, &cellsize);

    file.close();
}

bool GrainProcessor2D::PointInsideTriangle(Eigen::Vector2f point, Eigen::Vector2f t[3])
{
    constexpr float eps = 1e-3;
    Eigen::Matrix2f M,B;
    M << t[1]-t[0], t[2]-t[0];
    B = M.inverse();
    point -= t[0];
    Eigen::Vector2f b = B * point; // barycentric
    return (b[0]>-eps && b[1]>-eps && (b.sum() < 1+eps));
}


void GrainProcessor2D::IdentifyGrains()
{
    spdlog::info("identify grains");

    BVHN2D::BVHNFactory.release(leaves);
    leaves.reserve(tris2.size());

    for(int i=0;i<tris2.size();i++)
    {
        Triangle &t = tris2[i];
        BVHN2D* bvhn = BVHN2D::BVHNFactory.take();
        bvhn->isLeaf = true;
        bvhn->box.Reset();
        for(int j=0;j<3;j++) bvhn->Expand(t.nds[j]);
        bvhn->child1 = bvhn->child2 = nullptr;
        bvhn->elem = i;
        leaves.push_back(bvhn);
    }

    spdlog::info("building bvh");
    root.Build(leaves, 0);
    spdlog::info("finished building bvh");

    grainID.resize(buffer.size());
    spdlog::info("tetra2 {}; grainID {}", tris2.size(), grainID.size());

    // identify grains
    int unidentifiedPoints = 0;
    int problematicPoints = 0;

#pragma omp parallel for reduction(+:unidentifiedPoints,problematicPoints)
    for(int i=0;i<buffer.size();i++)
    {
        std::array<float, 2> &arr = buffer[i];
        Eigen::Vector2f v(arr[0],arr[1]);
        v *= (scale);
        v -= v.array().floor().matrix();

        BVHN2D bvhn;
        bvhn.isLeaf = true;
        bvhn.elem = -1;
        bvhn.box.Reset();
        bvhn.Expand(v);

        std::vector<std::pair<BVHN2D*,BVHN2D*>> broad_list;
        broad_list.reserve(10);
        bvhn.Collide(&root, broad_list);

        grainID[i] = -1;
        for(auto &pair : broad_list)
        {
            BVHN2D *bvhn2 = pair.second;
            int idx = bvhn2->elem;
            if(idx == -1) spdlog::critical("elem index -1");
            Triangle &t = tris2[idx];
            bool result = PointInsideTriangle(v, t.nds);
            if(result)
            {
                grainID[i] = (short)t.grain;
                break;
            }
        }

        if(grainID[i] == -1)
        {
            problematicPoints++;
            spdlog::warn("node {}; grain {}; broad_list {}",i,grainID[i],broad_list.size());
            if(broad_list.size()!=0)
            {
                auto &p = broad_list.front();
                short grain = (short)tris2[p.second->elem].grain;
                grainID[i] = grain;
                if(grain < 0) spdlog::critical("elem {} grain {}",p.second->elem,grain);
            }
            else
            {
                spdlog::critical("broad list 0");
            }
        }
        if(grainID[i] == -1) { grainID[i] = 0; unidentifiedPoints++; }
    }

    llGrainID.resize(grainID.size());
    for(int i=0;i<grainID.size();i++) llGrainID[i] = grainID[i];

    spdlog::info("finished processing points; problematic {}; unidentified {}",problematicPoints, unidentifiedPoints);
}


void GrainProcessor2D::LoadMSH()
{
    std::vector<Eigen::Vector2f> vertices;
    std::vector<std::array<int,4>> elems;
    std::vector<Triangle> tris1;

    spdlog::info("load {}", meshFileName);
    if(!std::filesystem::exists(meshFileName)) spdlog::critical("file does not exist");

    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::open(meshFileName);

    // get nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    std::unordered_map<std::size_t, int> mtags; // gmsh nodeTag -> node sequential number

    gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);
    vertices.reserve(nodeTags.size()*9);

    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in deformable");
        Eigen::Vector2f coords = Eigen::Vector2f(nodeCoords[i*3+0], nodeCoords[i*3+1]);
        mtags[tag] = vertices.size();
        vertices.push_back(coords);
    }
    spdlog::info("nodes read {}",vertices.size());
    int nVerticesOriginal = vertices.size();
    vertices2.resize(nVerticesOriginal*9);

    // create surrounding copies
    int count = 0;
    for(int i=-1;i<=1;i++)
        for(int j=-1;j<=1;j++)
        {
            for(int m=0;m<nVerticesOriginal;m++)
                vertices2[m + nVerticesOriginal*count] = vertices[m] + Eigen::Vector2f(i,j);
            count++;
        }

    std::vector<std::pair<int,int>> dimTagsGrains;
    gmsh::model::getEntities(dimTagsGrains,2);
    spdlog::info("dimTagsGrains size {}",dimTagsGrains.size());

    //    std::unordered_set<int> used_nodes;
    elems.reserve(dimTagsGrains.size()*50);
    tris1.clear();
    tris1.reserve(elems.size());

    std::vector<int> types;
    std::vector<std::vector<std::size_t>> elemtags, nodetags;
    gmsh::model::mesh::getElements(types, elemtags, nodetags);

    for(int i=0;i<types.size();i++)
        spdlog::info("type {}; elemtags {}", types[i], elemtags[i].size());


    for(std::size_t j=0;j<dimTagsGrains.size();j++)
    {
        std::vector<std::size_t> tetraTags, nodeTagsInTetra;
        int entityTag = dimTagsGrains[j].second;
        gmsh::model::mesh::getElementsByType(2, tetraTags, nodeTagsInTetra, entityTag);
        //        spdlog::info("grain {}; tetratags {}", j, tetraTags.size());

        for(std::size_t i=0;i<tetraTags.size();i++)
        {
            std::array<int,4> elem;
            Triangle t;

            for(int k=0;k<3;k++)
            {
                elem[k] = mtags.at(nodeTagsInTetra[i*3+k]);
                t.nds[k] = vertices[elem[k]];
            }
            t.grain = elem[3] = j;
            elems.push_back(elem);
            tris1.push_back(t);
        }
    }
    spdlog::info("tris read {}", tris1.size());

    int nTetraOriginal = tris1.size();
    elems2.resize(nTetraOriginal*9);
    tris2.resize(nTetraOriginal*9);

    count = 0;
    for(int i=-1;i<=1;i++)
        for(int j=-1;j<=1;j++)
        {
            for(int m=0;m<tris1.size();m++)
            {
                int idx = m + count*nTetraOriginal;
                elems2[idx] = elems[m];
                Triangle &t = tris2[idx];
                for(int n=0;n<3;n++)
                {
                    elems2[idx][n] += count*nVerticesOriginal;
                    t.nds[n] = vertices2[elems2[idx][n]];
                }
                t.grain = elems2[m][3];
            }
            count++;
        }

    gmsh::clear();
}




float GrainProcessor2D::pointToSegmentDistance(const Eigen::Vector3f& p, const Eigen::Vector3f& a, const Eigen::Vector3f& b)
{
    Eigen::Vector3f ab = b - a;
    Eigen::Vector3f ap = p - a;
    float t = ap.dot(ab) / ab.squaredNorm(); // Projection factor
    t = std::clamp(t, 0.0f, 1.0f);          // Clamp t to [0, 1]
    Eigen::Vector3f projection = a + t * ab;       // Projection point
    return (p - projection).squaredNorm();  // Squared distance
}


int GrainProcessor2D::categorizeColor(const float color[3])
{
    Eigen::Vector3f colorVec(color[0], color[1], color[2]);
    float minDistance = std::numeric_limits<float>::max();
    int category = -1;

    for (int i = 0; i < 7; ++i)
    { // Iterate through each segment
        Eigen::Vector3f a(colordata[i][0], colordata[i][1], colordata[i][2]);
        Eigen::Vector3f b(colordata[i + 1][0], colordata[i + 1][1], colordata[i + 1][2]);
        float distance = pointToSegmentDistance(colorVec, a, b);
        if (distance < minDistance) {
            minDistance = distance;
            category = i / 2; // Each category has 2 segments
        }
    }
    return category;
}
