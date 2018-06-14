// materials/flake.cpp*
#include "materials/flake.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"
#include "interaction.h"
#include <random>
#include <algorithm>
#define PI 3.14159265
namespace pbrt {

// FlakeMaterial Method Definitions
void FlakeMaterial::ComputeScatteringFunctions(
    SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
    bool allowMultipleLobes) const {
    int N = pow(10,7); // the total number of flakes
    int T = 16; // number of flakes in leaf
    int K = 10; // traversal depth
    Point2f uv = si->uv;
    double dxy = 0.01;  //epsilon

//     Calculate the four nodes of parallelogram A
    Point2f footprint_upperleft = Point2f(uv.x-si->dudx*dxy-si->dudy*dxy,uv.y-si->dvdy*dxy-si->dvdx*dxy);
    Point2f footprint_upperright = Point2f(uv.x+si->dudx*dxy-si->dudy*dxy,uv.y-si->dvdy*dxy+si->dvdx*dxy);
    Point2f footprint_lowerleft = Point2f(uv.x-si->dudx*dxy+si->dudy*dxy,uv.y+si->dvdy*dxy-si->dvdx*dxy);
    Point2f footprint_lowerright = Point2f(uv.x+si->dudx*dxy+si->dudy*dxy,uv.y+si->dvdy*dxy+si->dvdx*dxy);
    Float area = abs((footprint_lowerright.x-footprint_lowerleft.x)*(footprint_lowerleft.y-footprint_upperleft.y));

// Find the bounding box of A
    float x_min = std::min({footprint_upperleft.x,footprint_lowerleft.x, footprint_upperright.x,footprint_lowerright.x});
    float y_min = std::min({footprint_upperleft.y,footprint_upperright.y,footprint_lowerright.y,footprint_lowerleft.y});
    float x_max = std::max({footprint_upperleft.x,footprint_lowerleft.x, footprint_upperright.x,footprint_lowerright.x});
    float y_max = std::max({footprint_upperleft.y,footprint_upperright.y,footprint_lowerright.y,footprint_lowerleft.y});

    Float bounding_box_area = std::abs((y_max-y_min)*(x_max-x_min));
    Float area_ratio = area/bounding_box_area;
    // Traverse to the leaves
    Float leaf_length = 1/(pow(2,K));
    // calculate the index that the boundingbox overlapping
    Point2i upperleft_index = Point2i(int(x_min/leaf_length),int(y_min/leaf_length));
    Point2i lowerright_index = Point2i(int(x_max/leaf_length),int(y_max/leaf_length));
    
    int num_boxes = (lowerright_index.y-upperleft_index.y+1)*(lowerright_index.x-upperleft_index.x+1);
    int number_in_bounding_grid = SampleGaussian(K,num_boxes,N);
   
    int numberFlakes = int(abs(area_ratio*number_in_bounding_grid));

    RGBSpectrum kd;

    // Perform bump mapping with _bumpMap_, if present
    if (bumpMap) Bump(bumpMap, si);
    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);

    // Clamp the number of flakes to [0,300]
    int Maxflakes = 300;
    float flakecolor;
    if (numberFlakes>Maxflakes)
        numberFlakes = Maxflakes;
    else if (numberFlakes<=0){
        numberFlakes = 0;

    }

    kd = Kd->Evaluate(*si).Clamp();
    if (!kd.IsBlack()){
        si->bsdf->Add(ARENA_ALLOC(arena, LambertianReflection)(kd));
    }

    // Initialize specular component of flake material
    Spectrum ks = Ks->Evaluate(*si).Clamp();
    if (!ks.IsBlack()) {
    // Here use mirror's fresnel
        Fresnel *fresnel = ARENA_ALLOC(arena, FresnelNoOp)();
        // Create microfacet distribution _distrib_ for plastic material
        Float rough = roughness->Evaluate(*si);
        if (remapRoughness)
            rough = TrowbridgeReitzDistribution::RoughnessToAlpha(rough);
        MicrofacetDistribution *distrib =
            ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(rough, rough);
        // Pass to constructor of flakeBxDF
            BxDF *spec =
                ARENA_ALLOC(arena, flakeBxDF)(ks, fresnel, distrib, numberFlakes, gamma, area);
        si->bsdf->Add(spec);
    }
}


FlakeMaterial *CreateFlakeMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> Kd =
        mp.GetSpectrumTexture("Kd", Spectrum(0.25f));
    std::shared_ptr<Texture<Spectrum>> Ks =
        mp.GetSpectrumTexture("Ks", Spectrum(0.25f));
    std::shared_ptr<Texture<Float>> roughness =
        mp.GetFloatTexture("roughness", .1f);
    std::shared_ptr<Texture<Float>> bumpMap =
        mp.GetFloatTextureOrNull("bumpmap");
    bool remapRoughness = mp.FindBool("remaproughness", true);
    float gamma = mp.FindFloat("gamma",5.0f);
    gamma = gamma * PI / 180.0;
    return new FlakeMaterial(Kd, Ks, roughness, bumpMap, remapRoughness, gamma);
}

}  // namespace pbrt

// sample from a multivariate gaussian distribution with p = (1/4,1/4,1/4,1/4)
// Here a normal distribution is used to approximate
int SampleGaussian(int depth, int number_boxes, int N){
    double mu = 0.25;
    double sd = 0.1;
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(mu,sd);
    double total_number = 0;
    double number_in_grid = N;
    for (int j=1;j<=number_boxes;j++){
        for(int i=1;i<depth;i++){
            double x_grid = distribution(generator);
        if (x_grid>0){
            number_in_grid = number_in_grid*x_grid;
        }
        else{
            number_in_grid = number_in_grid*0.25;
        }
        // std::cout<<"number"<<std::endl;
        }
        total_number += number_in_grid;
    }
    // std::cout<<number_in_grid*number_boxes<<std::endl;
    return int(total_number);
    // return 16*number_boxes;
}
