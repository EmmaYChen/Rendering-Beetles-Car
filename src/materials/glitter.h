
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_flakeBXDF_H
#define PBRT_CORE_flakeBXDF_H

// core/reflection.h*
#include "pbrt.h"
#include "geometry.h"
#include "microfacet.h"
#include "shape.h"
#include "spectrum.h"
#include "core/reflection.h"

namespace pbrt {

class flakeBxDF : public BxDF {
  public:
    // MicrofacetReflection Public Methods
    flakeBxDF(const Spectrum &R, Fresnel *fresnel, MicrofacetDistribution *distribution, int flakenumber, float gamma, float area_A)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
          R(R),
          fresnel(fresnel),
          distribution(distribution),
          flakenumber(flakenumber),
          gamma(gamma),
          area_A(area_A){
      }

    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    // MicrofacetReflection Private Data
    const Spectrum R;
    const Fresnel *fresnel;
    const MicrofacetDistribution *distribution;
    // gamma is the cone radius
    float gamma;
    int flakenumber;
    float area_A;
};

}  // namespace pbrt

#endif  // PBRT_CORE_flakeBXDN_H
