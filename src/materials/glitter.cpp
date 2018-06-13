// core/reflection.cpp*
#include "spectrum.h"
#include "sampler.h"
#include "sampling.h"
#include "interpolation.h"
#include "scene.h"
#include "interaction.h"
#include "stats.h"
#include <stdarg.h>
#include <ctime>
#include "glitter.h"
#include "reflection.h"
#include "texture.h"
#include "paramset.h"
#include "math.h"
#include <vector>
#include "materials/sampleHSL.h"

#define PI 3.14159265
using namespace std;
namespace pbrt {

    float RandomFloat(float a, float b){
        // srand(time(0));
        float random = ((float) rand()) / (float) RAND_MAX;
        float diff = b - a;
        float r = random * diff;
        return a + r;
    }

    Vector3f RandomVector(){
        float phi = RandomFloat(0, PI * 2);
        float costheta = RandomFloat(-1, 1);
        float theta = acos(costheta);
        float x = sin( theta ) * cos( phi );
        float y = sin( theta ) * sin( phi );
        float z = cos( theta );
        return Vector3f(x,y,z);
    }

    void CreateReflectionCache(const Vector3f wo, vector<Vector3f> &normals,
                               vector<Vector3f> &rcache, vector<float> &weights, int n, vector<float> &ctable ){

        float c = 0;
        ctable.push_back(c);

        for(int i =0; i<n; i++){
            Vector3f m = Normalize(RandomVector());
            Vector3f r = Reflect(wo,m);
            if(!SameHemisphere(wo,r)) r *= -1;
            float w = std::abs(Dot(m,r));
            c += w;
            ctable.push_back(c);
            normals.push_back(m);
            rcache.push_back(r);
            weights.push_back(w);
        }
        //if(ctable.size()!= n+1) std::cout<<"something wrong" <<std::endl;
    }

    vector<double> CalculateDistribution(const Vector3f wi, vector<float> weights, vector<Vector3f> reflectioncache, int n, float gamma){
        double D = 0;
        int N = pow(10,8);
        for(int i = 0; i < n; i++){
            if ( Dot(wi,reflectioncache[i]) >= cos(gamma)){
                D += weights[i];
            }
        }
        vector<double> DK;
        DK.push_back(D/double(N));
        DK.push_back(D);
        return DK;
    }

    Vector3f BinarySearch(vector<Vector3f> cache, vector<float> ctable, double value, int n){
        int start = 0;
        int end = n;
        while(start <= end) {
            int middle = start + (end - start)/2;
            if (ctable[start] < value && ctable[start + 1] >= value) {
                return cache[start];
            }
            else if (ctable[middle] >= value) {
                end = middle;
            }
            else start = middle;
        }
        return Vector3f(0,0,0);
    }

    Vector3f GenerateBRDFi(Vector3f r, float gamma){
        Vector3f i;
        while(true){
            i = RandomVector();
            i = Normalize(i);
            if (Dot(i,r) >= cos(gamma)){
                return i;
            }
        }
        return Vector3f(0,0,0);
    }

    std::string flakeBxDF::ToString() const {
        return std::string("[ flakeBxDF R: ") + R.ToString() +
        std::string(" fresnel: ") + fresnel->ToString() + std::string(" ]");
    }

    // BxDF Method Definitions
    Spectrum flakeBxDF::f(const Vector3f &wo, const Vector3f &wi) const {

        Float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);
        if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0.);
        Vector3f wh = wi + wo;
        if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0.);
        wh = Normalize(wh);
        if(!flakenumber) return Spectrum(0.);

        int bmin = 1000, bmax = 1000;
        int E_x = area_A * pow(10,8);
        // cout <<"E_x" << E_x 《《

        vector<Vector3f> normals, reflectioncache;
        vector<float> weights, ctable;
        CreateReflectionCache(wo, normals, reflectioncache, weights, flakenumber, ctable);
        vector<double> DK = CalculateDistribution(wi, weights, reflectioncache, flakenumber, gamma);
        //if (DK[0] != DK[0]) return Spectrum(10000.);
        // check
        if(DK[0]!=DK[0]||!DK[0])return  Spectrum(0.);
        //bool flag = (DK[0] || DK[0] == DK[0]) && flakenumber;
        bool flag = (DK[0] == DK[0]);

        // discrete
        if((flakenumber < bmin) && flag ){
           // Handle degenerate cases for microfacet reflection
            wh = Normalize(wh);
            Spectrum F = fresnel->Evaluate(Dot(wi, wh));

            // color
            Point3f sampleColor = SampleFromHSL();
            float color[] = {sampleColor.x,sampleColor.y,sampleColor.z};
            F *= Spectrum::FromRGB(color);

            Float conearea = 2.0 * PI * (1.0 - cos(gamma));
            // Float ih = Dot(wh,wi);
            Float distri = DK[0] * (4.0 / (abs(area_A) * conearea));

            return R * F * distri / (4.0 * cosThetaI * cosThetaO * 200.0);
        }

        // continuous
        else if((flakenumber >= bmax) || !flag){

            Spectrum F = fresnel->Evaluate(Dot(wi, wh));

            //color
            Point3f sampleColor = AverageColor();
            float color[] = {sampleColor.x,sampleColor.y,sampleColor.z};
            F *= RGBSpectrum::FromRGB(color);
           // return Spectrum(0.);
            return R * distribution->D(wh) * distribution->G(wo, wi) * F /
               (4 * cosThetaI * cosThetaO);
       }
       return Spectrum(0.);

    }

    Spectrum flakeBxDF::Sample_f(const Vector3f &wo, Vector3f *wi,
                                 const Point2f &u, Float *pdf,
                                 BxDFType *sampledType) const {
        int bmin = 1000, bmax = 1000;
        if (wo.z == 0) return Spectrum(0.);
        if (flakenumber <=0 ) return Spectrum(0.);

        vector<Vector3f> normals, reflectioncache;
        vector<float> weights,ctable;
        CreateReflectionCache(wo, normals, reflectioncache, weights, flakenumber,ctable);

        // bool flag = flakenumber >0 ;
        int E_x = area_A*8*pow(10,10);

        // discrete
       if(flakenumber < bmin){

            float randomvalue = RandomFloat(0,ctable[flakenumber]);
            Vector3f r = BinarySearch(reflectioncache, ctable, randomvalue, flakenumber);
            *wi = GenerateBRDFi(r,gamma);
            if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);
            if(pdf) *pdf = Pdf(wo, *wi);
            return f(wo, *wi);
       }

        // continuous
        //  if(flakenumber >= bmax)
        else if(flakenumber >= bmax){

            Vector3f wh = distribution->Sample_wh(wo, u);
            *wi = Reflect(wo, wh);
            if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);
            if(pdf) *pdf = distribution->Pdf(wo, wh) / (4.0 * Dot(wo, wh));
            return f(wo, *wi);
        }

        return Spectrum(0.);
    }

    Float flakeBxDF:: Pdf(const Vector3f &wo, const Vector3f &wi) const {
        int bmin = 1000, bmax = 1000;
        if (!SameHemisphere(wo, wi)) return 0;
        if (flakenumber <= 0) return 0;
        vector<Vector3f> normals, reflectioncache;
        vector<float> weights, ctable;
        CreateReflectionCache(wo, normals, reflectioncache, weights, flakenumber, ctable);
        vector<double> DK = CalculateDistribution(wi, weights, reflectioncache, flakenumber, gamma);
        Vector3f wh = Normalize(wo + wi);

        //experiment
        if (DK[1] != DK[1]) return 0.0;
        //check
        // if(DK[1]!=DK[1]||!DK[1])return  0.;
        bool flag = (DK[1] && DK[1] == DK[1]) && flakenumber;
        int E_x = area_A * 8 * pow(10,10);

        // discrete
        if((flakenumber < bmin) && flag){
            return  DK[1] / (PI * (1.0 - cos(gamma)) * ctable[flakenumber]);
        }

        else if(flakenumber >= bmax || !flag){
            return distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));
        }

        return 0.0;
    }
}  // namespace pbrt
