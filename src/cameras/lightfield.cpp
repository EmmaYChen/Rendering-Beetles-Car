
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// cameras/perspective.cpp*
#include "cameras/lightfield.h" 
#include "paramset.h"
#include "sampler.h"
#include "sampling.h"
#include "light.h"
#include "stats.h"

namespace pbrt {

// PerspectiveCamera Method Definitions
LightFieldCamera::LightFieldCamera(const AnimatedTransform &CameraToWorld,
                                     const Bounds2f &screenWindow,
                                     Float shutterOpen, Float shutterClose,
                                     Float lensRadius, Float focalDistance,
                                     Float fov, Film *film,
                                     const Medium *medium,
                                     Vector2i camerasperdim, Bounds2f cameragridbounds)
    : ProjectiveCamera(CameraToWorld, Perspective(fov, 1e-2f, 1000.f),
                       screenWindow, shutterOpen, shutterClose, lensRadius,
                       focalDistance, film, medium), camerasperdim(camerasperdim),
                       cameragridbounds(cameragridbounds), CameraToScreen(Perspective(fov, 1e-2f, 1000.f)){
    // Compute differential changes in origin for perspective camera rays
    
    myScreenToRaster =
            Scale(film->fullResolution.x/camerasperdim.x, film->fullResolution.y/camerasperdim.y, 1) *
            Scale(1/(screenWindow.pMax.x- screenWindow.pMin.x),
                  1/(screenWindow.pMin.y - screenWindow.pMax.y), 1) *
            Translate(Vector3f(-screenWindow.pMin.x, -screenWindow.pMax.y, 0));
        myRasterToScreen = Inverse(myScreenToRaster);
        myRasterToCamera = Inverse(CameraToScreen) * myRasterToScreen;
    
    dxCamera =
        (myRasterToCamera(Point3f(1, 0, 0)) - myRasterToCamera(Point3f(0, 0, 0)));
    dyCamera =
        (myRasterToCamera(Point3f(0, 1, 0)) - myRasterToCamera(Point3f(0, 0, 0)));


    Point2i res = film->fullResolution;
    Point3f pMin = myRasterToCamera(Point3f(0, 0, 0));
    Point3f pMax = myRasterToCamera(Point3f(res.x, res.y, 0));
    pMin /= pMin.z;
    pMax /= pMax.z;
    A = std::abs((pMax.x - pMin.x) * (pMax.y - pMin.y));

    film_res_x = film->fullResolution.x/camerasperdim.x ;
    film_res_y = film->fullResolution.y/camerasperdim.y ;

    //Point3f CameraGrid[camerasperdim.x][camerasperdim.y];
    float delta_x = (cameragridbounds.pMax.x - cameragridbounds.pMin.x)/ (camerasperdim.x -1);
    float delta_y = (cameragridbounds.pMax.y - cameragridbounds.pMin.y)/ (camerasperdim.y -1);

    for (int i=0; i< camerasperdim.y; i++)
        for(int j=0; j< camerasperdim.x; j++){
            float x = cameragridbounds.pMin.x + i * delta_x;
            float y = cameragridbounds.pMax.y - j * delta_y;
            Point3f origin = Point3f(x,y,0.0);
            CameraGrid[i][j] = origin;
    }
}

Float LightFieldCamera::GenerateRay(const CameraSample &sample,
                                     Ray *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);
    // Compute raster and camera sample positions
    int index_x = sample.pFilm.x/film_res_x;
    int index_y = sample.pFilm.y/film_res_y;
    Point3f pFilm = Point3f(sample.pFilm.x - index_x * film_res_x, sample.pFilm.y - index_y * film_res_y, 0);
    Point3f pCamera = myRasterToCamera(pFilm);
    Point3f newOrigin = CameraGrid[index_x][index_y];
    *ray = Ray(newOrigin, Normalize(Vector3f(pCamera)));
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    ray->medium = medium;
    *ray = CameraToWorld(*ray);
    return 1;
}

Float LightFieldCamera::GenerateRayDifferential(const CameraSample &sample,
                                                 RayDifferential *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);
    // Compute raster and camera sample positions
    int index_x = sample.pFilm.x/film_res_x;
    int index_y = sample.pFilm.y/film_res_y;
    Point3f pFilm = Point3f(sample.pFilm.x - index_x * film_res_x, sample.pFilm.y - index_y * film_res_y, 0);
    Point3f pCamera = myRasterToCamera(pFilm);
    Point3f newOrigin = CameraGrid[index_x][index_y];
    *ray = Ray(newOrigin, Normalize(Vector3f(pCamera)));

    Vector3<Float> dxCamera2 =
    (RasterToCamera(Point3f(1, 0, 0) * camerasperdim.x) - RasterToCamera(Point3f(0, 0, 0)));
    Vector3<Float> dyCamera2 =
    (RasterToCamera(Point3f(0, 1, 0) * camerasperdim.y) - RasterToCamera(Point3f(0, 0, 0)));
    
    ray->rxOrigin = ray->ryOrigin = ray->o;
    ray->rxDirection = Normalize(Vector3f(pCamera) + dxCamera2);
    ray->ryDirection = Normalize(Vector3f(pCamera) + dyCamera2);
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    ray->medium = medium;
    *ray = CameraToWorld(*ray);
    ray->hasDifferentials = true;
    return 1;

    // // test generateray
    // Float wt = GenerateRay(sample, ray);
    // if (wt == 0) return 0;
    // ray->hasDifferentials = false;
    // return wt;
}


LightFieldCamera *CreateLightfieldCamera(const ParamSet &params,
                                           const AnimatedTransform &cam2world,
                                           Film *film, const Medium *medium) {
    // Extract common camera parameters from _ParamSet_
    Float shutteropen = params.FindOneFloat("shutteropen", 0.f);
    Float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    if (shutterclose < shutteropen) {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
                shutterclose, shutteropen);
        std::swap(shutterclose, shutteropen);
    }
    Float lensradius = params.FindOneFloat("lensradius", 0.f);
    Float focaldistance = params.FindOneFloat("focaldistance", 1e6);
    Float frame = params.FindOneFloat(
        "frameaspectratio",
        Float(film->fullResolution.x) / Float(film->fullResolution.y));
    Bounds2f screen;
    if (frame > 1.f) {
        screen.pMin.x = -frame;
        screen.pMax.x = frame;
        screen.pMin.y = -1.f;
        screen.pMax.y = 1.f;
    } else {
        screen.pMin.x = -1.f;
        screen.pMax.x = 1.f;
        screen.pMin.y = -1.f / frame;
        screen.pMax.y = 1.f / frame;
    }
    int swi;
    const Float *sw = params.FindFloat("screenwindow", &swi);
    if (sw) {
        if (swi == 4) {
            screen.pMin.x = sw[0];
            screen.pMax.x = sw[1];
            screen.pMin.y = sw[2];
            screen.pMax.y = sw[3];
        } else
            Error("\"screenwindow\" should have four values");
    }
    Float fov = params.FindOneFloat("fov", 90.);
    Float halffov = params.FindOneFloat("halffov", -1.f);
    if (halffov > 0.f)
        // hack for structure synth, which exports half of the full fov
        fov = 2.f * halffov;
    /////////////// camerasperdim, cameragridbounds
    int cpi;
    Vector2i camerasperdim;
    const int *cp = params.FindInt("camerasperdim", &cpi);
    if (cp){
        if(cpi == 2){
           camerasperdim.x = cp[0];
           camerasperdim.y = cp[1];
        }
        else
            Error("\"camerasperdim\" should have two values");
    }
    
    
    int cgi;
    Bounds2f cameragridbounds;
    const Float *cg = params.FindFloat("cameragridbounds", &cgi);
    if (cg) {
        if (cgi == 4) {
            cameragridbounds.pMin.x = cg[0];
            cameragridbounds.pMax.x = cg[1];
            cameragridbounds.pMin.y = cg[2];
            cameragridbounds.pMax.y = cg[3];
        } else
            Error("\"cameragridbounds\" should have four values");
    }

    return new LightFieldCamera(cam2world, screen, shutteropen, shutterclose,
                                 lensradius, focaldistance, fov, film, medium,
                                 camerasperdim, cameragridbounds);
}

}  // namespace pbrt
