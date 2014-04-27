#ifndef RIVER_FUNCTIONS_HPP
#define RIVER_FUNCTIONS_HPP

#include "river_functions.h"
#include <functional>
#include "palabos3D.h"
#include "palabos3D.hh"
#include "quadtree.h"
#include <algorithm>

namespace tree{
//function overload to allow triangles to be placed in the quadtree
template<typename T>
bool within(const plb::Array<plb::Array<T, 3>, 3>& shape,
    const tree::Rect<T>& r)
{
  for(int i = 0; i < 3; ++i) {
    if(shape[i][0] < r.xMin || shape[i][0] > r.xMax 
        || shape[i][1] < r.yMin || shape[i][1] > r.yMax) {
      return false;
    }
  }
  return true;
}
}

//Moeller-Trumbore intersection algorithm
//for rays and triangles
inline bool rayIntersects(const Triangle<double>& tri,
    const Vec3<double>& origin, const Vec3<double>& dir)
{
  static const double EPSILON = 0.00001;
  Vec3<double> edge1 = tri[1] - tri[0];
  Vec3<double> edge2 = tri[2] - tri[0];

  //calculate determinant
  Vec3<double> pVec = plb::crossProduct(dir, edge2);
  double det = dot(edge1, pVec);
  if(det > -EPSILON && det < EPSILON) {
    return false;
  }
  double inv_det = 1.0 / det;

  Vec3<double> tVec = origin - tri[0];
  double u = dot(tVec, pVec) * inv_det;
  //the intersection lies outside the triangle
  if(u < 0 || u > 1) {
    return false;
  }

  Vec3<double> qVec = plb::crossProduct(tVec, edge1);
  double v = dot(dir, qVec) * inv_det;
  if(v < 0 || (u + v) > 1) {
    return false;
  }

  double t = dot(edge2, qVec) * inv_det;

  if(t > EPSILON) {
    return true;
  }

  return false;
}

inline std::function<int (plb::plint, plb::plint, plb::plint)> wallFlagsFunction(
    const plb::TriangleSet<double>& mesh, Vec3<plb::plint> latticeSize,
    plb::Cuboid<double> meshBounds)
{
  using namespace plb;
  tree::QuadTree<double, Triangle<double>> tree(
      meshBounds.lowerLeftCorner[0], meshBounds.upperRightCorner[0],
      meshBounds.lowerLeftCorner[1], meshBounds.upperRightCorner[1]);
  Vec3<double> meshScale = (meshBounds.upperRightCorner - 
    meshBounds.lowerLeftCorner);
  for(int i = 0; i < 3; ++i) {
    meshScale[i] /= latticeSize[i];
  }
  for(auto& tri : mesh.getTriangles()){
    tree.insert(tri);
  }

  //copying the tree, but whatever we only have to do it once
  //and capturing by move is convoluted
  return [=](plint pX, plint pY, plint pZ){
    const Vec3<double> direction = { 0.0, 0.0, 1.0 };
    Vec3<double> origin = {static_cast<double>(pX),
      static_cast<double>(pY), static_cast<double>(pZ)};
    origin *= meshScale;
    origin += meshBounds.lowerLeftCorner;
    auto first = tree.beginAt(origin[0], origin[1]);
    auto last = tree.end();
    int hits = std::count_if(first, last, [&](const Triangle<double>& t){
        return rayIntersects(t, origin, direction);
        });
    return ((hits & 1) == 0) ? twoPhaseFlag::empty : twoPhaseFlag::wall;
  };
}
// ============================================
// CREATE LIQUID
// ============================================

template<typename T, template<typename U> class Descriptor>
void CreateLiquid3D<T,Descriptor>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks)
{
  typedef Descriptor<T> D;
  using namespace twoPhaseFlag;
  FreeSurfaceProcessorParam3D<T,Descriptor> param(atomicBlocks);

  for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
    for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
      for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
        if(iY & 0 == 0)
        {
          continue;
        }
        T iniRho = T(inflowSpeed);
        param.attributeDynamics (
            iX,iY,iZ, dynamicsTemplate->clone() );
        //iniCellAtEquilibrium(param.cell(iX,iY,iZ), iniRho, injectionVelocity);
        param.setDensity(iX,iY,iZ, iniRho);
        //param.setMomentum(iX,iY,iZ, iniRho*injectionVelocity);
        param.mass(iX,iY,iZ) = iniRho;
        param.volumeFraction(iX,iY,iZ) = (T)1;
        param.flag(iX,iY,iZ) = fluid;

      }
    }
  }
}

template<typename T, template<typename U> class Descriptor>
CreateLiquid3D<T,Descriptor>* CreateLiquid3D<T,Descriptor>::clone() const {
  return new CreateLiquid3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void DestroyLiquid3D<T,Descriptor>::processGenericBlocks(Box3D domain,std::vector<AtomicBlock3D*> atomicBlocks)
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam3D<T,Descriptor> param(atomicBlocks);
            
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                //param.attributeDynamics(iX,iY,iZ, new NoDynamics<T,Descriptor>());
                param.setDensity(iX,iY,iZ, suckPower);
                param.mass(iX,iY,iZ) = suckPower;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
DestroyLiquid3D<T,Descriptor>* DestroyLiquid3D<T,Descriptor>::clone() const {
    return new DestroyLiquid3D<T,Descriptor>(*this);
}

#endif
