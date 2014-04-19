/*
 *  Extended functionality for importing models to Palabos
 *  Brock Jackman
 *
 */
#ifndef RIVER_FUNCTIONS_H
#define RIVER_FUNCTIONS_H

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;

template <typename T>
using Vec3 = plb::Array<T, 3>;

template <typename T>
using Triangle = typename plb::TriangleSet<T>::Triangle;

//palabos is hiding its dot product someplace weird
double dot(const Vec3<double>& v1, const Vec3<double>& v2)
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

//returns true if the ray from origin in direction dir intersects tri
bool rayIntersects(const Triangle<double>& tri, 
    const Vec3<double>& origin, const Vec3<double>& dir);

typedef std::function<int (plb::plint, plb::plint, plb::plint)> FlagFunc;

//returns function to determine if a given point is within a mesh
FlagFunc wallFlagsFunction(
    const plb::TriangleSet<double>& mesh, Vec3<plb::plint> latticeSize,
    plb::Cuboid<double> meshBounds);


//functor to initialize both walls and water
struct InitialFlags
{
  InitialFlags(FlagFunc wall, FlagFunc fluid) : wall(wall), fluid(fluid) {}
  int operator()(plb::plint x, plb::plint y, plb::plint z)
  {
    if(wall(x, y, z) == plb::twoPhaseFlag::wall) {
      return plb::twoPhaseFlag::wall;
    }
    return fluid(x, y, z);
  }

  FlagFunc wall, fluid;
};

//Fluid creation data processor
template< typename T,template<typename U> class Descriptor>
class CreateLiquid3D : public BoxProcessingFunctional3D {
public:
    CreateLiquid3D(Dynamics<T,Descriptor>* dynamicsTemplate_, T inflowSpeed)
        : dynamicsTemplate(dynamicsTemplate_), inflowSpeed(inflowSpeed)
    { }
    CreateLiquid3D(CreateLiquid3D<T,Descriptor> const& rhs)
        : dynamicsTemplate(rhs.dynamicsTemplate->clone()),
          inflowSpeed(rhs.inflowSpeed)
    { }
    CreateLiquid3D<T,Descriptor>* operator=(CreateLiquid3D<T,Descriptor> const& rhs)
    { 
        CreateLiquid3D<T,Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(CreateLiquid3D<T,Descriptor>& rhs) {
        std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
        std::swap(inflowSpeed, rhs.inflowSpeed);
    }
    virtual ~CreateLiquid3D() {
        delete dynamicsTemplate;
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual CreateLiquid3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::dataStructure;    // Fluid
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass
        modified[4] = modif::staticVariables;  // Mass-fraction
        modified[5] = modif::staticVariables;  // Flag-status
        modified[6] = modif::nothing;          // Normal
        modified[7] = modif::nothing;          // Interface lists
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }
private:
    Dynamics<T,Descriptor>* dynamicsTemplate;
    T inflowSpeed;
};

#include "river_functions.hpp"
#endif
