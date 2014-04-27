/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2012 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* The breaking dam free surface problem. This code demonstrates the basic usage of the
 * free surface module in Palabos. Surface tension and contact angles are optional. 
 */

#include "palabos3D.h"
#include "palabos3D.hh"
#include "river_functions.h"

using namespace plb;
using namespace std;

#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
typedef double T;


// Smagorinsky constant for LES model.
const T cSmago = 0.14;

// Physical dimensions of the system (in meters).
const T lx = 2.0;
const T ly = 2.0;
const T lz = 1.0;

const T rhoEmpty = T(1);
Array<T,3> forceOrientation(T(),T(),(T)1);
    
plint writeImagesIter   = 100;
plint getStatisticsIter = 20;

plint maxIter;
plint N;
plint nx, ny, nz;
T delta_t, delta_x;
Array<T,3> externalForce;
T nuPhys, nuLB, tau, omega, Bo, surfaceTensionLB, contactAngle;

std::string outDir;

void setupParameters() {
    delta_x = lz / N;
    nx = util::roundToInt(lx / delta_x);
    ny = util::roundToInt(ly / delta_x);
    nz = util::roundToInt(lz / delta_x);

    // Gravity in lattice units.
    T gLB = 9.8 * delta_t * delta_t/delta_x;
    externalForce = Array<T,3>(0., 0., -gLB);
    tau            = (nuPhys*DESCRIPTOR<T>::invCs2*delta_t)/(delta_x*delta_x) + 0.5;
    omega          = 1./tau;    
    nuLB           = (tau-0.5)*DESCRIPTOR<T>::cs2; // Viscosity in lattice units.
    
    surfaceTensionLB = rhoEmpty * gLB * N * N / Bo;
}

//create zone that creates water, only works as long as it is touching fluid
void addSpout(Box3D inlet, T pressure,
    FreeSurfaceFields3D<T,DESCRIPTOR>& fields)
{
  setToConstant(fields.flag, inlet, (int)twoPhaseFlag::fluid);
  integrateProcessingFunctional(
      new CreateLiquid3D<T, DESCRIPTOR>(fields.dynamics->clone(), pressure),
      inlet.enlarge(-2), fields.twoPhaseArgs, 4);
}

//create zone that destroys incident fluid
void addDrain(Box3D drain, FreeSurfaceFields3D<T, DESCRIPTOR>& fields)
{
  integrateProcessingFunctional(
      new CreateLiquid3D<T, DESCRIPTOR>(fields.dynamics->clone(), 0.9),
        drain, fields.twoPhaseArgs, 4);
}

void readMesh(const std::string& meshFile, 
    FreeSurfaceFields3D<T,DESCRIPTOR>& fields)
{
  pcout << "Voxelizing mesh...\n";
  TriangleSet<double> mesh(meshFile);
  Cuboid<T> mBounds = mesh.getBoundingCuboid();
  T meshHeight = mBounds.upperRightCorner[2] - mBounds.lowerLeftCorner[2];
  mBounds.upperRightCorner[2] += meshHeight * 0.2;

  auto wallFunc = wallFlagsFunction(mesh, Vec3<plint>(nx, ny, nz), mBounds);
  auto fluidFunc = [](plint x, plint y, plint z){
    return twoPhaseFlag::empty; };
    /*
      return (z > (nz * 0.75) && x > (nx / 2)) 
        ? twoPhaseFlag::fluid : twoPhaseFlag::empty; };
        */

  setToConstant(fields.flag, fields.flag.getBoundingBox(), (int)twoPhaseFlag::wall);
  setToFunction(fields.flag, fields.flag.getBoundingBox().enlarge(-1),
      InitialFlags(std::move(wallFunc), std::move(fluidFunc)));
  pcout << "Voxelization complete!\n";
}


void writeResults(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, MultiScalarField3D<T>& volumeFraction, plint iT)
{
    static const plint nx = lattice.getNx();
    static const plint ny = lattice.getNy();
    static const plint nz = lattice.getNz();
    Box3D slice(0, nx-1, ny/2, ny/2, 0, nz-1);
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledPpm(createFileName("u", iT, 6),
                               *computeVelocityNorm(lattice, slice)); 

    imageWriter.writeScaledPpm(createFileName("rho", iT, 6),
                               *computeDensity(lattice, slice));
                   
    imageWriter.writeScaledPpm(createFileName("volumeFraction", iT, 6), *extractSubDomain(volumeFraction, slice));

    // Use a marching-cube algorithm to reconstruct the free surface and write an STL file.
    std::vector<T> isoLevels;
    isoLevels.push_back((T) 0.5);
    typedef TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> triangles;
    isoSurfaceMarchingCube(triangles, volumeFraction, isoLevels, volumeFraction.getBoundingBox());
    TriangleSet<T>(triangles).writeBinarySTL(createFileName(outDir+"/interface", iT, 6)+".stl");

    VtkImageOutput3D<T> vtkOut(createFileName("volumeFraction", iT, 6), 1.);
    vtkOut.writeData<float>(volumeFraction, "vf", 1.);
}

void writeStatistics(FreeSurfaceFields3D<T,DESCRIPTOR>& fields) {
    pcout << " -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- " << endl;
    T averageMass = freeSurfaceAverageMass<T,DESCRIPTOR>(fields.twoPhaseArgs, fields.lattice.getBoundingBox());
    pcout << "Average Mass: " << averageMass  << endl;
    T averageDensity = freeSurfaceAverageDensity<T,DESCRIPTOR>(fields.twoPhaseArgs, fields.lattice.getBoundingBox());
    pcout << "Average Density: " << setprecision(12) << averageDensity  << endl;

    T averageVolumeFraction = freeSurfaceAverageVolumeFraction<T,DESCRIPTOR>(fields.twoPhaseArgs, fields.lattice.getBoundingBox());
    pcout << "Average Volume-Fraction: " << setprecision(12) << averageVolumeFraction  << endl;

    pcout << " -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- " << endl;
}


int main(int argc, char **argv)
{
    plbInit(&argc, &argv);
    global::directories().setInputDir("./");
    std::string stlFile;
    T inPressure;
        
    if (global::argc() != 9) {
        pcout << "Error missing some input parameter\n";
    }

    try {
        global::argv(1).read(outDir);
        global::directories().setOutputDir(outDir+"/");
        global::argv(2).read(stlFile);
        global::argv(3).read(inPressure);
        global::argv(4).read(nuPhys);
        global::argv(5).read(Bo);
        global::argv(6).read(contactAngle);
        global::argv(7).read(N);
        global::argv(8).read(delta_t);
        global::argv(9).read(maxIter);
    }
    catch(PlbIOException& except) {
        pcout << except.what() << std::endl;
        pcout << "The parameters for this program are :\n";
        pcout << "1. Output directory name.\n";
        pcout << "2. Mesh input file.\n";
        pcout << "3. Input pressure factor.\n";
        pcout << "4. kinematic viscosity in physical Units (m^2/s) .\n";
        pcout << "5. Bond number (Bo = rho * g * L^2 / gamma).\n";
        pcout << "6. Contact angle (in degrees).\n";
        pcout << "7. number of lattice nodes for lz .\n";
        pcout << "8. delta_t .\n"; 
        pcout << "9. maxIter .\n";
        pcout << "Reasonable parameters on a desktop computer are: " << (std::string)global::argv(0) << " tmp blah.stl 1.02 1.e-5 100 80.0 40 1.e-3 80000\n";
        pcout << "Reasonable parameters on a parallel machine are: " << (std::string)global::argv(0) << " tmp blah.tel 1.05 1.e-6 100 80.0 100 1.e-4 80000\n";
        exit (EXIT_FAILURE);
    }
    
    setupParameters();
    
    pcout << "delta_t= " << delta_t << endl;
    pcout << "delta_x= " << delta_x << endl;
    pcout << "delta_t*delta_t/delta_x= " << delta_t*delta_t/delta_x << endl;
    pcout << "externalForce= " << externalForce[2] << endl;
    pcout << "relaxation time= " << tau << endl;
    pcout << "omega= " << omega << endl;
    pcout << "kinematic viscosity physical units = " << nuPhys << endl;
    pcout << "kinematic viscosity lattice units= " << nuLB << endl;
    
    global::timer("initialization").start();
    

    SparseBlockStructure3D blockStructure(createRegularDistribution3D(nx, ny, nz));

    Dynamics<T,DESCRIPTOR>* dynamics
        = new SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega, cSmago);

    // If surfaceTensionLB is 0, then the surface tension algorithm is deactivated.
    // If contactAngle is less than 0, then the contact angle algorithm is deactivated.
    FreeSurfaceFields3D<T,DESCRIPTOR> fields(blockStructure, dynamics->clone(), rhoEmpty,
                                             surfaceTensionLB, contactAngle, externalForce);
    //integrateProcessingFunctional(new ShortenBounceBack3D<T,DESCRIPTOR>, fields.lattice.getBoundingBox(), fields.twoPhaseArgs, 0);

    //read mesh file to set flags
    readMesh(stlFile, fields);
    Box3D inlet(nx * 0.85, nx - 1, ny * 0.85, ny - 1, 1, nz * 0.2);
    addSpout(inlet, inPressure, fields);

    Box3D drain(1, 2, 1, 2, 1, nz * 0.2);
    addDrain(drain, fields);
    
    fields.defaultInitialize();

    pcout << "Time spent for setting up lattices: "
          << global::timer("initialization").stop() << endl;
    T lastIterationTime = T();

    for (plint iT = 0; iT <= maxIter; ++iT) {
        global::timer("iteration").restart();
        
        T sum_of_mass_matrix = T();
        T lost_mass = T();
        if (iT % getStatisticsIter==0) {
            pcout << endl;
            pcout << "ITERATION = " << iT << endl;
            pcout << "Time of last iteration is " << lastIterationTime << " seconds" << endl;
            writeStatistics(fields);
            sum_of_mass_matrix = fields.lattice.getInternalStatistics().getSum(0);
            pcout << "Sum of mass matrix: " << sum_of_mass_matrix << std::endl;
            lost_mass = fields.lattice.getInternalStatistics().getSum(1);
            pcout << "Lost mass: " << lost_mass << std::endl;
            pcout << "Total mass: " << sum_of_mass_matrix + lost_mass << std::endl;
            pcout << "Interface cells: " << fields.lattice.getInternalStatistics().getIntSum(0) << std::endl;
        }

        if (iT % writeImagesIter == 0) {
            global::timer("images").start();
            writeResults(fields.lattice, fields.volumeFraction, iT);
            pcout << "Total time spent for writing images: "
                << global::timer("images").stop() << endl;
        }                           

        // This includes the collision-streaming cycle, plus all free-surface operations.
        fields.lattice.executeInternalProcessors();
        fields.lattice.evaluateStatistics();
        fields.lattice.incrementTime();

        lastIterationTime = global::timer("iteration").stop();
    }
}

