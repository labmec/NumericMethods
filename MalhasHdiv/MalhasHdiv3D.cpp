//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"
// #include "pzbiharmonic.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzcompel.h"
#include "TPZHybridizeHDiv.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include "TPZAnalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzgeopyramid.h"
#include "TPZGeoLinear.h"

#include <TPZRefPattern.h>

#include "TPZMaterial.h"
#include "Elasticity/TPZElasticity2D.h"
#include "Elasticity/TPZElasticity3D.h"
#include "pzlog.h"

#include "TPZGenGrid2D.h"
#include "TPZGenGrid3D.h"

#include <time.h>
#include <stdio.h>

#include <math.h>

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

#include <set>
#include <map>
#include <vector>
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZVTKGeoMesh.h"
#include "TPZExtendGridDimension.h"

#include "Poisson/TPZMatPoisson.h"
#include "TPZNullMaterial.h"
//#include "mixedpoisson.h"
#include "TPZBndCond.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "TPZLinearAnalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSolver.h"
#include "TPZPersistenceManager.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZCompMeshTools.h"
#include "TPZMixedDarcyWithFourSpaces.h"
#include "pzcmesh.h"
#include "pzbdstrmatrix.h"

#include "TPZHdivTransfer.h"
#include "pzseqsolver.h"
#include "pzmgsolver.h"
#include "TPZTimer.h"
#include "TPZPrintUtils.cpp"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "ScenarioConfig.h"

#ifdef _AUTODIFF
#include "tfad.h"
#include "fad.h"
#include "pzextractval.h"
#endif


ScenarioConfig myScenario;

std::ofstream rprint("Harmonic3D_Scenario8.txt",std::ofstream::out);

/**
 * @brief Generates the force function for the 1D case
 * @param pt: Points values
 * @param disp: Vector
 */
auto Ladoderecho_1D = [](const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    STATE x = pt[0];
    double fx= 4*M_PI*M_PI*sin(2*M_PI*x);        //Force function definition
    disp[0]=fx;
    
};

/**
 * @brief Generates the force function for the 2D case
 * @param pt: Points values
 * @param disp: Vector
 */

auto Ladoderecho_2D = [](const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    STATE x = pt[0];
    STATE y = pt[1];
    STATE z = pt[2];
    
    // double fx=  2*(x-1)*x*(y-1)*y + 2*(x-1)*x*(z-1)*z + 2*(y-1)*y*(z-1)*z; //Force function definition
    double fx = 0.;
    
    //    double fx =-4144.653167389283*pow(10,
    //                                      -pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    (-2*M_PI + 15.*x) + 4771.70829943056*
    //    pow(10,-pow(-2*M_PI + 15.*x,2) -
    //        pow(-2*M_PI + 15.*y,2))*pow(-2*M_PI + 15.*x,3) +
    //    4771.70829943056*pow(10,
    //                         -pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    (-2*M_PI + 15.*x)*pow(-2*M_PI + 15.*y,2);
    //
    //    double fx = -4144.653167389282*pow(2,2 - pow(-2*M_PI + 15.*x,2) -
    //                                       pow(-2*M_PI + 15.*y,2))*
    //    pow(5,-pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    (-2*M_PI + 15.*x) + 4771.708299430558*
    //    pow(2,2 - pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    pow(5,-pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    pow(-2*M_PI + 15.*x,3) +
    //    4771.708299430558*pow(2,
    //                          2 - pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    pow(5,-pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    (-2*M_PI + 15.*x)*pow(-2*M_PI + 15.*y,2);
    
    disp[0]=fx;
    
    
};

auto exactSol = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    // u[0]= std::sin(M_PI*x)*std::sin(M_PI*y);
    // u[0] = x*x*x*y - y*y*y*x;    

    // gradU(0,0) = -M_PI*cos(M_PI*x)*sin(M_PI*y);
    // gradU(1,0) = -M_PI*cos(M_PI*y)*sin(M_PI*x);
    // gradU(0,0) = (3.*x*x*y - y*y*y);
    // gradU(1,0) = (x*x*x - 3.*y*y*x);

    // u[0] = (x-1)*x*(y-1)*y*(z-1)*z;
    // gradU(0,0) = (x-1)*(y-1)*y*(z-1)*z + x*(y-1)*y*(z-1)*z;
    // gradU(1,0) = (x-1)*x*(y-1)*(z-1)*z + (x-1)*x*y*(z-1)*z;
    // gradU(2,0) = (x-1)*x*(y-1)*y*(z-1) + (x-1)*x*(y-1)*y*z;

    REAL aux = 1./sinh(sqrt(2.)*M_PI);
    u[0] = sin(M_PI*x)*sin(M_PI*y)*sinh(sqrt(2.)*M_PI*z)*aux;
    gradU(0,0) = M_PI*cos(M_PI*x)*sin(M_PI*y)*sinh(sqrt(2.)*M_PI*z)*aux;
    gradU(1,0) = M_PI*cos(M_PI*y)*sin(M_PI*x)*sinh(sqrt(2.)*M_PI*z)*aux;
    gradU(2,0) = sqrt(2.)*M_PI*cosh(sqrt(2.)*M_PI*z)*sin(M_PI*x)*sin(M_PI*y)*aux;


    //exact Divergent
    // gradU(2,0) = -2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
    

};



//Creating geometric 1D and 2D mesh
TPZGeoMesh * GenerateGmeshOne(int nx, double l);
TPZGeoMesh * GenerateGmesh(int nx, int ny, double l, double h);

//Stablish de force fuction for 1D and 2D mesh
// void Ladoderecho_1D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
// void Ladoderecho_2D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//Stablish an exact solution
void SolExact(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);

//Creating computational meshes
TPZCompMesh * GeneratePressureCmesh(TPZGeoMesh *Gmesh, int order_internal);
TPZCompMesh * GenerateConstantCmesh(TPZGeoMesh *Gmesh, bool third_LM);
TPZCompMesh * GenerateFluxCmesh(TPZGeoMesh *Gmesh, int order_internal, int order_border, bool coarse=false);
TPZMultiphysicsCompMesh * GenerateMixedCmesh(TPZVec<TPZCompMesh *> fvecmesh, int order, bool two_d_Q);

//Creates index vector
void IndexVectorCoFi(TPZMultiphysicsCompMesh *Coarse_sol, TPZMultiphysicsCompMesh *Fine_sol, TPZVec<int64_t> & indexvec);

//Hdiv Test
void HDiv(int nx, int order_small, int order_high, bool condense_equations_Q, bool two_d_Q);


//Transfer DOF from coarse mesh to fine mesh
void TransferDegreeOfFreedom(TPZFMatrix<STATE> & CoarseDoF, TPZFMatrix<STATE> & FineDoF, TPZVec<int64_t> & DoFIndexes);

//Analysis configuration
void ConfigurateAnalyses(TPZCompMesh * cmesh_c, TPZCompMesh * cmesh_f, bool must_opt_band_width_Q, int number_threads, TPZAnalysis &an_c,TPZAnalysis &an_f, bool UsePardiso_Q);

//Split connects from a certain mesh
void SplitConnects(TPZCompMesh *fluxmesh, TPZGeoEl *gel ,int j);

//Shows shape functions for a certain element
void ShowShape(TPZCompMesh * cmesh, int element, int funcion,std::string plotname);
void ShowShape(TPZCompMesh * cmesh, int element, int funcion,std::string plotname);

void HdiVSimple(int nx, int order_high, bool condense_equations_Q, bool two_d_Q);

void AssociateElements(TPZMultiphysicsCompMesh *cmesh, TPZVec<int64_t> &elementgroup, std::set<int> &matId);
void CondenseBCElements(TPZMultiphysicsCompMesh *cmesh, std::set<int> &matIdBC);


using namespace std;

int main(){
    
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    TPZTimer clock;
    clock.start();
    int k = 2;
    std::cout << "\nSolving problem 2, k = " << k << "\n\n";
    HDiv(2, 1, k, true, true);
    std::cout << "\nSolving problem 3, k = " << k << "\n\n";
    HDiv(3, 1, k, true, true);
    std::cout << "\nSolving problem 4, k = " << k << "\n\n";
    HDiv(4, 1, k, true, true);    
    std::cout << "\nSolving problem 5, k = " << k << "\n\n";
    HDiv(5, 1, k, true, true);    
    std::cout << "\nSolving problem 6, k = " << k << "\n\n";
    HDiv(6, 1, k, true, true);    
    std::cout << "\nSolving problem 7, k = " << k << "\n\n";
    HDiv(7, 1, k, true, true);    
    std::cout << "\nSolving problem 8, k = " << k << "\n\n";
    HDiv(8, 1, k, true, true);    
    std::cout << "\nSolving problem 9, k = " << k << "\n\n";
    HDiv(9, 1, k, true, true);    
    std::cout << "\nSolving problem 10, k = " << k << "\n\n";
    HDiv(10, 1, k, true, true);    
    std::cout << "\nSolving problem 11, k = " << k << "\n\n";
    HDiv(11, 1, k, true, true);    
    std::cout << "\nSolving problem 12, k = " << k << "\n\n";
    HDiv(12, 1, k, true, true);    
    std::cout << "\nSolving problem 13, k = " << k << "\n\n";
    HDiv(13, 1, k, true, true);    
    std::cout << "\nSolving problem 14, k = " << k << "\n\n";
    HDiv(14, 1, k, true, true);    
    std::cout << "\nSolving problem 15, k = " << k << "\n\n";
    HDiv(15, 1, k, true, true);    
    std::cout << "\nSolving problem 16, k = " << k << "\n\n";
    HDiv(16, 1, k, true, true);

    k = 3;
    std::cout << "\nSolving problem 2, k = " << k << "\n\n";
    HDiv(2, 1, k, true, true);
    std::cout << "\nSolving problem 3, k = " << k << "\n\n";
    HDiv(3, 1, k, true, true);
    std::cout << "\nSolving problem 4, k = " << k << "\n\n";
    HDiv(4, 1, k, true, true);    
    std::cout << "\nSolving problem 5, k = " << k << "\n\n";
    HDiv(5, 1, k, true, true);    
    std::cout << "\nSolving problem 6, k = " << k << "\n\n";
    HDiv(6, 1, k, true, true);    
    std::cout << "\nSolving problem 7, k = " << k << "\n\n";
    HDiv(7, 1, k, true, true);    
    std::cout << "\nSolving problem 8, k = " << k << "\n\n";
    HDiv(8, 1, k, true, true);    
    std::cout << "\nSolving problem 9, k = " << k << "\n\n";
    HDiv(9, 1, k, true, true);    
    std::cout << "\nSolving problem 10, k = " << k << "\n\n";
    HDiv(10, 1, k, true, true);    
    std::cout << "\nSolving problem 11, k = " << k << "\n\n";
    HDiv(11, 1, k, true, true);    
    std::cout << "\nSolving problem 12, k = " << k << "\n\n";
    HDiv(12, 1, k, true, true);    
    std::cout << "\nSolving problem 13, k = " << k << "\n\n";
    HDiv(13, 1, k, true, true);    
    std::cout << "\nSolving problem 14, k = " << k << "\n\n";
    HDiv(14, 1, k, true, true);    
    std::cout << "\nSolving problem 15, k = " << k << "\n\n";
    HDiv(15, 1, k, true, true);    
    std::cout << "\nSolving problem 16, k = " << k << "\n\n";
    HDiv(16, 1, k, true, true);

    k = 4;
    std::cout << "\nSolving problem 2, k = " << k << "\n\n";
    HDiv(2, 1, k, true, true);
    std::cout << "\nSolving problem 3, k = " << k << "\n\n";
    HDiv(3, 1, k, true, true);
    std::cout << "\nSolving problem 4, k = " << k << "\n\n";
    HDiv(4, 1, k, true, true);    
    std::cout << "\nSolving problem 5, k = " << k << "\n\n";
    HDiv(5, 1, k, true, true);    
    std::cout << "\nSolving problem 6, k = " << k << "\n\n";
    HDiv(6, 1, k, true, true);    
    std::cout << "\nSolving problem 7, k = " << k << "\n\n";
    HDiv(7, 1, k, true, true);    
    std::cout << "\nSolving problem 8, k = " << k << "\n\n";
    HDiv(8, 1, k, true, true);    
    std::cout << "\nSolving problem 9, k = " << k << "\n\n";
    HDiv(9, 1, k, true, true);    
    std::cout << "\nSolving problem 10, k = " << k << "\n\n";
    HDiv(10, 1, k, true, true);    
    std::cout << "\nSolving problem 11, k = " << k << "\n\n";
    HDiv(11, 1, k, true, true);    
    std::cout << "\nSolving problem 12, k = " << k << "\n\n";
    HDiv(12, 1, k, true, true);    
    std::cout << "\nSolving problem 13, k = " << k << "\n\n";
    HDiv(13, 1, k, true, true);    
    std::cout << "\nSolving problem 14, k = " << k << "\n\n";
    HDiv(14, 1, k, true, true);    
    std::cout << "\nSolving problem 15, k = " << k << "\n\n";
    HDiv(15, 1, k, true, true);    
    std::cout << "\nSolving problem 16, k = " << k << "\n\n";
    HDiv(16, 1, k, true, true);

    k = 5;
    std::cout << "\nSolving problem 2, k = " << k << "\n\n";
    HDiv(2, 1, k, true, true);
    std::cout << "\nSolving problem 3, k = " << k << "\n\n";
    HDiv(3, 1, k, true, true);
    std::cout << "\nSolving problem 4, k = " << k << "\n\n";
    HDiv(4, 1, k, true, true);    
    std::cout << "\nSolving problem 5, k = " << k << "\n\n";
    HDiv(5, 1, k, true, true);    
    std::cout << "\nSolving problem 6, k = " << k << "\n\n";
    HDiv(6, 1, k, true, true);    
    std::cout << "\nSolving problem 7, k = " << k << "\n\n";
    HDiv(7, 1, k, true, true);    
    std::cout << "\nSolving problem 8, k = " << k << "\n\n";
    HDiv(8, 1, k, true, true);    
    std::cout << "\nSolving problem 9, k = " << k << "\n\n";
    HDiv(9, 1, k, true, true);    
    std::cout << "\nSolving problem 10, k = " << k << "\n\n";
    HDiv(10, 1, k, true, true);    
    std::cout << "\nSolving problem 11, k = " << k << "\n\n";
    HDiv(11, 1, k, true, true);    
    std::cout << "\nSolving problem 12, k = " << k << "\n\n";
    HDiv(12, 1, k, true, true);    
    std::cout << "\nSolving problem 13, k = " << k << "\n\n";
    HDiv(13, 1, k, true, true);    
    std::cout << "\nSolving problem 14, k = " << k << "\n\n";
    HDiv(14, 1, k, true, true);    
    std::cout << "\nSolving problem 15, k = " << k << "\n\n";
    HDiv(15, 1, k, true, true);    
    std::cout << "\nSolving problem 16, k = " << k << "\n\n";
    HDiv(16, 1, k, true, true);     
    
    //    HdiVSimple(30, 2, true, true);
    clock.stop();
    std::ofstream Out("TotalTime.txt");
    operator<<(Out, clock );
}

/**
 * @brief Runs a HDiv problem with 4 spaces for 1D or 2D case
 * @param order_small: Low order for border elements
 * @param order_high: High order for internal elements
 * @param condense_equations_Q: Bool wether the problem is condensed or not
 * @param two_d_Q: Bool wether the problem es 1D (false) or 2D (true)
 */
void HDiv(int nx, int order_small, int order_high, bool condense_equations_Q, bool two_d_Q){
    
    //    TPZGeoMesh *gmesh_11 = GenerateGmesh(nx, nx, 2, 2);      // Generates a 2D geo mesh
    //    TPZCompMesh *flux = GenerateFluxCmesh(gmesh_11, 1, 1);
    //
    //    int el_index=3;
    //    int nfun=flux->Element(el_index)->NEquations();
    //    for (int i=0; i<nfun; i++) {
    //
    //        std::string filename("elementFunc.vtk");
    //        std::string file(filename+std::to_string(i)+".vtk");
    //        ShowShape(flux,el_index,i,file);
    //    };

    myScenario.BasePOrder = order_high;
    myScenario.Scenario = EScenario::Scenario8;
    myScenario.ConfigurateScenario();
    
    bool KeepOneLagrangian = true;
    bool KeepMatrix = false;
    bool render_shapes_Q = false;               //Prints a .VTK showing the render shapes
    bool must_opt_band_width_Q = true;
    int number_threads = 12;
    
    TPZGeoMesh *gmesh;
    
    // Generates a 2D geo mesh
    TPZGeoMesh *gmesh_1 = GenerateGmesh(nx, nx, 1, 1);
    
    // Generates a 1D geo mesh
    // TPZGeoMesh *gmesh_2 = GenerateGmeshOne(nx, 1);
    
    //Asks if the problem is 2D
    if (two_d_Q) {
        gmesh = gmesh_1;
    } else {
        DebugStop();
        // gmesh = gmesh_2;
    }
    {
        // TPZPrintUtils util;
        // util.PrintGeoMesh(gmesh);
    }


    TPZMultiphysicsCompMesh *MixedMesh_c = 0;
    TPZManVector<TPZCompMesh *> vecmesh_c(4);      //Vector for coarse mesh case (4 spaces)
    {
        TPZCompMesh *q_cmesh = GenerateFluxCmesh(gmesh, myScenario.CoarseBubbleFluxPOrder, order_small, true);
        TPZCompMesh *p_cmesh = GeneratePressureCmesh(gmesh, myScenario.CoarsePressurePOrder);
        TPZCompMesh *gavg_cmesh = GenerateConstantCmesh(gmesh,false);
        TPZCompMesh *pavg_cmesh = GenerateConstantCmesh(gmesh,true);
        vecmesh_c[0] = q_cmesh;              //Flux
        vecmesh_c[1] = p_cmesh;              //Pressure
        vecmesh_c[2] = gavg_cmesh;           //Average distribute flux
        vecmesh_c[3] = pavg_cmesh;           //Average pressure
        
        MixedMesh_c = GenerateMixedCmesh(vecmesh_c, 1, two_d_Q);       //1 Stands for the corse mesh order
        
        TPZPrintUtils util;
        util.PrintCompMesh(vecmesh_c[0],"q_cmeshC");
        util.PrintCompMesh(vecmesh_c[1],"p_cmeshC");
        util.PrintCompMesh(vecmesh_c[2],"gavg_cmeshC");
        util.PrintCompMesh(vecmesh_c[3],"pavg_cmeshC");
        util.PrintCompMesh(MixedMesh_c,"MixedMesh_C");
    }
    
    std::set<int> matIDBC = {-1,-2,-3,-4,-5,-6};
    if (condense_equations_Q) {             //Asks if you want to condesate the problem
        // Created condensed elements for the elements that have internal nodes
        // std::cout << "CNEQUATIONS1 = " << MixedMesh_c->NEquations() << std::endl;
        CondenseBCElements(MixedMesh_c,matIDBC);
        // std::cout << "CNEQUATIONS2 = " << MixedMesh_c->NEquations() << std::endl;
        TPZCompMeshTools::CondenseElements(MixedMesh_c, 3, KeepMatrix);
        // std::cout << "CNEQUATIONS3 = " << MixedMesh_c->NEquations() << std::endl;
        rprint << MixedMesh_c->NEquations() << " ";
    }
    
    TPZMultiphysicsCompMesh * MixedMesh_f = 0;
    TPZManVector<TPZCompMesh *> vecmesh_f(4);      //Vector for fine mesh case (4 spaces)
    {
        TPZCompMesh *q_cmesh = GenerateFluxCmesh(gmesh, myScenario.FineBubbleFluxPOrder, order_small, false);
        TPZCompMesh *p_cmesh = GeneratePressureCmesh(gmesh, myScenario.FinePressurePOrder);
        TPZCompMesh *gavg_cmesh = GenerateConstantCmesh(gmesh,false);
        TPZCompMesh *pavg_cmesh = GenerateConstantCmesh(gmesh,true);
        vecmesh_f[0] = q_cmesh;              //Flux
        vecmesh_f[1] = p_cmesh;              //Pressure
        vecmesh_f[2] = gavg_cmesh;           //Average distribute flux
        vecmesh_f[3] = pavg_cmesh;           //Average pressure
        
        MixedMesh_f = GenerateMixedCmesh(vecmesh_f, 2, two_d_Q); //2 Stands for the corse mesh order

        TPZPrintUtils util;
        util.PrintCompMesh(vecmesh_f[0],"q_cmeshF");
        util.PrintCompMesh(vecmesh_f[1],"p_cmeshF");
        util.PrintCompMesh(vecmesh_f[2],"gavg_cmeshF");
        util.PrintCompMesh(vecmesh_f[3],"pavg_cmeshF");
        util.PrintCompMesh(MixedMesh_f,"MixedMesh_F");
    }
    
    
    //Asks if you want to condesate the problem
    if (condense_equations_Q) {
        // Created condensed elements for the elements that have internal nodes
        // std::cout << "NEQUATIONS1 = " << MixedMesh_f->NEquations() << std::endl;
        rprint << MixedMesh_f->NEquations() << " ";
        CondenseBCElements(MixedMesh_f,matIDBC);
        rprint << MixedMesh_f->NEquations() << " ";
        // std::cout << "NEQUATIONS2 = " << MixedMesh_f->NEquations() << std::endl;
        TPZCompMeshTools::CondenseElements(MixedMesh_f, 3, KeepMatrix);
        // std::cout << "NEQUATIONS3 = " << MixedMesh_f->NEquations() << std::endl;

        
        // TPZPrintUtils util;
        // util.PrintCompMesh(MixedMesh_f,"MixedMesh_F_condensed");
    }
    
    
    //Solving the system:
    MixedMesh_c->InitializeBlock();    //Resequence the block object, remove unconnected connect objects
    MixedMesh_f->InitializeBlock();    //and reset the dimension of the solution vector
    TPZLinearAnalysis an_c(MixedMesh_c,must_opt_band_width_Q);// = new TPZLinearAnalysis();
    TPZLinearAnalysis an_f(MixedMesh_f,must_opt_band_width_Q);// = new TPZLinearAnalysis();
    //True to use pardiso
    ConfigurateAnalyses(MixedMesh_c, MixedMesh_f, must_opt_band_width_Q, number_threads, an_c, an_f, false);
    
    if(render_shapes_Q){
        TPZLinearAnalysis anloc(MixedMesh_f,false);
        std::string filename("Shape.vtk");
        TPZVec<int64_t> indices(18);
        for(int i=0; i<18; i++) indices[i] = i;
        const TPZVec<std::string> varname(1);
        varname[0]="Flux";
        anloc.ShowShape(filename, indices,1,varname);
    }
    
    
    TPZTimer AssemblyFine;
    AssemblyFine.start();
    
    an_c.SetExact(exactSol,3);
    an_f.SetExact(exactSol,3);

    // std::cout << "Matriz equations = " << an_f->StructMatrix()->Matrix()->Rows()<< std::endl;

    // Assembly fine operator
    std::cout << "Assembling fine model... \n"; 
    an_f.Assemble();
    std::cout << "Finish assembling fine model... \n";
    AssemblyFine.stop();
    std::ofstream Out_AssemblyFine("Assembly_Fine.txt");
    operator<<(Out_AssemblyFine, AssemblyFine );
    
    std::ofstream filemate("MatrixFine.txt");
    // an_f.MatrixSolver<STATE>().Matrix()->Print("K=",filemate,EMathematicaInput);
    // an_f.Rhs().Print("F=",filemate,EMathematicaInput);

    TPZTimer AssemblyAndSolvingCoarse;
    AssemblyAndSolvingCoarse.start();
    
    // Assembly for coarse operator
    std::cout << "Assembling coarse model... \n";
    an_c.Assemble();
    std::cout << "Finish Assembling coarse model... \n";
    //    an_c->Solve();      //Sin esto se pierde la solución
    
    AssemblyAndSolvingCoarse.stop();
    
    std::ofstream Out_AssemblyAndSolvingCoarse("Assembly_Coarse.txt");
    operator<<(Out_AssemblyAndSolvingCoarse, AssemblyAndSolvingCoarse );
    

    // An iterative solution
    {
        TPZTimer clock4;
        clock4.start();
        
        // Constructing block diagonal
        if(1){
            TPZBlockDiagonalStructMatrix<STATE> bdstr(MixedMesh_f);     //Give the fine
            
            
            TPZBlockDiagonal<STATE> * sp = new TPZBlockDiagonal<STATE>();
            
            bdstr.AssembleBlockDiagonal(*sp);
            
            TPZAutoPointer<TPZMatrix<STATE> > sp_auto(sp);
            int64_t n_con = MixedMesh_f->NConnects();
            for (int ic = 0; ic < n_con; ic++) {
                TPZConnect & con = MixedMesh_f->ConnectVec()[ic];
                bool check = con.IsCondensed() || con.HasDependency() || con.LagrangeMultiplier() == 0;
                if (check) {
                    continue;
                }
                
                int64_t seqnum = con.SequenceNumber();
                int block_size = MixedMesh_f->Block().Size(seqnum);
                if (block_size != 1) {
                    continue;
                }
                
                int64_t pos = MixedMesh_f->Block().Position(seqnum);
                (*sp).PutVal(pos, pos, 1.0);
            }
            
            clock4.stop();
            std::ofstream Out4("ConstructingDiagonalBlock.txt");
            operator<<(Out4, clock4 );
            
            TPZTimer clock5;
            clock5.start();
            
            TPZVec<int64_t> Indexes;
            IndexVectorCoFi(MixedMesh_c, MixedMesh_f, Indexes);
            int64_t neq_coarse = MixedMesh_c->NEquations();
            int64_t neq_fine = MixedMesh_f->NEquations();
            TPZHdivTransfer<STATE> *transfer = new TPZHdivTransfer<STATE>(neq_coarse, neq_fine, Indexes);
            TPZFMatrix<STATE> coarsesol(neq_coarse,1,1.), finesol(neq_fine,1,1.);
            transfer->Multiply(finesol, coarsesol,0); //It mutiplies itself by TPZMatrix<TVar>A adding the result in res
            transfer->Multiply(coarsesol, finesol,1); //z = beta * y(coarse) + alpha * opt(this)*x (fine)
            //             finesol.Print(std::cout);
            //             coarsesol.Print(std::cout);
            
            //Transfers the solution from coarse to fine mesh
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt);
            //            coarsesol.Print(std::cout);
            // an_c->Solver().Solve(coarsesol,coarsesol); //Force vector, solution
            an_c.Solve(); //Force vector, solution
            //            coarsesol.Print(std::cout);
            TPZMGSolver<STATE> mgsolve(transfer,an_c.MatrixSolver<STATE>(),1);
            // mgsolve.SetMatrix(an_f->Solver()->Matrix());
            mgsolve.SetMatrix(an_f.MatrixSolver<STATE>().Matrix());
            //            finesol.Print(std::cout);
            mgsolve.Solve(finesol, finesol);
            //            finesol.Print(std::cout);
            //End of transferation
            
            
            clock5.stop();
            std::ofstream Out5("Transferation.txt");
            operator<<(Out5, clock5 );
            
            TPZTimer clock6;
            clock6.start();
            
            //Iterative method process
            
            TPZFMatrix<STATE> rhscoarse = an_c.Rhs();
            //            rhscoarse.Print(std::cout);
            TPZFMatrix<STATE> rhsfine = an_f.Rhs();
            //            rhsfine.Print(std::cout);
            //            rhsfine.Print("rhsfine" , std::cout);
            //            rhscoarse.Print("rhscoarse" , std::cout);
            TPZStepSolver<STATE> BDSolve(sp_auto);
            BDSolve.SetDirect(ELU);
            TPZSequenceSolver<STATE> seqsolver;
            
            // seqsolver.SetMatrix(an_f->Solver()->Matrix());
            seqsolver.SetMatrix(an_f.MatrixSolver<STATE>().Matrix());
            seqsolver.AppendSolver(mgsolve); //Updates the values of the preconditioner based on the values of the matrix
            seqsolver.AppendSolver(BDSolve);
            seqsolver.AppendSolver(mgsolve);
            
            seqsolver.Solve(rhsfine, finesol);
            //            finesol.Print(std::cout);
            
            
            //            std::ofstream file("matblock.nb");
            //            sp->Print("k = ",file,EMathematicaInput);
            TPZStepSolver<STATE> cg_solve(an_f.MatrixSolver<STATE>().Matrix());
            int maxIt = 200;
            cg_solve.SetCG(maxIt, seqsolver, 1.e-10, 0);
            // TPZStepSolver<STATE> step2(an_f.MatrixSolver<STATE>().Matrix());
            // step2.SetDirect(ELDLt);
            
            // 

            //            finesol.Print(std::cout);
            finesol.Zero();
            //            finesol.Print(std::cout);
            cg_solve.Solve(rhsfine, finesol);
            // step2.Solve(rhsfine, finesol);
            //            finesol.Print(std::cout);

            rprint << cg_solve.NumIterations() << "\n";
            
            MixedMesh_f->LoadSolution(finesol);
            clock6.stop();
            std::ofstream Out6("Iterative_Method.txt");
            operator<<(Out6, clock6 );
        }
        
        
        if(0){
            
            TPZTimer clock7;
            clock7.start();
            std::cout << "Printing result in vtk file " << std::endl;
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(vecmesh_c, MixedMesh_c);
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(vecmesh_f, MixedMesh_f);
            
            //PostProcess
            TPZStack<std::string> scalar, vectors;
            TPZManVector<std::string,15> scalnames(8), vecnames(2);
            vecnames[0]  = "Flux";
            vecnames[1]  = "ExactFlux";
            scalnames[0] = "Pressure";
            scalnames[1] = "DivFlux";
            scalnames[2] = "g_average";
            scalnames[3] = "u_average";
            scalnames[4] = "Permeability";
            scalnames[5] = "ExactPressure";
            scalnames[6] = "Divergence";
            scalnames[7] = "ExactDiv";

            std::ofstream filePrint_coarse("MixedHdiv_coarse.txt");
            MixedMesh_c->Print(filePrint_coarse);
            std::string name_coarse = "MixedHdiv_coarse.vtk";
            
            
            std::ofstream filePrint_fine("MixedHdiv_fine.txt");
            MixedMesh_f->Print(filePrint_fine);
            std::string name_fine = "MixedHdiv_fine.vtk";
            
            int di;
            if (two_d_Q) {
                di = 3;                 //Dimension definition
            } else {
                di = 1;                 //Dimension definition
            }
            an_c.DefineGraphMesh(di, scalnames, vecnames, name_coarse);
            an_c.PostProcess(0,di);
            
            an_f.DefineGraphMesh(di, scalnames, vecnames, name_fine);
            an_f.PostProcess(0,di);
            
            std::cout << "Printing result in vtk file " << std::endl;

            clock7.stop();
            std::ofstream Out7("Postprocess.txt");
            operator<<(Out7, clock7 );
            
        }
        
    }
    
    
}

/**
 * @brief Generates the geometric mesh
 * @param nx: number of partions on x
 * @param ny: number of partions on y
 * @param l: lenght
 * @param h: height
 * @return Geometric mesh
 */
TPZGeoMesh * GenerateGmesh(int nx, int ny, double l, double h){
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    
    TPZVec<int> nels(3,0);
    nels[0]=nx;         //Elements over x
    nels[1]=ny;         //Elements over y
    nels[2]=ny;         //Elements over y
    
    TPZVec<REAL> x0(3,0.0);
    TPZVec<REAL> x1(3,l);
    // x1[1]=h;
    // x1[2]=0;
    
    //Setting boundary conditions (negative numbers to recognize them)
    TPZGenGrid3D gen(x0,x1,nels,MMeshType::EHexahedral);
    // TPZGenGrid3D gen(x0,x1,nels,MMeshType::ETetrahedral);
    gmesh = gen.BuildVolumetricElements(1);
    gmesh = gen.BuildBoundaryElements(-1,-2,-3,-4,-5,-6);
    // gen.SetElementType();
    // gen.Read(gmesh);
    // gen.SetBC(gmesh, 4, -1);
    // gen.SetBC(gmesh, 5, -2);
    // gen.SetBC(gmesh, 6, -3);
    // gen.SetBC(gmesh, 7, -4);
    // gen.SetBC(gmesh, 8, -5);
    // gen.SetBC(gmesh, 9, -6);
    return gmesh;
}

/**
 * @brief Generates the pressure computational mesh
 * @param Geometric mesh
 * @param Order
 * @return Pressure computational mesh
 */
TPZCompMesh * GeneratePressureCmesh(TPZGeoMesh *Gmesh, int order){
    
    TPZCompMesh *Cmesh= new TPZCompMesh (Gmesh);
    
    Cmesh->SetDimModel(Gmesh->Dimension());
    Cmesh->SetDefaultOrder(order);
    Cmesh->SetAllCreateFunctionsDiscontinuous();
    Cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    //Add material to the mesh
    int dimen = Gmesh->Dimension();
    int MaterialId = 1;
    STATE Permeability=1;
    
    TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(MaterialId, dimen);
    
    //No convection
    REAL conv=0;
    TPZVec<REAL> convdir(dimen, 0);
    // mat->SetConstantPermeability(Permeability);
    
    //Insert material to mesh
    Cmesh->InsertMaterialObject(mat);
    
    //Autobuild
    Cmesh->AutoBuild();
    
    int ncon = Cmesh->NConnects();
    
    //Set Lagrange multiplier
    for(int i=0; i<ncon; i++){
        TPZConnect &newnod = Cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    return Cmesh;
}

/**
 * @brief Generates the constant computational mesh
 * @param Gmesh: Geometric mesh
 * @param third_LM: Bool Third Lagrange multiplier
 * @return Constant computational mesh
 */
TPZCompMesh * GenerateConstantCmesh(TPZGeoMesh *Gmesh, bool third_LM)
{
    TPZCompMesh *Cmesh= new TPZCompMesh (Gmesh);
    
    Cmesh->SetDimModel(Gmesh->Dimension());
    Cmesh->SetDefaultOrder(0);
    Cmesh->SetAllCreateFunctionsDiscontinuous();
    
    //Add material to the mesh
    int dimen = Gmesh->Dimension();
    int MaterialId = 1;
    
    TPZNullMaterial<> *mat =new TPZNullMaterial<>(MaterialId);
    mat->SetDimension(dimen);
    mat->SetNStateVariables(1);
    
    //Insert material to mesh
    Cmesh->InsertMaterialObject(mat);
    
    //Autobuild
    Cmesh->AutoBuild();
    
    int ncon = Cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = Cmesh->ConnectVec()[i];
        if (third_LM) {
            newnod.SetLagrangeMultiplier(3);
        }else{
            newnod.SetLagrangeMultiplier(2);
        }
    }
    return Cmesh;
}


/**
 * @brief Generates the flux computational mesh
 * @param mesh: Geometric mesh
 * @param order_internal: Order used for internal elements
 * @param order_border: Order used for border elements
 * @return Flux computational mesh
 */
TPZCompMesh * GenerateFluxCmesh(TPZGeoMesh *mesh, int order_internal, int order_border, bool coarse){
    
    int dimen = mesh->Dimension();
    TPZCompMesh *Cmesh = new TPZCompMesh(mesh);
    Cmesh->SetDimModel(dimen);
    Cmesh->SetDefaultOrder(order_border);
    
    //Definition of the approximation space
    int perm=1;
    REAL conv=0;
    REAL perme=1;
    TPZVec<REAL> convdir(dimen , 0.0);
    TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(perm , dimen);
    // mat->SetParameters(perme, conv, convdir);
    
    //Inserting volumetric materials objects
    Cmesh->InsertMaterialObject(mat);
    
    //Create H(div) functions
    
    Cmesh->ApproxSpace().SetHDivFamily(myScenario.HDivFam);
    Cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dimen);
    
    //Insert boundary conditions
    int D=0;
    int N=1;
    int BC1=-1;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZManVector<STATE> val2(1,0.0);
    TPZBndCondT<STATE> *bc1 = mat->CreateBC(mat, BC1, D, val1, val2);
    Cmesh->InsertMaterialObject(bc1);
    
    int BC2=-2;
    val2[0]=00;
    TPZBndCondT<STATE> *bc2 = mat->CreateBC(mat, BC2, D, val1, val2);
    Cmesh->InsertMaterialObject(bc2);
    
    int BC3=-3;
    val2[0]=0;
    TPZBndCondT<STATE> *bc3 = mat->CreateBC(mat, BC3, D, val1, val2);
    Cmesh->InsertMaterialObject(bc3);
    
    int BC4=-4;
    val2[0]=0;
    TPZBndCondT<STATE> *bc4 = mat->CreateBC(mat, BC4, D, val1, val2);
    Cmesh->InsertMaterialObject(bc4);
    
    int BC5=-5;
    val2[0]=0;
    TPZBndCondT<STATE> *bc5 = mat->CreateBC(mat, BC5, D, val1, val2);
    Cmesh->InsertMaterialObject(bc5);

    int BC6=-6;
    val2[0]=0;
    TPZBndCondT<STATE> *bc6 = mat->CreateBC(mat, BC6, D, val1, val2);
    Cmesh->InsertMaterialObject(bc6);
    
    Cmesh->AutoBuild();
    
    int64_t nel = Cmesh->NElements();
    for (int el=0; el<nel; el++) {
        TPZCompEl *cel = Cmesh->Element(el);
        if(!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        if (coarse && gel->Dimension() != Cmesh->Dimension()) {
            intel->SetSideOrder(gel->NSides()-1, myScenario.CoarseBoundaryFluxPOrder);
        } else if (coarse) {
            intel->SetSideOrder(gel->NSides()-1, myScenario.CoarseBubbleFluxPOrder);
        } else if (gel->Dimension() != Cmesh->Dimension()) {
            intel->SetSideOrder(gel->NSides()-1, myScenario.FineBoundaryFluxPOrder);
        } else {
            intel->SetSideOrder(gel->NSides()-1, myScenario.FineBubbleFluxPOrder);
        }
        
    }
    Cmesh->ExpandSolution();
    
    return Cmesh;
}

/**
 * @brief Generates the mixed computational mesh
 * @param fvecmesh: Vector thats contains flux and pressure computational mesh
 * @param order: Ordertwo_d_Q
 * @param two_d_Q: Bool wether is 1D (false) or 2D (true)
 * @return Mixed computational mesh
 */

TPZMultiphysicsCompMesh * GenerateMixedCmesh(TPZVec<TPZCompMesh *> fvecmesh, int order, bool two_d_Q){
    TPZGeoMesh *gmesh = fvecmesh[1]->Reference();
    TPZMultiphysicsCompMesh *MixedMesh = new TPZMultiphysicsCompMesh(gmesh);
    
    //Definition of the approximation space
    int dimen= gmesh->Dimension();
    int matnum=1;
    REAL perm=1;
    
    //Inserting material
    // TPZMixedDarcyWithFourSpaces * mat = new TPZMixedDarcyWithFourSpaces(matnum, dimen);
    TPZMixedDarcyFlow * mat = new TPZMixedDarcyFlow(matnum, dimen);
    mat->SetConstantPermeability(perm);
    
    if (two_d_Q) {
        // TPZFunction<STATE> sourceterm = new TPZDummyFunction<STATE>(Ladoderecho_2D, 10);
        mat->SetForcingFunction(Ladoderecho_2D,3);
    } else {
        // TPZFunction<STATE> sourceterm = new TPZDummyFunction<STATE>(Ladoderecho_1D, 10);
        mat->SetForcingFunction(Ladoderecho_1D,3);
    }
    
    //Inserting volumetric materials objects
    MixedMesh->InsertMaterialObject(mat);
    
    //Boundary conditions
    int D=0;
    int N=1;
    int BC1=-1;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZManVector<STATE> val2(1,0.0);
    
    val2[0]=0.0;
    TPZBndCondT<STATE> *bc1 = mat->CreateBC(mat, BC1, D, val1, val2);
    bc1->SetForcingFunctionBC(exactSol,3);
    MixedMesh->InsertMaterialObject(bc1);
    
    int BC2=-2;
    val2[0]=0;
    TPZBndCondT<STATE> *bc2 = mat->CreateBC(mat, BC2, D, val1, val2);
    bc2->SetForcingFunctionBC(exactSol,3);
    MixedMesh->InsertMaterialObject(bc2);
    
    int BC3=-3;
    val2[0]=0;
    TPZBndCondT<STATE> *bc3 = mat->CreateBC(mat, BC3, D, val1, val2);
    bc3->SetForcingFunctionBC(exactSol,3);
    MixedMesh->InsertMaterialObject(bc3);
    
    int BC4=-4;
    val2[0]=0;
    TPZBndCondT<STATE> *bc4 = mat->CreateBC(mat, BC4, D, val1, val2);
    bc4->SetForcingFunctionBC(exactSol,3);
    MixedMesh->InsertMaterialObject(bc4);

    int BC5=-5;
    val2[0]=0;
    TPZBndCondT<STATE> *bc5 = mat->CreateBC(mat, BC5, D, val1, val2);
    bc5->SetForcingFunctionBC(exactSol,3);
    MixedMesh->InsertMaterialObject(bc5);

    int BC6=-6;
    val2[0]=0;
    TPZBndCondT<STATE> *bc6 = mat->CreateBC(mat, BC6, D, val1, val2);
    bc6->SetForcingFunctionBC(exactSol,3);
    MixedMesh->InsertMaterialObject(bc6);
    
    MixedMesh->SetAllCreateFunctionsMultiphysicElem();
    MixedMesh->SetDimModel(dimen);
    
    //Autobuild
    TPZManVector<int,5> active_approx_spaces(4); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    active_approx_spaces[2] = 1;
    active_approx_spaces[3] = 1;
    MixedMesh->BuildMultiphysicsSpace(active_approx_spaces,fvecmesh);
    
    TPZBuildMultiphysicsMesh::AddElements(fvecmesh, MixedMesh);
    TPZBuildMultiphysicsMesh::AddConnects(fvecmesh,MixedMesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fvecmesh, MixedMesh);
    
    std::cout<<"n equ Mixed: "<<MixedMesh->NEquations()<<std::endl;
    std::cout<<"n equ Flux: "<<fvecmesh[0]->NEquations()<<std::endl;
    std::cout<<"n equ Pressure: "<<fvecmesh[1]->NEquations()<<std::endl;
    std::cout<<"n equ Constant: "<<fvecmesh[2]->NEquations()<<std::endl;
    std::cout<<"n equ Constant: "<<fvecmesh[3]->NEquations()<<std::endl;
    std::cout<<"--------------------"<<std::endl;
    return MixedMesh;
};


/**
 * @brief Generates a index vector which relates coarse and fine mesh
 * @param Coarse_sol: Coarse mesh
 * @param Fine_sol: Fine mesh
 * @param indexvec:
 */
void IndexVectorCoFi(TPZMultiphysicsCompMesh *Coarse_sol, TPZMultiphysicsCompMesh *Fine_sol, TPZVec<int64_t> & indexvec)
{
    int64_t maxcone = Coarse_sol->NConnects();
    
    int64_t indexvecsize = 0;
    for (int j=0; j<maxcone; j++) {
        bool is_condensed = Coarse_sol->ConnectVec()[j].IsCondensed();
        if (is_condensed == true ) continue;
        
        int64_t sequence_coarse = Coarse_sol->ConnectVec()[j].SequenceNumber();
        int blocksize_coarse = Coarse_sol->Block().Size(sequence_coarse);
        
        indexvecsize += blocksize_coarse;
    }
    indexvec.Resize(indexvecsize,-1);
    
    for (int j=0; j<maxcone; j++) {
        bool is_condensed = Coarse_sol->ConnectVec()[j].IsCondensed();
        if (is_condensed == true ) continue;
        
        int64_t sequence_coarse = Coarse_sol->ConnectVec()[j].SequenceNumber();
        int64_t sequence_fine = Fine_sol->ConnectVec()[j].SequenceNumber();
        int blocksize_coarse = Coarse_sol->Block().Size(sequence_coarse);
        int blocksize_fine = Fine_sol->Block().Size(sequence_fine);
        
        if (blocksize_coarse > blocksize_fine) {
            DebugStop();
            
        }
        
        for(int i=0; i<blocksize_coarse; i++){
            int64_t pos_coarse = Coarse_sol->Block().Position(sequence_coarse);
            int64_t pos_fine = Fine_sol->Block().Position(sequence_fine);
            indexvec[pos_coarse+i] = pos_fine+i;
        }
    }
}

/**
 * @brief Transfer DOF from coarse to fine mesh
 * @param CoarseDoF: Solution coarse matrix
 * @param FineDoF: Solution fine matrix
 * @param DoFIndexes: DOF index vector
 */
void TransferDegreeOfFreedom(TPZFMatrix<STATE> & CoarseDoF, TPZFMatrix<STATE> & FineDoF, TPZVec<int64_t> & DoFIndexes){
    
    int64_t n_data = DoFIndexes.size();
    for (int64_t i = 0 ; i < n_data; i++) {
        FineDoF(DoFIndexes[i],0) = CoarseDoF(i,0);
    }
    
}

/**
 * @brief Generates a geometric 1D mesh
 * @param nx: number of partions on x
 * @param l: lenght
 * @return Geometric 1D mesh
 */
TPZGeoMesh * GenerateGmeshOne(int nx, double l){
    //Creates vector nodes
    double h = l/nx;
    int Domain_Mat_Id = 1;
    int Inlet_bc_Id = -1;
    int Outletbc_Id = -2;
    TPZVec<REAL> xp(3,0.0);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nx+1);
    for (int64_t i=0; i<nx+1; i++) {
        xp[0] =(i)*h;
        gmesh->NodeVec()[i]= TPZGeoNode(i, xp, *gmesh);
    }
    
    
    //Creates elements
    TPZVec<int64_t> cornerindexes(2);
    for (int64_t iel=0; iel<nx; iel++) {
        cornerindexes[0]=iel;
        cornerindexes[1]=iel+1;
        gmesh->CreateGeoElement(EOned, cornerindexes, Domain_Mat_Id, iel);
    }
    //Set BC
    gmesh->Element(0)->CreateBCGeoEl(0, Inlet_bc_Id);
    gmesh->Element(nx-1)->CreateBCGeoEl(1, Outletbc_Id);
    gmesh->SetDimension(1);
    gmesh->BuildConnectivity();
    
    return gmesh;
}

/**
 * @brief Configurate the matrix for analysis
 * @param cmesh_c: Computational coarse mesh
 * @param cmesh_f: Computational fine mesh
 * @param must_opt_band_width_Q: Wether the band width is optimized or not (it is neccesaty to rearrange the matrix in order to be less sparse)
 * @param number_threads:
 * @param an_c: Coarse analysis mesh
 * @param an_f: Fine analysis mesh
 * @param UsePardiso_Q: Wether using Pardiso or not
 */
void ConfigurateAnalyses(TPZCompMesh * cmesh_c, TPZCompMesh * cmesh_f, bool must_opt_band_width_Q, int number_threads, TPZAnalysis &an_c,TPZAnalysis &an_f, bool UsePardiso_Q){
    
    // an_c->SetCompMesh(cmesh_c,must_opt_band_width_Q);
    // an_f->SetCompMesh(cmesh_f,must_opt_band_width_Q);
    TPZStepSolver<STATE> step;
    if (UsePardiso_Q) {
        
        TPZSSpStructMatrix<STATE> sparse_matrix_coarse(cmesh_c);
        TPZSSpStructMatrix<STATE> sparse_matrix_fine(cmesh_f);
        sparse_matrix_coarse.SetNumThreads(number_threads);
        sparse_matrix_fine.SetNumThreads(number_threads);
        an_c.SetStructuralMatrix(sparse_matrix_coarse);
        an_f.SetStructuralMatrix(sparse_matrix_fine);
        
    }else{
        
        TPZSkylineStructMatrix<STATE> sparse_matrix_coarse(cmesh_c);
        TPZSkylineStructMatrix<STATE> sparse_matrix_fine(cmesh_f);
        sparse_matrix_coarse.SetNumThreads(number_threads);
        sparse_matrix_fine.SetNumThreads(number_threads);
        an_c.SetStructuralMatrix(sparse_matrix_coarse);
        an_f.SetStructuralMatrix(sparse_matrix_fine);
        
    }
    step.SetDirect(ELDLt);
    an_c.SetSolver(step);
    an_f.SetSolver(step);

}

/**
 * @brief Shows shape functions of a certain element
 * @param cmesh: Computational mesh
 * @param element: Element number to show it shapes functions
 * @param funcion: Shape function number to show
 * @return plotfile: Name for the plot to print
 */
void ShowShape(TPZCompMesh * cmesh, int element, int funcion, std::string plotfile){
    TPZGeoMesh *gmesh = cmesh->Reference();
    gmesh->ResetReference();
    
    int nels = cmesh->NElements();
    for (int iel =0 ; iel<nels; iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        TPZGeoEl *gel = cel->Reference();
        if (gel->Index() == element) {
            gel->SetMaterialId(100);
        }
    }
    
    TPZMaterial  * mat_2(cmesh->MaterialVec()[1]);
    std::map<int, TPZMaterial *> matvec;
    mat_2->Clone(matvec);
    cmesh->CleanUp();
    
    TPZMixedDarcyFlow *mat = new TPZMixedDarcyFlow(100 , 2);
    TPZMaterial *perf_mat( matvec[1]);
    TPZMixedDarcyFlow *aux_mat = dynamic_cast<TPZMixedDarcyFlow *>(perf_mat);
    
    cmesh->InsertMaterialObject(mat);
    cmesh->AutoBuild();
    
    int gdl = cmesh->Solution().Rows();
    int nfuncols = cmesh->Solution().Cols();
    int nfunrows = cmesh->Solution().Rows();
    TPZFMatrix<STATE> solu(nfunrows,nfuncols,1.0);
    TPZFMatrix<STATE> sol(gdl,1,0.0);
    
    sol.Resize(gdl, 1);
    int nel = cmesh->NElements();
    int val =1;
    for (int i=0; i<nel; i++) {
        TPZCompEl *cel = cmesh->Element(i);
        TPZGeoEl *gel =cel->Reference();
        int index = gel->Index();
        if (index != element) {
            continue;
        }
        
        int nconnects = cel->NConnects();
        int acum =0;
        
        for(int j=0; j<nconnects; j++){
            int sec_number = cel->Connect(j).fSequenceNumber;
            if(sec_number > -1)
            {
                acum = acum + cmesh->Block().Size(sec_number);
                if (acum > funcion) {
                    
                    //                    if (val) {
                    
                    //SplitConnects(cmesh, gel ,j);
                    //                        cmesh->CleanUp();
                    val=0;
                    
                    int delta = funcion + cmesh->Block().Size(sec_number) - acum;
                    int64_t pos = cmesh->Block().Position(sec_number);
                    pos = pos +  delta;
                    sol.Zero();
                    sol(pos,0)=1.0;
                    cmesh->LoadSolution(sol);
                    
                    TPZLinearAnalysis *an = new TPZLinearAnalysis(cmesh,false);
                    {
                        const int dim = an->Mesh()->Dimension();
                        int div = 5;
                        // std::string plotfile = "SHAPES.vtk";
                        TPZStack<std::string> scalar_names, vec_names;
                        
                        scalar_names.push_back("Solution");
                        an->DefineGraphMesh(dim,scalar_names,vec_names,plotfile);
                        an->PostProcess(div,dim);
                        std::cout << "The function :"<<funcion<<" of element "<<element<<" has been printed on .vtk file" << std::endl;
                        
                    }
                    int nelss = gmesh->NElements();
                    for (int iel =0 ; iel<nelss; iel++) {
                        TPZGeoEl *gel = gmesh->Element(iel);
                        if (gel->Index() == element) {
                            gel->SetMaterialId(1);
                        }
                        cmesh->CleanUp();
                        cmesh->InsertMaterialObject(aux_mat);
                        cmesh->AutoBuild();
                        
                        return;
                        
                    }
                }
            }
        }
    }
    
    
}
void SplitConnects(TPZCompMesh *fluxmesh, TPZGeoEl *gel , int j){
    TPZCompEl *cel = gel->Reference();
    int nelconnect = cel->Connect(j).NElConnected();
    if (nelconnect>1) {
        
    }
    
    int iside = j + 4;
    if (iside >= gel->NSides()-1) {
        return;
    }
    TPZGeoElSide neighb(gel, iside);
    TPZStack<TPZGeoElSide> allneigh;
    TPZGeoElSide neig = neighb.Neighbour();
    if (neig.Element()->Dimension() !=2) {
        return;
    }
    
    TPZCompElSide left(cel,iside);
    TPZGeoElSide gelr = gel->Neighbour(iside);
    TPZCompElSide right(gelr.Element()->Reference(), gelr.Side());
    
    TPZGeoElSide gleft(gel,iside);
    TPZGeoElSide gright(gleft.Neighbour());
    TPZInterpolatedElement *intelleft = dynamic_cast<TPZInterpolatedElement *> (left.Element());
    TPZInterpolatedElement *intelright = dynamic_cast<TPZInterpolatedElement *> (right.Element());
    intelleft->SetSideOrient(left.Side(), 1);
    intelright->SetSideOrient(right.Side(), 1);
    TPZStack<TPZCompElSide> equalright;
    TPZConnect &cleft = intelleft->SideConnect(0, left.Side());
    int conr = right.Element()->Connect(right.Side()-4).SequenceNumber();
    int conl = left.Element()->Connect(left.Side()-4).SequenceNumber();
    
    if (conr != conl) {
        return;
    }
    
    
    int64_t index = fluxmesh->AllocateNewConnect(cleft);
    TPZConnect &newcon = fluxmesh->ConnectVec()[index];
    cleft.DecrementElConnected();
    newcon.ResetElConnected();
    newcon.IncrementElConnected();
    newcon.SetSequenceNumber(fluxmesh->NConnects() - 1);
    
    int rightlocindex = intelright->SideConnectLocId(0, right.Side());
    intelright->SetConnectIndex(rightlocindex, index);
    
    int sideorder = cleft.Order();
    fluxmesh->SetDefaultOrder(sideorder);
    
    int gdl = fluxmesh->Solution().Rows();
    fluxmesh->Solution().Resize(gdl+2, 1);
    TPZFMatrix<STATE> sol(gdl+2,1,0.0);
    fluxmesh->LoadSolution(sol);
    
    
    return;
}

void HdiVSimple(int nx, int order_high, bool condense_equations_Q, bool two_d_Q){
    
    bool KeepOneLagrangian = false;
    bool KeepMatrix = false;
    
    TPZGeoMesh *gmesh;
    
    // Generates a 2D geo mesh
    TPZGeoMesh *gmesh_1 = GenerateGmesh(nx, nx, 1, 1);
    
    // Generates a 1D geo mesh
    TPZGeoMesh *gmesh_2 = GenerateGmeshOne(nx, 1);
    
    //Asks if the problem is 2D
    if (two_d_Q) {
        gmesh = gmesh_1;
    } else {
        gmesh = gmesh_2;
    }
    
    TPZMultiphysicsCompMesh *MixedMesh_c = 0;
    TPZManVector<TPZCompMesh *> vecmesh_c(2);      //Vector for coarse mesh case (4 spaces)
    {
        TPZCompMesh *q_cmesh = GenerateFluxCmesh(gmesh, order_high, order_high);
        TPZCompMesh *p_cmesh = GeneratePressureCmesh(gmesh, order_high);
        
        vecmesh_c[0] = q_cmesh;              //Flux
        vecmesh_c[1] = p_cmesh;              //Pressure
        
        MixedMesh_c = GenerateMixedCmesh(vecmesh_c, 1, two_d_Q);       //1 Stands for the corse mesh order
    }
    
    if (condense_equations_Q) {             //Asks if you want to condesate the problem
        MixedMesh_c->ComputeNodElCon();
        int dim = MixedMesh_c->Dimension();
        int64_t nel = MixedMesh_c->NElements();
        for (int64_t el =0; el<nel; el++) {
            TPZCompEl *cel = MixedMesh_c->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            if(gel->Dimension() != dim) continue;
            int nc = cel->NConnects();
            cel->Connect(nc-1).IncrementElConnected();
        }
        // Created condensed elements for the elements that have internal nodes
        TPZCompMeshTools::CreatedCondensedElements(MixedMesh_c, KeepOneLagrangian, KeepMatrix);
    }
    
    
    //Solving the system:
    MixedMesh_c->InitializeBlock();    //Resequence the block object, remove unconnected connect objects
    TPZLinearAnalysis *an_c = new TPZLinearAnalysis;
    
    TPZTimer AssemblyFine;
    AssemblyFine.start();
    
    // Assembly fine operator
    an_c->Assemble();
    AssemblyFine.stop();
    std::ofstream Out_AssemblyFine("Assembly.txt");
    operator<<(Out_AssemblyFine, AssemblyFine );
    
    
    
    
    if(1){
        
        TPZTimer clock7;
        clock7.start();
        
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(vecmesh_c, MixedMesh_c);
        
        //PostProcess
        TPZStack<std::string> scalar, vectors;
        TPZManVector<std::string,10> scalnames(4), vecnames(1);
        vecnames[0]  = "q";
        scalnames[0] = "p";
        scalnames[1] = "kappa";
        scalnames[1] = "div_q";
        scalnames[2] = "g_average";
        scalnames[3] = "u_average";
        
        std::ofstream filePrint_coarse("MixedHdiv_coarse.txt");
        MixedMesh_c->Print(filePrint_coarse);
        std::string name_coarse = "MixedHdiv_coarse.vtk";
        
        
        int di;
        if (two_d_Q) {
            di = 2;                 //Dimension definition
        } else {
            di = 1;                 //Dimension definition
        }
        an_c->DefineGraphMesh(di, scalnames, vecnames, name_coarse);
        an_c->PostProcess(0,di);
        
        
        clock7.stop();
        std::ofstream Out7("Postprocess.txt");
        operator<<(Out7, clock7 );
        
    }
    
}


void AssociateElements(TPZMultiphysicsCompMesh *cmesh, TPZVec<int64_t> &elementgroup, std::set<int> &matId)
{
    // for (auto i:matId)
    // {
    //     std::cout << " " << i << "\n";
    // }   
    
    int64_t nel = cmesh->NElements();
    elementgroup.Resize(nel, -1);
    elementgroup.Fill(-1);
    
    int64_t nconnects = cmesh->NConnects();
    TPZVec<int64_t> groupindex(nconnects, -1);
    int dim = cmesh->Dimension();
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        cel->LoadElementReference();
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != dim) {
            continue;
        }
        
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);

        int nfacets = cel->Reference()->NSides(cel->Dimension()-1);
        
        int k = -1;
        auto nconnects = connectlist.size();

        for (int i=0; i<nfacets; i++) {
            int cindex = connectlist[i];

            k++;
            auto gel = cel->Reference();
            auto nnodes = gel->NCornerNodes();
            auto nsides = gel->NSides();
            TPZGeoElSide geoside(gel,nsides-nfacets-1+k);
            auto neig = geoside.Neighbour(); 
            if (!neig.Element())continue;
            int neigMatId = neig.Element()->MaterialId();

            if (matId.find(neigMatId) == matId.end()) {
                continue;
            }
            
            if (groupindex[cindex] == -1) {
                groupindex[cindex] = cel->Index();
            }
        }
    }

    // std::cout << "Groups of connects " << groupindex << std::endl;
    // std::cout << "Groups of connects " << groupindex2 << std::endl;
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);

        int64_t groupfound = -1;
        int k = -1;
        for (auto cindex : connectlist) {           
            if (groupindex[cindex] != -1) {
                // assign the element to the group
                if(groupfound != -1 && groupfound != groupindex[cindex])
                {
                    //Do nothing
                }else{
                    elementgroup[cel->Index()] = groupindex[cindex];
                    groupfound = groupindex[cindex];
                }
            }
        }
    }
    // std::cout << "Element group = " << elementgroup << std::endl;
}


void CondenseBCElements(TPZMultiphysicsCompMesh *cmesh, std::set<int> &matIdBC){

    int64_t nel = cmesh->NElements();
    TPZVec<int64_t> groupnumber(nel,-1);

    auto aux=matIdBC;
    /// compute a groupnumber associated with each element
    AssociateElements(cmesh,groupnumber,aux);

    std::map<int64_t, TPZElementGroup *> groupmap;
    //    std::cout << "Groups of connects " << groupnumber << std::endl;
    for (int64_t el = 0; el<nel; el++) {
        int64_t groupnum = groupnumber[el];

        if(groupnum == -1) continue;
        auto iter = groupmap.find(groupnum);
        if (groupmap.find(groupnum) == groupmap.end()) {
            TPZElementGroup *elgr = new TPZElementGroup(*cmesh);
            groupmap[groupnum] = elgr;
            elgr->AddElement(cmesh->Element(el));
        }
        else
        {
            iter->second->AddElement(cmesh->Element(el));
        }
    }

    cmesh->ComputeNodElCon();
    // increment elconnect em um elemento de pressao 
    // Compmesh tools condense equations
    nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
        if (elgr) {
            TPZCondensedCompEl *cond = new TPZCondensedCompEl(elgr);
            cond->SetKeepMatrix(false);
        }
    }

    cmesh->InitializeBlock();
    cmesh->ComputeNodElCon();

}