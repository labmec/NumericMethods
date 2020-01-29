//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "ConfigurateCase.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"

SimulationCase SimulationCase2DMHM();
SimulationCase SimulationCase2D();

int main(){
    
    ConfigurateCase confi;
    confi.SetSimulationCase(SimulationCase2DMHM());
//    confi.SetFineOrder(1);
//    confi.SetCoarseOrder(1);
    
    TPZAutoPointer<TPZMHMixedMesh4SpacesControl> mhmcontrol = confi.CreateMHMMixedMesh4Spaces();
    
    TPZCompMesh* MHMIxed= mhmcontrol->CMesh().operator->();
    SimulationCase sim;
   
    if (0) {
        std::ofstream filefinal("MHM_After.txt");
        MHMIxed->Print(filefinal);
        std::cout<<MHMIxed->NEquations()<<std::endl;
        std::ofstream fileflux("fluxmhm_After.txt");
        std::ofstream filepressure("pressuremhm_After.txt");
        std::ofstream filefluxavg("fluxmhmavg_After.txt");
        std::ofstream filepressureavg("pressuremhmavg_After.txt");
        mhmcontrol->GetMeshes()[0]->Print(fileflux);
        mhmcontrol->GetMeshes()[1]->Print(filepressure);
        mhmcontrol->GetMeshes()[2]->Print(filefluxavg);
        mhmcontrol->GetMeshes()[3]->Print(filepressureavg);
    }

    bool shouldrenumber = true;
    TPZAnalysis an_coarse(MHMIxed,shouldrenumber);
    
    TPZSymetricSpStructMatrix strmat(MHMIxed);
    strmat.SetNumThreads(confi.GetSimulationCase().n_threads);
    an_coarse.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an_coarse.SetSolver(step);
    std::cout << "Assembling\n";
    an_coarse.Assemble();
    std::ofstream filemate("MatrixCoarse.txt");
    an_coarse.Solver().Matrix()->Print("EkRs",filemate,EMathematicaInput);
    
    std::cout << "Solving\n";
    an_coarse.Solve();
    std::cout << "Finished\n";
    an_coarse.LoadSolution(); // Compute internal dofs
    
    
    // PostProcess
    TPZStack<std::string> scalar, vectors;
    TPZManVector<std::string,10> scalnames(4), vecnames(1);
    vecnames[0]  = "q";             //Flux
    scalnames[0] = "p";             //Pressure
    scalnames[1] = "kappa";         //Permeability
    scalnames[1] = "div_q";         //Flux Divergent
    scalnames[2] = "g_average";     //Distributed flux
    scalnames[3] = "u_average";     //Average pressure
    
    
    
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(MHMIxed->Reference(), out, true);
    }
    
    std::ofstream filefinal2("mhmfinal.txt");
    MHMIxed->Print(filefinal2);
    std::cout<<MHMIxed->NEquations()<<std::endl;
    
    
    std::string name_coarse("Simulation_Results.vtk");
    
    an_coarse.DefineGraphMesh(2, scalnames, vecnames, name_coarse);
    an_coarse.PostProcess(0,2);
    
    return 0;
};

SimulationCase SimulationCase2DMHM(){
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.KeepMatrixQ =false;
    sim.KeepOneLagrangianQ = false;
    sim.IsCondensedQ = true;
    sim.n_threads = 0;
    sim.omega_ids.push_back(1);
    sim.omega_dim.push_back(2);
    sim.omega_ids.push_back(2);
    sim.omega_dim.push_back(2);
    
    int bc_non_flux = -1;
    int bc_non_flux2 = -3;
    int bc_non_flux3 = -2;
    int bc_non_flux4 = -4;
    int bc_inlet  = -5;
    int bc_outlet = -6;
    
    sim.gamma_ids.push_back(bc_non_flux);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_non_flux2);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_non_flux3);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_non_flux4);
    sim.gamma_dim.push_back(1);
    
    sim.gamma_ids.push_back(bc_inlet);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_outlet);
    sim.gamma_dim.push_back(1);
    
    int bc_type_D = 0;
    int bc_type_N = 1;
    REAL p_inlet  = 5000.0;
    REAL p_outlet = 1.0;
    REAL qn       = 0.0;
    
    sim.type.push_back(bc_type_N);
    sim.type.push_back(bc_type_N);
    sim.type.push_back(bc_type_N);
    sim.type.push_back(bc_type_N);
    
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_D);
    
    sim.vals.push_back(qn);
    sim.vals.push_back(qn);
    sim.vals.push_back(qn);
    sim.vals.push_back(qn);
    sim.vals.push_back(p_inlet);
    sim.vals.push_back(p_outlet);
    
    return sim;
};


SimulationCase SimulationCase2D(){
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.KeepMatrixQ =false;
    sim.KeepOneLagrangianQ = false;
    sim.IsCondensedQ = true;
    sim.n_threads = 24;
    sim.omega_ids.push_back(1);
    sim.omega_dim.push_back(2);
    sim.permeabilities.push_back(1.0);
    
    
    int bc_non_flux = -1;
    int bc_inlet  = -2;
    int bc_non_flux2 = -3;
    int bc_outlet = -4;
    
    
    sim.gamma_ids.push_back(bc_non_flux);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_inlet);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_non_flux2);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_outlet);
    sim.gamma_dim.push_back(1);
    
    int bc_type_D = 0;    //    D = 0;
    int bc_type_N = 1;    //    N = 1;
    REAL p_inlet  = 0.0;
    REAL p_outlet = 0.0;
    REAL qn       = 0.0;
    
    sim.type.push_back(bc_type_N);
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_N);
    sim.type.push_back(bc_type_D);
    
    sim.vals.push_back(qn);
    sim.vals.push_back(p_inlet);
    sim.vals.push_back(qn);
    sim.vals.push_back(p_outlet);
    
    return sim;
};
