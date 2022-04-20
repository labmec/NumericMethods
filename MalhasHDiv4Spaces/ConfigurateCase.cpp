//
//  ConfigurateCase.h
//  Numeric Methods
//
//  Created by Jorge Paúl Ordóñez Andrade on 10/11/19.

#include "ConfigurateCase.h"
#include "TPZRefPattern.h"
#include "TPZRefPatternTools.h"
#include "pzrefquad.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMHMixedMesh4SpacesControl.h"
#include "TPZMHMeshControl.h"
#include "pzvec.h"
#include "TPZPrintUtils.h"
#include <string>

void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);
void DivideMesh(TPZGeoMesh * gmesh);

ConfigurateCase::ConfigurateCase(){
    
}
void PrintPoint(TPZVec<double> point){
    std::cout<<point[0]<<std::endl;
    std::cout<<point[1]<<std::endl;
    std::cout<<point[2]<<std::endl;
}

/**
 * @brief Creates the computatinal Hdiv mesh
 * @param gmesh: geometric mesh
 * @param orderfine: order for the fine mesh
 * @param ordercoarse: order fot the coarse mesh
 * @return computatinal Hdiv mesh
 */
TPZCompMesh *ConfigurateCase::HDivMesh(TPZGeoMesh * gmesh, int orderfine, int ordercoarse){
    
    int dimension = gmesh->Dimension();         //Gets mesh dimension
    int nvols = fsim_case.omega_ids.size();     //Gets the number of materials id
    int nbound = fsim_case.gamma_ids.size();    //Gets the number of BC
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh); //Creates a computational mesh
    
    cmesh->SetDefaultOrder(orderfine);        //Sets a default order
    TPZFMatrix<STATE> val1(dimension,dimension,0.0);
    TPZManVector<STATE> val2(dimension,0.0);
    
    for (int ivol=0; ivol<nvols; ivol++) {
        //Sets data needed to materials
        TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(fsim_case.omega_ids[ivol],fsim_case.omega_dim[ivol]);
        volume->SetConstantPermeability(fsim_case.permeabilities[ivol]);
        
        cmesh->InsertMaterialObject(volume);
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {
                val2[0]=fsim_case.vals[ibound];
                int condType=fsim_case.type[ibound];
                //Creates BC conditions and insert them to the mesh
                TPZBndCondT<STATE> * face = volume->CreateBC(volume,fsim_case.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
            }
        }
    }

    
    cmesh->SetAllCreateFunctionsHDiv();

    cmesh->AutoBuild();
    std::stringstream file_name;
    file_name << "q_cmesh_DELETE" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
    
        int64_t nel = cmesh->NElements();
        for (int el=0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            if(!cel) continue;
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            if(!intel){
                std::cout<<"Interpolated element not found"<<std::endl;
                DebugStop();
            };
            TPZGeoEl *gel = intel->Reference();
            intel->SetSideOrder(gel->NSides()-1, orderfine);
        }
        cmesh->ExpandSolution();
    
    cmesh->InitializeBlock();
    
#ifdef PZDEBUG
//    std::stringstream file_name;
//    file_name << "q_cmesh_raw" << ".txt";
//    std::ofstream sout(file_name.str().c_str());
//    cmesh->Print(sout);
#endif
    
    return cmesh;
}

/**
 * @brief Creates the discontinous mesh
 * @param gmesh: geometric mesh
 * @param order: order for the mesh
 * @param lagrangemult: lagrange multiplier number
 * @return computatinal discontinuos mesh
 */
TPZCompMesh *ConfigurateCase::DiscontinuousMesh(TPZGeoMesh * gmesh, int order, int lagrangemult){
    int dim = gmesh->Dimension();
    int nvols = fsim_case.omega_ids.size();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZFMatrix<STATE> val1(dim,dim,0.0),val2(dim,1,0.0);
    
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    //Creating approximation space for flux mesh
    if (lagrangemult == 1) {
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    
    //Setting material for avg pressure and distributed flux
     std::set<int> matids;
    if (order == 0) {
        for (int ivol=0; ivol < nvols; ivol++) {
            TPZNullMaterial<> * volume = new TPZNullMaterial<>(fsim_case.omega_ids[ivol]);
            cmesh->InsertMaterialObject(volume);
            if (fsim_case.omega_dim[ivol] == dim) {
                matids.insert(fsim_case.omega_ids[ivol]);
            }
        }
    }
    //Setting material for pressure

    if (order != 0) {
        
        for (int ivol=0; ivol < nvols; ivol++) {
            TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(fsim_case.omega_ids[ivol],dim);
            volume->SetConstantPermeability(fsim_case.permeabilities[ivol]);
            cmesh->InsertMaterialObject(volume);
            if (fsim_case.omega_dim[ivol] == dim) {
                matids.insert(fsim_case.omega_ids[ivol]);
            }
        }
    }
    
    cmesh->AutoBuild(matids);
    cmesh->InitializeBlock();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(lagrangemult);
    }
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name << "p_cmesh_raw" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
};

/**
 * @brief Creates the multiphysics computational mesh
 * @param gmesh: geometric mesh
 * @param order: order for the mesh
 * @param lagrangemult: lagrange multiplier number
 * @return computatinal discontinuos mesh
 */
TPZMultiphysicsCompMesh *ConfigurateCase::CreateMultCompMesh(){
    TPZVec<TPZCompMesh *> meshvec(4);
    
    
    m_gmesh=CreateGeowithRefPattern();
    int orderfine = m_fineorder;
    int ordercoarse = m_coarseorder;
    int dimension = m_gmesh->Dimension();
    int nvols = fsim_case.omega_ids.size();     //Materials
    int nbound= fsim_case.gamma_ids.size();     //BC
    bool useSubstructure_Q = true;
    
    if (nvols<1) {
        std::cout<<"Error: Omega is not defined."<<std::endl;
        DebugStop();
    }
    if (nbound<1) {
        std::cout<<"Error: Gamma is not defined."<<std::endl;
        DebugStop();
    }
    
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(m_gmesh);
    TPZFNMatrix<9,STATE> val1(dimension,dimension,0.0);
    TPZManVector<STATE> val2(dimension,0.0);
    
    for (int ivol=0; ivol<nvols; ivol++) {
        
        TPZMixedDarcyWithFourSpaces * volume = new TPZMixedDarcyWithFourSpaces(fsim_case.omega_ids[ivol],dimension);
        volume->SetConstantPermeability(fsim_case.permeabilities[ivol]);
        cmesh->InsertMaterialObject(volume);

        
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {
                val2[0]=fsim_case.vals[ibound];
                int condType=fsim_case.type[ibound];
                TPZBndCondT<STATE> * face = volume->CreateBC(volume,fsim_case.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
            }
        }
    }
    
    cmesh->SetDimModel(dimension);
    
    meshvec[0] = HDivMesh(m_gmesh, orderfine, ordercoarse);     //Flux
    meshvec[1] = DiscontinuousMesh(m_gmesh, ordercoarse, 1);    //Pressure
    meshvec[2] = DiscontinuousMesh(m_gmesh, 0, 2);              //Avg Pressure
    meshvec[3] = DiscontinuousMesh(m_gmesh, 0, 3);              //Distributed flux
    TPZManVector<int,5> active_approx_spaces(4,1); 
  
    TPZPrintUtils utils;
    utils.PrintCompMesh(meshvec[0],"mesh0");
    utils.PrintCompMesh(meshvec[1],"mesh1");
    utils.PrintCompMesh(meshvec[2],"mesh2");
    utils.PrintCompMesh(meshvec[3],"mesh3");

    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->BuildMultiphysicsSpace(active_approx_spaces,meshvec);

    TPZMHMixedMesh4SpacesControl control(m_gmesh);
    control.SubStructure();
    if (useSubstructure_Q) {
        HideTheElements(cmesh);
    }
    
    if (fsim_case.IsCondensedQ){      //Asks if you want to condesate the problem
        cmesh->ComputeNodElCon();
        int dim = cmesh->Dimension();
        int64_t nel = cmesh->NElements();
        for (int64_t el =0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            if(gel->Dimension() != dim) continue;
            int nc = cel->NConnects();
            cel->Connect(nc-1).IncrementElConnected();  //Increment avg pressure connect in order to not condense
        }
        
        TPZCompMeshTools::CreatedCondensedElements(cmesh, fsim_case.KeepOneLagrangianQ, fsim_case.KeepMatrixQ);
        
        if (1) {
            std::ofstream filePrint("MixedHdiv.txt");
            cmesh->Print(filePrint);
        }
    }
    
    GroupAndCondense(cmesh);
    cmesh->InitializeBlock();
    
    return cmesh;
    
};

/**
 * @brief Creates the analysis for the problem
 * @param mcmesh: multiphysics computational mesh
 * @return The analysis
 */
TPZAnalysis *ConfigurateCase::CreateAnalysis(TPZMultiphysicsCompMesh *mcmesh){
    
    TPZLinearAnalysis *an = new TPZLinearAnalysis(mcmesh);
    
    TPZStepSolver<STATE> step;
    if (fsim_case.UsePardisoQ) {
        
        TPZSSpStructMatrix<STATE> sparse_matrix(mcmesh);
        sparse_matrix.SetNumThreads(fsim_case.n_threads);
        an->SetStructuralMatrix(sparse_matrix);
        
    }else{
        
        TPZSkylineStructMatrix<STATE> sparse_matrix(mcmesh);
        sparse_matrix.SetNumThreads(fsim_case.n_threads);
        an->SetStructuralMatrix(sparse_matrix);
        
    }
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    
    return an;
}

/**
 * @brief Creates a uniform mesh
 * @param nx: partitions number on x
 * @param L: length
 * @param ny: partitions number on y
 * @param h: height
 * @param nz:partitions number on z
 * @param w: width
 * @return The analysis
 */
TPZGeoMesh * ConfigurateCase::CreateUniformMesh(int nx, REAL L, int ny, REAL h, int nz, REAL w){
    
    TPZVec<int> nels;   //Creating the elements object
    nels.Resize(3);
    nels[0]=nx;         //Elements over x
    nels[1]=ny;         //Elements over y
    
    TPZVec<REAL> x0(3,0.0);     //Creating coordinates point x0 (0,0,0)
    TPZVec<REAL> x1(3,0.0);     //Creating coordinates point x1 (0,0,0)
    x1[0]=L;                    //Changing the x coordinate of x1 (L,0,0) to later create the box
    
    if (ny!=0) {        //If ny!=0 means that is a 2DGeoMesh then x1 (L,h,0)
        x1[1]=h;
    }
    TPZGeoMesh *gmesh = new TPZGeoMesh;         //Creates the geometric mesh
    TPZGenGrid2D gen(nels,x0,x1);             
    gen.SetRefpatternElements(true);
    gen.Read(gmesh);                            //Add elements and nodes, created with defaut matid=1
    
    if (nz!=0 ) {        //If nz!=0 means that is a 3DGeoMesh then x1 (L,h,0 to var)
        double var = w/nz;
        TPZExtendGridDimension extend(gmesh,var);
        extend.SetElType(1);            // 1 for RefPattern and 0 for UniformRefinement
        gmesh = extend.ExtendedMesh(nz);
    }
    
    if (nz!=0) {                                //Setting boundary conditions for a 3D Case
        for (auto gel:gmesh->ElementVec()) {
            TPZFMatrix<REAL> coordinates;
            gel->NodesCoordinates(coordinates);
            if(coordinates(2,0)==0){
                gel->CreateBCGeoEl(20, -1);
            }
            if(coordinates(2,4)==w){
                gel->CreateBCGeoEl(25, -2);
            }
            
            if(coordinates(0,0)==0.0 ){
                gel->CreateBCGeoEl(24, -3);
            }
            if(coordinates(1,0)==0.0 ){
                gel->CreateBCGeoEl(21, -4);
            }
            
            if(coordinates(0,1)== L ){
                gel->CreateBCGeoEl(22, -5);
            }
            if(coordinates(1,3)==h){
                gel->CreateBCGeoEl(23, -6);
            }
        };
        gmesh->SetDimension(3);             //After creating the BC sets a logical dim=3
    }
    
    if (ny!=0 && nz==0) {                   //If 2D sets the BC (commented because the need of other refinement)
//        gen.SetBC(gmesh, 4, -1);
//        gen.SetBC(gmesh, 5, -2);
//        gen.SetBC(gmesh, 6, -3);
//        gen.SetBC(gmesh, 7, -4);
        gmesh->SetDimension(2);             //After creating the BC sets a logical dim=2
    }
  
    if (ny==0 && nz==0) {               //If 2D sets the matid of the BC
        double dh = L/nx;
        int Domain_Mat_Id = 1;
        int Inlet_bc_Id = -1;
        int Outletbc_Id = -2;
        TPZVec<REAL> xp(3,0.0);         //Creates vector of nodes
        
        gmesh->NodeVec().Resize(nx+1);
        for (int64_t i=0; i<nx+1; i++) {
            xp[0] =(i)*dh;
            gmesh->NodeVec()[i]= TPZGeoNode(i, xp, *gmesh);
        }
        
        TPZVec<int64_t> cornerindexes(2);   //Creates elements
        for (int64_t iel=0; iel<nx; iel++) {
            cornerindexes[0]=iel;
            cornerindexes[1]=iel+1;
            gmesh->CreateGeoElement(EOned, cornerindexes, Domain_Mat_Id, iel);
        }
        gmesh->Element(0)->CreateBCGeoEl(0, Inlet_bc_Id);     //Sets BC
        gmesh->Element(nx-1)->CreateBCGeoEl(1, Outletbc_Id);
        gmesh->SetDimension(1);
        gmesh->BuildConnectivity();
    }
    
    gmesh->BuildConnectivity();
    m_gmesh = gmesh;
    return gmesh;
    
}

/**
 * @brief Creates a MHM Mixed Mesh with 2 Spaces
 * @return The analysis
 */
TPZAutoPointer<TPZMHMixedMeshControl> ConfigurateCase::CreateMHMMixedMesh2Spaces(){
    TPZGeoMesh *gmeshcoarse = CreateGeowithRefPattern();
    TPZAutoPointer<TPZGeoMesh> gmeshauto = gmeshcoarse;
    
    
    TPZMHMixedMeshControl *mhm = new TPZMHMixedMeshControl(gmeshauto);
    
    TPZVec<int64_t> coarseindices;
    ComputeCoarseIndices(gmeshauto.operator->(), coarseindices);
    
    
    mhm->DefinePartitionbyCoarseIndices(coarseindices);
    // Create geometric elements
    
    {
        std::set<int> matids;
        for (auto omId:fsim_case.omega_ids ) {
            matids.insert(omId);
        } ;
        mhm->fMaterialIds = matids;
        matids.clear();
        
        for (auto omId:fsim_case.gamma_ids ) {
            matids.insert(omId);
        } ;
        mhm->fMaterialBCIds = matids;
    }
    
    
    InsertMaterialObjects(*mhm);
    mhm->SetInternalPOrder(1);
    mhm->SetSkeletonPOrder(1);
    
    if(0){
        std::ofstream filemesh("After_Interfaces.txt");
        mhm->Print(filemesh);
    }

    mhm->DivideSkeletonElements(0);
    mhm->DivideBoundarySkeletonElements();
    if (0) {
        std::ofstream file_geo("Geometry.txt");
        mhm->CMesh()->Reference()->Print(file_geo);
    }
    
    bool substructure = true;
    
    mhm->BuildComputationalMesh(substructure);
    
    
    if(0)
    {
        std::ofstream file("GMeshControl.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(mhm->CMesh()->Reference(), file);
    }
    
    
#ifdef PZDEBUG2
    if(0)
    {
        std::ofstream out("MixedMeshControlHDiv.txt");
        meshcontrol.Print(out);
    }
#endif
    
    std::cout << "MHM Hdiv Computational meshes created\n";
#ifdef PZDEBUG2
    if(1)
    {
        std::ofstream gfile("geometryMHMHdiv.txt");
        gmeshauto->Print(gfile);
        std::ofstream out_mhm("MHM_hdiv.txt");
        meshcontrol.CMesh()->Print(out_mhm);
        
    }
#endif
    
    std::cout << "Number of equations MHMixed " << mhm->CMesh()->NEquations() << std::endl;
    
    if (0) {
        TPZCompMesh *MixedMesh = mhm->CMesh().operator->();
        std::ofstream multcmesh("multcompmesh.txt");
        MixedMesh->Print(multcmesh);
    }
    
    
    return mhm;
    
}

/**
 * @brief Creates a MHM Mixed mesh
 * @return MHM mesh
 */
TPZAutoPointer<TPZMHMixedMesh4SpacesControl> ConfigurateCase::CreateMHMMixedMesh4Spaces(){
    TPZGeoMesh *gmeshcoarse = CreateGeowithRefPattern();
    TPZMHMixedMesh4SpacesControl *mhm = new TPZMHMixedMesh4SpacesControl(gmeshcoarse);
    
    TPZVec<int64_t> coarseindices;
    ComputeCoarseIndices(gmeshcoarse, coarseindices);
    
    mhm->DefinePartitionbyCoarseIndices(coarseindices);
    // Create geometric elements
    
    {
        std::set<int> matids;
        for (auto omId:fsim_case.omega_ids ) {      //Materials
            matids.insert(omId);
        } ;
        mhm->fMaterialIds = matids;
        matids.clear();
        
        for (auto omId:fsim_case.gamma_ids ) {      //BC
            matids.insert(omId);
        } ;
        mhm->fMaterialBCIds = matids;
    }
    
    InsertMaterialObjects(*mhm);
    
    mhm->SetInternalPOrder(1);
    mhm->SetSkeletonPOrder(1);
    
    if (0) {
        std::ofstream filemesh("AfterInterface.txt");
        mhm->Print(filemesh);
    }

    mhm->DivideSkeletonElements(0);
    mhm->DivideBoundarySkeletonElements();
    if (0) {
        std::ofstream file_geo("geometry.txt");
        mhm->CMesh()->Reference()->Print(file_geo);
    }
    
    bool substructure = true;
    
    mhm->BuildComputationalMesh(substructure);
    std::cout << "Number of equations MHMixed " << mhm->CMesh()->NEquations() << std::endl;
    
    if (0) {
        TPZCompMesh *MixedMesh = mhm->CMesh().operator->();
        std::ofstream multcmesh("multcompmesh.txt");
        MixedMesh->Print(multcmesh);
    }
    
    return mhm;
}

/**
 * @brief Creates a MHM Mixed mesh
 * @param is4Spaces: Wether is 4 or 2 spaces
 * @return MHM mesh
 */
TPZAutoPointer<TPZMHMixedMesh4SpacesControl> ConfigurateCase::CreateMHMMixedMesh(bool is4Spaces){
    DebugStop();
    
//    if (is4Spaces==true) {
//        return CreateMHMMixedMesh4Spaces();
//    }
//    else{
//        return CreateMHMMixedMesh2Spaces();
//    }
}

/**
 * @brief Compute the Coarse Indixes
 * @param gmesh: geometric mesh
 * @param coarseindices: vector with coarse indexes
 * @return MHM mesh
 */
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices)
    {
      
        coarseindices.Resize(gmesh->NElements());
        int count = 0;
        for (int64_t el=0; el<gmesh->NElements(); el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel || gel->Dimension() != gmesh->Dimension()) continue;
            if(gel->Father()) continue;
            coarseindices[count] = el;
            count++;
        }
        coarseindices.Resize(count);
    }

/**
 * @brief Create Geometric mesh with the selected RefPattern
 * @return Geometric mesh
 */
TPZGeoMesh* ConfigurateCase::CreateGeowithRefPattern(){
    TPZGeoMesh *gmeshrefpattern = CreateUniformMesh(1,1,1,1);
    
    TPZGeoMesh *gmesgline =new TPZGeoMesh;
    
    TPZGeoEl *gel = gmeshrefpattern->Element(0);
    
    TPZFMatrix<REAL> coordinates;
    gel->NodesCoordinates(coordinates);             //Gets the element coordinates
    TPZVec<int64_t> nodeindex(4);                   //Quadrilateral coordinates
    TPZVec<int64_t> nodeindexline(2);               //Line coordinates
    int nnodes = gmeshrefpattern->NNodes();         //Gets the number of nodes
    
    //Creating a RefPattern
    gmeshrefpattern->NodeVec().Resize(nnodes +2);   //Adds two new nodes to the system (nodes to create new elements)
    //Gets the middle nodes between the "atomic" nodes
    TPZVec<REAL> x(3,0);
    x[0]= 0.0;
    x[1]= 0.5*(coordinates(1,0)+coordinates(1,2));
    gmeshrefpattern->NodeVec()[nnodes].Initialize(nnodes, x, *gmeshrefpattern);
    x[0]=coordinates(0,1);
    x[1]=0.5*(coordinates(1,1)+coordinates(1,3));
    gmeshrefpattern->NodeVec()[nnodes+1].Initialize(nnodes+1, x, *gmeshrefpattern);
    
    //Sets the partition definition for lines
    gmesgline->NodeVec().Resize(2);
    x[0]=0.0;
    x[1]=1.0;
    gmesgline->NodeVec()[0].Initialize(0, x, *gmesgline);
    x[0]=0.5;
    x[1]=1.0;
    gmesgline->NodeVec()[1].Initialize(1, x, *gmesgline);
    nodeindexline[0]=0;
    nodeindexline[1]=1;
    int64_t indexline=0;
    gmesgline->CreateGeoElement(EOned, nodeindexline, 1,indexline );    //type, corner index, matid, index
//    gmesgline->Element(1)->SetFather(gmesgline->Element(0));
    indexline=1;
    gmesgline->CreateGeoElement(EOned, nodeindexline, 1,indexline );    //type, corner index, matid, index
    gmesgline->Element(1)->SetFather(gmesgline->Element(0));
    gmesgline->BuildConnectivity();
    
    TPZAutoPointer<TPZRefPattern> refpatline = new TPZRefPattern(*gmesgline);
    
    //Sets the partition definition for faces
    int64_t index = gmeshrefpattern->NElements();
    nodeindex[0]=0;
    nodeindex[1]=1;
    nodeindex[2]=5;
    nodeindex[3]=4;
    gmeshrefpattern->CreateGeoElement(EQuadrilateral, nodeindex, 1,index );
    gmeshrefpattern->Element(index)->SetFather(gmeshrefpattern->Element(0));
    
    index = gmeshrefpattern->NElements();
    nodeindex[0]=4;
    nodeindex[1]=5;
    nodeindex[2]=3;
    nodeindex[3]=2;
    gmeshrefpattern->CreateGeoElement(EQuadrilateral, nodeindex, 1,index );
    gmeshrefpattern->Element(index)->SetFather(gmeshrefpattern->Element(0));
    
    gmeshrefpattern->BuildConnectivity();
    
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(*gmeshrefpattern);
    
    TPZGeoMesh *gmesh2 = CreateUniformMesh(2,1,1,1);
    
    //Sets the new refpattern to divide the mesh as decided
    gmesh2->Element(0)->SetRefPattern(refpat);
    gmesh2->Element(1)->SetRefPattern(refpat);
    TPZVec<TPZGeoEl *> sons;
    
    //Divides father elements and creates children
    gmesh2->Element(0)->Divide(sons);
    gmesh2->Element(1)->Divide(sons);
    gmesh2->BuildConnectivity();

    if (1) {
        std::ofstream file2("REFPatternBeforeBC.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh2, file2); //Prints a VTK before adding the BC elements
    }

    
    gmesh2->Element(0)->CreateBCGeoEl(4, -1);
    gmesh2->Element(0)->CreateBCGeoEl(6, -3);
    gmesh2->Element(0)->CreateBCGeoEl(7, -4);
    
    gmesh2->Element(1)->CreateBCGeoEl(4, -1);
    gmesh2->Element(1)->CreateBCGeoEl(5, -2);
    gmesh2->Element(1)->CreateBCGeoEl(6, -3);
    
    gmesh2->Element(8)->Divide(sons);
    gmesh2->Element(10)->Divide(sons);
    
    
    std::set<int> matids;
    matids.insert(-1);
    
    if (1) {
        std::ofstream fileafter("REFPatternAfterBC.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh2, fileafter);  //Prints a VTK after adding the BC elements
    }

    //Creates new lines to be children of the originals in order to generate the MHM mesh
    gmesh2->Element(6)->SetRefPattern(refpatline);
    gmesh2->Element(7)->SetRefPattern(refpatline);
    gmesh2->Element(9)->SetRefPattern(refpatline);
    gmesh2->Element(11)->SetRefPattern(refpatline);
    
    TPZVec<TPZGeoEl *> sons2;
    gmesh2->Element(6)->Divide(sons2);
    gmesh2->BuildConnectivity();
    gmesh2->Element(7)->Divide(sons2);
    gmesh2->BuildConnectivity();
    gmesh2->Element(9)->Divide(sons2);
    gmesh2->BuildConnectivity();
    gmesh2->Element(11)->Divide(sons2);
    gmesh2->BuildConnectivity();
    
    gmesh2->Element(18)->SetMaterialId(-5);
    gmesh2->Element(19)->SetMaterialId(-6);
    
    gmesh2->Element(4)->SetMaterialId(2);
    gmesh2->Element(5)->SetMaterialId(2);
    
    if (1) {
        std::ofstream file("WithREFPattern.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh2, file);   //Prints a final VTK of the geometric mesh after all the addings
    }
    
    return gmesh2;
}

/**
 * @brief Inserts the materials objects to the MHM mesh
 * @param control: Pointer to the mesh
 */
void ConfigurateCase::InsertMaterialObjects(TPZMHMeshControl &control)
{
    TPZCompMesh &cmesh = control.CMesh();
    TPZGeoMesh &gmesh = control.GMesh();
    
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);
    
    int dim = gmesh.Dimension();
    cmesh.SetDimModel(dim);
    
    TPZCompMesh *MixedFluxPressureCmesh = &cmesh;

    //Setting four spaces material
    int iter=0;
    // Materials
    for (auto matindex:fsim_case.omega_ids) {
        int dim = fsim_case.omega_dim[iter];
        
        TPZMixedDarcyWithFourSpaces * mat = new TPZMixedDarcyWithFourSpaces(matindex,dim);
        
        mat->SetConstantPermeability(1. + iter*1000);   //Problem set with 2 different K
        MixedFluxPressureCmesh->InsertMaterialObject(mat);
        if (iter==0) {
            int iterbc=0;
            for (auto bcIndex:fsim_case.gamma_ids) {
                int type = fsim_case.type[iterbc];
                val2[0] = fsim_case.vals[iterbc];
                TPZBndCondT<STATE> * bc = mat->CreateBC(mat, bcIndex, type, val1, val2);
                MixedFluxPressureCmesh->InsertMaterialObject(bc);
                iterbc++;
            }
        }
        iter++;
    }
}

/**
 * @brief Inserts the materials objects to the MHM mesh
 * @param control: Pointer to the mesh
 */
void ConfigurateCase::GroupAndCondense(TPZMultiphysicsCompMesh *cmesh)
{

    std::ofstream out("Malha_noMHM_cmeshBC.txt");
    cmesh->Print(out);
    int nels = cmesh->NElements();
//for (std::map<int64_t,int64_t>::iterator it=fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
    for (int el=0; el<nels; el++) {
    TPZCompEl *cel = cmesh->Element(el);
    TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
    if (!subcmesh) {
        //            DebugStop();
        continue;
    }

    TPZCompMeshTools::GroupElements(subcmesh);
    subcmesh->ComputeNodElCon();

    //Prints the submeshes (txt, vtk)
    if (0) {
        std::stringstream sout;
        std::stringstream sout2;
        sout << "submesh_" << el << ".txt";
        std::ofstream filehide(sout.str());
        sout2 << "submesh_" << el << ".vtk";
        std::ofstream file(sout2.str());
        TPZVTKGeoMesh::PrintCMeshVTK(subcmesh, file,true);
        subcmesh->Print(filehide);
        el++;
        }
        //End
        
        //        TPZAutoPointer<TPZGuiInterface> guiInter = new TPZGuiInterface;
        //        bool shouldrenumber = true;
        //        subcmesh->SetAnalysisSkyline(0,1, guiInter);
        //        subcmesh->MakeAllInternal();
        //        subcmesh->Assemble();
        //        TPZCompMesh *csubmesh = dynamic_cast<TPZCompMesh *>(subcmesh);
        //        TPZAnalysis an_sub(subcmesh,shouldrenumber);
        //        TPZSymetricSpStructMatrix strmat(subcmesh);
        //        strmat.SetNumThreads(0);
        //
        //        an_sub.SetStructuralMatrix(strmat);
        //        TPZStepSolver<STATE> step;
        //        step.SetDirect(ELDLt);
        //        an_sub.SetSolver(step);
        //        an_sub.Assemble();
        //        std::ofstream file_sub("MatrixSubmesh.txt");
        //        an_sub.Solver().Matrix()->Print("EkSm",file_sub,EMathematicaInput);
        
        //        subcmesh->Analysis();
        //        subcmesh->CalcStiff(ek, ef);
        
        //        // AQUI NACE EL PROBLEMA
        //        bool shouldrenumber = false;
        //        TPZAnalysis an_sub(subcmesh,shouldrenumber);
        //        TPZSymetricSpStructMatrix strmat(subcmesh);
        //        strmat.SetNumThreads(0);
        //        an_sub.SetStructuralMatrix(strmat);
        //        TPZStepSolver<STATE> step;
        //        step.SetDirect(ELDLt);
        //        an_sub.SetSolver(step);
        //        std::cout << "Assembling submesh\n";
        //        an_sub.Assemble();
        //        std::ofstream filemate("MatrixSubmesh.txt");
        //        an_sub.Solver().Matrix()->Print("EkSb",filemate,EMathematicaInput);
        //        TPZElementMatrix ek;
        //        TPZElementMatrix ef;
        //        cel->CalcStiff(ek, ef);
        //        // AQUI MUERE EL PROBLEMA
        
#ifdef LOG4CXX2
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            subcmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        // Increment nelconnected of exterior connects
        
        int nel = subcmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = subcmesh->Element(el);
            if (!cel) {
                continue;
            }
            int nconnects = cel->NConnects();
            for (int icon=0; icon<nconnects; icon++) {
                TPZConnect &connect = cel->Connect(icon);
                
                int lagrangemult = connect.LagrangeMultiplier();
                //Increment the number of connected elements for the avg pressure in order to not condense them
                if (lagrangemult==3) {
                    connect.IncrementElConnected();
                }
            }
        }
        
        TPZCompMeshTools::CreatedCondensedElements(subcmesh, false);
        subcmesh->CleanUpUnconnectedNodes();
        
        int numthreads = 0;
        int preconditioned = 0;
#ifdef LOG4CXX2
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            subcmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        TPZAutoPointer<TPZGuiInterface> guiInterface;
        subcmesh->SetAnalysisSkyline(numthreads, preconditioned, guiInterface);
        std::ofstream filehide2("subcmeshAfter.txt");
        subcmesh->Print(filehide2);
        }
    }
void ConfigurateCase::HideTheElements(TPZMultiphysicsCompMesh *cmesh)
{
    
    bool KeepOneLagrangian = true;
    TPZMHMixedMesh4SpacesControl control(cmesh->Reference());

    
    typedef std::set<int64_t> TCompIndexes;
    std::map<int64_t, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = cmesh->Reference();
    gmesh->ResetReference();
    cmesh->LoadReferences();
    int64_t nel = cmesh->NElements();

    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        int64_t domain = control.WhichSubdomain(cel);
        if (domain == -1) {
            continue;
        }
        ElementGroups[domain].insert(el);
    }
    if (ElementGroups.size() <= 5)
    {
        std::cout << "Number of element groups " << ElementGroups.size() << std::endl;
        std::map<int64_t,TCompIndexes>::iterator it;
        for (it=ElementGroups.begin(); it != ElementGroups.end(); it++) {
            std::cout << "Group " << it->first << " group size " << it->second.size() << std::endl;
            std::cout << " elements ";
            std::set<int64_t>::iterator its;
            for (its = it->second.begin(); its != it->second.end(); its++) {
                std::cout << *its << "|" << cmesh->Element(*its)->Reference()->Index() << " ";
            }
            std::cout << std::endl;
        }
    }
    
    std::map<int64_t,int64_t> submeshindices;
    TPZCompMeshTools::PutinSubmeshes(cmesh, ElementGroups, submeshindices, KeepOneLagrangian);
    std::cout << "After putting in substructures\n";
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
    
    GroupAndCondense(cmesh);
    
    std::cout << "Finished substructuring\n";
}

TPZCompMesh ConfigurateCase::CreateSubStructure(){
    
    m_gmesh = CreateGeowithRefPattern();
    std::ofstream outdelete("Gmesh_test_1.txt");
    m_gmesh->Print(outdelete);
    TPZCompMesh *cmesh = HDivMesh(m_gmesh, 2, 2);
    std::ofstream out("Hdiv_test_1.txt");
    cmesh->Print(out);
    int64_t nels = cmesh->NElements();
    TPZVec<int64_t> submeshes_index(nels,-1);
    int64_t index;
    int dimen;
    int index_BC;
    int index_BC2;
    
    for (int64_t el=0; el<nels; el++) {
        TPZCompEl *element = cmesh->Element(el);
        TPZGeoMesh *gmsh = cmesh->Reference();
        TPZGeoEl *geoel = element->Reference();
        dimen = element->Dimension();
        if (!element) {
            DebugStop();}
        index = geoel->FatherIndex();

        if (dimen == 2) {
            submeshes_index[el] = index;
            std::cout<<index<<std::endl;
        }
        if (dimen == 1) {
            index_BC = geoel->NeighbourIndex(2);
            index_BC2 = gmsh->Element(index_BC)->FatherIndex();
            std::cout<<index_BC2<<std::endl;
        }
    }

    TPZSubCompMesh *subcmesh0 = new TPZSubCompMesh(*cmesh);
    TPZSubCompMesh *subcmesh1 = new TPZSubCompMesh(*cmesh);
    int nd;
    
    for (int ind=0; ind<nels; ind++) {
        nd =submeshes_index[ind];
        if (nd==0) {
            subcmesh0->TransferElement(cmesh, ind);
        }
        if (nd==1) {
            subcmesh1->TransferElement(cmesh, ind);
        }
    }
//    subcmesh0->CleanUpUnconnectedNodes();
//    subcmesh1->CleanUpUnconnectedNodes();
    subcmesh0->ComputeNodElCon();
    subcmesh1->ComputeNodElCon();
    subcmesh0->MakeAllInternal();
    subcmesh1->MakeAllInternal();
    
    
    
    std::ofstream outdel0("submesh_0.txt");
    subcmesh0->Print(outdel0);
    std::ofstream outdel1("submesh_1.txt");
    subcmesh1->Print(outdel1);
    
    
    std::ofstream outdel_fat("Father_mesh_0.txt");
    subcmesh0->FatherMesh()->Print(outdel_fat);
    
    return 0;
}
