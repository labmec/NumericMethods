
#include "TPZMHMixedMesh4SpacesControl.h"
#include "pzelementgroup.h"


void TPZMHMixedMesh4SpacesControl::BuildComputationalMesh(bool usersubstructure){
    
    if (fpOrderInternal == 0 || fpOrderSkeleton == 0) {
        std::cout<<"Wrong polynomial order set"<<std::endl;
        DebugStop();
    }
    InsertPeriferalMaterialObjects();
    CreateHDivMHMMesh();
    
    InsertPeriferalPressureMaterialObjects();
    if(fNState > 1){
        fRotationMesh = new TPZCompMesh(fGMesh);
        InsertPeriferalRotationMaterialObjects();
    }
#ifdef PZDEBUG
    if (fFluxMesh->Dimension() != fGMesh->Dimension()) {
        DebugStop();
    }
#endif
    CreatePressureMHMMesh();
    CreateAverageFlux();
    CreateAveragePressure();
    
   
    
   
    if(fNState > 1)
    {
        CreateRotationMesh();
    }
    
    CreateHDivPressureMHMMesh();
    
    
    std::cout << "Total number of equations " << fCMesh->Solution().Rows() << std::endl;
    fGlobalSystemWithLocalCondensationSize = fCMesh->NEquations();
    fGlobalSystemSize = fCMesh->Solution().Rows();
    
//    fCMesh->ComputeNodElCon();
#ifdef PZDEBUG2
    {
        std::ofstream out("Friendly.txt");
        PrintFriendly(out);
    }
    CheckMeshConsistency();
#endif
    
    
    
    if (usersubstructure) {
        HideTheElements();
    }
  
    
    fNumeq = fCMesh->NEquations();
    
#ifdef PZDEBUG2
    {
        int64_t nel = fCMesh->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fCMesh->Element(el);
            TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
            if(sub)
            {
                std::stringstream sout;
                sout << "submesh_" << el << ".vtk";
                std::ofstream file(sout.str());
                TPZVTKGeoMesh::PrintCMeshVTK(sub, file,true);
            }
        }
        
    }
#endif
}
void TPZMHMixedMesh4SpacesControl::CreateHDivPressureMHMMesh()
{
    TPZManVector<TPZCompMesh *,4 > cmeshes(4);
    cmeshes[0] = fFluxMesh.operator->();
    cmeshes[1] = fPressureFineMesh.operator->();
    cmeshes[2] = fcmeshFluxAverg.operator->();
    cmeshes[3] = fcmeshPressureAverg.operator->();
    if (0){
        std::ofstream out("PressureMesh_MultiPhis.txt");
        cmeshes[1] ->Print(out);
        
        std::ofstream out2("FluxMesh_MultiPhis.txt");
        cmeshes[0] ->Print(out2);
        
        std::ofstream out3("FluxAverage_MultiPhis.txt");
        cmeshes[2] ->Print(out3);
        
        std::ofstream out4("PressureAverage_MultiPhis.txt");
        cmeshes[3] ->Print(out4);
    };
    
    
    if(fNState > 1) {
        cmeshes.Resize(3);
        cmeshes[2] = fRotationMesh.operator->();
    }
    TPZGeoMesh *gmesh = cmeshes[0]->Reference();
    if(!gmesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    // Multiphysics mesh
    TPZCompMesh * MixedFluxPressureCmesh = fCMesh.operator->();
    MixedFluxPressureCmesh->SetDimModel(dim);
    MixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElem();
    
    BuildMultiPhysicsMesh();
    TPZManVector<TPZCompMesh * ,4> meshvector;
    
#ifdef PZDEBUG
    if(1)
    {
        std::ofstream out2("gmesh.txt");
        gmesh->Print(out2);
        std::ofstream out3("HDivMesh.txt");
        fFluxMesh->Print(out3);
        std::ofstream out4("PressureMesh.txt");
        fPressureFineMesh->Print(out4);
    }
#endif
    
    meshvector = cmeshes;
  
    
    // populate the connect to subdomain data structure for the multiphysics mesh
    JoinSubdomains(meshvector, MixedFluxPressureCmesh);
    
    // Transferindo para a multifisica
    //    TPZBuildMultiphysicsMesh::AddElements(meshvector, MixedFluxPressureCmesh);
    //    TPZBuildMultiphysicsMesh::AddConnects(meshvector, MixedFluxPressureCmesh);
    // copy the solution of the atomic meshes to the multiphysics mesh
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, MixedFluxPressureCmesh);
    
   
    
    std::pair<int,int> skelmatid(fSkeletonMatId,fSecondSkeletonMatId);
//    CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-1);
//    CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-2);
#ifdef PZDEBUG
    if(1)
    {
        MixedFluxPressureCmesh->ComputeNodElCon();
        std::ofstream file("cmeshmphys.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(MixedFluxPressureCmesh, file,true);
        std::ofstream out("cmeshmphys.txt");
        MixedFluxPressureCmesh->Print(out);
    }
#endif
    MixedFluxPressureCmesh->CleanUpUnconnectedNodes();
    
#ifdef PZDEBUG
    if(0)
    {
        std::ofstream out("multiphysics.txt");
        MixedFluxPressureCmesh->Print(out);
    }
#endif
    
    return;
    
}

void TPZMHMixedMesh4SpacesControl::CreateAverageFlux()
{
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    TPZCompMesh * cmeshtemp = new TPZCompMesh(gmesh);
    fcmeshFluxAverg = cmeshtemp;
    
    //The pressure mesh should be empty when calling this method
    int64_t nskeletonconnects = fcmeshFluxAverg->NConnects();
    if(nskeletonconnects != 0){
        DebugStop();
    }

    // create and organize the pressure mesh
    // the pressure mesh is composed of discontinuous H1 elements
    TPZCompMesh * cmeshfluxavg = fcmeshFluxAverg.operator->();
    gmesh->ResetReference();
    cmeshfluxavg->SetName("PressureMeshAverage");
    cmeshfluxavg->SetDimModel(gmesh->Dimension());
    cmeshfluxavg->SetAllCreateFunctionsDiscontinuous(); //AQUI
    cmeshfluxavg->ApproxSpace().CreateDisconnectedElements(true);
    cmeshfluxavg->SetDefaultOrder(0);
    
//    generate elements for all material ids of meshdim
    std::set<int> matids;

    TPZNullMaterial * volume = new TPZNullMaterial(1);
    cmeshfluxavg->InsertMaterialObject(volume);
    TPZNullMaterial * volume2 = new TPZNullMaterial(2);
    cmeshfluxavg->InsertMaterialObject(volume2);
    matids.insert(1);
    matids.insert(2);
    cmeshfluxavg->AutoBuild(matids);
    fcmeshFluxAverg->ExpandSolution();
    
    if(0)
    {
        std::ofstream out("PressureAVERAGEMesh.txt");
        cmeshfluxavg->Print(out);
    }
    
    
#ifdef PZDEBUG
    // a very strange check!! Why does material id 1 need to be volumetric?
    {
        int64_t nel = fGMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if (gel && gel->MaterialId() == 1) {
                if (gel->Dimension() != fGMesh->Dimension()) {
                    DebugStop();
                }
            }
        }
    }
#endif
    
    int64_t nc = cmeshfluxavg->NConnects();
    for (int64_t ic=nskeletonconnects; ic<nc; ic++) {
        cmeshfluxavg->ConnectVec()[ic].SetLagrangeMultiplier(2); //Distributed flux
    }
    // associate the connects with the proper subdomain
    gmesh->ResetReference();
    int64_t nel = cmeshfluxavg->NElements();
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshfluxavg->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        // if a computational element was created outside the range of material ids
        // something very strange happened...
        TPZGeoEl *gel = cel->Reference();
        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
        {
            DebugStop();
        }
        int domain = fGeoToMHMDomain[gel->Index()];
#ifdef PZDEBUG
        if (domain == -1) {
            DebugStop();
        }
#endif
        
        SetSubdomain(cel, domain);
    }
    
    return;
}
void TPZMHMixedMesh4SpacesControl::CreateAveragePressure()
{
    TPZGeoMesh * gmesh = fGMesh.operator->();

    gmesh->ResetReference();
    TPZCompMesh * cmeshtemp = new TPZCompMesh(gmesh);
    fcmeshPressureAverg = cmeshtemp;
    // the pressure mesh should be empty when calling this method
    int64_t nskeletonconnects = fcmeshPressureAverg->NConnects();
    if(nskeletonconnects != 0){
        DebugStop();
    }
    
    // create and organize the pressure mesh
    // the pressure mesh is composed of discontinuous H1 elements
    TPZCompMesh * cmeshpressureavr = fcmeshPressureAverg.operator->();
    gmesh->ResetReference();
    cmeshpressureavr->SetName("PressureMeshAverage");
    cmeshpressureavr->SetDimModel(gmesh->Dimension());
    cmeshpressureavr->SetAllCreateFunctionsDiscontinuous(); //AQUI
    cmeshpressureavr->ApproxSpace().CreateDisconnectedElements(true);
    cmeshpressureavr->SetDefaultOrder(0);
    int meshdim = cmeshpressureavr->Dimension();
    // generate elements for all material ids of meshdim
    std::set<int> matids;
    //    for (auto it:fMaterialIds) {
    //        TPZMaterial *mat = fPressureFineMesh->FindMaterial(it);
    //        if (mat && mat->Dimension() == meshdim) {
    //            matids.insert(it);
    //            cmeshfluxavg->InsertMaterialObject(mat);
    //        }
    //    }
    TPZNullMaterial * volume = new TPZNullMaterial(1);
    cmeshpressureavr->InsertMaterialObject(volume);
    TPZNullMaterial * volume2 = new TPZNullMaterial(2);
    cmeshpressureavr->InsertMaterialObject(volume2);
    matids.insert(1);
    matids.insert(2);
    cmeshpressureavr->AutoBuild(matids);
    fcmeshPressureAverg->ExpandSolution();
    
    if(0)
    {
        std::ofstream out("PressureAVERAGEMesh.txt");
        cmeshpressureavr->Print(out);
    }
    
    
#ifdef PZDEBUG
    // a very strange check!! Why does material id 1 needs to be volumetric?
    {
        int64_t nel = fGMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if (gel && gel->MaterialId() == 1) {
                if (gel->Dimension() != fGMesh->Dimension()) {
                    DebugStop();
                }
            }
        }
    }
#endif
    
    int64_t nc = cmeshpressureavr->NConnects();
    for (int64_t ic=nskeletonconnects; ic<nc; ic++) {
        cmeshpressureavr->ConnectVec()[ic].SetLagrangeMultiplier(3);
    }
    // associate the connects with the proper subdomain
    gmesh->ResetReference();
    int64_t nel = cmeshpressureavr->NElements();
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshpressureavr->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        // if a computational element was created outside the range of material ids
        // something very strange happened...
        TPZGeoEl *gel = cel->Reference();
        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
        {
            DebugStop();
        }
        int domain = fGeoToMHMDomain[gel->Index()];
#ifdef PZDEBUG
        if (domain == -1) {
            DebugStop();
        }
#endif//
       
        
        
        //
        
        SetSubdomain(cel, domain);
    }
    
  
    
    return;
}
void TPZMHMixedMesh4SpacesControl::BuildMultiPhysicsMesh()
{
    if (fCMesh->NElements() != 0) {
        std::cout<<"The mesh has to be empty to build it at this stage"<<std::endl;
        DebugStop();
    }
    fCMesh->SetAllCreateFunctionsMultiphysicElem();
    TPZMultiphysicsCompMesh *mphysics = dynamic_cast<TPZMultiphysicsCompMesh *>(fCMesh.operator->());
   
    int vecsize = 4;
    TPZManVector<TPZCompMesh *> meshvec(vecsize);
    meshvec[0] = fFluxMesh.operator->();
    meshvec[1] = fPressureFineMesh.operator->();
    meshvec[2] = fcmeshFluxAverg.operator->();
    meshvec[3] = fcmeshPressureAverg.operator->();
//    if(fNState > 1)
//    {
//        meshvec[2] = this->fRotationMesh.operator->();
//    }
    TPZManVector<int64_t> shouldcreate(fGMesh->NElements(),0);
    std::set<int> matids;
    for (auto it : fCMesh->MaterialVec()) {
        matids.insert(it.first);
    }
    int64_t nel = fFluxMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fFluxMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        // this means that all geometric elements associated with flux elements will generate a computational element
        if (matids.find(gel->MaterialId()) == matids.end()) {
            DebugStop();
        }
        if (cel) {
            shouldcreate[cel->Reference()->Index()] = 1;
        }
    }
    nel = fcmeshPressureAverg->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fcmeshPressureAverg->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if (matids.find(gel->MaterialId()) == matids.end()) {
            DebugStop();
        }
        if (cel) {
            shouldcreate[cel->Reference()->Index()] = 1;
        }
    }
    // define the intersection of the finest references
    nel = fGMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (shouldcreate[el])
        {
            TPZGeoEl *fat = gel->Father();
            while(fat)
            {
                if(shouldcreate[fat->Index()] == 1)
                {
                    shouldcreate[fat->Index()] = 0;
                }
                fat = fat->Father();
            }
        }
    }
    TPZStack<int64_t> gelindexes;
    for (int64_t el=0; el<nel; el++) {
        if (shouldcreate[el])
        {
            gelindexes.Push(el);
        }
    }
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Geometric indices for which we will create multiphysics elements" << std::endl;
        sout << gelindexes;
        //        std::cout << sout.str() << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    

    mphysics->BuildMultiphysicsSpace(meshvec,gelindexes);
    
  
    
}

void TPZMHMixedMesh4SpacesControl::HideTheElements()
{
   
    bool KeepOneLagrangian = true;
    if (fHybridize) {
        KeepOneLagrangian = false;
    }
    //
   

    //
    
    typedef std::set<int64_t> TCompIndexes;
    std::map<int64_t, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = fCMesh->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    fCMesh->LoadReferences();
    int64_t nel = fCMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->Element(el);
        int64_t domain = WhichSubdomain(cel);
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
                std::cout << *its << "|" << fCMesh->Element(*its)->Reference()->Index() << " ";
            }
            std::cout << std::endl;
        }
    }
    
    std::map<int64_t,int64_t> submeshindices;
    TPZCompMeshTools::PutinSubmeshes(fCMesh.operator->(), ElementGroups, submeshindices, KeepOneLagrangian);
    std::cout << "After putting in substructures\n";
    fMHMtoSubCMesh = submeshindices;
    fCMesh->ComputeNodElCon();
    fCMesh->CleanUpUnconnectedNodes();
   
    GroupandCondenseElements();
   
    
    std::cout << "Finished substructuring\n";
}

void TPZMHMixedMesh4SpacesControl::GroupandCondenseElements()
{
    
    for (std::map<int64_t,int64_t>::iterator it=fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
        TPZCompEl *cel = fCMesh->Element(it->second);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            DebugStop();
        }
      
        TPZCompMeshTools::GroupElements(subcmesh);
        subcmesh->ComputeNodElCon();
        std::ofstream filehide("subcmesh1.txt");
        subcmesh->Print(filehide);

#ifdef LOG4CXX2
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            subcmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        // TODO: increment nelconnected of exterior connects
        bool keeplagrange = true;
        
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
            std::cout<<lagrangemult<<std::endl;
            if (lagrangemult==3) {
                connect.IncrementElConnected();
            }
        }
        }
        //
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
    }
}
/// Put the element set into a subcompmesh and make the connects internal
void TPZMHMixedMesh4SpacesControl::PutinSubmeshes(TPZCompMesh *cmesh, std::map<int64_t,std::set<int64_t> >&elindices, std::map<int64_t,int64_t> &indices, bool KeepOneLagrangian)

{
    for (std::map<int64_t,std::set<int64_t> >::iterator it = elindices.begin(); it != elindices.end(); it++) {
        int64_t index;
        TPZSubCompMesh *subcmesh = new TPZSubCompMesh(*cmesh,index);
        indices[it->first] = index;
        for (std::set<int64_t>::iterator itloc = it->second.begin(); itloc != it->second.end(); itloc++) {
            subcmesh->TransferElement(cmesh, *itloc);
        }
    }
    cmesh->ComputeNodElCon();
    for (std::map<int64_t,int64_t>::iterator it = indices.begin(); it != indices.end(); it++) {
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cmesh->Element(it->second));
        if (!subcmesh) {
            DebugStop();
        }
        int count = 0;
        if (KeepOneLagrangian)
        {
            int64_t nconnects = subcmesh->NConnects();
            for (int64_t ic=0; ic<nconnects; ic++) {
                TPZConnect &c = subcmesh->Connect(ic);
                if (c.LagrangeMultiplier() > 0) {
                    c.IncrementElConnected();
                    count++;
                    if(count == 1 && c.NState() == 1)
                    {
                        break;
                    }
                    else if(count == 2 && c.NState() == 2)
                    {
                        break;
                    }
                    else if(count == 3 && c.NState() == 3)
                    {
                        break;
                    }
                }
            }
        }
        subcmesh->MakeAllInternal();
    }
}
