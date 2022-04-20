//
// Created by Jeferson Fernandes on 11/08/21.
//

#ifndef PRINT_UTILS_H
#define PRINT_UTILS_H

// TODO add doc

#include <iostream>
#include <string>
#include <map>
#include "pzstack.h"
#include <TPZVTKGeoMesh.h>

template<class T>
class TPZVec;

class TPZCompMesh;
class TPZGeoMesh;
class TPZMaterial;
class TPZCompElSide;
class TPZMultiphysicsCompMesh;


class TPZPrintUtils {

public:
    /**
     * @brief Prints the computational mesh information (specially the connects)
     * 
     * @param cmesh 
     */
    void PrintCMeshConnects(TPZCompMesh *cmesh);

    /**
     * @brief Prints geometric mesh information (specially neighbours)
     * 
     * @param geomesh 
     */
    void PrintGeoMesh(TPZGeoMesh *geomesh);

    /**
     * @brief Prints computational mesh information to a file
     * 
     * @param cmesh 
     * @param file_name 
     */
    void PrintCompMesh(TPZCompMesh *cmesh, std::string file_name);

    
};
#endif