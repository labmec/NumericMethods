########### next target ###############

file(GLOB HEADERS *.h)

file(GLOB SOURCE *.cpp)


SET(Commonfiles_SRC "TPZCondensedMeshAlgebra.cpp" "TPZDistributedMeshAlgebra.cpp" "TPZHdivTransfer.cpp" "TPZMixedDarcyWithFourSpaces.cpp" "ConfigurateCase.cpp" "TPZMHMixedMesh4SpacesControl.cpp")

#add_executable(MalhasHdiv MalhasHdiv.cpp ${Commonfiles_SRC} ${HEADERS})
#add_executable(MalhasHdivRefactor MalhasHdivRefactor.cpp ${TPZCondensedMeshAlgebra.cpp TPZDistributedMeshAlgebra.cpp  TPZHdivTransfer.cpp TPZMixedDarcyWithFourSpaces.cpp }  ${HEADERS})
add_executable(MalhasHdiv4Spaces MalhasHdivRefactor.cpp ${Commonfiles_SRC}  ${HEADERS})


#target_include_directories(MalhasHdiv PUBLIC ${lib_dir})
target_include_directories(MalhasHdiv4Spaces PUBLIC ${lib_dir})


#target_link_libraries(MalhasHdiv pz ${OpenCV_LIBS})
target_link_libraries(MalhasHdiv4Spaces pz )

file (COPY PlotGenerator.nb/ DESTINATION Debug/PlotGenerator.nb)

