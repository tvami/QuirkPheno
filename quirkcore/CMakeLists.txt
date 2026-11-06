# Quirk Geant4 extension: ported from Athena Simulation/G4Extensions/Quirks
file(GLOB QUIRKCORE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
add_library(quirkcore SHARED ${QUIRKCORE_SRC})
target_include_directories(quirkcore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(quirkcore PUBLIC ${G4_CFLAGS} -Wno-deprecated-declarations)
target_link_libraries(quirkcore PUBLIC ${G4_LIBS})
