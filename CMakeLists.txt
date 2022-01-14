cmake_minimum_required (VERSION 3.10)
project (identity)

# Credit: https://stackoverflow.com/questions/52180281/cmake-cxx-compiler-version-is-pointing-to-the-old-gcc-version
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.5)
        message(FATAL_ERROR "GNU g++ 7.5.0 or later is required. Your current version is: " ${CMAKE_CXX_COMPILER_VERSION})
    endif()  
else()
    message(FATAL_ERROR "Your compiler is currently unsupported: " ${CMAKE_CXX_COMPILER_ID})
endif()

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS " -O3 -ffast-math -fopenmp -fmessage-length=0")

add_library(main SHARED		
${CMAKE_SOURCE_DIR}/src/FastaReader.cpp		
${CMAKE_SOURCE_DIR}/src/Mutator.cpp			
${CMAKE_SOURCE_DIR}/src/ReaderAlignerCoordinator.cpp			
${CMAKE_SOURCE_DIR}/src/Parameters.cpp			
${CMAKE_SOURCE_DIR}/src/Matrix.cpp
${CMAKE_SOURCE_DIR}/src/DataGenerator.cpp
${CMAKE_SOURCE_DIR}/src/SynDataGenerator.cpp
${CMAKE_SOURCE_DIR}/src/GLMClassifier.cpp
${CMAKE_SOURCE_DIR}/src/GLMRegressor.cpp
${CMAKE_SOURCE_DIR}/src/Normalizer.cpp
${CMAKE_SOURCE_DIR}/src/SimConverter.cpp
${CMAKE_SOURCE_DIR}/src/FeatureExpander.cpp
${CMAKE_SOURCE_DIR}/src/GLM.cpp
${CMAKE_SOURCE_DIR}/src/Node.cpp
${CMAKE_SOURCE_DIR}/src/Evaluator.cpp
${CMAKE_SOURCE_DIR}/src/Feature.cpp
${CMAKE_SOURCE_DIR}/src/FeatureSquared.cpp
${CMAKE_SOURCE_DIR}/src/FeaturePaired.cpp
${CMAKE_SOURCE_DIR}/src/StatisticInfo.cpp
${CMAKE_SOURCE_DIR}/src/ITransformer.cpp
${CMAKE_SOURCE_DIR}/src/GLMPredictor.cpp
${CMAKE_SOURCE_DIR}/src/Serializer.cpp
)

add_library(clustering SHARED
${CMAKE_SOURCE_DIR}/src/meshclust/MeShClust.cpp
${CMAKE_SOURCE_DIR}/src/meshclust/ClusterInfo.cpp
${CMAKE_SOURCE_DIR}/src/meshclust/ClusterEvaluator.cpp
)

target_link_libraries(clustering main)

set(HEADER_FILES 
${CMAKE_SOURCE_DIR}/src/Util.h
${CMAKE_SOURCE_DIR}/src/Aligner.h	
${CMAKE_SOURCE_DIR}/src/KmerHistogram.h
${CMAKE_SOURCE_DIR}/src/Statistician.h
${CMAKE_SOURCE_DIR}/src/BestFirst.h
${CMAKE_SOURCE_DIR}/src/LockFreeQueue.h
${CMAKE_SOURCE_DIR}/src/AlignerParallel.h
${CMAKE_SOURCE_DIR}/src/IdentityCalculator.h
${CMAKE_SOURCE_DIR}/src/IdentityCalculator1.h
)

set(CLUSTER_HEADER_FILES 
${CMAKE_SOURCE_DIR}/src/meshclust/ClusteringUtil.h
${CMAKE_SOURCE_DIR}/src/meshclust/MeanShift.h
${CMAKE_SOURCE_DIR}/src/meshclust/MeanShiftLarge.h
${CMAKE_SOURCE_DIR}/src/meshclust/Cluster.h
${CMAKE_SOURCE_DIR}/src/meshclust/Reservoir.h
)

add_executable(identity ${CMAKE_SOURCE_DIR}/src/Identity.cpp ${HEADER_FILES})
target_link_libraries(identity main)

add_executable(meshclust ${CMAKE_SOURCE_DIR}/src/meshclust/MeShClust.cpp ${HEADER_FILES} ${CLUSTER_HEADER_FILES})
target_link_libraries(meshclust clustering)
