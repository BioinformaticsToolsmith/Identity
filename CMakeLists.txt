cmake_minimum_required (VERSION 3.10)

set(CMAKE_CXX_COMPILER g++-7)

project (identity)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set(SOURCES 
${CMAKE_SOURCE_DIR}/src/Aligner.cpp			
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
${CMAKE_SOURCE_DIR}/src/Identity.cpp
${CMAKE_SOURCE_DIR}/src/GLMPredictor.cpp
)

set(HEADER_FILES 
${CMAKE_SOURCE_DIR}/src/KmerHistogram.h
${CMAKE_SOURCE_DIR}/src/Statistician.h
${CMAKE_SOURCE_DIR}/src/BestFirst.h
${CMAKE_SOURCE_DIR}/src/LockFreeQueue.h
${CMAKE_SOURCE_DIR}/src/AlignerParallel.h
)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS " -O3 -ffast-math -fopenmp -fmessage-length=0")

add_executable(identity ${CMAKE_SOURCE_DIR}/src/Identity.cpp ${SOURCES} ${HEADER_FILES})

