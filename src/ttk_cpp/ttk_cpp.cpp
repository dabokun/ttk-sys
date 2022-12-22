/// \ingroup examples
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2017.
///
/// \brief Minimalist C++-only TTK example pipeline, including:
///  -# The computation of a persistence curve
///  -# The computation of a persistence diagram
///  -# The selection of the most persistent pairs of the diagram
///  -# The pre-simplification of the data according to this selection
///  -# The computation of the Morse-Smale complex on this simplified data
///  -# The storage of the output of this pipeline to disk.
///
/// This example reproduces the Figure 1 of the TTK companion paper:
/// "The Topology ToolKit", J. Tierny, G. Favelier, J. Levine, C. Gueunet, M.
/// Michaux., IEEE Transactions on Visualization and Computer Graphics, Proc.
/// of IEEE VIS 2017.
///
/// See the individual VTK wrappers (core/vtk/) to see how to use each ttk::base
/// (C++-only) TTK component.

// include the local headers

#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <CommandLineParser.h>
#include <MorseSmaleComplex.h>
#include <PersistenceCurve.h>
#include <PersistenceDiagram.h>
#include <ContourTree.h>
#include <TopologicalSimplification.h>
#include <HelloWorld.h>
#include <FTRGraph.h>
#include <Triangulation.h>
#include "ttk_cpp.hpp"

Ttk_rs::Ttk_rs()
{
}

void Ttk_rs::test_ftr(
    float *data_3d,       // 一次元ベクタに変換後の面分光データ(三次元スカラー場)
    unsigned int datalen, // data_3dのサイズ
    unsigned int xdim,    // もとの三次元データのそれぞれの次元のサイズ
    unsigned int ydim,    // xdim,ydimは赤経赤緯方向、zdimは波長方向
    unsigned int zdim)
{
  ttk::globalDebugLevel_ = 3;
  ttk::Triangulation triangulation;
  triangulation.setInputGrid(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, xdim, ydim, zdim);

  ttk::ftr::FTRGraph<float, ttk::Triangulation> myftr;
  std::vector<float> pointSet(data_3d, data_3d + datalen);
  std::vector<ttk::SimplexId> order(pointSet.size());
  ttk::preconditionOrderArray(pointSet.size(), pointSet.data(), order.data());
  myftr.preconditionTriangulation(&triangulation);
  myftr.setScalars(pointSet.data());
  myftr.setVertexSoSoffsets(order.data());
  myftr.build();

  ttk::ftr::Graph mygraph;
  mygraph = myftr.extractOutputGraph();
  std::string str = mygraph.print(1);
  std::cout << str << std::endl;
  return;
}

void Ttk_rs::ttk_helloworld()
{
  ttk::HelloWorld myhello;
  ttk::globalDebugLevel_ = 3;
  myhello.printMsg("Hello World!");
  return;
}

void Ttk_rs::compute_persistence_pairs(float *data_2d, unsigned int datalen, unsigned int xdim, unsigned int ydim, int *birth, int *death, int *len)
{
  ttk::globalDebugLevel_ = 3;
  // validate
  if (datalen != xdim * ydim)
  {
    return;
  }
  std::vector<float> pointSet(data_2d, data_2d + datalen);
  ttk::Triangulation triangulation;
  triangulation.setInputGrid(0.0, 0.0, 0.0, 1.0, 1.0, 0.0, xdim, ydim, 1);

  std::vector<ttk::SimplexId> order(pointSet.size());
  ttk::preconditionOrderArray(pointSet.size(), pointSet.data(), order.data(), 1);

  ttk::PersistenceDiagram diagram;
  std::vector<ttk::PersistencePair> diagramOutput(5000);
  diagram.preconditionTriangulation(&triangulation);
  diagram.execute(diagramOutput, pointSet.data(), 0, order.data(), &triangulation);
  for (int i = 0; i < diagramOutput.size(); ++i)
  {
    birth[i] = diagramOutput[i].birth.id;
    death[i] = diagramOutput[i].death.id;
  }
  *len = diagramOutput.size();
  /*std::cout << "PersistentDiagram======================================" << diagramOutput.size() << std::endl;
  for (int i = 0; i < diagramOutput.size(); ++i)
  {
    std::cout << "birth: " << diagramOutput[i].birth << ", death: " << diagramOutput[i].death << std::endl;
    std::cout << "birthType: " << (int)diagramOutput[i].birthType << ", deathType: " << (int)diagramOutput[i].deathType << std::endl;
    std::cout << "persistence: " << diagramOutput[i].persistence << std::endl;
    std::cout << "birth_scalar: " << pointSet[diagramOutput[i].birth] << ", death_scalar: " << pointSet[diagramOutput[i].death] << std::endl;
  }*/
  std::cout << "compute CC End" << std::endl;
  return;
}

void Ttk_rs::compute_persistence_pairs_3d(float *data_3d, unsigned int datalen, unsigned int xdim, unsigned int ydim, unsigned int zdim, int *birth, int *death, int *len)
{
  ttk::globalDebugLevel_ = 3;
  // validate
  if (datalen != xdim * ydim * zdim)
  {
    return;
  }
  std::vector<float> pointSet(data_3d, data_3d + datalen);
  ttk::Triangulation triangulation;
  triangulation.setInputGrid(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, xdim, ydim, zdim);

  std::vector<ttk::SimplexId> order(pointSet.size());
  ttk::preconditionOrderArray(pointSet.size(), pointSet.data(), order.data(), 1);

  ttk::PersistenceDiagram diagram;
  std::vector<ttk::PersistencePair> diagramOutput(50000);
  diagram.preconditionTriangulation(&triangulation);
  diagram.execute(diagramOutput, pointSet.data(), 0, order.data(), &triangulation);
  std::cout << "Number of Pairs: " << diagramOutput.size();
  for (int i = 0; i < diagramOutput.size(); ++i)
  {
    birth[i] = diagramOutput[i].birth.id;
    death[i] = diagramOutput[i].death.id;
  }
  *len = diagramOutput.size();
  /*std::cout << "PersistentDiagram======================================" << diagramOutput.size() << std::endl;
  for (int i = 0; i < diagramOutput.size(); ++i)
  {
    std::cout << "birth: " << diagramOutput[i].birth << ", death: " << diagramOutput[i].death << std::endl;
    std::cout << "birthType: " << (int)diagramOutput[i].birthType << ", deathType: " << (int)diagramOutput[i].deathType << std::endl;
    std::cout << "persistence: " << diagramOutput[i].persistence << std::endl;
    std::cout << "birth_scalar: " << pointSet[diagramOutput[i].birth] << ", death_scalar: " << pointSet[diagramOutput[i].death] << std::endl;
  }*/
  std::cout << "compute CC End" << std::endl;
  return;
}

void Ttk_rs::simplification(float *data_2d, unsigned int datalen, unsigned int xdim, unsigned int ydim, int *authorized_birth_ids, int *authorized_death_ids, unsigned int authorized_datalen,
                            unsigned int *cp_point_type, float *cp_coordx, float *cp_coordy, float *cp_value, unsigned int *cp_cellid, unsigned int *cp_pl_vertex_identifier, unsigned int *cp_manifold_size,
                            unsigned int *cp_len, unsigned int *sp_id, float *sp_coordx, float *sp_coordy, unsigned int *sp_point_type, unsigned int *sp_cellid, unsigned int *sp_len,
                            unsigned int *sc_id, unsigned int *sc_source, unsigned int *sc_dest, unsigned int *sc_connectivity_s, unsigned int *sc_connectivity_d, unsigned int *sc_separatrix_id,
                            unsigned int *sc_separatrix_type, unsigned int *sc_f_maxima, unsigned int *sc_f_minima, float *sc_f_diff, unsigned int *sc_len)
{
  ttk::globalDebugLevel_ = 3;
  ttk::Triangulation triangulation;
  std::vector<float> pointSet(data_2d, data_2d + datalen);
  triangulation.setInputGrid(0.0, 0.0, 0.0, 1.0, 1.0, 0.0, xdim, ydim, 1);
  std::vector<ttk::SimplexId> order(pointSet.size());
  ttk::preconditionOrderArray(pointSet.size(), pointSet.data(), order.data());
  ttk::TopologicalSimplification simplification;
  simplification.preconditionTriangulation(&triangulation);
  std::vector<float> simplifiedHeight = pointSet;
  std::vector<ttk::SimplexId> authorizedCriticalPoints, simplifiedOrder = order;
  for (int i = 0; i < authorized_datalen; ++i)
  {
    authorizedCriticalPoints.push_back(authorized_birth_ids[i]);
    authorizedCriticalPoints.push_back(authorized_death_ids[i]);
  }
  simplification.execute<float>(pointSet.data(), simplifiedHeight.data(),
                                authorizedCriticalPoints.data(), order.data(),
                                simplifiedOrder.data(),
                                authorizedCriticalPoints.size(), triangulation);
  // 7. computing the Morse-Smale complex
  ttk::MorseSmaleComplex morseSmaleComplex;
  // critical points
  ttk::MorseSmaleComplex::OutputCriticalPoints outCriticalPoints{};
  // 1-separatrices
  ttk::MorseSmaleComplex::Output1Separatrices out1Separatrices{};
  // 2-separatrices
  ttk::MorseSmaleComplex::Output2Separatrices out2Separatrices{};
  std::vector<ttk::SimplexId> ascendingSegmentation(
      triangulation.getNumberOfVertices(), -1),
      descendingSegmentation(triangulation.getNumberOfVertices(), -1),
      mscSegmentation(triangulation.getNumberOfVertices(), -1);
  ttk::MorseSmaleComplex::OutputManifold outSegmentation{
      ascendingSegmentation.data(), descendingSegmentation.data(),
      mscSegmentation.data()};
  morseSmaleComplex.preconditionTriangulation(&triangulation);
  morseSmaleComplex.execute(
      outCriticalPoints, out1Separatrices, out2Separatrices, outSegmentation,
      simplifiedHeight.data(), 0, simplifiedOrder.data(), triangulation);
  // segmentation
  /*// critical points
  ttk::SimplexId criticalPoints_numberOfPoints{};
  std::vector<std::array<float, 3>> criticalPoints_points;
  std::vector<char> criticalPoints_points_cellDimensions;
  std::vector<ttk::SimplexId> criticalPoints_points_cellIds;
  std::vector<char> criticalPoints_points_isOnBoundary;
  std::vector<ttk::SimplexId> criticalPoints_points_PLVertexIdentifiers;
  std::vector<ttk::SimplexId> criticalPoints_points_manifoldSize;
  // 1-separatrices
  ttk::SimplexId separatrices1_numberOfPoints{};
  std::vector<float> separatrices1_points;
  std::vector<char> separatrices1_points_smoothingMask;
  std::vector<char> separatrices1_points_cellDimensions;
  std::vector<ttk::SimplexId> separatrices1_points_cellIds;
  ttk::SimplexId separatrices1_numberOfCells{};
  std::vector<ttk::SimplexId> separatrices1_cells_connectivity;
  std::vector<ttk::SimplexId> separatrices1_cells_sourceIds;
  std::vector<ttk::SimplexId> separatrices1_cells_destinationIds;
  std::vector<ttk::SimplexId> separatrices1_cells_separatrixIds;
  std::vector<char> separatrices1_cells_separatrixTypes;
  std::vector<char> separatrices1_cells_isOnBoundary;
  std::vector<ttk::SimplexId> separatrices1_cells_separatrixFunctionMaximaId;
  std::vector<ttk::SimplexId> separatrices1_cells_separatrixFunctionMinimaId;
  std::vector<double> separatrices1_cells_separatrixFunctionDiffs;
  // segmentation
  std::vector<ttk::SimplexId> ascendingSegmentation(
      triangulation.getNumberOfVertices(), -1),
      descendingSegmentation(triangulation.getNumberOfVertices(), -1),
      mscSegmentation(triangulation.getNumberOfVertices(), -1);
  morseSmaleComplex.preconditionTriangulation(&triangulation);
  morseSmaleComplex.setInputScalarField(simplifiedHeight.data());
  morseSmaleComplex.setInputOffsets(simplifiedOrder.data());
  morseSmaleComplex.setOutputMorseComplexes(ascendingSegmentation.data(),
                                            descendingSegmentation.data(),
                                            mscSegmentation.data());
  morseSmaleComplex.setOutputCriticalPoints(
      &criticalPoints_points,
      &criticalPoints_points_cellDimensions, &criticalPoints_points_cellIds,
      &criticalPoints_points_isOnBoundary,
      &criticalPoints_points_PLVertexIdentifiers,
      &criticalPoints_points_manifoldSize);
  morseSmaleComplex.setOutputSeparatrices1(
      &separatrices1_numberOfPoints, &separatrices1_points,
      &separatrices1_points_smoothingMask, &separatrices1_points_cellDimensions,
      &separatrices1_points_cellIds, &separatrices1_numberOfCells,
      &separatrices1_cells_connectivity, &separatrices1_cells_sourceIds,
      &separatrices1_cells_destinationIds, &separatrices1_cells_separatrixIds,
      &separatrices1_cells_separatrixTypes,
      &separatrices1_cells_separatrixFunctionMaximaId,
      &separatrices1_cells_separatrixFunctionMinimaId,
      &separatrices1_cells_isOnBoundary);

  morseSmaleComplex.execute<float>(triangulation);*/
  std::cout << "Number Of CP: " << (unsigned int)outCriticalPoints.points_.size() << std::endl;
  std::cout << "Number Of SP: " << (unsigned int)out1Separatrices.pt.numberOfPoints_ << std::endl;
  std::cout << "Number Of SC: " << (unsigned int)out1Separatrices.cl.numberOfCells_ << std::endl;

  *cp_len = outCriticalPoints.cellDimensions_.size();
  std::cout << "Input CP" << std::endl;
  for (int i = 0; i < outCriticalPoints.cellDimensions_.size(); ++i)
  {
    cp_point_type[i] = outCriticalPoints.cellDimensions_[i];
    cp_coordx[i] = outCriticalPoints.points_[i][0];
    cp_coordy[i] = ydim - outCriticalPoints.points_[i][1];
    cp_value[i] = outCriticalPoints.points_[i][2];
    cp_cellid[i] = outCriticalPoints.cellIds_[i];
    cp_pl_vertex_identifier[i] = outCriticalPoints.PLVertexIdentifiers_[i];
    cp_manifold_size[i] = outCriticalPoints.manifoldSize_[i];
  }
  std::cout << "Input SP" << std::endl;
  *sp_len = (unsigned int)out1Separatrices.pt.numberOfPoints_;
  for (int i = 0; i < out1Separatrices.pt.numberOfPoints_; ++i)
  {
    sp_id[i] = i;
    sp_coordx[i] = out1Separatrices.pt.points_[i * 3];
    sp_coordy[i] = ydim - out1Separatrices.pt.points_[i * 3 + 1];
    sp_point_type[i] = (unsigned int)out1Separatrices.pt.cellDimensions_[i];
    sp_cellid[i] = out1Separatrices.pt.cellIds_[i];
  }
  std::cout << "Input SC" << std::endl;
  *sc_len = (unsigned int)out1Separatrices.cl.numberOfCells_;
  for (int i = 0; i < out1Separatrices.cl.numberOfCells_; ++i)
  {
    sc_id[i] = (unsigned int)i;

    sc_source[i] = (unsigned int)out1Separatrices.cl.sourceIds_[i];

    sc_dest[i] = (unsigned int)out1Separatrices.cl.destinationIds_[i];

    sc_connectivity_s[i] = (unsigned int)out1Separatrices.cl.connectivity_[i * 2];

    sc_connectivity_d[i] = (unsigned int)out1Separatrices.cl.connectivity_[i * 2 + 1];

    sc_separatrix_id[i] = (unsigned int)out1Separatrices.cl.separatrixIds_[i];

    sc_separatrix_type[i] = (unsigned int)out1Separatrices.cl.separatrixTypes_[i];

    /*std::cout << "sc_f_maxima[" << i << "]" << std::endl;
    std::cout << "sc_size: " << separatrices1_cells_separatrixFunctionMaximaId.size() << std::endl;
    sc_f_maxima[i] = (unsigned int)separatrices1_cells_separatrixFunctionMaximaId[i];
    std::cout << "OK" << std::endl;
    sc_f_minima[i] = (unsigned int)separatrices1_cells_separatrixFunctionMinimaId[i];*/
    //  sc_f_diff[i] = separatrices1_cells_separatrixFunctionDiffs[i];
  }
  /*std::cout << "Separatrices1================================== Points" << std::endl;
  std::cout << separatrices1_numberOfPoints << std::endl;
  std::cout << separatrices1_points.size() << std::endl;
  std::cout << separatrices1_points_cellDimensions.size() << std::endl;
  std::cout << separatrices1_points_cellIds.size() << std::endl;*/
  /*for (int i = 0; i < separatrices1_numberOfPoints; ++i)
  {
    std::cout << "ID:" << i << std::endl;
    std::cout << "Coord: " << separatrices1_points[i * 3] << ", " << ydim - separatrices1_points[i * 3 + 1] << std::endl;
    std::cout << "Kind: " << (int)separatrices1_points_cellDimensions[i] << std::endl;
    std::cout << "Cell ID:" << separatrices1_points_cellIds[i] << std::endl;
    std::cout << std::endl;
  }*/

  /*std::cout << separatrices1_numberOfCells << std::endl;
  std::cout << separatrices1_cells_connectivity.size() << std::endl;
  std::cout << separatrices1_cells_sourceIds.size() << ", " << separatrices1_cells_destinationIds.size() << std::endl;
  std::cout << separatrices1_cells_separatrixIds.size() << std::endl;
  std::cout << separatrices1_cells_separatrixTypes.size() << std::endl;
  std::cout << separatrices1_cells_separatrixFunctionMaxima.size() << ", " << separatrices1_cells_separatrixFunctionMinima.size() << ", " << separatrices1_cells_separatrixFunctionDiffs.size() << std::endl;*/
  /*std::cout << "Separatrices1================================= Cells" << std::endl;
  for (int i = 0; i < separatrices1_numberOfCells; ++i)
  {
    std::cout << "Source: " << separatrices1_cells_sourceIds[i] << ", Dest: " << separatrices1_cells_destinationIds[i] << std::endl;
    std::cout << "Connectivity?: " << separatrices1_cells_connectivity[i * 2] << ", " << separatrices1_cells_connectivity[i * 2 + 1] << std::endl;
    std::cout << "Separatrix Id: " << separatrices1_cells_separatrixIds[i] << std::endl;
    std::cout << "Separatrix Type: " << (int)separatrices1_cells_separatrixTypes[i] << std::endl;
    std::cout << "F Maxima: " << separatrices1_cells_separatrixFunctionMaxima[i] << ", "
              << "F Minima: " << separatrices1_cells_separatrixFunctionMinima[i] << ", "
              << "F Diff: " << separatrices1_cells_separatrixFunctionDiffs[i] << std::endl
              << std::endl;
  }*/
  std::cout << "compute Simplification End" << std::endl;
  return;
}

void Ttk_rs::simplification_3d(float *data_3d, unsigned int datalen, unsigned int xdim, unsigned int ydim, unsigned int zdim, int *authorized_birth_ids, int *authorized_death_ids, unsigned int authorized_datalen,
                               unsigned int *cp_point_type, float *cp_coordx, float *cp_coordy, float *cp_coordz, float *cp_value, unsigned int *cp_cellid, unsigned int *cp_pl_vertex_identifier, unsigned int *cp_manifold_size,
                               unsigned int *cp_len, unsigned int *sp_id, float *sp_coordx, float *sp_coordy, float *sp_coordz, unsigned int *sp_point_type, unsigned int *sp_cellid, unsigned int *sp_len,
                               unsigned int *sc_id, unsigned int *sc_source, unsigned int *sc_dest, unsigned int *sc_connectivity_s, unsigned int *sc_connectivity_d, unsigned int *sc_separatrix_id,
                               unsigned int *sc_separatrix_type, unsigned int *sc_f_maxima, unsigned int *sc_f_minima, float *sc_f_diff, unsigned int *sc_len)
{
  ttk::globalDebugLevel_ = 3;
  ttk::Triangulation triangulation;
  std::vector<float> pointSet(data_3d, data_3d + datalen);
  triangulation.setInputGrid(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, xdim, ydim, zdim);
  std::vector<ttk::SimplexId> order(pointSet.size());
  ttk::preconditionOrderArray(pointSet.size(), pointSet.data(), order.data());
  ttk::TopologicalSimplification simplification;
  simplification.preconditionTriangulation(&triangulation);
  std::vector<float> simplifiedHeight = pointSet;
  std::vector<ttk::SimplexId> authorizedCriticalPoints, simplifiedOrder = order;
  for (int i = 0; i < authorized_datalen; ++i)
  {
    authorizedCriticalPoints.push_back(authorized_birth_ids[i]);
    authorizedCriticalPoints.push_back(authorized_death_ids[i]);
  }
  simplification.execute<float>(pointSet.data(), simplifiedHeight.data(),
                                authorizedCriticalPoints.data(), order.data(),
                                simplifiedOrder.data(),
                                authorizedCriticalPoints.size(), triangulation);
  /*// 7. computing the Morse-Smale complex
  ttk::MorseSmaleComplex morseSmaleComplex;
  ttk::MorseSmaleComplex::OutputCriticalPoints outCriticalPoints{};
  // critical points
  ttk::SimplexId criticalPoints_numberOfPoints{};
  std::vector<std::array<float, 3>> criticalPoints_points;
  std::vector<char> criticalPoints_points_cellDimensions;
  std::vector<ttk::SimplexId> criticalPoints_points_cellIds;
  std::vector<char> criticalPoints_points_isOnBoundary;
  std::vector<ttk::SimplexId> criticalPoints_points_PLVertexIdentifiers;
  std::vector<ttk::SimplexId> criticalPoints_points_manifoldSize;
  // 1-separatrices
  ttk::SimplexId separatrices1_numberOfPoints{};
  std::vector<float> separatrices1_points;
  std::vector<char> separatrices1_points_smoothingMask;
  std::vector<char> separatrices1_points_cellDimensions;
  std::vector<ttk::SimplexId> separatrices1_points_cellIds;
  ttk::SimplexId separatrices1_numberOfCells{};
  std::vector<ttk::SimplexId> separatrices1_cells_connectivity;
  std::vector<ttk::SimplexId> separatrices1_cells_sourceIds;
  std::vector<ttk::SimplexId> separatrices1_cells_destinationIds;
  std::vector<ttk::SimplexId> separatrices1_cells_separatrixIds;
  std::vector<char> separatrices1_cells_separatrixTypes;
  std::vector<char> separatrices1_cells_isOnBoundary;
  std::vector<ttk::SimplexId> separatrices1_cells_separatrixFunctionMaximaId;
  std::vector<ttk::SimplexId> separatrices1_cells_separatrixFunctionMinimaId;
  std::vector<double> separatrices1_cells_separatrixFunctionDiffs;
  // segmentation
  std::vector<ttk::SimplexId> ascendingSegmentation(
      triangulation.getNumberOfVertices(), -1),
      descendingSegmentation(triangulation.getNumberOfVertices(), -1),
      mscSegmentation(triangulation.getNumberOfVertices(), -1);
  morseSmaleComplex.preconditionTriangulation(&triangulation);
  morseSmaleComplex.setInputScalarField(simplifiedHeight.data());
  morseSmaleComplex.setInputOffsets(simplifiedOrder.data());
  morseSmaleComplex.setOutputMorseComplexes(ascendingSegmentation.data(),
                                            descendingSegmentation.data(),
                                            mscSegmentation.data());
  morseSmaleComplex.setOutputCriticalPoints(
      &criticalPoints_points,
      &criticalPoints_points_cellDimensions, &criticalPoints_points_cellIds,
      &criticalPoints_points_isOnBoundary,
      &criticalPoints_points_PLVertexIdentifiers,
      &criticalPoints_points_manifoldSize);
  morseSmaleComplex.setOutputSeparatrices1(
      &separatrices1_numberOfPoints, &separatrices1_points,
      &separatrices1_points_smoothingMask, &separatrices1_points_cellDimensions,
      &separatrices1_points_cellIds, &separatrices1_numberOfCells,
      &separatrices1_cells_connectivity, &separatrices1_cells_sourceIds,
      &separatrices1_cells_destinationIds, &separatrices1_cells_separatrixIds,
      &separatrices1_cells_separatrixTypes,
      &separatrices1_cells_separatrixFunctionMaximaId,
      &separatrices1_cells_separatrixFunctionMinimaId,
      &separatrices1_cells_isOnBoundary);

  morseSmaleComplex.execute<float>(triangulation);*/
  // 7. computing the Morse-Smale complex
  ttk::MorseSmaleComplex morseSmaleComplex;
  // critical points
  ttk::MorseSmaleComplex::OutputCriticalPoints outCriticalPoints{};
  // 1-separatrices
  ttk::MorseSmaleComplex::Output1Separatrices out1Separatrices{};
  // 2-separatrices
  ttk::MorseSmaleComplex::Output2Separatrices out2Separatrices{};
  std::vector<ttk::SimplexId> ascendingSegmentation(
      triangulation.getNumberOfVertices(), -1),
      descendingSegmentation(triangulation.getNumberOfVertices(), -1),
      mscSegmentation(triangulation.getNumberOfVertices(), -1);
  ttk::MorseSmaleComplex::OutputManifold outSegmentation{
      ascendingSegmentation.data(), descendingSegmentation.data(),
      mscSegmentation.data()};
  morseSmaleComplex.preconditionTriangulation(&triangulation);
  morseSmaleComplex.execute(
      outCriticalPoints, out1Separatrices, out2Separatrices, outSegmentation,
      simplifiedHeight.data(), 0, simplifiedOrder.data(), triangulation);
  std::cout << "Number Of CP: " << (unsigned int)outCriticalPoints.points_.size() << std::endl;
  std::cout << "Number Of SP: " << (unsigned int)out1Separatrices.pt.numberOfPoints_ << std::endl;
  std::cout << "Number Of SC: " << (unsigned int)out1Separatrices.cl.numberOfCells_ << std::endl;

  *cp_len = outCriticalPoints.cellDimensions_.size();
  std::cout << "Input CP" << std::endl;
  for (int i = 0; i < outCriticalPoints.cellDimensions_.size(); ++i)
  {
    cp_point_type[i] = outCriticalPoints.cellDimensions_[i];
    cp_coordx[i] = outCriticalPoints.points_[i][0];
    cp_coordy[i] = ydim - outCriticalPoints.points_[i][1];
    cp_coordz[i] = outCriticalPoints.points_[i][2];
    cp_cellid[i] = outCriticalPoints.cellIds_[i];
    cp_pl_vertex_identifier[i] = outCriticalPoints.PLVertexIdentifiers_[i];
    cp_manifold_size[i] = outCriticalPoints.manifoldSize_[i];
  }
  std::cout << "Input SP" << std::endl;
  *sp_len = (unsigned int)out1Separatrices.pt.numberOfPoints_;
  for (int i = 0; i < out1Separatrices.pt.numberOfPoints_; ++i)
  {
    sp_id[i] = i;
    sp_coordx[i] = out1Separatrices.pt.points_[i * 3];
    sp_coordy[i] = ydim - out1Separatrices.pt.points_[i * 3 + 1];
    sp_coordz[i] = out1Separatrices.pt.points_[i * 3 + 2];
    sp_point_type[i] = (unsigned int)out1Separatrices.pt.cellDimensions_[i];
    sp_cellid[i] = out1Separatrices.pt.cellIds_[i];
  }
  std::cout << "Input SC" << std::endl;
  *sc_len = (unsigned int)out1Separatrices.cl.numberOfCells_;
  for (int i = 0; i < out1Separatrices.cl.numberOfCells_; ++i)
  {
    sc_id[i] = (unsigned int)i;

    sc_source[i] = (unsigned int)out1Separatrices.cl.sourceIds_[i];

    sc_dest[i] = (unsigned int)out1Separatrices.cl.destinationIds_[i];

    sc_connectivity_s[i] = (unsigned int)out1Separatrices.cl.connectivity_[i * 2];

    sc_connectivity_d[i] = (unsigned int)out1Separatrices.cl.connectivity_[i * 2 + 1];

    sc_separatrix_id[i] = (unsigned int)out1Separatrices.cl.separatrixIds_[i];

    sc_separatrix_type[i] = (unsigned int)out1Separatrices.cl.separatrixTypes_[i];

    /*std::cout << "sc_f_maxima[" << i << "]" << std::endl;
    std::cout << "sc_size: " << separatrices1_cells_separatrixFunctionMaximaId.size() << std::endl;
    sc_f_maxima[i] = (unsigned int)separatrices1_cells_separatrixFunctionMaximaId[i];
    std::cout << "OK" << std::endl;
    sc_f_minima[i] = (unsigned int)separatrices1_cells_separatrixFunctionMinimaId[i];*/
    //  sc_f_diff[i] = separatrices1_cells_separatrixFunctionDiffs[i];
  }
  /*std::cout << "Separatrices1================================== Points" << std::endl;
  std::cout << separatrices1_numberOfPoints << std::endl;
  std::cout << separatrices1_points.size() << std::endl;
  std::cout << separatrices1_points_cellDimensions.size() << std::endl;
  std::cout << separatrices1_points_cellIds.size() << std::endl;*/
  /*for (int i = 0; i < separatrices1_numberOfPoints; ++i)
  {
    std::cout << "ID:" << i << std::endl;
    std::cout << "Coord: " << separatrices1_points[i * 3] << ", " << ydim - separatrices1_points[i * 3 + 1] << std::endl;
    std::cout << "Kind: " << (int)separatrices1_points_cellDimensions[i] << std::endl;
    std::cout << "Cell ID:" << separatrices1_points_cellIds[i] << std::endl;
    std::cout << std::endl;
  }*/

  /*std::cout << separatrices1_numberOfCells << std::endl;
  std::cout << separatrices1_cells_connectivity.size() << std::endl;
  std::cout << separatrices1_cells_sourceIds.size() << ", " << separatrices1_cells_destinationIds.size() << std::endl;
  std::cout << separatrices1_cells_separatrixIds.size() << std::endl;
  std::cout << separatrices1_cells_separatrixTypes.size() << std::endl;
  std::cout << separatrices1_cells_separatrixFunctionMaxima.size() << ", " << separatrices1_cells_separatrixFunctionMinima.size() << ", " << separatrices1_cells_separatrixFunctionDiffs.size() << std::endl;*/
  /*std::cout << "Separatrices1================================= Cells" << std::endl;
  for (int i = 0; i < separatrices1_numberOfCells; ++i)
  {
    std::cout << "Source: " << separatrices1_cells_sourceIds[i] << ", Dest: " << separatrices1_cells_destinationIds[i] << std::endl;
    std::cout << "Connectivity?: " << separatrices1_cells_connectivity[i * 2] << ", " << separatrices1_cells_connectivity[i * 2 + 1] << std::endl;
    std::cout << "Separatrix Id: " << separatrices1_cells_separatrixIds[i] << std::endl;
    std::cout << "Separatrix Type: " << (int)separatrices1_cells_separatrixTypes[i] << std::endl;
    std::cout << "F Maxima: " << separatrices1_cells_separatrixFunctionMaxima[i] << ", "
              << "F Minima: " << separatrices1_cells_separatrixFunctionMinima[i] << ", "
              << "F Diff: " << separatrices1_cells_separatrixFunctionDiffs[i] << std::endl
              << std::endl;
  }*/
  std::cout << "compute Simplification End" << std::endl;
  return;
}

/*void Ttk_rs::compute_critical_points(float *data_2d, unsigned int datalen, unsigned int xdim, unsigned int ydim)
{
  ttk::globalDebugLevel_ = 3;
  //validate
  if (datalen != xdim * ydim)
  {
    return;
  }
  std::vector<float> pointSet(data_2d, data_2d + datalen);
  ttk::Triangulation triangulation;
  triangulation.setInputGrid(0.0, 0.0, 0.0, 1.0, 1.0, 0.0, xdim, ydim, 1);

  std::vector<ttk::SimplexId> order(pointSet.size());
  ttk::preconditionOrderArray(pointSet.size(), pointSet.data(), order.data());
  */

/*ttk::PersistenceCurve curve;
  std::vector<std::pair<float, ttk::SimplexId>> outputCurve;
  curve.preconditionTriangulation(&triangulation);
  curve.setOutputCTPlot(&outputCurve);
  curve.execute<float>(pointSet.data(), order.data(), &triangulation);
  std::cout << "PersistentCurve======================================" << outputCurve.size() << std::endl;
  for (int i = 0; i < outputCurve.size(); ++i)
  {
    std::cout << outputCurve[i].first << "," << outputCurve[i].second << std::endl;
  }

  ttk::PersistenceDiagram diagram;
  std::vector<ttk::PersistencePair> diagramOutput;
  diagram.preconditionTriangulation(&triangulation);
  diagram.execute(diagramOutput, pointSet.data(), order.data(), &triangulation);

  std::cout << "PersistentDiagram======================================" << diagramOutput.size() << std::endl;
  for (int i = 0; i < diagramOutput.size(); ++i)
  {
    std::cout << "birth: " << diagramOutput[i].birth << ", death: " << diagramOutput[i].death << std::endl;
    std::cout << "birthType: " << (int)diagramOutput[i].birthType << ", deathType: " << (int)diagramOutput[i].deathType << std::endl;
    std::cout << "persistence: " << diagramOutput[i].persistence << std::endl;
    std::cout << "birth_scalar: " << pointSet[diagramOutput[i].birth] << ", death_scalar: " << pointSet[diagramOutput[i].death] << std::endl;
  }*/

// std::vector<float> simplifiedHeight = pointSet;
// std::vector<ttk::SimplexId> authorizedCriticalPoints, simplifiedOrder = order;
/*for (int i = 0; i < (int)diagramOutput.size(); i++)
  {
    if (diagramOutput[i].persistence > filter)
    {
      // 5. selecting the most persistent pairs
      authorizedCriticalPoints.push_back(diagramOutput[i].birth);
      authorizedCriticalPoints.push_back(diagramOutput[i].death);
    }
  }

  std::cout << "authorizedCriticalPoints======================================" << authorizedCriticalPoints.size() << std::endl;
  for (int i = 0; i < authorizedCriticalPoints.size(); i += 2)
  {
    std::cout << "birth: " << authorizedCriticalPoints[i] << ", death: " << authorizedCriticalPoints[i + 1] << std::endl;
  }*/
// 6. simplifying the input data to remove non-persistent pairs
/*ttk::TopologicalSimplification simplification;
simplification.preconditionTriangulation(&triangulation);
simplification.execute<float>(pointSet.data(), simplifiedHeight.data(),
                              authorizedCriticalPoints.data(), order.data(),
                              simplifiedOrder.data(),
                              authorizedCriticalPoints.size(), triangulation);
// 7. computing the Morse-Smale complex
ttk::MorseSmaleComplex morseSmaleComplex;
// critical points
ttk::SimplexId criticalPoints_numberOfPoints{};
std::vector<float> criticalPoints_points;
std::vector<char> criticalPoints_points_cellDimensions;
std::vector<ttk::SimplexId> criticalPoints_points_cellIds;
std::vector<char> criticalPoints_points_isOnBoundary;
std::vector<float> criticalPoints_points_cellScalars;
std::vector<ttk::SimplexId> criticalPoints_points_PLVertexIdentifiers;
std::vector<ttk::SimplexId> criticalPoints_points_manifoldSize;
// 1-separatrices
ttk::SimplexId separatrices1_numberOfPoints{};
std::vector<float> separatrices1_points;
std::vector<char> separatrices1_points_smoothingMask;
std::vector<char> separatrices1_points_cellDimensions;
std::vector<ttk::SimplexId> separatrices1_points_cellIds;
ttk::SimplexId separatrices1_numberOfCells{};
std::vector<long long> separatrices1_cells_connectivity;
std::vector<ttk::SimplexId> separatrices1_cells_sourceIds;
std::vector<ttk::SimplexId> separatrices1_cells_destinationIds;
std::vector<ttk::SimplexId> separatrices1_cells_separatrixIds;
std::vector<char> separatrices1_cells_separatrixTypes;
std::vector<char> separatrices1_cells_isOnBoundary;
std::vector<double> separatrices1_cells_separatrixFunctionMaxima;
std::vector<double> separatrices1_cells_separatrixFunctionMinima;
std::vector<double> separatrices1_cells_separatrixFunctionDiffs;
// segmentation
std::vector<ttk::SimplexId> ascendingSegmentation(
    triangulation.getNumberOfVertices(), -1),
    descendingSegmentation(triangulation.getNumberOfVertices(), -1),
    mscSegmentation(triangulation.getNumberOfVertices(), -1);
morseSmaleComplex.preconditionTriangulation(&triangulation);
morseSmaleComplex.setInputScalarField(simplifiedHeight.data());
morseSmaleComplex.setInputOffsets(simplifiedOrder.data());
morseSmaleComplex.setOutputMorseComplexes(ascendingSegmentation.data(),
                                          descendingSegmentation.data(),
                                          mscSegmentation.data());
morseSmaleComplex.setOutputCriticalPoints(
    &criticalPoints_numberOfPoints, &criticalPoints_points,
    &criticalPoints_points_cellDimensions, &criticalPoints_points_cellIds,
    &criticalPoints_points_cellScalars, &criticalPoints_points_isOnBoundary,
    &criticalPoints_points_PLVertexIdentifiers,
    &criticalPoints_points_manifoldSize);
morseSmaleComplex.setOutputSeparatrices1(
    &separatrices1_numberOfPoints, &separatrices1_points,
    &separatrices1_points_smoothingMask, &separatrices1_points_cellDimensions,
    &separatrices1_points_cellIds, &separatrices1_numberOfCells,
    &separatrices1_cells_connectivity, &separatrices1_cells_sourceIds,
    &separatrices1_cells_destinationIds, &separatrices1_cells_separatrixIds,
    &separatrices1_cells_separatrixTypes,
    &separatrices1_cells_separatrixFunctionMaxima,
    &separatrices1_cells_separatrixFunctionMinima,
    &separatrices1_cells_separatrixFunctionDiffs,
    &separatrices1_cells_isOnBoundary);

morseSmaleComplex.execute<float>(triangulation);*/
/*for (int i = 0; i < criticalPoints_points_cellDimensions.size(); ++i)
  {
    std::cout << "Critical Points #" << i << "===================" << std::endl;
    std::cout << "Kind: ";
    switch (criticalPoints_points_cellDimensions[i])
    {
    case 0:
      std::cout << "Minimum" << std::endl;
      break;
    case 1:
      std::cout << "Saddle" << std::endl;
      break;
    case 2:
      std::cout << "Maximum" << std::endl;
      break;
    default:
      std::cout << "Unreachable!" << std::endl;
      break;
    }
    std::cout << "Coord: " << criticalPoints_points[i * 3] << ", " << ydim - criticalPoints_points[i * 3 + 1] << std::endl;
    std::cout << "Value: " << criticalPoints_points_cellScalars[i] << std::endl;
    std::cout << "CellId: " << criticalPoints_points_cellIds[i] << std::endl;
    std::cout << "PLVertexIdentifier: " << criticalPoints_points_PLVertexIdentifiers[i] << std::endl;
    std::cout << "manifoldSize: " << criticalPoints_points_manifoldSize[i] << std::endl;
  }*/
/*std::cout << "Separatrices1================================== Points" << std::endl;
  std::cout << separatrices1_numberOfPoints << std::endl;
  std::cout << separatrices1_points.size() << std::endl;
  std::cout << separatrices1_points_cellDimensions.size() << std::endl;
  std::cout << separatrices1_points_cellIds.size() << std::endl;*/
/*for (int i = 0; i < separatrices1_numberOfPoints; ++i)
  {
    std::cout << "ID:" << i << std::endl;
    std::cout << "Coord: " << separatrices1_points[i * 3] << ", " << ydim - separatrices1_points[i * 3 + 1] << std::endl;
    std::cout << "Kind: " << (int)separatrices1_points_cellDimensions[i] << std::endl;
    std::cout << "Cell ID:" << separatrices1_points_cellIds[i] << std::endl;
    std::cout << std::endl;
  }*/

/*std::cout << separatrices1_numberOfCells << std::endl;
  std::cout << separatrices1_cells_connectivity.size() << std::endl;
  std::cout << separatrices1_cells_sourceIds.size() << ", " << separatrices1_cells_destinationIds.size() << std::endl;
  std::cout << separatrices1_cells_separatrixIds.size() << std::endl;
  std::cout << separatrices1_cells_separatrixTypes.size() << std::endl;
  std::cout << separatrices1_cells_separatrixFunctionMaxima.size() << ", " << separatrices1_cells_separatrixFunctionMinima.size() << ", " << separatrices1_cells_separatrixFunctionDiffs.size() << std::endl;*/
/*std::cout << "Separatrices1================================= Cells" << std::endl;
  for (int i = 0; i < separatrices1_numberOfCells; ++i)
  {
    std::cout << "Source: " << separatrices1_cells_sourceIds[i] << ", Dest: " << separatrices1_cells_destinationIds[i] << std::endl;
    std::cout << "Connectivity?: " << separatrices1_cells_connectivity[i * 2] << ", " << separatrices1_cells_connectivity[i * 2 + 1] << std::endl;
    std::cout << "Separatrix Id: " << separatrices1_cells_separatrixIds[i] << std::endl;
    std::cout << "Separatrix Type: " << (int)separatrices1_cells_separatrixTypes[i] << std::endl;
    std::cout << "F Maxima: " << separatrices1_cells_separatrixFunctionMaxima[i] << ", "
              << "F Minima: " << separatrices1_cells_separatrixFunctionMinima[i] << ", "
              << "F Diff: " << separatrices1_cells_separatrixFunctionDiffs[i] << std::endl
              << std::endl;
  }*/

/*ttk::ContourTree contour_tree(pointSet.data(), criticalPoints_points_manifoldSize.data(),
                                criticalPoints_points_cellIds.data(),
                                separatrices1_numberOfPoints, separatrices1_numberOfCells);*/
//}

void Ttk_rs::data_handling(float *data1, char *data2, float *data3, unsigned int data3_len)
{
  std::string cpp_string = data2;
  std::vector<float> cpp_vector(data3, data3 + data3_len);
  std::cout << "in C++, data1 = " << *data1 << std::endl;
  std::cout << "in C++, data2 = " << cpp_string << std::endl;
  std::cout << "in C++, data3 = [";
  for (int i = 0; i < cpp_vector.size(); ++i)
  {
    std::cout << std::fixed << std::setprecision(1) << cpp_vector[i];
    if (i == cpp_vector.size() - 1)
    {
      std::cout << "]" << std::endl;
    }
    else
    {
      std::cout << ", ";
    }
  }
  // Rewriting data
  cpp_vector.push_back(3.0);
  cpp_string = u8"abcdefghijk";
  std::cout << cpp_string << std::endl;
  *data1 = 2.5f;
  for (int i = 0; i < cpp_string.size(); ++i)
  {
    data2[i] = cpp_string[i];
  }
  for (int i = 0; i < cpp_vector.size(); ++i)
  {
    data3[i] = cpp_vector[i];
  }
  return;
}

/*myvec<float, long long int> Ttk_rs::load_off_file(const std::string &inputPath,
                                                  float *pointSet, int pointSetSize,
                                                  long long int *trianglesetCo, int trianglesetCoSize,
                                                  long long int *trianglesetOff, int trianglesetOffSize)
{
  ttk::globalDebugLevel_ = 3;
  // load some terrain from some OFF file.
  if (inputPath.empty())
  {
    Ttk_rs::rett.state = -1;
    return Ttk_rs::rett;
  }

  ttk::Debug dbg;
  dbg.setDebugLevel(ttk::globalDebugLevel_);
  dbg.setDebugMsgPrefix("main::load");
  dbg.printMsg("Reading input mesh...");
  int vertexNumber = 0, triangleNumber = 0;
  std::string keyword;
  std::ifstream f;
  f.open(inputPath.data(), std::ios_base::in);

  if (!f)
  {
    dbg.printErr("Cannot read file `" + inputPath + "'!");
    Ttk_rs::rett.state = -1;
    return Ttk_rs::rett;
  }
  f >> keyword;
  std::cout << "keyword:" << keyword << std::endl;
  if (keyword != "OFF")
  {
    dbg.printErr("Input OFF file `" + inputPath + "' seems invalid :(");
    Ttk_rs::rett.state = -2;
    return Ttk_rs::rett;
  }
  f >> vertexNumber;
  f >> triangleNumber;
  f >> keyword;

  std::vector<float> actpointSet(pointSet, pointSet + pointSetSize);
  std::vector<long long int> acttrianglesetCo(trianglesetCo, trianglesetCo + trianglesetCoSize);
  std::vector<long long int> acttrianglesetOff(trianglesetOff, trianglesetOff + trianglesetOffSize);

  actpointSet.resize(3 * vertexNumber);
  acttrianglesetCo.resize(3 * triangleNumber);
  acttrianglesetOff.resize(triangleNumber + 1);
  std::cout << "vertexNumber:" << vertexNumber << std::endl;
  std::cout << "triangleNumber:" << triangleNumber << std::endl;
  std::cout << "keyword:" << keyword << std::endl;
  for (int i = 0; i < 3 * vertexNumber; i++)
  {
    f >> actpointSet[i];
  }
  int offId = 0;
  int coId = 0;
  for (int i = 0; i < triangleNumber; i++)
  {
    int cellSize;
    f >> cellSize;
    if (cellSize != 3)
    {
      std::cerr << "cell size " << cellSize << " != 3" << std::endl;
      Ttk_rs::rett.state = -3;
      return Ttk_rs::rett;
    }
    acttrianglesetOff[offId++] = coId;
    for (int j = 0; j < 3; j++)
    {
      int cellId;
      f >> cellId;
      acttrianglesetCo[coId++] = cellId;
    }
  }
  acttrianglesetOff[offId] = coId; // the last one

  f.close();

  for (int i = 0; i < actpointSet.size(); ++i)
  {
    pointSet[i] = actpointSet[i];
  }
  for (int i = 0; i < acttrianglesetCo.size(); ++i)
  {
    trianglesetCo[i] = acttrianglesetCo[i];
  }
  for (int i = 0; i < acttrianglesetOff.size(); ++i)
  {
    trianglesetOff[i] = acttrianglesetOff[i];
  }
  Ttk_rs::rett.ptr = pointSet;
  Ttk_rs::rett.length = actpointSet.size();
  Ttk_rs::rett.ptr2 = trianglesetCo;
  Ttk_rs::rett.length2 = acttrianglesetCo.size();
  Ttk_rs::rett.ptr3 = trianglesetOff;
  Ttk_rs::rett.length3 = acttrianglesetOff.size();
  Ttk_rs::rett.state = 0;
  dbg.printMsg("... done! (read " + std::to_string(vertexNumber) + " vertices, " + std::to_string(triangleNumber) + " triangles)");
  return Ttk_rs::rett;
}

void Ttk_rs::test_ttk_processing1(float *pointSet, int pointSetSize,
                                  long long int *trianglesetCo, int trianglesetCoSize,
                                  long long int *trianglesetOff, int trianglesetOffSize)
{
  std::cout << "process1" << std::endl;
  ttk::Debug dbg;
  dbg.setDebugLevel(ttk::globalDebugLevel_);
  std::vector<float> actpointSet(pointSet, pointSet + pointSetSize);
  actpointSet.resize(pointSetSize);
  std::vector<long long int> acttrianglesetCo(trianglesetCo, trianglesetCo + trianglesetCoSize);
  std::vector<long long int> acttrianglesetOff(trianglesetOff, trianglesetOff + trianglesetOffSize);
  acttrianglesetCo.resize(trianglesetCoSize);
  acttrianglesetOff.resize(trianglesetOffSize);
  std::cout << actpointSet.size() / 3 << ", " << actpointSet.data() << std::endl;
  triangulation.setInputPoints(actpointSet.size() / 3, actpointSet.data());
  long long int triangleNumber = acttrianglesetOff.size() - 1;
#ifdef TTK_CELL_ARRAY_NEW
  triangulation.setInputCells(
      triangleNumber, trianglesetCo, trianglesetOff);
#else
  LongSimplexId *triangleSet;
  CellArray::TranslateToFlatLayout(acttrianglesetCo, acttrianglesetOff, triangleset);
  triangulation.setInputCells(triangleNumber, triangleSet);
#endif
}

void Ttk_rs::test_ttk_processing2(float *pointSet, int pointSetSize)
{
  std::cout << "process2" << std::endl;
  std::vector<float> actpointSet(pointSet, pointSet + pointSetSize);
  actpointSet.resize(pointSetSize);

  std::vector<float> height(actpointSet.size() / 3);
  int vertexId = 0;
  // use the z-coordinate here
  for (int i = 2; i < (int)actpointSet.size(); i += 3)
  {
    height[vertexId] = actpointSet[i];
    vertexId++;
  }
  // order array: every vertex sorted according to the elevation field
  std::vector<ttk::SimplexId> order(height.size());
  // precondition/fill in the order array according to the elevation field
  ttk::preconditionOrderArray(height.size(), height.data(), order.data());

  // 2. computing the persistence curve
  ttk::PersistenceCurve curve;
  std::vector<std::pair<float, ttk::SimplexId>> outputCurve;
  curve.preconditionTriangulation(&triangulation);
  curve.setOutputCTPlot(&outputCurve);
  curve.execute<float>(height.data(), order.data(), &triangulation);
  return;
}

int save(const std::vector<float> &pointSet,
         const std::vector<long long int> &triangleSetCo,
         const std::vector<long long int> &triangleSetOff,
         const std::string &outputPath)
{

  // save the simplified terrain in some OFF file
  std::string fileName(outputPath);

  std::ofstream f(fileName.data(), std::ios::out);

  if (!f)
  {
    ttk::Debug dbg;
    dbg.setDebugLevel(ttk::globalDebugLevel_);
    dbg.setDebugMsgPrefix("main::save");
    dbg.printErr("Could not write output file `" + fileName + "'!");
    return -1;
  }

  const int nbTriangles = triangleSetOff.size() - 1;

  f << "OFF" << std::endl;
  f << pointSet.size() / 3 << " " << nbTriangles << " 0" << std::endl;

  for (int i = 0; i < (int)pointSet.size() / 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      f << pointSet[3 * i + j];
      f << " ";
    }
    f << std::endl;
  }

  for (int i = 0; i < nbTriangles; i++)
  {
    int cellSize = triangleSetOff[i + 1] - triangleSetOff[i];
    assert(cellSize == 3);
    f << cellSize << " ";
    for (int j = triangleSetOff[i]; j < triangleSetOff[i + 1]; j++)
    {
      f << triangleSetCo[j];
      f << " ";
    }
    f << std::endl;
  }

  f.close();

  return 0;
}

int hello(int argc, char **argv)
{

  std::string inputFilePath;
  ttk::CommandLineParser parser;

  ttk::globalDebugLevel_ = 3;

  // register the arguments to the command line parser
  parser.setArgument("i", &inputFilePath, "Path to input OFF file");
  // parse
  parser.parse(argc, argv);

  std::vector<float> pointSet;
  std::vector<long long int> triangleSetCo, triangleSetOff;
  ttk::Triangulation triangulation;

  // load the input
  load(inputFilePath, pointSet, triangleSetCo, triangleSetOff);
  triangulation.setInputPoints(pointSet.size() / 3, pointSet.data());
  long long int triangleNumber = triangleSetOff.size() - 1;
#ifdef TTK_CELL_ARRAY_NEW
  triangulation.setInputCells(
      triangleNumber, triangleSetCo.data(), triangleSetOff.data());
#else
  //LongSimplexId *triangleSet;
  //CellArray::TranslateToFlatLayout(triangleSetCo, triangleSetOff, triangleSet);
  //triangulation.setInputCells(triangleNumber, triangleSet);
#endif

  // NOW, do the TTK processing

  // computing some elevation
  std::vector<float> height(pointSet.size() / 3);
  int vertexId = 0;
  // use the z-coordinate here
  for (int i = 2; i < (int)pointSet.size(); i += 3)
  {
    height[vertexId] = pointSet[i];
    vertexId++;
  }

  // order array: every vertex sorted according to the elevation field
  std::vector<ttk::SimplexId> order(height.size());
  // precondition/fill in the order array according to the elevation field
  ttk::preconditionOrderArray(height.size(), height.data(), order.data());

  // 2. computing the persistence curve
  ttk::PersistenceCurve curve;
  std::vector<std::pair<float, ttk::SimplexId>> outputCurve;
  curve.preconditionTriangulation(&triangulation);
  curve.setOutputCTPlot(&outputCurve);
  curve.execute<float>(height.data(), order.data(), &triangulation);

  // 3. computing the persitence diagram
  ttk::PersistenceDiagram diagram;
  std::vector<ttk::PersistencePair> diagramOutput;
  diagram.preconditionTriangulation(&triangulation);
  diagram.execute(diagramOutput, height.data(), order.data(), &triangulation);

  // 4. selecting the critical point pairs
  std::vector<float> simplifiedHeight = height;
  std::vector<ttk::SimplexId> authorizedCriticalPoints, simplifiedOrder = order;
  for (int i = 0; i < (int)diagramOutput.size(); i++)
  {
    if (diagramOutput[i].persistence > 0.05)
    {
      // 5. selecting the most persistent pairs
      authorizedCriticalPoints.push_back(diagramOutput[i].birth);
      authorizedCriticalPoints.push_back(diagramOutput[i].death);
    }
  }

  // 6. simplifying the input data to remove non-persistent pairs
  ttk::TopologicalSimplification simplification;
  simplification.preconditionTriangulation(&triangulation);
  simplification.execute<float>(height.data(), simplifiedHeight.data(),
                                authorizedCriticalPoints.data(), order.data(),
                                simplifiedOrder.data(),
                                authorizedCriticalPoints.size(), triangulation);

  // assign the simplified values to the input mesh
  for (int i = 0; i < (int)simplifiedHeight.size(); i++)
  {
    pointSet[3 * i + 2] = simplifiedHeight[i];
  }

  // 7. computing the Morse-Smale complex
  ttk::MorseSmaleComplex morseSmaleComplex;
  // critical points
  ttk::SimplexId criticalPoints_numberOfPoints{};
  std::vector<float> criticalPoints_points;
  std::vector<char> criticalPoints_points_cellDimensions;
  std::vector<ttk::SimplexId> criticalPoints_points_cellIds;
  std::vector<char> criticalPoints_points_isOnBoundary;
  std::vector<float> criticalPoints_points_cellScalars;
  std::vector<ttk::SimplexId> criticalPoints_points_PLVertexIdentifiers;
  std::vector<ttk::SimplexId> criticalPoints_points_manifoldSize;
  // 1-separatrices
  ttk::SimplexId separatrices1_numberOfPoints{};
  std::vector<float> separatrices1_points;
  std::vector<char> separatrices1_points_smoothingMask;
  std::vector<char> separatrices1_points_cellDimensions;
  std::vector<ttk::SimplexId> separatrices1_points_cellIds;
  ttk::SimplexId separatrices1_numberOfCells{};
  std::vector<long long> separatrices1_cells_connectivity;
  std::vector<ttk::SimplexId> separatrices1_cells_sourceIds;
  std::vector<ttk::SimplexId> separatrices1_cells_destinationIds;
  std::vector<ttk::SimplexId> separatrices1_cells_separatrixIds;
  std::vector<char> separatrices1_cells_separatrixTypes;
  std::vector<char> separatrices1_cells_isOnBoundary;
  std::vector<double> separatrices1_cells_separatrixFunctionMaxima;
  std::vector<double> separatrices1_cells_separatrixFunctionMinima;
  std::vector<double> separatrices1_cells_separatrixFunctionDiffs;
  // segmentation
  std::vector<ttk::SimplexId> ascendingSegmentation(
      triangulation.getNumberOfVertices(), -1),
      descendingSegmentation(triangulation.getNumberOfVertices(), -1),
      mscSegmentation(triangulation.getNumberOfVertices(), -1);
  morseSmaleComplex.preconditionTriangulation(&triangulation);
  morseSmaleComplex.setInputScalarField(simplifiedHeight.data());
  morseSmaleComplex.setInputOffsets(simplifiedOrder.data());
  morseSmaleComplex.setOutputMorseComplexes(ascendingSegmentation.data(),
                                            descendingSegmentation.data(),
                                            mscSegmentation.data());
  morseSmaleComplex.setOutputCriticalPoints(
      &criticalPoints_numberOfPoints, &criticalPoints_points,
      &criticalPoints_points_cellDimensions, &criticalPoints_points_cellIds,
      &criticalPoints_points_cellScalars, &criticalPoints_points_isOnBoundary,
      &criticalPoints_points_PLVertexIdentifiers,
      &criticalPoints_points_manifoldSize);
  morseSmaleComplex.setOutputSeparatrices1(
      &separatrices1_numberOfPoints, &separatrices1_points,
      &separatrices1_points_smoothingMask, &separatrices1_points_cellDimensions,
      &separatrices1_points_cellIds, &separatrices1_numberOfCells,
      &separatrices1_cells_connectivity, &separatrices1_cells_sourceIds,
      &separatrices1_cells_destinationIds, &separatrices1_cells_separatrixIds,
      &separatrices1_cells_separatrixTypes,
      &separatrices1_cells_separatrixFunctionMaxima,
      &separatrices1_cells_separatrixFunctionMinima,
      &separatrices1_cells_separatrixFunctionDiffs,
      &separatrices1_cells_isOnBoundary);

  morseSmaleComplex.execute<float>(triangulation);

  // save the output
  save(pointSet, triangleSetCo, triangleSetOff, "output.off");

  return 0;
}*/
