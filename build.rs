extern crate bindgen;

use cmake;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("cargo:rerun-if-changed=src/ttk_cpp/ttk_cpp.cpp");
    println!("cargo:rerun-if-changed=src/ttk_cpp/ttk_cpp.hpp");
    println!("cargo:rerun-if-changed=src/ttk_cpp/CMakeList.txt");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    let out_path = out_path.to_str().unwrap();
    let mut ttk_libraries_path = String::new();
    let mut libomp_path = String::new();

    for result in BufReader::new(File::open("Settings.txt")?).lines() {
        let l = result?;
        let v: Vec<&str> = l.split('=').collect();
        if v[0] == "TTK_LIBRARIES" {
            ttk_libraries_path = v[1].to_string();
        } else if v[0] == "LIBOMP" {
            libomp_path = v[1].to_string();
        }
    }

    println!("cargo:rustc-link-search={}/build/Release/", out_path); //Search for ttk_cpp.lib
    println!("cargo:rustc-link-search=native={}", ttk_libraries_path); //Search for ttk*.lib
    println!("cargo:rustc-link-search=native={}", libomp_path); //Search for libomp.lib
                                                                //println!("cargo:rustc-link-search=native=C:/Users/kstrb/repo/PCL/lib/windows/x64");
    println!("cargo:rustc-link-lib=static=ttk_cpp");
    /*println!("cargo:rustc-link-lib=static=ttkAlgorithm");
    println!("cargo:rustc-link-lib=static=ttkAlgorithmCS");
    println!("cargo:rustc-link-lib=static=ttkArrayEditor");
    println!("cargo:rustc-link-lib=static=ttkArrayEditorCS");
    println!("cargo:rustc-link-lib=static=ttkArrayPreconditioning");
    println!("cargo:rustc-link-lib=static=ttkArrayPreconditioningCS");
    println!("cargo:rustc-link-lib=static=ttkBarycentricSubdivision");
    println!("cargo:rustc-link-lib=static=ttkBarycentricSubdivisionCS");
    println!("cargo:rustc-link-lib=static=ttkBaseabstractMorseSmaleComplex");
    println!("cargo:rustc-link-lib=static=ttkBaseabstractTriangulation");
    println!("cargo:rustc-link-lib=static=ttkBaseauction");
    println!("cargo:rustc-link-lib=static=ttkBasebarycentricSubdivision");
    println!("cargo:rustc-link-lib=static=ttkBasebottleneckDistance");
    println!("cargo:rustc-link-lib=static=ttkBaseboundingVolumeHierarchy");
    println!("cargo:rustc-link-lib=static=ttkBasecinemaImaging");
    println!("cargo:rustc-link-lib=static=ttkBasecinemaQuery");
    println!("cargo:rustc-link-lib=static=ttkBasecommon");
    println!("cargo:rustc-link-lib=static=ttkBasecontinuousScatterPlot");
    println!("cargo:rustc-link-lib=static=ttkBasecontourAroundPoint");
    println!("cargo:rustc-link-lib=static=ttkBasecontourForests");
    println!("cargo:rustc-link-lib=static=ttkBasecontourTree");
    println!("cargo:rustc-link-lib=static=ttkBasecontourTreeAlignment");
    println!("cargo:rustc-link-lib=static=ttkBasedepthImageBasedGeometryApproximation");
    println!("cargo:rustc-link-lib=static=ttkBasedimensionReduction");
    println!("cargo:rustc-link-lib=static=ttkBasediscreteGradient");
    println!("cargo:rustc-link-lib=static=ttkBasedistanceField");
    println!("cargo:rustc-link-lib=static=ttkBaseeigenField");
    println!("cargo:rustc-link-lib=static=ttkBaseexplicitTriangulation");
    println!("cargo:rustc-link-lib=static=ttkBasefiberSurface");
    println!("cargo:rustc-link-lib=static=ttkBaseftmTree");
    println!("cargo:rustc-link-lib=static=ttkBaseftmTreePP");
    println!("cargo:rustc-link-lib=static=ttkBaseftrGraph");
    println!("cargo:rustc-link-lib=static=ttkBasegaussianPointCloud");
    println!("cargo:rustc-link-lib=static=ttkBasegeometry");
    println!("cargo:rustc-link-lib=static=ttkBaseharmonicField");
    println!("cargo:rustc-link-lib=static=ttkBasehelloWorld");
    println!("cargo:rustc-link-lib=static=ttkBaseicosphere");
    println!("cargo:rustc-link-lib=static=ttkBaseimplicitTriangulation");
    println!("cargo:rustc-link-lib=static=ttkBaseintegralLines");
    println!("cargo:rustc-link-lib=static=ttkBasejacobiSet");
    println!("cargo:rustc-link-lib=static=ttkBasekdTree");
    println!("cargo:rustc-link-lib=static=ttkBaselaplacian");
    println!("cargo:rustc-link-lib=static=ttkBaselDistance");
    println!("cargo:rustc-link-lib=static=ttkBaselDistanceMatrix");
    println!("cargo:rustc-link-lib=static=ttkBaselocalizedTopologicalSimplification");
    println!("cargo:rustc-link-lib=static=ttkBaselowestCommonAncestor");
    println!("cargo:rustc-link-lib=static=ttkBasemandatoryCriticalPoints");
    println!("cargo:rustc-link-lib=static=ttkBasemanifoldCheck");
    println!("cargo:rustc-link-lib=static=ttkBasemeshGraph");
    println!("cargo:rustc-link-lib=static=ttkBasemorseSmaleComplex");
    println!("cargo:rustc-link-lib=static=ttkBasemorseSmaleComplex2D");
    println!("cargo:rustc-link-lib=static=ttkBasemorseSmaleComplex3D");
    println!("cargo:rustc-link-lib=static=ttkBasemorseSmaleQuadrangulation");
    println!("cargo:rustc-link-lib=static=ttkBaseperiodicImplicitTriangulation");
    println!("cargo:rustc-link-lib=static=ttkBasepersistenceCurve");
    println!("cargo:rustc-link-lib=static=ttkBasepersistenceDiagram");
    println!("cargo:rustc-link-lib=static=ttkBasepersistenceDiagramClustering");
    println!("cargo:rustc-link-lib=static=ttkBasepersistenceDiagramDistanceMatrix");
    println!("cargo:rustc-link-lib=static=ttkBaseplanarGraphLayout");
    println!("cargo:rustc-link-lib=static=ttkBasequadrangulationSubdivision");
    println!("cargo:rustc-link-lib=static=ttkBaserangeDrivenOctree");
    println!("cargo:rustc-link-lib=static=ttkBasereebSpace");
    println!("cargo:rustc-link-lib=static=ttkBasescalarFieldCriticalPoints");
    println!("cargo:rustc-link-lib=static=ttkBasescalarFieldSmoother");
    println!("cargo:rustc-link-lib=static=ttkBaseskeleton");
    println!("cargo:rustc-link-lib=static=ttkBasetopologicalCompression");
    println!("cargo:rustc-link-lib=static=ttkBasetopologicalSimplification");
    println!("cargo:rustc-link-lib=static=ttkBasetrackingFromFields");
    println!("cargo:rustc-link-lib=static=ttkBasetrackingFromOverlap");
    println!("cargo:rustc-link-lib=static=ttkBasetrackingFromPersistenceDiagrams");
    println!("cargo:rustc-link-lib=static=ttkBasetriangulation");
    println!("cargo:rustc-link-lib=static=ttkBaseuncertainDataEstimator");
    println!("cargo:rustc-link-lib=static=ttkBaseunionFind");
    println!("cargo:rustc-link-lib=static=ttkBlank");
    println!("cargo:rustc-link-lib=static=ttkBlankCS");
    println!("cargo:rustc-link-lib=static=ttkBlockAggregator");
    println!("cargo:rustc-link-lib=static=ttkBlockAggregatorCS");
    println!("cargo:rustc-link-lib=static=ttkBottleneckDistance");
    println!("cargo:rustc-link-lib=static=ttkBottleneckDistanceCS");
    println!("cargo:rustc-link-lib=static=ttkCinemaDarkroom");
    println!("cargo:rustc-link-lib=static=ttkCinemaDarkroomCS");
    println!("cargo:rustc-link-lib=static=ttkCinemaImaging");
    println!("cargo:rustc-link-lib=static=ttkCinemaImagingCS");
    println!("cargo:rustc-link-lib=static=ttkCinemaProductReader");
    println!("cargo:rustc-link-lib=static=ttkCinemaProductReaderCS");
    println!("cargo:rustc-link-lib=static=ttkCinemaQuery");
    println!("cargo:rustc-link-lib=static=ttkCinemaQueryCS");
    println!("cargo:rustc-link-lib=static=ttkCinemaReader");
    println!("cargo:rustc-link-lib=static=ttkCinemaReaderCS");
    println!("cargo:rustc-link-lib=static=ttkCinemaWriter");
    println!("cargo:rustc-link-lib=static=ttkCinemaWriterCS");
    println!("cargo:rustc-link-lib=static=ttkComponentSize");
    println!("cargo:rustc-link-lib=static=ttkComponentSizeCS");
    println!("cargo:rustc-link-lib=static=ttkContinuousScatterPlot");
    println!("cargo:rustc-link-lib=static=ttkContinuousScatterPlotCS");
    println!("cargo:rustc-link-lib=static=ttkContourAroundPoint");
    println!("cargo:rustc-link-lib=static=ttkContourAroundPointCS");
    println!("cargo:rustc-link-lib=static=ttkContourForests");
    println!("cargo:rustc-link-lib=static=ttkContourForestsCS");
    println!("cargo:rustc-link-lib=static=ttkContourTreeAlignment");
    println!("cargo:rustc-link-lib=static=ttkContourTreeAlignmentCS");
    println!("cargo:rustc-link-lib=static=ttkDataSetInterpolator");
    println!("cargo:rustc-link-lib=static=ttkDataSetInterpolatorCS");
    println!("cargo:rustc-link-lib=static=ttkDataSetToTable");
    println!("cargo:rustc-link-lib=static=ttkDataSetToTableCS");
    println!("cargo:rustc-link-lib=static=ttkDepthImageBasedGeometryApproximation");
    println!("cargo:rustc-link-lib=static=ttkDepthImageBasedGeometryApproximationCS");
    println!("cargo:rustc-link-lib=static=ttkDimensionReduction");
    println!("cargo:rustc-link-lib=static=ttkDimensionReductionCS");
    println!("cargo:rustc-link-lib=static=ttkDiscreteGradient");
    println!("cargo:rustc-link-lib=static=ttkDiscreteGradientCS");
    println!("cargo:rustc-link-lib=static=ttkDistanceField");
    println!("cargo:rustc-link-lib=static=ttkDistanceFieldCS");
    println!("cargo:rustc-link-lib=static=ttkEigenField");
    println!("cargo:rustc-link-lib=static=ttkEigenFieldCS");
    println!("cargo:rustc-link-lib=static=ttkEndFor");
    println!("cargo:rustc-link-lib=static=ttkEndForCS");
    println!("cargo:rustc-link-lib=static=ttkExtract");
    println!("cargo:rustc-link-lib=static=ttkExtractCS");
    println!("cargo:rustc-link-lib=static=ttkFiber");
    println!("cargo:rustc-link-lib=static=ttkFiberCS");
    println!("cargo:rustc-link-lib=static=ttkFiberSurface");
    println!("cargo:rustc-link-lib=static=ttkFiberSurfaceCS");
    println!("cargo:rustc-link-lib=static=ttkForEach");
    println!("cargo:rustc-link-lib=static=ttkForEachCS");
    println!("cargo:rustc-link-lib=static=ttkFTMTree");
    println!("cargo:rustc-link-lib=static=ttkFTMTreeCS");
    println!("cargo:rustc-link-lib=static=ttkFTRGraph");
    println!("cargo:rustc-link-lib=static=ttkFTRGraphCS");
    println!("cargo:rustc-link-lib=static=ttkGaussianPointCloud");
    println!("cargo:rustc-link-lib=static=ttkGaussianPointCloudCS");
    println!("cargo:rustc-link-lib=static=ttkGeometrySmoother");
    println!("cargo:rustc-link-lib=static=ttkGeometrySmootherCS");
    println!("cargo:rustc-link-lib=static=ttkGridLayout");
    println!("cargo:rustc-link-lib=static=ttkGridLayoutCS");
    println!("cargo:rustc-link-lib=static=ttkHarmonicField");
    println!("cargo:rustc-link-lib=static=ttkHarmonicFieldCS");
    println!("cargo:rustc-link-lib=static=ttkHelloWorld");
    println!("cargo:rustc-link-lib=static=ttkHelloWorldCS");
    println!("cargo:rustc-link-lib=static=ttkIcosphere");
    println!("cargo:rustc-link-lib=static=ttkIcosphereCS");
    println!("cargo:rustc-link-lib=static=ttkIcosphereFromObject");
    println!("cargo:rustc-link-lib=static=ttkIcosphereFromObjectCS");
    println!("cargo:rustc-link-lib=static=ttkIcospheresFromPoints");
    println!("cargo:rustc-link-lib=static=ttkIcospheresFromPointsCS");
    println!("cargo:rustc-link-lib=static=ttkIdentifierRandomizer");
    println!("cargo:rustc-link-lib=static=ttkIdentifierRandomizerCS");
    println!("cargo:rustc-link-lib=static=ttkIdentifiers");
    println!("cargo:rustc-link-lib=static=ttkIdentifiersCS");
    println!("cargo:rustc-link-lib=static=ttkIdentifyByScalarField");
    println!("cargo:rustc-link-lib=static=ttkIdentifyByScalarFieldCS");
    println!("cargo:rustc-link-lib=static=ttkImportEmbeddingFromTable");
    println!("cargo:rustc-link-lib=static=ttkImportEmbeddingFromTableCS");
    println!("cargo:rustc-link-lib=static=ttkIntegralLines");
    println!("cargo:rustc-link-lib=static=ttkIntegralLinesCS");
    println!("cargo:rustc-link-lib=static=ttkJacobiSet");
    println!("cargo:rustc-link-lib=static=ttkJacobiSetCS");
    println!("cargo:rustc-link-lib=static=ttkLDistance");
    println!("cargo:rustc-link-lib=static=ttkLDistanceCS");
    println!("cargo:rustc-link-lib=static=ttkLDistanceMatrix");
    println!("cargo:rustc-link-lib=static=ttkLDistanceMatrixCS");
    println!("cargo:rustc-link-lib=static=ttkMandatoryCriticalPoints");
    println!("cargo:rustc-link-lib=static=ttkMandatoryCriticalPointsCS");
    println!("cargo:rustc-link-lib=static=ttkManifoldCheck");
    println!("cargo:rustc-link-lib=static=ttkManifoldCheckCS");
    println!("cargo:rustc-link-lib=static=ttkMatrixToHeatMap");
    println!("cargo:rustc-link-lib=static=ttkMatrixToHeatMapCS");
    println!("cargo:rustc-link-lib=static=ttkMeshGraph");
    println!("cargo:rustc-link-lib=static=ttkMeshGraphCS");
    println!("cargo:rustc-link-lib=static=ttkMeshSubdivision");
    println!("cargo:rustc-link-lib=static=ttkMeshSubdivisionCS");
    println!("cargo:rustc-link-lib=static=ttkMorseSmaleComplex");
    println!("cargo:rustc-link-lib=static=ttkMorseSmaleComplexCS");
    println!("cargo:rustc-link-lib=static=ttkMorseSmaleQuadrangulation");
    println!("cargo:rustc-link-lib=static=ttkMorseSmaleQuadrangulationCS");
    println!("cargo:rustc-link-lib=static=ttkOBJWriter");
    println!("cargo:rustc-link-lib=static=ttkOBJWriterCS");
    println!("cargo:rustc-link-lib=static=ttkOFFReader");
    println!("cargo:rustc-link-lib=static=ttkOFFReaderCS");
    println!("cargo:rustc-link-lib=static=ttkOFFWriter");
    println!("cargo:rustc-link-lib=static=ttkOFFWriterCS");
    println!("cargo:rustc-link-lib=static=ttkPeriodicGrid");
    println!("cargo:rustc-link-lib=static=ttkPeriodicGridCS");
    println!("cargo:rustc-link-lib=static=ttkPersistenceCurve");
    println!("cargo:rustc-link-lib=static=ttkPersistenceCurveCS");
    println!("cargo:rustc-link-lib=static=ttkPersistenceDiagram");
    println!("cargo:rustc-link-lib=static=ttkPersistenceDiagramClustering");
    println!("cargo:rustc-link-lib=static=ttkPersistenceDiagramClusteringCS");
    println!("cargo:rustc-link-lib=static=ttkPersistenceDiagramCS");
    println!("cargo:rustc-link-lib=static=ttkPersistenceDiagramDistanceMatrix");
    println!("cargo:rustc-link-lib=static=ttkPersistenceDiagramDistanceMatrixCS");
    println!("cargo:rustc-link-lib=static=ttkPlanarGraphLayout");
    println!("cargo:rustc-link-lib=static=ttkPlanarGraphLayoutCS");
    println!("cargo:rustc-link-lib=static=ttkPointDataConverter");
    println!("cargo:rustc-link-lib=static=ttkPointDataConverterCS");
    println!("cargo:rustc-link-lib=static=ttkPointDataSelector");
    println!("cargo:rustc-link-lib=static=ttkPointDataSelectorCS");
    println!("cargo:rustc-link-lib=static=ttkPointMerger");
    println!("cargo:rustc-link-lib=static=ttkPointMergerCS");
    println!("cargo:rustc-link-lib=static=ttkPointSetToCurve");
    println!("cargo:rustc-link-lib=static=ttkPointSetToCurveCS");
    println!("cargo:rustc-link-lib=static=ttkProgramBase");
    println!("cargo:rustc-link-lib=static=ttkProgramBaseCS");
    println!("cargo:rustc-link-lib=static=ttkProjectionFromField");
    println!("cargo:rustc-link-lib=static=ttkProjectionFromFieldCS");
    println!("cargo:rustc-link-lib=static=ttkQuadrangulationSubdivision");
    println!("cargo:rustc-link-lib=static=ttkQuadrangulationSubdivisionCS");
    println!("cargo:rustc-link-lib=static=ttkRangePolygon");
    println!("cargo:rustc-link-lib=static=ttkRangePolygonCS");
    println!("cargo:rustc-link-lib=static=ttkReebSpace");
    println!("cargo:rustc-link-lib=static=ttkReebSpaceCS");
    println!("cargo:rustc-link-lib=static=ttkScalarFieldCriticalPoints");
    println!("cargo:rustc-link-lib=static=ttkScalarFieldCriticalPointsCS");
    println!("cargo:rustc-link-lib=static=ttkScalarFieldNormalizer");
    println!("cargo:rustc-link-lib=static=ttkScalarFieldNormalizerCS");
    println!("cargo:rustc-link-lib=static=ttkScalarFieldSmoother");
    println!("cargo:rustc-link-lib=static=ttkScalarFieldSmootherCS");
    println!("cargo:rustc-link-lib=static=ttkSphereFromPoint");
    println!("cargo:rustc-link-lib=static=ttkSphereFromPointCS");
    println!("cargo:rustc-link-lib=static=ttkStringArrayConverter");
    println!("cargo:rustc-link-lib=static=ttkStringArrayConverterCS");
    println!("cargo:rustc-link-lib=static=ttkTableDataSelector");
    println!("cargo:rustc-link-lib=static=ttkTableDataSelectorCS");
    println!("cargo:rustc-link-lib=static=ttkTextureMapFromField");
    println!("cargo:rustc-link-lib=static=ttkTextureMapFromFieldCS");
    println!("cargo:rustc-link-lib=static=ttkTopologicalCompression");
    println!("cargo:rustc-link-lib=static=ttkTopologicalCompressionCS");
    println!("cargo:rustc-link-lib=static=ttkTopologicalCompressionReader");
    println!("cargo:rustc-link-lib=static=ttkTopologicalCompressionReaderCS");
    println!("cargo:rustc-link-lib=static=ttkTopologicalCompressionWriter");
    println!("cargo:rustc-link-lib=static=ttkTopologicalCompressionWriterCS");
    println!("cargo:rustc-link-lib=static=ttkTopologicalSimplification");
    println!("cargo:rustc-link-lib=static=ttkTopologicalSimplificationCS");
    println!("cargo:rustc-link-lib=static=ttkTrackingFromFields");
    println!("cargo:rustc-link-lib=static=ttkTrackingFromFieldsCS");
    println!("cargo:rustc-link-lib=static=ttkTrackingFromOverlap");
    println!("cargo:rustc-link-lib=static=ttkTrackingFromOverlapCS");
    println!("cargo:rustc-link-lib=static=ttkTrackingFromPersistenceDiagrams");
    println!("cargo:rustc-link-lib=static=ttkTrackingFromPersistenceDiagramsCS");
    println!("cargo:rustc-link-lib=static=ttkTriangulationAlgorithm");
    println!("cargo:rustc-link-lib=static=ttkTriangulationAlgorithmCS");
    println!("cargo:rustc-link-lib=static=ttkTriangulationRequest");
    println!("cargo:rustc-link-lib=static=ttkTriangulationRequestCS");
    println!("cargo:rustc-link-lib=static=ttkUncertainDataEstimator");
    println!("cargo:rustc-link-lib=static=ttkUncertainDataEstimatorCS");
    println!("cargo:rustc-link-lib=static=ttkUserInterfaceBase");
    println!("cargo:rustc-link-lib=static=ttkUserInterfaceBaseCS");
    println!("cargo:rustc-link-lib=static=ttkWRLExporter");
    println!("cargo:rustc-link-lib=static=ttkWRLExporterCS");*/
    println!("cargo:rustc-link-lib=static=ttkAlgorithm");
    println!("cargo:rustc-link-lib=static=ttkAlgorithmCS");
    println!("cargo:rustc-link-lib=static=ttkArrayEditor");
    println!("cargo:rustc-link-lib=static=ttkArrayEditorCS");
    println!("cargo:rustc-link-lib=static=ttkArrayPreconditioning");
    println!("cargo:rustc-link-lib=static=ttkArrayPreconditioningCS");
    println!("cargo:rustc-link-lib=static=ttkBarycentricSubdivision");
    println!("cargo:rustc-link-lib=static=ttkBarycentricSubdivisionCS");
    println!("cargo:rustc-link-lib=static=ttkBaseabstractMorseSmaleComplex");
    println!("cargo:rustc-link-lib=static=ttkBaseabstractTriangulation");
    println!("cargo:rustc-link-lib=static=ttkBaseassignmentSolver");
    println!("cargo:rustc-link-lib=static=ttkBasebarycentricSubdivision");
    println!("cargo:rustc-link-lib=static=ttkBasebottleneckDistance");
    println!("cargo:rustc-link-lib=static=ttkBaseboundingVolumeHierarchy");
    println!("cargo:rustc-link-lib=static=ttkBasecinemaImaging");
    println!("cargo:rustc-link-lib=static=ttkBasecinemaQuery");
    println!("cargo:rustc-link-lib=static=ttkBasecommon");
    println!("cargo:rustc-link-lib=static=ttkBasecontinuousScatterPlot");
    println!("cargo:rustc-link-lib=static=ttkBasecontourAroundPoint");
    println!("cargo:rustc-link-lib=static=ttkBasecontourForests");
    println!("cargo:rustc-link-lib=static=ttkBasecontourForestsTree");
    println!("cargo:rustc-link-lib=static=ttkBasecontourTree");
    println!("cargo:rustc-link-lib=static=ttkBasecontourTreeAlignment");
    println!("cargo:rustc-link-lib=static=ttkBasedimensionReduction");
    println!("cargo:rustc-link-lib=static=ttkBasediscreteGradient");
    println!("cargo:rustc-link-lib=static=ttkBasedistanceField");
    println!("cargo:rustc-link-lib=static=ttkBasedynamicTree");
    println!("cargo:rustc-link-lib=static=ttkBaseeigenField");
    println!("cargo:rustc-link-lib=static=ttkBaseexplicitTriangulation");
    println!("cargo:rustc-link-lib=static=ttkBasefiberSurface");
    println!("cargo:rustc-link-lib=static=ttkBaseftmTree");
    println!("cargo:rustc-link-lib=static=ttkBaseftmTreePP");
    println!("cargo:rustc-link-lib=static=ttkBaseftrGraph");
    println!("cargo:rustc-link-lib=static=ttkBasegaussianPointCloud");
    println!("cargo:rustc-link-lib=static=ttkBasegeometry");
    println!("cargo:rustc-link-lib=static=ttkBaseharmonicField");
    println!("cargo:rustc-link-lib=static=ttkBasehelloWorld");
    println!("cargo:rustc-link-lib=static=ttkBaseimplicitTriangulation");
    println!("cargo:rustc-link-lib=static=ttkBaseintegralLines");
    println!("cargo:rustc-link-lib=static=ttkBasejacobiSet");
    println!("cargo:rustc-link-lib=static=ttkBasekdTree");
    println!("cargo:rustc-link-lib=static=ttkBaselaplacian");
    println!("cargo:rustc-link-lib=static=ttkBaselDistance");
    println!("cargo:rustc-link-lib=static=ttkBaselDistanceMatrix");
    println!("cargo:rustc-link-lib=static=ttkBaselowestCommonAncestor");
    println!("cargo:rustc-link-lib=static=ttkBasemandatoryCriticalPoints");
    println!("cargo:rustc-link-lib=static=ttkBasemanifoldCheck");
    println!("cargo:rustc-link-lib=static=ttkBasemergeTreeClustering");
    println!("cargo:rustc-link-lib=static=ttkBasemergeTreeDistanceMatrix");
    println!("cargo:rustc-link-lib=static=ttkBasemergeTreeTemporalReductionDecoding");
    println!("cargo:rustc-link-lib=static=ttkBasemergeTreeTemporalReductionEncoding");
    println!("cargo:rustc-link-lib=static=ttkBasemorseSmaleComplex");
    println!("cargo:rustc-link-lib=static=ttkBasemorseSmaleComplex2D");
    println!("cargo:rustc-link-lib=static=ttkBasemorseSmaleComplex3D");
    println!("cargo:rustc-link-lib=static=ttkBasemorseSmaleQuadrangulation");
    println!("cargo:rustc-link-lib=static=ttkBasemultiresTriangulation");
    println!("cargo:rustc-link-lib=static=ttkBaseperiodicImplicitTriangulation");
    println!("cargo:rustc-link-lib=static=ttkBasepersistenceCurve");
    println!("cargo:rustc-link-lib=static=ttkBasepersistenceDiagram");
    println!("cargo:rustc-link-lib=static=ttkBasepersistenceDiagramAuction");
    println!("cargo:rustc-link-lib=static=ttkBasepersistenceDiagramClustering");
    println!("cargo:rustc-link-lib=static=ttkBasepersistenceDiagramDistanceMatrix");
    println!("cargo:rustc-link-lib=static=ttkBaseplanarGraphLayout");
    println!("cargo:rustc-link-lib=static=ttkBaseprogressiveTopology");
    println!("cargo:rustc-link-lib=static=ttkBasequadrangulationSubdivision");
    println!("cargo:rustc-link-lib=static=ttkBaserangeDrivenOctree");
    println!("cargo:rustc-link-lib=static=ttkBasereebSpace");
    println!("cargo:rustc-link-lib=static=ttkBasescalarFieldCriticalPoints");
    println!("cargo:rustc-link-lib=static=ttkBasescalarFieldSmoother");
    println!("cargo:rustc-link-lib=static=ttkBaseskeleton");
    println!("cargo:rustc-link-lib=static=ttkBasetopologicalCompression");
    println!("cargo:rustc-link-lib=static=ttkBasetopologicalSimplification");
    println!("cargo:rustc-link-lib=static=ttkBasetrackingFromFields");
    println!("cargo:rustc-link-lib=static=ttkBasetrackingFromPersistenceDiagrams");
    println!("cargo:rustc-link-lib=static=ttkBasetriangulation");
    println!("cargo:rustc-link-lib=static=ttkBaseuncertainDataEstimator");
    println!("cargo:rustc-link-lib=static=ttkBaseunionFind");
    println!("cargo:rustc-link-lib=static=ttkBasewebSocketIO");
    println!("cargo:rustc-link-lib=static=ttkBlockAggregator");
    println!("cargo:rustc-link-lib=static=ttkBlockAggregatorCS");
    println!("cargo:rustc-link-lib=static=ttkBottleneckDistance");
    println!("cargo:rustc-link-lib=static=ttkBottleneckDistanceCS");
    println!("cargo:rustc-link-lib=static=ttkCinemaDarkroom");
    println!("cargo:rustc-link-lib=static=ttkCinemaDarkroomCS");
    println!("cargo:rustc-link-lib=static=ttkCinemaImaging");
    println!("cargo:rustc-link-lib=static=ttkCinemaImagingCS");
    println!("cargo:rustc-link-lib=static=ttkCinemaProductReader");
    println!("cargo:rustc-link-lib=static=ttkCinemaProductReaderCS");
    println!("cargo:rustc-link-lib=static=ttkCinemaQuery");
    println!("cargo:rustc-link-lib=static=ttkCinemaQueryCS");
    println!("cargo:rustc-link-lib=static=ttkCinemaReader");
    println!("cargo:rustc-link-lib=static=ttkCinemaReaderCS");
    println!("cargo:rustc-link-lib=static=ttkCinemaWriter");
    println!("cargo:rustc-link-lib=static=ttkCinemaWriterCS");
    println!("cargo:rustc-link-lib=static=ttkComponentSize");
    println!("cargo:rustc-link-lib=static=ttkComponentSizeCS");
    println!("cargo:rustc-link-lib=static=ttkContinuousScatterPlot");
    println!("cargo:rustc-link-lib=static=ttkContinuousScatterPlotCS");
    println!("cargo:rustc-link-lib=static=ttkContourAroundPoint");
    println!("cargo:rustc-link-lib=static=ttkContourAroundPointCS");
    println!("cargo:rustc-link-lib=static=ttkContourForests");
    println!("cargo:rustc-link-lib=static=ttkContourForestsCS");
    println!("cargo:rustc-link-lib=static=ttkContourTreeAlignment");
    println!("cargo:rustc-link-lib=static=ttkContourTreeAlignmentCS");
    println!("cargo:rustc-link-lib=static=ttkDataSetInterpolator");
    println!("cargo:rustc-link-lib=static=ttkDataSetInterpolatorCS");
    println!("cargo:rustc-link-lib=static=ttkDataSetToTable");
    println!("cargo:rustc-link-lib=static=ttkDataSetToTableCS");
    println!("cargo:rustc-link-lib=static=ttkDepthImageBasedGeometryApproximation");
    println!("cargo:rustc-link-lib=static=ttkDepthImageBasedGeometryApproximationCS");
    println!("cargo:rustc-link-lib=static=ttkDimensionReduction");
    println!("cargo:rustc-link-lib=static=ttkDimensionReductionCS");
    println!("cargo:rustc-link-lib=static=ttkDiscreteGradient");
    println!("cargo:rustc-link-lib=static=ttkDiscreteGradientCS");
    println!("cargo:rustc-link-lib=static=ttkDistanceField");
    println!("cargo:rustc-link-lib=static=ttkDistanceFieldCS");
    println!("cargo:rustc-link-lib=static=ttkEigenField");
    println!("cargo:rustc-link-lib=static=ttkEigenFieldCS");
    println!("cargo:rustc-link-lib=static=ttkEndFor");
    println!("cargo:rustc-link-lib=static=ttkEndForCS");
    println!("cargo:rustc-link-lib=static=ttkExtract");
    println!("cargo:rustc-link-lib=static=ttkExtractCS");
    println!("cargo:rustc-link-lib=static=ttkFiber");
    println!("cargo:rustc-link-lib=static=ttkFiberCS");
    println!("cargo:rustc-link-lib=static=ttkFiberSurface");
    println!("cargo:rustc-link-lib=static=ttkFiberSurfaceCS");
    println!("cargo:rustc-link-lib=static=ttkFlattenMultiBlock");
    println!("cargo:rustc-link-lib=static=ttkFlattenMultiBlockCS");
    println!("cargo:rustc-link-lib=static=ttkForEach");
    println!("cargo:rustc-link-lib=static=ttkForEachCS");
    println!("cargo:rustc-link-lib=static=ttkFTMTree");
    println!("cargo:rustc-link-lib=static=ttkFTMTreeCS");
    println!("cargo:rustc-link-lib=static=ttkFTRGraph");
    println!("cargo:rustc-link-lib=static=ttkFTRGraphCS");
    println!("cargo:rustc-link-lib=static=ttkGaussianPointCloud");
    println!("cargo:rustc-link-lib=static=ttkGaussianPointCloudCS");
    println!("cargo:rustc-link-lib=static=ttkGeometrySmoother");
    println!("cargo:rustc-link-lib=static=ttkGeometrySmootherCS");
    println!("cargo:rustc-link-lib=static=ttkGridLayout");
    println!("cargo:rustc-link-lib=static=ttkGridLayoutCS");
    println!("cargo:rustc-link-lib=static=ttkHarmonicField");
    println!("cargo:rustc-link-lib=static=ttkHarmonicFieldCS");
    println!("cargo:rustc-link-lib=static=ttkHelloWorld");
    println!("cargo:rustc-link-lib=static=ttkHelloWorldCS");
    println!("cargo:rustc-link-lib=static=ttkIcosphere");
    println!("cargo:rustc-link-lib=static=ttkIcosphereCS");
    println!("cargo:rustc-link-lib=static=ttkIcosphereFromObject");
    println!("cargo:rustc-link-lib=static=ttkIcosphereFromObjectCS");
    println!("cargo:rustc-link-lib=static=ttkIcospheresFromPoints");
    println!("cargo:rustc-link-lib=static=ttkIcospheresFromPointsCS");
    println!("cargo:rustc-link-lib=static=ttkIdentifierRandomizer");
    println!("cargo:rustc-link-lib=static=ttkIdentifierRandomizerCS");
    println!("cargo:rustc-link-lib=static=ttkIdentifiers");
    println!("cargo:rustc-link-lib=static=ttkIdentifiersCS");
    println!("cargo:rustc-link-lib=static=ttkIdentifyByScalarField");
    println!("cargo:rustc-link-lib=static=ttkIdentifyByScalarFieldCS");
    println!("cargo:rustc-link-lib=static=ttkImportEmbeddingFromTable");
    println!("cargo:rustc-link-lib=static=ttkImportEmbeddingFromTableCS");
    println!("cargo:rustc-link-lib=static=ttkIntegralLines");
    println!("cargo:rustc-link-lib=static=ttkIntegralLinesCS");
    println!("cargo:rustc-link-lib=static=ttkJacobiSet");
    println!("cargo:rustc-link-lib=static=ttkJacobiSetCS");
    println!("cargo:rustc-link-lib=static=ttkLDistance");
    println!("cargo:rustc-link-lib=static=ttkLDistanceCS");
    println!("cargo:rustc-link-lib=static=ttkLDistanceMatrix");
    println!("cargo:rustc-link-lib=static=ttkLDistanceMatrixCS");
    println!("cargo:rustc-link-lib=static=ttkMandatoryCriticalPoints");
    println!("cargo:rustc-link-lib=static=ttkMandatoryCriticalPointsCS");
    println!("cargo:rustc-link-lib=static=ttkManifoldCheck");
    println!("cargo:rustc-link-lib=static=ttkManifoldCheckCS");
    println!("cargo:rustc-link-lib=static=ttkMatrixToHeatMap");
    println!("cargo:rustc-link-lib=static=ttkMatrixToHeatMapCS");
    println!("cargo:rustc-link-lib=static=ttkMergeBlockTables");
    println!("cargo:rustc-link-lib=static=ttkMergeBlockTablesCS");
    println!("cargo:rustc-link-lib=static=ttkMergeTreeClustering");
    println!("cargo:rustc-link-lib=static=ttkMergeTreeClusteringCS");
    println!("cargo:rustc-link-lib=static=ttkMergeTreeDistanceMatrix");
    println!("cargo:rustc-link-lib=static=ttkMergeTreeDistanceMatrixCS");
    println!("cargo:rustc-link-lib=static=ttkMergeTreeTemporalReductionDecoding");
    println!("cargo:rustc-link-lib=static=ttkMergeTreeTemporalReductionDecodingCS");
    println!("cargo:rustc-link-lib=static=ttkMergeTreeTemporalReductionEncoding");
    println!("cargo:rustc-link-lib=static=ttkMergeTreeTemporalReductionEncodingCS");
    println!("cargo:rustc-link-lib=static=ttkMeshGraph");
    println!("cargo:rustc-link-lib=static=ttkMeshGraphCS");
    println!("cargo:rustc-link-lib=static=ttkMeshSubdivision");
    println!("cargo:rustc-link-lib=static=ttkMeshSubdivisionCS");
    println!("cargo:rustc-link-lib=static=ttkMorphologicalOperators");
    println!("cargo:rustc-link-lib=static=ttkMorphologicalOperatorsCS");
    println!("cargo:rustc-link-lib=static=ttkMorseSmaleComplex");
    println!("cargo:rustc-link-lib=static=ttkMorseSmaleComplexCS");
    println!("cargo:rustc-link-lib=static=ttkMorseSmaleQuadrangulation");
    println!("cargo:rustc-link-lib=static=ttkMorseSmaleQuadrangulationCS");
    println!("cargo:rustc-link-lib=static=ttkOBJWriter");
    println!("cargo:rustc-link-lib=static=ttkOBJWriterCS");
    println!("cargo:rustc-link-lib=static=ttkOFFReader");
    println!("cargo:rustc-link-lib=static=ttkOFFReaderCS");
    println!("cargo:rustc-link-lib=static=ttkOFFWriter");
    println!("cargo:rustc-link-lib=static=ttkOFFWriterCS");
    println!("cargo:rustc-link-lib=static=ttkPeriodicGrid");
    println!("cargo:rustc-link-lib=static=ttkPeriodicGridCS");
    println!("cargo:rustc-link-lib=static=ttkPersistenceCurve");
    println!("cargo:rustc-link-lib=static=ttkPersistenceCurveCS");
    println!("cargo:rustc-link-lib=static=ttkPersistenceDiagram");
    println!("cargo:rustc-link-lib=static=ttkPersistenceDiagramClustering");
    println!("cargo:rustc-link-lib=static=ttkPersistenceDiagramClusteringCS");
    println!("cargo:rustc-link-lib=static=ttkPersistenceDiagramCS");
    println!("cargo:rustc-link-lib=static=ttkPersistenceDiagramDistanceMatrix");
    println!("cargo:rustc-link-lib=static=ttkPersistenceDiagramDistanceMatrixCS");
    println!("cargo:rustc-link-lib=static=ttkPlanarGraphLayout");
    println!("cargo:rustc-link-lib=static=ttkPlanarGraphLayoutCS");
    println!("cargo:rustc-link-lib=static=ttkPointDataConverter");
    println!("cargo:rustc-link-lib=static=ttkPointDataConverterCS");
    println!("cargo:rustc-link-lib=static=ttkPointDataSelector");
    println!("cargo:rustc-link-lib=static=ttkPointDataSelectorCS");
    println!("cargo:rustc-link-lib=static=ttkPointMerger");
    println!("cargo:rustc-link-lib=static=ttkPointMergerCS");
    println!("cargo:rustc-link-lib=static=ttkPointSetToCurve");
    println!("cargo:rustc-link-lib=static=ttkPointSetToCurveCS");
    println!("cargo:rustc-link-lib=static=ttkProgramBase");
    println!("cargo:rustc-link-lib=static=ttkProjectionFromField");
    println!("cargo:rustc-link-lib=static=ttkProjectionFromFieldCS");
    println!("cargo:rustc-link-lib=static=ttkQuadrangulationSubdivision");
    println!("cargo:rustc-link-lib=static=ttkQuadrangulationSubdivisionCS");
    println!("cargo:rustc-link-lib=static=ttkRangePolygon");
    println!("cargo:rustc-link-lib=static=ttkRangePolygonCS");
    println!("cargo:rustc-link-lib=static=ttkReebSpace");
    println!("cargo:rustc-link-lib=static=ttkReebSpaceCS");
    println!("cargo:rustc-link-lib=static=ttkScalarFieldCriticalPoints");
    println!("cargo:rustc-link-lib=static=ttkScalarFieldCriticalPointsCS");
    println!("cargo:rustc-link-lib=static=ttkScalarFieldNormalizer");
    println!("cargo:rustc-link-lib=static=ttkScalarFieldNormalizerCS");
    println!("cargo:rustc-link-lib=static=ttkScalarFieldSmoother");
    println!("cargo:rustc-link-lib=static=ttkScalarFieldSmootherCS");
    println!("cargo:rustc-link-lib=static=ttkSphereFromPoint");
    println!("cargo:rustc-link-lib=static=ttkSphereFromPointCS");
    println!("cargo:rustc-link-lib=static=ttkStableManifoldPersistence");
    println!("cargo:rustc-link-lib=static=ttkStableManifoldPersistenceCS");
    println!("cargo:rustc-link-lib=static=ttkStringArrayConverter");
    println!("cargo:rustc-link-lib=static=ttkStringArrayConverterCS");
    println!("cargo:rustc-link-lib=static=ttkTableDataSelector");
    println!("cargo:rustc-link-lib=static=ttkTableDataSelectorCS");
    println!("cargo:rustc-link-lib=static=ttkTextureMapFromField");
    println!("cargo:rustc-link-lib=static=ttkTextureMapFromFieldCS");
    println!("cargo:rustc-link-lib=static=ttkTopologicalCompression");
    println!("cargo:rustc-link-lib=static=ttkTopologicalCompressionCS");
    println!("cargo:rustc-link-lib=static=ttkTopologicalCompressionReader");
    println!("cargo:rustc-link-lib=static=ttkTopologicalCompressionReaderCS");
    println!("cargo:rustc-link-lib=static=ttkTopologicalCompressionWriter");
    println!("cargo:rustc-link-lib=static=ttkTopologicalCompressionWriterCS");
    println!("cargo:rustc-link-lib=static=ttkTopologicalSimplification");
    println!("cargo:rustc-link-lib=static=ttkTopologicalSimplificationByPersistence");
    println!("cargo:rustc-link-lib=static=ttkTopologicalSimplificationByPersistenceCS");
    println!("cargo:rustc-link-lib=static=ttkTopologicalSimplificationCS");
    println!("cargo:rustc-link-lib=static=ttkTrackingFromFields");
    println!("cargo:rustc-link-lib=static=ttkTrackingFromFieldsCS");
    println!("cargo:rustc-link-lib=static=ttkTrackingFromOverlap");
    println!("cargo:rustc-link-lib=static=ttkTrackingFromOverlapCS");
    println!("cargo:rustc-link-lib=static=ttkTrackingFromPersistenceDiagrams");
    println!("cargo:rustc-link-lib=static=ttkTrackingFromPersistenceDiagramsCS");
    println!("cargo:rustc-link-lib=static=ttkTriangulationReader");
    println!("cargo:rustc-link-lib=static=ttkTriangulationReaderCS");
    println!("cargo:rustc-link-lib=static=ttkTriangulationRequest");
    println!("cargo:rustc-link-lib=static=ttkTriangulationRequestCS");
    println!("cargo:rustc-link-lib=static=ttkTriangulationWriter");
    println!("cargo:rustc-link-lib=static=ttkTriangulationWriterCS");
    println!("cargo:rustc-link-lib=static=ttkUncertainDataEstimator");
    println!("cargo:rustc-link-lib=static=ttkUncertainDataEstimatorCS");
    println!("cargo:rustc-link-lib=static=ttkUserInterfaceBase");
    println!("cargo:rustc-link-lib=static=ttkWebSocketIO");
    println!("cargo:rustc-link-lib=static=ttkWebSocketIOCS");
    println!("cargo:rustc-link-lib=static=ttkWRLExporter");
    println!("cargo:rustc-link-lib=static=ttkWRLExporterCS");
    //println!("cargo:rustc-link-lib=static=libiomp5md");
    println!("cargo:rustc-link-lib=static=libomp");
    //println!("cargo:rustc-link-lib=static=PCL-pxi");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());

    //generate ttk_cpp.lib from ttk_cpp.cpp
    let _dst = cmake::Config::new("src/ttk_cpp")
        .generator("Visual Studio 16 2019")
        .build();

    //generate ttk_cpp.rs from ttk_cpp.hpp
    let bindings = bindgen::Builder::default()
        .clang_arg("-x")
        .clang_arg("c++")
        .clang_arg("-fopenmp")
        .header("src/ttk_cpp/ttk_cpp.hpp")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        .generate()
        .expect("Unable to generate bindings");

    bindings
        .write_to_file(out_path.join("ttk_cpp.rs"))
        .expect("Couldn't write bindings!");
    Ok(())
}
