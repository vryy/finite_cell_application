//  see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Feb 7, 2017$
//   Revision:            $Revision: 1.0 $
//
//


// System includes


// External includes


// Project includes
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "custom_geometries/finite_cell_geometry.h"
#include "finite_cell_application.h"
#include "finite_cell_application_variables.h"

#ifdef SD_APP_FORWARD_COMPATIBILITY
#define FINITE_CELL_APP_CREATE_ELEMENT(element_type, geometry_type, number_of_nodes) \
    element_type( 0, Element::GeometryType::Pointer( new geometry_type <Node<3> >( Element::GeometryType::PointsArrayType( number_of_nodes ) ) ) )
#define FINITE_CELL_APP_CREATE_CONDITION(condition_type, geometry_type, number_of_nodes) \
    condition_type( 0, Condition::GeometryType::Pointer( new geometry_type <Node<3> >( Condition::GeometryType::PointsArrayType( number_of_nodes ) ) ) )
#else
#define FINITE_CELL_APP_CREATE_ELEMENT(element_type, geometry_type, number_of_nodes) \
    element_type( 0, Element::GeometryType::Pointer( new geometry_type <Node<3> >( Element::GeometryType::PointsArrayType( number_of_nodes, Node<3>() ) ) ) )
#define FINITE_CELL_APP_CREATE_CONDITION(condition_type, geometry_type, number_of_nodes) \
    condition_type( 0, Condition::GeometryType::Pointer( new geometry_type <Node<3> >( Condition::GeometryType::PointsArrayType( number_of_nodes, Node<3>() ) ) ) )
#endif

namespace Kratos
{

    KratosFiniteCellApplication::KratosFiniteCellApplication()
    #ifdef SD_APP_FORWARD_COMPATIBILITY
    : KratosApplication("KratosFiniteCellApplication")
    #else
    : KratosApplication()
    #endif
    , FINITE_CELL_APP_CREATE_CONDITION( mDummySurfaceCondition2D3N, Triangle2D3, 3 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummySurfaceCondition2D6N, Triangle2D6, 6 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummySurfaceCondition2D4N, Quadrilateral2D4, 4 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummySurfaceCondition2D8N, Quadrilateral2D8, 8 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummySurfaceCondition2D9N, Quadrilateral2D9, 9 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummySurfaceCondition3D3N, Triangle3D3, 3 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummySurfaceCondition3D6N, Triangle3D6, 6 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummySurfaceCondition3D4N, Quadrilateral3D4, 4 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummySurfaceCondition3D8N, Quadrilateral3D8, 8 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummySurfaceCondition3D9N, Quadrilateral3D9, 9 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyConditionPoint2D, Point2D, 1 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyConditionPoint3D, Point3D, 1 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyConditionLine2N, Line3D2, 2 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyConditionLine3N, Line3D3, 3 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyCondition2D3N, Triangle2D3, 3 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyCondition2D4N, Quadrilateral2D4, 4 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyCondition2D6N, Triangle2D6, 6 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyCondition2D8N, Quadrilateral2D8, 8 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyCondition2D9N, Quadrilateral2D9, 9 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyCondition3D4N, Tetrahedra3D4, 4 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyCondition3D10N, Tetrahedra3D10, 10 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyCondition3D8N, Hexahedra3D8, 8 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyCondition3D20N, Hexahedra3D20, 20 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyCondition3D27N, Hexahedra3D27, 27 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyCondition3D6N, Prism3D6, 6 )
    , FINITE_CELL_APP_CREATE_CONDITION( mDummyCondition3D15N, Prism3D15, 15 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummySurfaceElement2D3N, Triangle2D3, 3 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummySurfaceElement2D6N, Triangle2D6, 6 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummySurfaceElement2D4N, Quadrilateral2D4, 4 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummySurfaceElement2D8N, Quadrilateral2D8, 8 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummySurfaceElement2D9N, Quadrilateral2D9, 9 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummySurfaceElement3D3N, Triangle3D3, 3 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummySurfaceElement3D6N, Triangle3D6, 6 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummySurfaceElement3D4N, Quadrilateral3D4, 4 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummySurfaceElement3D8N, Quadrilateral3D8, 8 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummySurfaceElement3D9N, Quadrilateral3D9, 9 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummyVolumeElement3D4N, Tetrahedra3D4, 4 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummyVolumeElement3D10N, Tetrahedra3D10, 10 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummyVolumeElement3D8N, Hexahedra3D8, 8 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummyVolumeElement3D20N, Hexahedra3D20, 20 )
    , FINITE_CELL_APP_CREATE_ELEMENT( mDummyVolumeElement3D27N, Hexahedra3D27, 27 )
    {}

    void KratosFiniteCellApplication::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosFiniteCellApplication... " << std::endl;

        // register variables to Kratos kernel
        KRATOS_REGISTER_VARIABLE( SUBCELL_WEIGHTS )
        KRATOS_REGISTER_VARIABLE( SUBCELL_DOMAIN_SIZE )
        KRATOS_REGISTER_VARIABLE( SUBCELL_DOMAIN_SIZES )
        KRATOS_REGISTER_VARIABLE( PHYSICAL_INTEGRATION_POINT_THREED_STRESSES )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( PHYSICAL_INTEGRATION_POINT_LOCAL )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( PHYSICAL_INTEGRATION_POINT_GLOBAL )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( PHYSICAL_INTEGRATION_POINT_DISPLACEMENT )
        KRATOS_REGISTER_VARIABLE( OTHER_NODE_ID )
        KRATOS_REGISTER_VARIABLE( OTHER_ID )
        KRATOS_REGISTER_VARIABLE( NUMBER_OF_PHYSICAL_POINTS )

        // register conditions to Kratos kernel
        KRATOS_REGISTER_CONDITION("DummyConditionPoint2D", mDummyConditionPoint2D)
        KRATOS_REGISTER_CONDITION("DummyConditionPoint3D", mDummyConditionPoint3D)
        KRATOS_REGISTER_CONDITION( "DummyConditionLine2N", mDummyConditionLine2N )
        KRATOS_REGISTER_CONDITION( "DummyConditionLine3N", mDummyConditionLine3N )
        KRATOS_REGISTER_CONDITION( "DummyCondition2D3N", mDummyCondition2D3N )
        KRATOS_REGISTER_CONDITION( "DummyCondition2D4N", mDummyCondition2D4N )
        KRATOS_REGISTER_CONDITION( "DummyCondition2D6N", mDummyCondition2D6N )
        KRATOS_REGISTER_CONDITION( "DummyCondition2D8N", mDummyCondition2D8N )
        KRATOS_REGISTER_CONDITION( "DummyCondition2D9N", mDummyCondition2D9N )
        KRATOS_REGISTER_CONDITION( "DummyCondition3D4N", mDummyCondition3D4N )
        KRATOS_REGISTER_CONDITION( "DummyCondition3D10N", mDummyCondition3D10N )
        KRATOS_REGISTER_CONDITION( "DummyCondition3D8N", mDummyCondition3D8N )
        KRATOS_REGISTER_CONDITION( "DummyCondition3D20N", mDummyCondition3D20N )
        KRATOS_REGISTER_CONDITION( "DummyCondition3D27N", mDummyCondition3D27N )
        KRATOS_REGISTER_CONDITION( "DummyCondition3D6N", mDummyCondition3D6N )
        KRATOS_REGISTER_CONDITION( "DummyCondition3D15N", mDummyCondition3D15N )

        KRATOS_REGISTER_CONDITION( "DummySurfaceCondition2D3N", mDummySurfaceCondition2D3N )
        KRATOS_REGISTER_CONDITION( "DummySurfaceCondition2D6N", mDummySurfaceCondition2D6N )
        KRATOS_REGISTER_CONDITION( "DummySurfaceCondition2D4N", mDummySurfaceCondition2D4N )
        KRATOS_REGISTER_CONDITION( "DummySurfaceCondition2D8N", mDummySurfaceCondition2D8N )
        KRATOS_REGISTER_CONDITION( "DummySurfaceCondition2D9N", mDummySurfaceCondition2D9N )
        KRATOS_REGISTER_CONDITION( "DummySurfaceCondition3D3N", mDummySurfaceCondition3D3N )
        KRATOS_REGISTER_CONDITION( "DummySurfaceCondition3D6N", mDummySurfaceCondition3D6N )
        KRATOS_REGISTER_CONDITION( "DummySurfaceCondition3D4N", mDummySurfaceCondition3D4N )
        KRATOS_REGISTER_CONDITION( "DummySurfaceCondition3D8N", mDummySurfaceCondition3D8N )
        KRATOS_REGISTER_CONDITION( "DummySurfaceCondition3D9N", mDummySurfaceCondition3D9N )

        // register elements to Kratos kernel
        KRATOS_REGISTER_ELEMENT( "DummySurfaceElement2D3N", mDummySurfaceElement2D3N )
        KRATOS_REGISTER_ELEMENT( "DummySurfaceElement2D6N", mDummySurfaceElement2D6N )
        KRATOS_REGISTER_ELEMENT( "DummySurfaceElement2D4N", mDummySurfaceElement2D4N )
        KRATOS_REGISTER_ELEMENT( "DummySurfaceElement2D8N", mDummySurfaceElement2D8N )
        KRATOS_REGISTER_ELEMENT( "DummySurfaceElement2D9N", mDummySurfaceElement2D9N )
        KRATOS_REGISTER_ELEMENT( "DummySurfaceElement3D3N", mDummySurfaceElement3D3N )
        KRATOS_REGISTER_ELEMENT( "DummySurfaceElement3D6N", mDummySurfaceElement3D6N )
        KRATOS_REGISTER_ELEMENT( "DummySurfaceElement3D4N", mDummySurfaceElement3D4N )
        KRATOS_REGISTER_ELEMENT( "DummySurfaceElement3D8N", mDummySurfaceElement3D8N )
        KRATOS_REGISTER_ELEMENT( "DummySurfaceElement3D9N", mDummySurfaceElement3D9N )

        KRATOS_REGISTER_ELEMENT( "DummyVolumeElement3D4N", mDummyVolumeElement3D4N )
        KRATOS_REGISTER_ELEMENT( "DummyVolumeElement3D10N", mDummyVolumeElement3D10N )
        KRATOS_REGISTER_ELEMENT( "DummyVolumeElement3D8N", mDummyVolumeElement3D8N )
        KRATOS_REGISTER_ELEMENT( "DummyVolumeElement3D20N", mDummyVolumeElement3D20N )
        KRATOS_REGISTER_ELEMENT( "DummyVolumeElement3D27N", mDummyVolumeElement3D27N )
    }

} // namespace Kratos

