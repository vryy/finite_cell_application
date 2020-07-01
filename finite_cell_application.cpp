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
#include "finite_cell_application/finite_cell_application.h"
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
#include "geometries/line_2d.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "custom_geometries/finite_cell_geometry.h"


namespace Kratos
{

    KRATOS_CREATE_VARIABLE( Matrix, SUBCELL_WEIGHTS )
    KRATOS_CREATE_VARIABLE( double, SUBCELL_DOMAIN_SIZE )
    KRATOS_CREATE_VARIABLE( Vector, SUBCELL_DOMAIN_SIZES )
    KRATOS_CREATE_VARIABLE( Vector, PHYSICAL_INTEGRATION_POINT_THREED_STRESSES )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( PHYSICAL_INTEGRATION_POINT_LOCAL )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( PHYSICAL_INTEGRATION_POINT_GLOBAL )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( PHYSICAL_INTEGRATION_POINT_DISPLACEMENT )
    KRATOS_CREATE_VARIABLE( int, OTHER_NODE_ID )
    KRATOS_CREATE_VARIABLE( int, OTHER_ID )
    KRATOS_CREATE_VARIABLE( int, NUMBER_OF_PHYSICAL_POINTS )

    KratosFiniteCellApplication::KratosFiniteCellApplication()
    : KratosApplication()
    , mDummySurfaceCondition2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mDummySurfaceCondition2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mDummySurfaceCondition2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mDummySurfaceCondition2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mDummySurfaceCondition2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mDummySurfaceCondition3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mDummySurfaceCondition3D6N( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mDummySurfaceCondition3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mDummySurfaceCondition3D8N( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mDummySurfaceCondition3D9N( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mDummyConditionPoint2D( 0, Element::GeometryType::Pointer( new Point2D <Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) )
    , mDummyConditionPoint3D( 0, Element::GeometryType::Pointer( new Point3D <Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) )
    , mDummyConditionLine2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
    , mDummyConditionLine3N( 0, Element::GeometryType::Pointer( new Line3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mDummyCondition2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mDummyCondition2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mDummyCondition2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mDummyCondition2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mDummyCondition2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mDummyCondition3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mDummyCondition3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) )
    , mDummyCondition3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mDummyCondition3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) )
    , mDummyCondition3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )
    , mDummyCondition3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mDummyCondition3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15, Node<3>() ) ) ) )
    , mDummySurfaceElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mDummySurfaceElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mDummySurfaceElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mDummySurfaceElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mDummySurfaceElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mDummySurfaceElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) )
    , mDummySurfaceElement3D6N( 0, Element::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) )
    , mDummySurfaceElement3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mDummySurfaceElement3D8N( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mDummySurfaceElement3D9N( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )
    , mDummyVolumeElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    , mDummyVolumeElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) )
    , mDummyVolumeElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) )
    , mDummyVolumeElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) )
    , mDummyVolumeElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )
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

