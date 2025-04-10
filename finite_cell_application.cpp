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
    element_type( 0, Element::GeometryType::Pointer( new geometry_type <Node>( Element::GeometryType::PointsArrayType( number_of_nodes ) ) ) )
#define FINITE_CELL_APP_CREATE_CONDITION(condition_type, geometry_type, number_of_nodes) \
    condition_type( 0, Condition::GeometryType::Pointer( new geometry_type <Node>( Condition::GeometryType::PointsArrayType( number_of_nodes ) ) ) )
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
{}

void KratosFiniteCellApplication::Register()
{
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
}

} // namespace Kratos
