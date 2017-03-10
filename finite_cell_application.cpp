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
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "finite_cell_application.h"


namespace Kratos
{

    KRATOS_CREATE_VARIABLE( boost::python::object, LOAD_FUNCTION )

    KratosFiniteCellApplication::KratosFiniteCellApplication()
    : mDummyPointCondition2D( 0, Element::GeometryType::Pointer( new Point2D <Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) )
    , mDummyPointCondition3D( 0, Element::GeometryType::Pointer( new Point3D <Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) )
    {}

    void KratosFiniteCellApplication::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosFiniteCellApplication... " << std::endl;

        // register variables to Kratos kernel
        KRATOS_REGISTER_VARIABLE( LOAD_FUNCTION )

        // register conditions to Kratos kernel
        KRATOS_REGISTER_CONDITION("DummyPointCondition2D", mDummyPointCondition2D)
        KRATOS_REGISTER_CONDITION("DummyPointCondition3D", mDummyPointCondition3D)
    }

} // namespace Kratos

