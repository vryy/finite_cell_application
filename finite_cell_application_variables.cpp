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
#include "finite_cell_application_variables.h"

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

} // namespace Kratos

