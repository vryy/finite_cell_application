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
#include "finite_cell_application.h"


namespace Kratos
{

    KRATOS_CREATE_VARIABLE( boost::python::object, LOAD_FUNCTION )

    KratosFiniteCellApplication::KratosFiniteCellApplication()
    {}

    void KratosFiniteCellApplication::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosFiniteCellApplication... " << std::endl;

        // register variables to Kratos kernel
        KRATOS_REGISTER_VARIABLE( LOAD_FUNCTION )

    }

} // namespace Kratos

