// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "custom_python/add_div_free_basis_utility_to_python.h"
#include "custom_utilities/quadrature_utility.h"
#include "custom_utilities/div_free_basis_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

Matrix ComputeDivFreeBasis(DivFreeBasisUtility& dummy, const std::size_t& Dim, const std::size_t& Degree, const DivFreeBasisUtility::CoordinatesArrayType& rPoint)
{
    if (Dim == 2)
    {
        if (Degree == 0)
        {
            return dummy.GetValues<2, 0>(rPoint);
        }
        else if (Degree == 1)
        {
            return dummy.GetValues<2, 1>(rPoint);
        }
        else if (Degree == 2)
        {
            return dummy.GetValues<2, 2>(rPoint);
        }
        else if (Degree == 3)
        {
            return dummy.GetValues<2, 3>(rPoint);
        }
        else if (Degree == 4)
        {
            return dummy.GetValues<2, 4>(rPoint);
        }
    }
    else
        KRATOS_ERROR << "Invalid dimension " << Dim;
}

void DivFreeBasisUtility_AssignQuadrature2D(DivFreeBasisUtility& rDummy, Element::Pointer p_elem,
        LevelSet::Pointer p_level_set, const unsigned int& integration_order, const unsigned int& div_free_order)
{
    rDummy.AssignQuadrature2D(p_elem->GetGeometry(), *p_level_set, integration_order, div_free_order);
}

void FiniteCellApplication_AddDivFreeBasisUtilityToPython()
{

    class_<DivFreeBasisUtility, DivFreeBasisUtility::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("DivFreeBasisUtility", init<>() )
    .def("GetValues", &ComputeDivFreeBasis)
    .def("AssignQuadrature2D", &DivFreeBasisUtility_AssignQuadrature2D)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

