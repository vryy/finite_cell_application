// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//



// Project includes
#include "includes/element.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_algebra/brep.h"
#include "custom_utilities/quadrature_utility.h"
#include "custom_utilities/quad_tree.h"
#include "custom_utilities/div_free_basis_utility.h"
#include "custom_utilities/moment_fitting_utility.h"
#include "custom_utilities/moment_fitted_quad_tree_subcell.h"
#include "custom_utilities/finite_cell_auxilliary_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

Matrix ComputeDivFreeBasis(DivFreeBasisUtility& dummy, const std::size_t& Dim, const std::size_t& Degree, const DivFreeBasisUtility::CoordinatesArrayType& rPoint)
{
    if(Dim == 2)
    {
        if(Degree == 0)
            return dummy.GetValues<2, 0>(rPoint);
        else if(Degree == 1)
            return dummy.GetValues<2, 1>(rPoint);
        else if(Degree == 2)
            return dummy.GetValues<2, 2>(rPoint);
        else if(Degree == 3)
            return dummy.GetValues<2, 3>(rPoint);
        else if(Degree == 4)
            return dummy.GetValues<2, 4>(rPoint);
    }
}

std::size_t FiniteCellAuxilliaryUtility_GetLastNodeId(FiniteCellAuxilliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastNodeId(r_model_part);
}

std::size_t FiniteCellAuxilliaryUtility_GetLastElementId(FiniteCellAuxilliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastElementId(r_model_part);
}

std::size_t FiniteCellAuxilliaryUtility_GetLastConditionId(FiniteCellAuxilliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastConditionId(r_model_part);
}

std::size_t FiniteCellAuxilliaryUtility_GetLastPropertiesId(FiniteCellAuxilliaryUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastPropertiesId(r_model_part);
}

void FiniteCellAuxilliaryUtility_AddElement(FiniteCellAuxilliaryUtility& rDummy, ModelPart::ElementsContainerType& rpElements,
        Element::Pointer pElement)
{
    rDummy.AddElement(rpElements, pElement);
}

int QuadratureUtility_GetQuadratureType(QuadratureUtility& rDummy, const int& integration_method)
{
    return rDummy.GetQuadratureType(integration_method);
}

int QuadratureUtility_GetQuadratureOrder(QuadratureUtility& rDummy, const int& integration_method)
{
    return rDummy.GetQuadratureOrder(integration_method);
}

template<class TTreeType, class TBRepType>
void FiniteCellAuxilliaryUtility_MultithreadedRefineBy(FiniteCellAuxilliaryUtility& rDummy, boost::python::list& r_trees,
        const TBRepType& r_brep)
{
    typedef typename TTreeType::Pointer TTreePointerType;
    std::vector<TTreePointerType> trees;
    typedef boost::python::stl_input_iterator<TTreePointerType> iterator_tree_type;
    BOOST_FOREACH(const typename iterator_tree_type::value_type& t,
                  std::make_pair(iterator_tree_type(r_trees), // begin
                  iterator_tree_type() ) ) // end
    {
        trees.push_back(t);
    }

    rDummy.MultithreadedRefineBy<TTreeType, TBRepType>(trees, r_brep);
}

void FiniteCellApplication_AddCustomUtilitiesToPython()
{
    void(QuadratureUtility::*pointer_to_CreateConditionFromQuadraturePoint)(ModelPart&, boost::python::list&,
            const std::string&, const double&, const double&) const = &QuadratureUtility::PyCreateConditionFromQuadraturePoint;
    ModelPart::ConditionsContainerType(QuadratureUtility::*pointer_to_CreateConditionFromPoint)(ModelPart&,
            boost::python::list&, const std::string&) const = &QuadratureUtility::PyCreateConditionFromPoint;

    class_<QuadratureUtility, QuadratureUtility::Pointer, boost::noncopyable>
    ("QuadratureUtility", init<>())
    .def("GetDefaultIntegrationMethod", &QuadratureUtility::GetDefaultIntegrationMethod<Element>)
    .def("GetDefaultIntegrationMethod", &QuadratureUtility::GetDefaultIntegrationMethod<Condition>)
    .def("GetQuadratureType", QuadratureUtility_GetQuadratureType)
    .def("GetQuadratureOrder", QuadratureUtility_GetQuadratureOrder)
    .def("ScaleQuadrature", &QuadratureUtility::PyScaleQuadrature)
    .def("SaveQuadrature", &QuadratureUtility::PySaveQuadrature)
    .def("SaveQuadrature", &QuadratureUtility::PySaveQuadratureAdvanced)
    .def("SaveQuadratureSubCell", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<1> >)
    .def("SaveQuadratureSubCell2", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<2> >)
    .def("SaveQuadratureSubCell3", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<3> >)
    .def("SaveQuadratureSubCell4", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<4> >)
    .def("SaveQuadratureSubCell5", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<5> >)
    .def("SetQuadrature", &QuadratureUtility::PySetQuadrature)
    .def("CreateConditionFromQuadraturePoint", pointer_to_CreateConditionFromQuadraturePoint)
    .def("CreateConditionFromPoint", pointer_to_CreateConditionFromPoint)
    .def("CreatePoint", &QuadratureUtility::CreatePoint)
    ;

    void(DivFreeBasisUtility::*pointer_to_AssignQuadrature2D)(Element::Pointer&, const LevelSet&, const unsigned int&, const unsigned int&) = &DivFreeBasisUtility::AssignQuadrature2D;

    class_<DivFreeBasisUtility, DivFreeBasisUtility::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("DivFreeBasisUtility", init<>() )
    .def("GetValues", &ComputeDivFreeBasis)
    .def("AssignQuadrature2D", pointer_to_AssignQuadrature2D)
    ;

    class_<MomentFittingUtility, MomentFittingUtility::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("MomentFittingUtility", init<>())
    .def("FitQuadrature", &MomentFittingUtility::PyFitQuadrature<QuadTree<1> >)
    .def("MultithreadedFitQuadrature", &MomentFittingUtility::PyMultithreadedFitQuadrature<QuadTree<1> >)
    .def("MultithreadedFitQuadratureSubCell", &MomentFittingUtility::PyMultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<1> >)
    .def("FitQuadrature2", &MomentFittingUtility::PyFitQuadrature<QuadTree<2> >)
    .def("MultithreadedFitQuadrature2", &MomentFittingUtility::PyMultithreadedFitQuadrature<QuadTree<2> >)
    .def("MultithreadedFitQuadratureSubCell2", &MomentFittingUtility::PyMultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<2> >)
    .def("FitQuadrature3", &MomentFittingUtility::PyFitQuadrature<QuadTree<3> >)
    .def("MultithreadedFitQuadrature3", &MomentFittingUtility::PyMultithreadedFitQuadrature<QuadTree<3> >)
    .def("MultithreadedFitQuadratureSubCell3", &MomentFittingUtility::PyMultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<3> >)
    .def("FitQuadrature4", &MomentFittingUtility::PyFitQuadrature<QuadTree<4> >)
    .def("MultithreadedFitQuadrature4", &MomentFittingUtility::PyMultithreadedFitQuadrature<QuadTree<4> >)
    .def("MultithreadedFitQuadratureSubCell4", &MomentFittingUtility::PyMultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<4> >)
    .def("FitQuadrature5", &MomentFittingUtility::PyFitQuadrature<QuadTree<5> >)
    .def("MultithreadedFitQuadrature5", &MomentFittingUtility::PyMultithreadedFitQuadrature<QuadTree<5> >)
    .def("MultithreadedFitQuadratureSubCell5", &MomentFittingUtility::PyMultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell<5> >)
    ;

    Condition::Pointer(FiniteCellAuxilliaryUtility::*pointer_to_PyCreateCondition)(ModelPart&, const std::string&,
            const std::size_t&, Properties::Pointer, boost::python::list&) const = &FiniteCellAuxilliaryUtility::PyCreateCondition;

    ModelPart::ElementsContainerType(FiniteCellAuxilliaryUtility::*pointer_to_PyGetElements)(ModelPart&,
            boost::python::list&) const = &FiniteCellAuxilliaryUtility::PyGetElements;

    void(FiniteCellAuxilliaryUtility::*pointer_to_PyGetElements2)(ModelPart::ElementsContainerType&, ModelPart&,
            boost::python::list&) const = &FiniteCellAuxilliaryUtility::PyGetElements;

    void(FiniteCellAuxilliaryUtility::*pointer_to_Clean)(ModelPart&,
            ModelPart::ConditionsContainerType&, const int&) const = &FiniteCellAuxilliaryUtility::Clean;

    void(FiniteCellAuxilliaryUtility::*pointer_to_PrintGeometry)(Element::GeometryType::Pointer) const = &FiniteCellAuxilliaryUtility::Print;

    class_<FiniteCellAuxilliaryUtility, FiniteCellAuxilliaryUtility::Pointer, boost::noncopyable>
    ("FiniteCellAuxilliaryUtility", init<>())
    .def("CreateCondition", pointer_to_PyCreateCondition)
    .def("GetElements", pointer_to_PyGetElements)
    .def("GetElements", pointer_to_PyGetElements2)
    .def("Clean", pointer_to_Clean)
    .def("GetLastNodeId", &FiniteCellAuxilliaryUtility_GetLastNodeId)
    .def("GetLastElementId", &FiniteCellAuxilliaryUtility_GetLastElementId)
    .def("GetLastConditionId", &FiniteCellAuxilliaryUtility_GetLastConditionId)
    .def("GetLastPropertiesId", &FiniteCellAuxilliaryUtility_GetLastPropertiesId)
    .def("AddElement", &FiniteCellAuxilliaryUtility_AddElement)
    .def("MultithreadedRefineBy", &FiniteCellAuxilliaryUtility_MultithreadedRefineBy<RefinableTree, BRep>)
    .def("Print", pointer_to_PrintGeometry)
    ;
}

}  // namespace Python.

}  // namespace Kratos.

