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
//#include "custom_utilities/binary_tree.h"
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

template<class TBRepType>
void FiniteCellAuxilliaryUtility_MultithreadedQuadTreeRefineBy(FiniteCellAuxilliaryUtility& rDummy, boost::python::list& r_trees,
        const TBRepType& r_brep)
{
    FiniteCellAuxilliaryUtility_MultithreadedRefineBy<QuadTree, TBRepType>(rDummy, r_trees, r_brep);
}

template<class TBRepType>
void FiniteCellAuxilliaryUtility_MultithreadedQuadTreeSubCellRefineBy(FiniteCellAuxilliaryUtility& rDummy, boost::python::list& r_trees,
        const TBRepType& r_brep)
{
    FiniteCellAuxilliaryUtility_MultithreadedRefineBy<QuadTreeSubCell, TBRepType>(rDummy, r_trees, r_brep);
}

boost::python::list MomentFittedQuadTreeSubCell_CreateSubCellElements(MomentFittedQuadTreeSubCell& rDummy,
        ModelPart& r_model_part,
        const std::string& subcell_element_type,
        const int& cut_cell_quadrature_order,
        boost::python::list& cut_cell_full_quadrature,
        boost::python::list& subcell_weights,
        std::size_t lastElementId,
        std::size_t lastCondId)
{
    typedef Element::GeometryType GeometryType;

    ModelPart::ElementsContainerType NewElements;

//    if(!boost::python::list::is_empty(subcell_weights))
//    {
    // TODO find a way to check if a list is empty
        GeometryType::IntegrationPointsArrayType integration_points;
        Matrix Weights;

        std::size_t num_physical_point = boost::python::len(subcell_weights);
        std::size_t weight_length = boost::python::len(subcell_weights[0]);
        Weights.resize(num_physical_point, weight_length, false);

        for(std::size_t i = 0; i < num_physical_point; ++i)
        {
            boost::python::list weights = boost::python::extract<boost::python::list>(subcell_weights[i]);
            for(std::size_t j = 0; j < weight_length; ++j)
            {
                Weights(i, j) = boost::python::extract<double>(weights[j]);
            }
        }
//        KRATOS_WATCH(Weights)

        for(std::size_t i = 0; i < boost::python::len(cut_cell_full_quadrature); ++i)
        {
            boost::python::list point = boost::python::extract<boost::python::list>(cut_cell_full_quadrature[i]);
            GeometryType::IntegrationPointType integration_point;
            integration_point.X() = boost::python::extract<double>(point[0]);
            integration_point.Y() = boost::python::extract<double>(point[1]);
            integration_point.Z() = boost::python::extract<double>(point[2]);
            integration_point.Weight() = boost::python::extract<double>(point[3]);
//            KRATOS_WATCH(integration_point)
            integration_points.push_back(integration_point);
        }

        NewElements = rDummy.CreateSubCellElements(r_model_part,
            subcell_element_type,
            cut_cell_quadrature_order,
            integration_points,
            Weights,
            lastElementId,
            lastCondId);
//        std::cout << "----------------------" << std::endl;
//    }

    boost::python::list Output;
    Output.append(lastElementId);
    Output.append(lastCondId);
    Output.append(NewElements);

    return Output;
}

void FiniteCellApplication_AddCustomUtilitiesToPython()
{
    typedef Element::GeometryType::PointType NodeType;

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
    .def("SaveQuadratureSubCell", &QuadratureUtility::PySaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell>)
    .def("SetQuadrature", &QuadratureUtility::PySetQuadrature)
    .def("CreateConditionFromQuadraturePoint", pointer_to_CreateConditionFromQuadraturePoint)
    .def("CreateConditionFromPoint", pointer_to_CreateConditionFromPoint)
    .def("CreatePoint", &QuadratureUtility::CreatePoint)
    ;

//    class_<BinaryTree<1>, BinaryTree<1>::Pointer, boost::noncopyable, bases<QuadratureUtility> >
//    ("BinaryTree", init<Element::Pointer&>())
//    .def(init<NodeType::Pointer&, NodeType::Pointer&>())
//    .def("Refine", &BinaryTree<1>::Refine)
//    .def(self_ns::str(self))
//    ;

//    typedef BinaryTree<2> QuadTreeType;
//    double(QuadTreeType::*pointer_to_Integrate_double_quadtree)(const FunctionR3R1&, const int) const = &QuadTreeType::Integrate<double>;
//    Vector(QuadTreeType::*pointer_to_Integrate_Vector_quadtree)(const FunctionR3Rn&, const int) const = &QuadTreeType::Integrate<Vector>;
//    void(QuadTreeType::*pointer_to_ConstructQuadrature_quadtree)(const LevelSet&, const int) const = &QuadTreeType::ConstructQuadrature;

//    class_<QuadTreeType, QuadTreeType::Pointer, boost::noncopyable, bases<QuadratureUtility> >
//    ("QuadTree", init<Element::Pointer&>())
//    .def(init<NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&>())
//    .def("Refine", &QuadTreeType::Refine)
//    .def("RefineBy", &QuadTreeType::RefineBy)
//    .def("ResetId", &QuadTreeType::ResetId)
//    .def("Renumber", &QuadTreeType::PyRenumber)
//    .def("AddToModelPart", &QuadTreeType::PyAddToModelPart)
//    .def("Integrate", pointer_to_Integrate_double_quadtree)
//    .def("Integrate", pointer_to_Integrate_Vector_quadtree)
//    .def("ConstructQuadrature", pointer_to_ConstructQuadrature_quadtree)
//    .def(self_ns::str(self))
//    ;

//    typedef BinaryTree<3> OctTreeType;
//    double(OctTreeType::*pointer_to_Integrate_double_octtree)(const FunctionR3R1&, const int) const = &OctTreeType::Integrate<double>;

//    class_<OctTreeType, OctTreeType::Pointer, boost::noncopyable, bases<QuadratureUtility> >
//    ("OctTree", init<Element::Pointer&>())
//    .def(init<NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&, NodeType::Pointer&>())
//    .def("Refine", &OctTreeType::Refine)
//    .def("RefineBy", &OctTreeType::RefineBy)
//    .def("ResetId", &OctTreeType::ResetId)
//    .def("Renumber", &OctTreeType::PyRenumber)
//    .def("AddToModelPart", &OctTreeType::PyAddToModelPart)
//    .def("Integrate", pointer_to_Integrate_double_octtree)
//    .def(self_ns::str(self))
//    ;

    double(QuadTree::*pointer_to_Integrate_double_quadtree_local)(const FunctionR3R1&, const int&) const = &QuadTree::Integrate<double>;

    class_<QuadTree, QuadTree::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("QuadTreeLocal", init<Element::Pointer&>())
    .def(init<Condition::Pointer&>())
    .def("pGetGeometry", &QuadTree::pGetGeometry)
    .def("pCreateGeometry", &QuadTree::pCreateGeometry)
    .def("DomainSize", &QuadTree::DomainSize)
    .def("CenterOfGravity", &QuadTree::CenterOfGravity)
    .def("Refine", &QuadTree::Refine)
    .def("RefineBy", &QuadTree::RefineBy)
    .def("AddToModelPart", &QuadTree::PyAddToModelPart)
    .def("Integrate", pointer_to_Integrate_double_quadtree_local)
    .def("ConstructQuadrature", &QuadTree::ConstructQuadrature)
    .def(self_ns::str(self))
    ;

    class_<QuadTreeSubCell, QuadTreeSubCell::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("QuadTreeSubCell", init<Element::Pointer&>())
    .def(init<Condition::Pointer&>())
    .def("NumberOfSubCells", &QuadTreeSubCell::NumberOfSubCells)
    .def("DomainSize", &QuadTreeSubCell::DomainSize)
    .def("Refine", &QuadTreeSubCell::Refine)
    .def("RefineBy", &QuadTreeSubCell::RefineBy)
    .def("CreateQuadTree", &QuadTreeSubCell::CreateQuadTree)
    .def("ShallowAddToModelPart", &MomentFittedQuadTreeSubCell::PyAddToModelPart<true>) // only add the sub-cell
    .def("DeepAddToModelPart", &MomentFittedQuadTreeSubCell::PyAddToModelPart<false>) // add the sub-cell and all the quad-trees
    ;

    class_<MomentFittedQuadTreeSubCell, MomentFittedQuadTreeSubCell::Pointer, bases<QuadTreeSubCell>, boost::noncopyable>
    ("MomentFittedQuadTreeSubCell", init<Element::Pointer>())
    .def(init<Element::Pointer, const std::string&, const int&>())
    .def(init<Element::Pointer, const std::string&, const std::size_t&, const std::size_t&>())
    .def(init<Element::Pointer, const std::string&, const std::size_t&, const std::size_t&, const std::size_t&>())
    .def("GetElement", &MomentFittedQuadTreeSubCell::pGetElement)
    .def("FitQuadraturePhysicalPoints", &MomentFittedQuadTreeSubCell::PyFitQuadraturePhysicalPoints)
    .def("FitAndCreateSubCellElements", &MomentFittedQuadTreeSubCell::PyFitAndCreateSubCellElements)
    .def("CreateSubCellElements", MomentFittedQuadTreeSubCell_CreateSubCellElements)
    .def("CreateParasiteElement", &MomentFittedQuadTreeSubCell::CreateParasiteElement)
    ;

    void(DivFreeBasisUtility::*pointer_to_AssignQuadrature2D)(Element::Pointer&, const LevelSet&, const unsigned int&, const unsigned int&) = &DivFreeBasisUtility::AssignQuadrature2D;

    class_<DivFreeBasisUtility, DivFreeBasisUtility::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("DivFreeBasisUtility", init<>() )
    .def("GetValues", &ComputeDivFreeBasis)
    .def("AssignQuadrature2D", pointer_to_AssignQuadrature2D)
    ;

    class_<MomentFittingUtility, MomentFittingUtility::Pointer, boost::noncopyable, bases<QuadratureUtility> >
    ("MomentFittingUtility", init<>())
//    .def("FitQuadrature", &MomentFittingUtility::PyFitQuadrature<BinaryTree<2> >)
//    .def("FitQuadrature", &MomentFittingUtility::PyFitQuadrature<BinaryTree<3> >)
    .def("FitQuadrature", &MomentFittingUtility::PyFitQuadrature<QuadTree>)
    .def("MultithreadedFitQuadrature", &MomentFittingUtility::PyMultithreadedFitQuadrature<QuadTree>)
    .def("MultithreadedFitQuadratureSubCell", &MomentFittingUtility::PyMultithreadedFitQuadratureSubCell<MomentFittedQuadTreeSubCell>)
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
    .def("MultithreadedQuadTreeRefineBy", &FiniteCellAuxilliaryUtility_MultithreadedQuadTreeRefineBy<BRep>)
    .def("MultithreadedQuadTreeSubCellRefineBy", &FiniteCellAuxilliaryUtility_MultithreadedQuadTreeSubCellRefineBy<BRep>)
    .def("Print", pointer_to_PrintGeometry)
    ;

}
}  // namespace Python.
}  // namespace Kratos.

