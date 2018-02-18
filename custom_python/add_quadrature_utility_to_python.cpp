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
#include "custom_python/add_quadrature_utility_to_python.h"
#include "custom_utilities/quadrature_utility.h"
#include "custom_utilities/quad_tree.h"
#include "custom_utilities/moment_fitted_quad_tree_subcell.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

int QuadratureUtility_GetQuadratureType(QuadratureUtility& rDummy, const int& integration_method)
{
    return rDummy.GetQuadratureType(integration_method);
}

int QuadratureUtility_GetQuadratureOrder(QuadratureUtility& rDummy, const int& integration_method)
{
    return rDummy.GetQuadratureOrder(integration_method);
}

void QuadratureUtility_ScaleQuadrature(QuadratureUtility& rDummy,
        Element::Pointer p_elem, const int quadrature_order, const double ScaleFactor)
{
    GeometryData::IntegrationMethod ElementalIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(quadrature_order);

    rDummy.ScaleQuadrature(p_elem->GetGeometry(), ElementalIntegrationMethod, ScaleFactor);
}

void QuadratureUtility_SetQuadrature(QuadratureUtility& rDummy,
        Element::Pointer p_elem,
        const int& integration_order,
        boost::python::list& quadrature_data)
{
    Element::GeometryType::IntegrationPointsArrayType integration_points;

    for(std::size_t i = 0; i < boost::python::len(quadrature_data); ++i)
    {
        boost::python::list point = boost::python::extract<boost::python::list>(quadrature_data[i]);

        Element::GeometryType::IntegrationPointType integration_point;
        integration_point.X() = boost::python::extract<double>(point[0]);
        integration_point.Y() = boost::python::extract<double>(point[1]);
        integration_point.Z() = boost::python::extract<double>(point[2]);
        integration_point.Weight() = boost::python::extract<double>(point[3]);

        integration_points.push_back(integration_point);
    }

    Element::GeometryType::IntegrationMethod ThisIntegrationMethod = Function<double, double>::GetIntegrationMethod(integration_order);
    FiniteCellGeometryUtility::AssignGeometryData(p_elem->GetGeometry(), ThisIntegrationMethod, integration_points);
//        std::cout << "set quadrature for element " << p_elem->Id() << " completed" << std::endl;
}

void QuadratureUtility_SaveQuadrature(QuadratureUtility& rDummy,
    boost::python::list& pyElemList, const std::string& fileName,
    const std::string& writeMode)
{
    std::ofstream myFile;

    int mode;
    if(writeMode == "binary")
    {
        myFile.open(fileName.c_str(), std::ios::out | std::ios::binary);
        mode = 0;
        myFile.write((char*)&mode, sizeof(int));
    }
    else if(writeMode == "binary-full")
    {
        myFile.open(fileName.c_str(), std::ios::out | std::ios::binary);
        mode = 1;
        myFile.write((char*)&mode, sizeof(int));
    }
    else if(writeMode == "ascii")
    {
        myFile.open(fileName.c_str(), std::ios::out);
        mode = 2;
        myFile << mode << std::endl;
    }
    else if(writeMode == "ascii-full")
    {
        myFile.open(fileName.c_str(), std::ios::out);
        mode = 3;
        myFile << mode << std::endl;
    }
    else
        KRATOS_THROW_ERROR(std::logic_error, "Unknown write mode", writeMode)

    myFile.precision(16);
    myFile << std::scientific;

    typedef boost::python::stl_input_iterator<Element::Pointer> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& p_elem,
                    std::make_pair(iterator_value_type(pyElemList), // begin
                    iterator_value_type() ) ) // end
    {
        switch(mode)
        {
            case 0:
                rDummy.SaveQuadrature<0>(myFile, p_elem, p_elem->GetGeometry().GetDefaultIntegrationMethod());
                break;
            case 1:
                rDummy.SaveQuadrature<1>(myFile, p_elem, p_elem->GetGeometry().GetDefaultIntegrationMethod());
                break;
            case 2:
                rDummy.SaveQuadrature<2>(myFile, p_elem, p_elem->GetGeometry().GetDefaultIntegrationMethod());
                break;
            case 3:
                rDummy.SaveQuadrature<3>(myFile, p_elem, p_elem->GetGeometry().GetDefaultIntegrationMethod());
                break;
        }
    }

    myFile.close();
    std::cout << "Write quadrature to " << fileName << " successfully" << std::endl;
}

void QuadratureUtility_SaveQuadratureAdvanced(QuadratureUtility& rDummy,
        const std::string& fileName,
        const std::string& fileType,
        boost::python::list& pyCutElems,
        boost::python::list& pyExcludeElems,
        const int& accuracy)
{
    typedef boost::python::stl_input_iterator<Element::Pointer> iterator_value_type;
    std::ofstream myFile;

    if(fileType == std::string("python"))
    {
        myFile.open(fileName.c_str(), std::ios::out);
        myFile.precision(accuracy);
        myFile << std::scientific;

        myFile << "def GetCutElements():\n";
        myFile << "\tcut_elems = [";

        int number_of_element_per_row = 20;

        int cnt = 0;
        BOOST_FOREACH(const iterator_value_type::value_type& p_elem,
                std::make_pair(iterator_value_type(pyCutElems), // begin
                iterator_value_type() ) ) // end
        {
            if(cnt < number_of_element_per_row)
            {
                ++cnt;
            }
            else
            {
                cnt = 0;
                myFile << "\n\t";
            }
            myFile << p_elem->Id() << ", ";
        }

        myFile << "\t]\n";
        myFile << "\treturn cut_elems\n";

        /////////////////////////////////

        myFile << "\ndef GetExcludeElements():\n";
        myFile << "\texclude_elems = [";

        cnt = 0;
        BOOST_FOREACH(const iterator_value_type::value_type& p_elem,
                std::make_pair(iterator_value_type(pyExcludeElems), // begin
                iterator_value_type() ) ) // end
        {
            if(cnt < number_of_element_per_row)
            {
                ++cnt;
            }
            else
            {
                cnt = 0;
                myFile << "\n\t";
            }
            myFile << p_elem->Id() << ", ";
        }

        myFile << "\t]\n";
        myFile << "\treturn exclude_elems\n";

        ////////////////////////////////

        myFile << "\ndef GetCutCellQuadrature():\n";
        myFile << "\tquad_data = {}\n";
        BOOST_FOREACH(const iterator_value_type::value_type& p_elem,
                std::make_pair(iterator_value_type(pyCutElems), // begin
                iterator_value_type() ) ) // end
        {
            const Element::GeometryType::IntegrationPointsArrayType& integration_points
                = p_elem->GetGeometry().IntegrationPoints( p_elem->GetGeometry().GetDefaultIntegrationMethod() );

            myFile << "\tquad_data[" << p_elem->Id() << "] = [";
            for(std::size_t i = 0; i < integration_points.size(); ++i)
            {
                myFile << "\n\t\t[" << integration_points[i].X()
                             << ", " << integration_points[i].Y()
                             << ", " << integration_points[i].Z()
                             << ", " << integration_points[i].Weight() << "],";
            }
            myFile << "] # " << integration_points.size() << "\n";
        }
        myFile << "\treturn quad_data\n";

        myFile.close();
        std::cout << "Write quadrature to " << fileName << " successfully" << std::endl;
    }
    else
        KRATOS_THROW_ERROR(std::logic_error, "Unknown file type", fileType)
}

template<class TTreeType>
void QuadratureUtility_SaveQuadratureAdvancedSubCell(QuadratureUtility& rDummy,
        const std::string& fileName,
        const std::string& fileType,
        boost::python::list& pyCutTrees,
        boost::python::list& pyExcludeTrees,
        boost::python::list& pyQuadTreeTrees,
        const int& accuracy)
{
    typedef typename TTreeType::Pointer TTreePointerType;
    typedef boost::python::stl_input_iterator<TTreePointerType> iterator_value_type;
    std::ofstream myFile;

    if(fileType == std::string("python"))
    {
        myFile.open(fileName.c_str(), std::ios::out);
        myFile.precision(accuracy);
        myFile << std::scientific;

        myFile << "def GetCutElements():\n";
        myFile << "\tcut_elems = [";

        int number_of_element_per_row = 20;

        int cnt = 0;
        std::size_t number_of_cut_elems = 0;
        BOOST_FOREACH(const typename iterator_value_type::value_type& p_tree,
                std::make_pair(iterator_value_type(pyCutTrees), // begin
                iterator_value_type() ) ) // end
        {
            if(cnt < number_of_element_per_row)
            {
                ++cnt;
            }
            else
            {
                cnt = 0;
                myFile << "\n\t";
            }
            ++number_of_cut_elems;
            myFile << p_tree->pGetElement()->Id() << ", ";
        }

        myFile << "\t]\n";
        myFile << "##number of cut elements: " << number_of_cut_elems << "\n";
        myFile << "\treturn cut_elems\n";

        /////////////////////////////////

        myFile << "\ndef GetExcludeElements():\n";
        myFile << "\texclude_elems = [";

        cnt = 0;
        std::size_t number_of_exclude_elems = 0;
        BOOST_FOREACH(const typename iterator_value_type::value_type& p_tree,
                std::make_pair(iterator_value_type(pyExcludeTrees), // begin
                iterator_value_type() ) ) // end
        {
            if(cnt < number_of_element_per_row)
            {
                ++cnt;
            }
            else
            {
                cnt = 0;
                myFile << "\n\t";
            }
            ++number_of_exclude_elems;
            myFile << p_tree->pGetElement()->Id() << ", ";
        }

        myFile << "\t]\n";
        myFile << "##number of exclude elements: " << number_of_exclude_elems << "\n";
        myFile << "\treturn exclude_elems\n";

        ////////////////////////////////

        myFile << "\ndef GetQuadTreeElements():\n";
        myFile << "\tquadtree_elems = [";

        cnt = 0;
        std::size_t number_of_quadtree_elems = 0;
        BOOST_FOREACH(const typename iterator_value_type::value_type& p_tree,
                std::make_pair(iterator_value_type(pyQuadTreeTrees), // begin
                iterator_value_type() ) ) // end
        {
            if(cnt < number_of_element_per_row)
            {
                ++cnt;
            }
            else
            {
                cnt = 0;
                myFile << "\n\t";
            }
            ++number_of_quadtree_elems;
            myFile << p_tree->pGetElement()->Id() << ", ";
        }

        myFile << "\t]\n";
        myFile << "##number of quadtree elements: " << number_of_quadtree_elems << "\n";
        myFile << "\treturn quadtree_elems\n";

        ////////////////////////////////

        myFile << "\ndef GetCutCellQuadrature():\n";
        myFile << "\tcq = {}\n";
        BOOST_FOREACH(const typename iterator_value_type::value_type& p_tree,
                std::make_pair(iterator_value_type(pyCutTrees), // begin
                iterator_value_type() ) ) // end
        {
            const Element::GeometryType::IntegrationPointsArrayType& integration_points
                = p_tree->pGetElement()->GetGeometry().IntegrationPoints( p_tree->pGetElement()->GetGeometry().GetDefaultIntegrationMethod() );

            myFile << "\tcq[" << p_tree->pGetElement()->Id() << "] = [";
            for(std::size_t i = 0; i < integration_points.size(); ++i)
            {
                myFile << "\n\t\t[" << integration_points[i].X()
                             << ", " << integration_points[i].Y()
                             << ", " << integration_points[i].Z()
                             << ", " << integration_points[i].Weight() << "],";
            }
            myFile << "] # " << integration_points.size() << "\n";
        }

        myFile << "\treturn cq\n";

        ////////////////////////////////

        myFile << "\ndef GetCutCellFictitiousQuadrature():\n";
        myFile << "\tfiq = {}\n";
        BOOST_FOREACH(const typename iterator_value_type::value_type& p_tree,
                std::make_pair(iterator_value_type(pyCutTrees), // begin
                iterator_value_type() ) ) // end
        {
            const Element::GeometryType::IntegrationPointsArrayType& integration_points
                = p_tree->GetFictitiousIntegrationPoints();

            myFile << "\tfiq[" << p_tree->pGetElement()->Id() << "] = [";
            for(std::size_t i = 0; i < integration_points.size(); ++i)
            {
                myFile << "\n\t\t[" << integration_points[i].X()
                             << ", " << integration_points[i].Y()
                             << ", " << integration_points[i].Z()
                             << ", " << integration_points[i].Weight() << "],";
            }
            myFile << "] # " << integration_points.size() << "\n";
        }

        myFile << "\treturn fiq\n";

        ////////////////////////////////

        myFile << "\ndef GetCutCellFullQuadrature():\n";
        myFile << "\tfq = {}\n";
        BOOST_FOREACH(const typename iterator_value_type::value_type& p_tree,
                std::make_pair(iterator_value_type(pyCutTrees), // begin
                iterator_value_type() ) ) // end
        {
            const Element::GeometryType::IntegrationPointsArrayType& integration_points
                = p_tree->GetRepresentativeIntegrationPoints();

            myFile << "\tfq[" << p_tree->pGetElement()->Id() << "] = [";
            for(std::size_t i = 0; i < integration_points.size(); ++i)
            {
                myFile << "\n\t\t[" << integration_points[i].X()
                             << ", " << integration_points[i].Y()
                             << ", " << integration_points[i].Z()
                             << ", " << integration_points[i].Weight() << "],";
            }
            myFile << "]\n";
        }

        myFile << "\treturn fq\n";

        ////////////////////////////////

        myFile << "\ndef GetCutCellSubCellWeights():\n";
        myFile << "\tsw = {}\n";

        BOOST_FOREACH(const typename iterator_value_type::value_type& p_tree,
                std::make_pair(iterator_value_type(pyCutTrees), // begin
                iterator_value_type() ) ) // end
        {
            const Matrix& Weights = p_tree->pGetElement()->GetValue(SUBCELL_WEIGHTS);

            myFile << "\tsw[" << p_tree->pGetElement()->Id() << "] = [";
            for(std::size_t i = 0; i < Weights.size1(); ++i)
            {
                myFile << "\n\t\t[";
                for(std::size_t j = 0; j < Weights.size2(); ++j)
                    myFile << Weights(i, j) << ", ";
                myFile << "],";
            }
            myFile << "]\n";
        }

        myFile << "\treturn sw\n";

        ////////////////////////////////

        myFile << "\ndef GetCutCellSubCellDomainSizes():\n";
        myFile << "\tds = {}\n";

        BOOST_FOREACH(const typename iterator_value_type::value_type& p_tree,
                std::make_pair(iterator_value_type(pyCutTrees), // begin
                iterator_value_type() ) ) // end
        {
            const Vector& DomainSizes = p_tree->pGetElement()->GetValue(SUBCELL_DOMAIN_SIZES);

            myFile << "\tds[" << p_tree->pGetElement()->Id() << "] = [";
            for(std::size_t i = 0; i < DomainSizes.size(); ++i)
            {
                myFile << DomainSizes(i) << ", ";
            }
            myFile << "]\n";
        }

        myFile << "\treturn ds\n";

        ////////////////////////////////

        myFile << "\ndef GetQuadTreeQuadrature():\n";
        myFile << "\tqq = {}\n";

        BOOST_FOREACH(const typename iterator_value_type::value_type& p_tree,
                std::make_pair(iterator_value_type(pyQuadTreeTrees), // begin
                iterator_value_type() ) ) // end
        {
            const Element::GeometryType::IntegrationPointsArrayType& integration_points
                = p_tree->pGetElement()->GetGeometry().IntegrationPoints( p_tree->pGetElement()->GetGeometry().GetDefaultIntegrationMethod() );

            myFile << "\tqq[" << p_tree->pGetElement()->Id() << "] = [";
            for(std::size_t i = 0; i < integration_points.size(); ++i)
            {
                myFile << "\n\t\t[" << integration_points[i].X()
                             << ", " << integration_points[i].Y()
                             << ", " << integration_points[i].Z()
                             << ", " << integration_points[i].Weight() << "],";
            }
            myFile << "]\n";
        }

        myFile << "\treturn qq\n";

        ////////////////////////////////

        myFile.close();
        std::cout << "Write quadrature to " << fileName << " successfully" << std::endl;
    }
    else
        KRATOS_THROW_ERROR(std::logic_error, "Unknown file type", fileType)
}

/// Create a new condition from point
ModelPart::ConditionsContainerType QuadratureUtility_CreateConditionFromPoint(QuadratureUtility& rDummy,
        ModelPart& r_model_part,
        boost::python::list& pyPointList,
        const std::string& sample_cond_name,
        Properties::Pointer pProperties)
{
    // find the maximum node Id
    std::size_t lastNodeId = FiniteCellAuxilliaryUtility::GetLastNodeId(r_model_part);

    // find the maximum condition Id
    std::size_t lastCondId = FiniteCellAuxilliaryUtility::GetLastConditionId(r_model_part);

    // find the maximum properties Id
    std::size_t lastPropId = FiniteCellAuxilliaryUtility::GetLastPropertiesId(r_model_part);

    // get the sample condition
    if(!KratosComponents<Condition>::Has(sample_cond_name))
        KRATOS_THROW_ERROR(std::logic_error, sample_cond_name, "is not registered to the KRATOS kernel")
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_cond_name);

    ModelPart::ConditionsContainerType NewConditions;

    typedef Element::GeometryType::PointType::PointType PointType;
    typedef boost::python::stl_input_iterator<PointType> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& point,
                    std::make_pair(iterator_value_type(pyPointList), // begin
                    iterator_value_type() ) ) // end
    {
        Condition::Pointer pNewCond = rDummy.CreateConditionFromPoint(r_model_part, point, r_clone_condition, pProperties, lastNodeId, lastCondId);
        NewConditions.push_back(pNewCond);
    }

    for(typename ModelPart::ConditionsContainerType::ptr_iterator it = NewConditions.ptr_begin();
            it != NewConditions.ptr_end(); ++it)
        r_model_part.Conditions().push_back(*it);

    return NewConditions;
}

/// Create a new condition from point
ModelPart::ConditionsContainerType QuadratureUtility_CreateConditionFromPoint2(QuadratureUtility& rDummy,
        ModelPart& r_model_part,
        boost::python::list& pyPointList,
        const std::string& sample_cond_name)
{
    // find the maximum properties Id
    std::size_t lastPropId = FiniteCellAuxilliaryUtility::GetLastPropertiesId(r_model_part);

    // create new properties
    Properties::Pointer pProperties = Properties::Pointer(new Properties(++lastPropId));
    r_model_part.AddProperties(pProperties);

    return QuadratureUtility_CreateConditionFromPoint(rDummy, r_model_part, pyPointList, sample_cond_name, pProperties);
}

void QuadratureUtility_CreateConditionFromQuadraturePoint(QuadratureUtility& rDummy,
        ModelPart& r_model_part,
        boost::python::list& pyElemList,
        const std::string& sample_cond_name,
        const double& min_weight,
        const double& max_weight)
{
    // find the maximum node Id
    std::size_t lastNodeId = FiniteCellAuxilliaryUtility::GetLastNodeId(r_model_part);

    // find the maximum condition Id
    std::size_t lastCondId = FiniteCellAuxilliaryUtility::GetLastConditionId(r_model_part);

    // find the maximum properties Id
    std::size_t lastPropId = FiniteCellAuxilliaryUtility::GetLastPropertiesId(r_model_part);

    // create new properties
    Properties::Pointer pProperties = Properties::Pointer(new Properties(++lastPropId));
    r_model_part.AddProperties(pProperties);

    // get the sample condition
    if(!KratosComponents<Condition>::Has(sample_cond_name))
        KRATOS_THROW_ERROR(std::logic_error, sample_cond_name, "is not registered to the KRATOS kernel")
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_cond_name);

    std::size_t num_conds = 0;
    typedef boost::python::stl_input_iterator<Element::Pointer> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& p_elem,
                    std::make_pair(iterator_value_type(pyElemList), // begin
                    iterator_value_type() ) ) // end
    {
        GeometryData::IntegrationMethod ThisIntegrationMethod = p_elem->GetGeometry().GetDefaultIntegrationMethod();
        num_conds += rDummy.CreateConditionFromQuadraturePoint(r_model_part, p_elem, ThisIntegrationMethod,
                r_clone_condition, pProperties, lastNodeId, lastCondId, min_weight, max_weight);
    }

    std::cout << num_conds << " conditions of type " << sample_cond_name << " is added to the model_part" << std::endl;
}

void FiniteCellApplication_AddQuadratureUtilityToPython()
{

    class_<QuadratureUtility, QuadratureUtility::Pointer, boost::noncopyable>
    ("QuadratureUtility", init<>())
    .def("GetDefaultIntegrationMethod", &QuadratureUtility::GetDefaultIntegrationMethod<Element>)
    .def("GetDefaultIntegrationMethod", &QuadratureUtility::GetDefaultIntegrationMethod<Condition>)
    .def("GetQuadratureType", QuadratureUtility_GetQuadratureType)
    .def("GetQuadratureOrder", QuadratureUtility_GetQuadratureOrder)
    .def("ScaleQuadrature", &QuadratureUtility_ScaleQuadrature)
    .def("SaveQuadrature", &QuadratureUtility_SaveQuadrature)
    .def("SaveQuadrature", &QuadratureUtility_SaveQuadratureAdvanced)
    .def("SaveQuadratureSubCell", &QuadratureUtility_SaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<1> >)
    .def("SaveQuadratureSubCell2", &QuadratureUtility_SaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<2> >)
    .def("SaveQuadratureSubCell3", &QuadratureUtility_SaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<3> >)
    .def("SaveQuadratureSubCell4", &QuadratureUtility_SaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<4> >)
    .def("SaveQuadratureSubCell5", &QuadratureUtility_SaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<5> >)
    .def("SaveQuadratureSubCell6", &QuadratureUtility_SaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<6> >)
    .def("SaveQuadratureSubCell7", &QuadratureUtility_SaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<7> >)
    .def("SaveQuadratureSubCell8", &QuadratureUtility_SaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<8> >)
    .def("SaveQuadratureSubCell9", &QuadratureUtility_SaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<9> >)
    .def("SaveQuadratureSubCell10", &QuadratureUtility_SaveQuadratureAdvancedSubCell<MomentFittedQuadTreeSubCell<10> >)
    .def("SetQuadrature", &QuadratureUtility_SetQuadrature)
    .def("CreateConditionFromQuadraturePoint", &QuadratureUtility_CreateConditionFromQuadraturePoint)
    .def("CreateConditionFromPoint", &QuadratureUtility_CreateConditionFromPoint)
    .def("CreateConditionFromPoint", &QuadratureUtility_CreateConditionFromPoint2)
    .def("CreatePoint", &QuadratureUtility::CreatePoint)
    ;

}

}  // namespace Python.

}  // namespace Kratos.

