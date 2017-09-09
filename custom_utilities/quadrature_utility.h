//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         finite_cell_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            23 Feb 2017
//


#if !defined(KRATOS_QUADRATURE_UTILITY_H_INCLUDED )
#define  KRATOS_QUADRATURE_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <fstream>


// External includes
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/numeric/ublas/lu.hpp> 


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_flags.h"
#include "geometries/point_3d.h"
#include "custom_geometries/finite_cell_geometry.h"
#include "custom_algebra/function/function.h"
#include "custom_utilities/finite_cell_geometry_utility.h"
#include "custom_utilities/finite_cell_auxilliary_utility.h"
#include "finite_cell_application/finite_cell_application.h"


namespace Kratos
{
///@addtogroup FiniteCellApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Abstract class to compute/assign the quadrature weights for elements/conditions
*/
class QuadratureUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of QuadratureUtility
    KRATOS_CLASS_POINTER_DEFINITION(QuadratureUtility);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    QuadratureUtility() {}

    /// Destructor.
    virtual ~QuadratureUtility() {}


    template<class TEntityType>
    int GetDefaultIntegrationMethod(typename TEntityType::Pointer p_elem) const
    {
        return p_elem->GetIntegrationMethod();
//        return p_elem->GetGeometry().GetDefaultIntegrationMethod();
    }


    /// The integration method is decoded as follow:
    /// It is a hex number, i.e 0x02
    /// The first hex digit is the quadrature method: i.e
    ///     0: Kratos built-in quadrature (it's technically Gauss-Legendre)
    ///     1: extended Gauss-Legendre
    ///     2: Gauss-Lobatto quadrature
    ///     3: Gauss-Radau quadrature
    /// The second hex digit is the order of the quadrature
    /// By default, there're 16 types of quadrature and order upto 15 is supported. Probably it's enough for conventional problem.
    static inline int GetQuadratureType(const int& integration_method)
    {
        return (integration_method & 0xF0) >> 4;
    }

    static inline int GetQuadratureOrder(const int& integration_method)
    {
        return (integration_method & 0x0F);
    }


    void PyScaleQuadrature(Element::Pointer& p_elem,
            const int quadrature_order, const double ScaleFactor) const
    {
        GeometryData::IntegrationMethod ElementalIntegrationMethod
                = Function<double, double>::GetIntegrationMethod(quadrature_order);

        this->ScaleQuadrature(p_elem->GetGeometry(), ElementalIntegrationMethod, ScaleFactor);
    }


    void ScaleQuadrature(GeometryType& r_geom,
            const GeometryData::IntegrationMethod& ElementalIntegrationMethod,
            const double ScaleFactor) const
    {
        const GeometryType::IntegrationPointsArrayType& integration_points
                = r_geom.IntegrationPoints( ElementalIntegrationMethod );

        Vector Mw(integration_points.size());
        for(std::size_t i = 0; i < integration_points.size(); ++i)
        {
            Mw(i) = integration_points[i].Weight() * ScaleFactor;
        }

        /* create new quadrature and assign to the geometry */
        FiniteCellGeometryUtility::AssignGeometryData(r_geom, ElementalIntegrationMethod, Mw);
    }


    void PySetQuadrature(Element::Pointer& p_elem,
            const int& integration_order,
            boost::python::list& quadrature_data) const
    {
        GeometryType::IntegrationPointsArrayType integration_points;

        for(std::size_t i = 0; i < boost::python::len(quadrature_data); ++i)
        {
            boost::python::list point = boost::python::extract<boost::python::list>(quadrature_data[i]);

            GeometryType::IntegrationPointType integration_point;
            integration_point.X() = boost::python::extract<double>(point[0]);
            integration_point.Y() = boost::python::extract<double>(point[1]);
            integration_point.Z() = boost::python::extract<double>(point[2]);
            integration_point.Weight() = boost::python::extract<double>(point[3]);

            integration_points.push_back(integration_point);
        }

        GeometryType::IntegrationMethod ThisIntegrationMethod = Function<double, double>::GetIntegrationMethod(integration_order);
        FiniteCellGeometryUtility::AssignGeometryData(p_elem->GetGeometry(), ThisIntegrationMethod, integration_points);
//        std::cout << "set quadrature for element " << p_elem->Id() << " completed" << std::endl;
    }


    void PySaveQuadratureAdvanced(const std::string& fileName,
            const std::string& fileType,
            boost::python::list& pyCutElems,
            boost::python::list& pyExcludeElems,
            const int& accuracy) const
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
                const GeometryType::IntegrationPointsArrayType& integration_points
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
    void PySaveQuadratureAdvancedSubCell(const std::string& fileName,
            const std::string& fileType,
            boost::python::list& pyCutTrees,
            boost::python::list& pyExcludeTrees,
            boost::python::list& pyQuadTreeTrees,
            const int& accuracy) const
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
                const GeometryType::IntegrationPointsArrayType& integration_points
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

            myFile << "\ndef GetCutCellFullQuadrature():\n";
            myFile << "\tfq = {}\n";
            BOOST_FOREACH(const typename iterator_value_type::value_type& p_tree,
                    std::make_pair(iterator_value_type(pyCutTrees), // begin
                    iterator_value_type() ) ) // end
            {
                const GeometryType::IntegrationPointsArrayType& integration_points
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
                const GeometryType::IntegrationPointsArrayType& integration_points
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


    void PySaveQuadrature(boost::python::list& pyElemList, const std::string& fileName,
            const std::string& writeMode) const
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
                    SaveQuadrature<0>(myFile, p_elem, p_elem->GetGeometry().GetDefaultIntegrationMethod());
                    break;
                case 1:
                    SaveQuadrature<1>(myFile, p_elem, p_elem->GetGeometry().GetDefaultIntegrationMethod());
                    break;
                case 2:
                    SaveQuadrature<2>(myFile, p_elem, p_elem->GetGeometry().GetDefaultIntegrationMethod());
                    break;
                case 3:
                    SaveQuadrature<3>(myFile, p_elem, p_elem->GetGeometry().GetDefaultIntegrationMethod());
                    break;
            }
        }

        myFile.close();
        std::cout << "Write quadrature to " << fileName << " successfully" << std::endl;
    }


    template<int TMode>
    void SaveQuadrature(std::ofstream& rOStream, const Element::Pointer& p_elem,
            const GeometryData::IntegrationMethod& ElementalIntegrationMethod) const
    {
        const GeometryType::IntegrationPointsArrayType& integration_points
                = p_elem->GetGeometry().IntegrationPoints( ElementalIntegrationMethod );

        if(TMode == 0 || TMode == 1) // binary
        {
            std::size_t elem_id = p_elem->Id();
            std::size_t num_points = integration_points.size();
            rOStream.write((char*)&elem_id, sizeof(std::size_t));
            rOStream.write((char*)&ElementalIntegrationMethod, sizeof(int));
            rOStream.write((char*)&num_points, sizeof(std::size_t));

            double aux;
            if(TMode == 1)
            {
                for(std::size_t i = 0; i < integration_points.size(); ++i)
                {
                    aux = integration_points[i].X();
                    rOStream.write((char*)&aux, sizeof(double));

                    aux = integration_points[i].Y();
                    rOStream.write((char*)&aux, sizeof(double));

                    aux = integration_points[i].Z();
                    rOStream.write((char*)&aux, sizeof(double));
                }
            }

            for(std::size_t i = 0; i < integration_points.size(); ++i)
            {
                aux = integration_points[i].Weight();
                rOStream.write((char*)&aux, sizeof(double));
            }
        }
        else if(TMode == 2 || TMode == 3) // ascii
        {
            rOStream << p_elem->Id() << " " << ElementalIntegrationMethod
                     << " " << integration_points.size() << std::endl;

            for(std::size_t i = 0; i < integration_points.size(); ++i)
            {
                if(TMode == 3)
                {
                    rOStream << "  " << integration_points[i].X()
                             << " \t" << integration_points[i].Y()
                             << " \t" << integration_points[i].Z();
                }
                rOStream << " \t" << integration_points[i].Weight() << std::endl;
            }
        }
    }


    void PyCreateConditionFromQuadraturePoint(ModelPart& r_model_part,
        boost::python::list& pyElemList,
        const std::string& sample_cond_name,
        const double& min_weight,
        const double& max_weight) const
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
            num_conds += this->CreateConditionFromQuadraturePoint(r_model_part, p_elem, ThisIntegrationMethod,
                    r_clone_condition, pProperties, lastNodeId, lastCondId, min_weight, max_weight);
        }

        std::cout << num_conds << " conditions of type " << sample_cond_name << " is added to the model_part" << std::endl;
    }


    /// Create an array of condition (technically point condition) from quadrature point of the element
    /// r_model_part    the conditions after creation will be added to the model_part
    /// ElementalIntegrationMethod  the integration method of the element, which the quadrature points are read
    /// r_sample_cond   the condition name in Kratos database
    /// pProperties     the Properties used to create new condition
    /// lastNodeId      starting node id to assign to new nodes
    /// lastCondId      starting condition id to assign to new conditions
    /// min_weight      minimum weight to create new condition
    /// max_weight      maximum weight to create new condition
    std::size_t CreateConditionFromQuadraturePoint(ModelPart& r_model_part,
            const Element::Pointer& p_elem,
            const GeometryData::IntegrationMethod& ElementalIntegrationMethod,
            Condition const& r_sample_cond,
            Properties::Pointer pProperties,
            std::size_t& lastNodeId,
            std::size_t& lastCondId,
            const double& min_weight,
            const double& max_weight
            ) const
    {
        const GeometryType::IntegrationPointsArrayType& integration_points
                = p_elem->GetGeometry().IntegrationPoints( ElementalIntegrationMethod );

        ModelPart::ConditionsContainerType NewConditions;

        PointType GlobalCoords;
        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            if(integration_points[point].Weight() > min_weight && integration_points[point].Weight() < max_weight)
            {
                p_elem->GetGeometry().GlobalCoordinates(GlobalCoords, integration_points[point]);
                Condition::Pointer pNewCond = CreateConditionFromPoint(r_model_part, GlobalCoords, r_sample_cond, pProperties, lastNodeId, lastCondId);
                NewConditions.push_back(pNewCond);
            }
        }

        for(typename ModelPart::ConditionsContainerType::ptr_iterator it = NewConditions.ptr_begin();
                it != NewConditions.ptr_end(); ++it)
            r_model_part.Conditions().push_back(*it);

        return NewConditions.size();
    }

    /// Create a new condition from point
    ModelPart::ConditionsContainerType PyCreateConditionFromPoint(ModelPart& r_model_part,
            boost::python::list& pyPointList,
            const std::string& sample_cond_name
            ) const
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

        ModelPart::ConditionsContainerType NewConditions;

        typedef boost::python::stl_input_iterator<PointType> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& point,
                        std::make_pair(iterator_value_type(pyPointList), // begin
                        iterator_value_type() ) ) // end
        {
            Condition::Pointer pNewCond = CreateConditionFromPoint(r_model_part, point, r_clone_condition, pProperties, lastNodeId, lastCondId);
            NewConditions.push_back(pNewCond);
        }

        for(typename ModelPart::ConditionsContainerType::ptr_iterator it = NewConditions.ptr_begin();
                it != NewConditions.ptr_end(); ++it)
            r_model_part.Conditions().push_back(*it);

        return NewConditions;
    }

    /// Create a new condition from point
    Condition::Pointer CreateConditionFromPoint(ModelPart& r_model_part,
            const PointType& r_point,
            Condition const& r_sample_cond,
            Properties::Pointer pProperties,
            std::size_t& lastNodeId,
            std::size_t& lastCondId
            ) const
    {
        NodeType::Pointer pNewNode(new NodeType(++lastNodeId, r_point));
        pNewNode->SetSolutionStepVariablesList(&r_model_part.GetNodalSolutionStepVariablesList());
        pNewNode->SetBufferSize(r_model_part.GetBufferSize());
        r_model_part.AddNode(pNewNode);
//            NodeType::Pointer pNewNode = r_model_part.CreateNewNode(++lastNodeId, GlobalCoords[0], GlobalCoords[1], GlobalCoords[2]);
//            KRATOS_WATCH(*pNewNode)
        GeometryType::Pointer pTempGeometry = GeometryType::Pointer( new Point3D<NodeType>(pNewNode) );
        Condition::Pointer pNewCond = r_sample_cond.Create(++lastCondId, pTempGeometry, pProperties);
        pNewCond->Set(ACTIVE, true);
        return pNewCond;
    }

    /// Create a point from coordinates
    PointType CreatePoint(const double& rX, const double& rY, const double& rZ) const
    {
        PointType P(rX, rY, rZ);
        return P;
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Quadrature Utility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    QuadratureUtility& operator=(QuadratureUtility const& rOther);

    /// Copy constructor.
    QuadratureUtility(QuadratureUtility const& rOther);


    ///@}

}; // Class QuadratureUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, QuadratureUtility& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const QuadratureUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_QUADRATURE_UTILITY_H_INCLUDED  defined
