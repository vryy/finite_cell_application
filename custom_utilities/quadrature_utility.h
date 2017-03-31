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
        FiniteCellGeometry<GeometryType>::AssignGeometryData(r_geom, ElementalIntegrationMethod, Mw);
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
        FiniteCellGeometry<GeometryType>::AssignGeometryData(p_elem->GetGeometry(), ThisIntegrationMethod, integration_points);
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
                myFile << "]\n";
            }
            myFile << "\treturn quad_data\n";

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
        std::size_t lastNodeId = 0;
        for(typename ModelPart::NodesContainerType::ptr_iterator it = r_model_part.Nodes().ptr_begin();
                it != r_model_part.Nodes().ptr_end(); ++it)
        {
            if((*it)->Id() > lastNodeId)
                lastNodeId = (*it)->Id();
        }

        // find the maximum condition Id
        std::size_t lastCondId = 0;
        for(typename ModelPart::ConditionsContainerType::ptr_iterator it = r_model_part.Conditions().ptr_begin();
                it != r_model_part.Conditions().ptr_end(); ++it)
        {
            if((*it)->Id() > lastCondId)
                lastCondId = (*it)->Id();
        }

        // find the maximum properties Id
        std::size_t lastPropId = 0;
        for(typename ModelPart::PropertiesContainerType::ptr_iterator it = r_model_part.rProperties().ptr_begin();
                it != r_model_part.rProperties().ptr_end(); ++it)
        {
            if((*it)->Id() > lastPropId)
                lastPropId = (*it)->Id();
        }

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
        GeometryType::Pointer pTempGeometry;
        for(std::size_t point = 0; point < integration_points.size(); ++point)
        {
            if(integration_points[point].Weight() > min_weight && integration_points[point].Weight() < max_weight)
            {
                p_elem->GetGeometry().GlobalCoordinates(GlobalCoords, integration_points[point]);
                NodeType::Pointer pNewNode(new NodeType(++lastNodeId, GlobalCoords));
                pNewNode->SetSolutionStepVariablesList(&r_model_part.GetNodalSolutionStepVariablesList());
                pNewNode->SetBufferSize(r_model_part.GetBufferSize());
                r_model_part.AddNode(pNewNode);
    //            NodeType::Pointer pNewNode = r_model_part.CreateNewNode(++lastNodeId, GlobalCoords[0], GlobalCoords[1], GlobalCoords[2]);
    //            KRATOS_WATCH(*pNewNode)
                pTempGeometry = GeometryType::Pointer( new Point3D<NodeType>(pNewNode) );
                Condition::Pointer pNewCond = r_sample_cond.Create(++lastCondId, pTempGeometry, pProperties);
                pNewCond->Set(ACTIVE, true);
                NewConditions.push_back(pNewCond);
            }
        }

        for(typename ModelPart::ConditionsContainerType::ptr_iterator it = NewConditions.ptr_begin();
                it != NewConditions.ptr_end(); ++it)
            r_model_part.Conditions().push_back(*it);

        return NewConditions.size();
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
