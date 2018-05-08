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
//  Date:            3 Feb 2018
//


#if !defined(KRATOS_FINITE_CELL_MESH_UTILITY_H_INCLUDED )
#define  KRATOS_FINITE_CELL_MESH_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <fstream>


// External includes
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/foreach.hpp>
#include <boost/progress.hpp>


// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "custom_utilities/finite_cell_auxilliary_utility.h"
#include "custom_utilities/finite_cell_geometry_utility.h"


namespace Kratos
{
///@addtogroup FiniteCellApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

template<int TDim, int TOrder>
struct GenerateStructuredMesh_Helper
{
    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    GenerateStructuredMesh_Helper(const PointType& StartPoint,
            const PointType& EndPoint,
            const std::vector<std::size_t>& nsampling
    ) : mStartPoint(StartPoint), mEndPoint(EndPoint), mnsampling(nsampling)
    {}

    void Execute(std::vector<std::vector<PointType> >& rPoints)
    {
        std::stringstream ss;
        ss << "Error calling unimplemented " << __FUNCTION__ << "_" << TDim << "_" << TOrder;
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    void Execute(std::vector<std::vector<std::vector<PointType> > >& rPoints)
    {
        std::stringstream ss;
        ss << "Error calling unimplemented " << __FUNCTION__ << "_" << TDim << "_" << TOrder;
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    const PointType& mStartPoint;
    const PointType& mEndPoint;
    const std::vector<std::size_t>& mnsampling;
};

template<>
struct GenerateStructuredMesh_Helper<2, 1> : public GenerateStructuredMesh_Helper<0, 0>
{
    typedef GenerateStructuredMesh_Helper<0, 0> BaseType;

    GenerateStructuredMesh_Helper(const PointType& StartPoint,
            const PointType& EndPoint,
            const std::vector<std::size_t>& nsampling)
    : BaseType(StartPoint, EndPoint, nsampling)
    {}

    void Execute(std::vector<std::vector<PointType> >& rPoints)
    {
        double x, y;
        double dx = (mEndPoint[0] - mStartPoint[0]) / mnsampling[0];
        double dy = (mEndPoint[1] - mStartPoint[1]) / mnsampling[1];

        rPoints.resize(mnsampling[1]+1);
        for (std::size_t j = 0; j < mnsampling[1]+1; ++j)
        {
            rPoints[j].resize(mnsampling[0]+1);
            y = mStartPoint[1] + j*dy;

            for (std::size_t i = 0; i < mnsampling[0]+1; ++i)
            {
                x = mStartPoint[0] + i*dx;
                rPoints[j][i] = PointType(x, y, 0.0);
            }
        }
    }
};

template<>
struct GenerateStructuredMesh_Helper<3, 1> : public GenerateStructuredMesh_Helper<0, 0>
{
    typedef GenerateStructuredMesh_Helper<0, 0> BaseType;

    GenerateStructuredMesh_Helper(const PointType& StartPoint,
            const PointType& EndPoint,
            const std::vector<std::size_t>& nsampling)
    : BaseType(StartPoint, EndPoint, nsampling)
    {}

    void Execute(std::vector<std::vector<std::vector<PointType> > >& rPoints)
    {
        double x, y, z;
        double dx = (mEndPoint[0] - mStartPoint[0]) / mnsampling[0];
        double dy = (mEndPoint[1] - mStartPoint[1]) / mnsampling[1];
        double dz = (mEndPoint[2] - mStartPoint[2]) / mnsampling[2];

        rPoints.resize(mnsampling[2]+1);
        for (std::size_t k = 0; k < mnsampling[2]+1; ++k)
        {
            rPoints[k].resize(mnsampling[1]+1);
            z = mStartPoint[2] + k*dz;

            for (std::size_t j = 0; j < mnsampling[1]+1; ++j)
            {
                rPoints[k][j].resize(mnsampling[0]+1);
                y = mStartPoint[1] + j*dy;

                for (std::size_t i = 0; i < mnsampling[0]+1; ++i)
                {
                    x = mStartPoint[0] + i*dx;
                    rPoints[k][j][i] = PointType(x, y, z);
                }
            }
        }
    }
};

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
/** class for auxiliary mesh routines
*/
class FiniteCellMeshUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FiniteCellMeshUtility
    KRATOS_CLASS_POINTER_DEFINITION(FiniteCellMeshUtility);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FiniteCellMeshUtility() {}

    /// Destructor.
    virtual ~FiniteCellMeshUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Generate the background structure mesh
    static void GenerateStructuredMesh2D(std::vector<std::vector<PointType> >& sampling_points,
        const int& type,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::size_t>& nsampling)
    {
        if (type == 1)
        {
            GenerateStructuredMesh_Helper<2, 1>(StartPoint, EndPoint, nsampling).Execute(sampling_points);
        }
        else if (type == 2)
        {
            GenerateStructuredMesh_Helper<2, 2>(StartPoint, EndPoint, nsampling).Execute(sampling_points);
        }
        else if (type == 3)
        {
            GenerateStructuredMesh_Helper<2, 3>(StartPoint, EndPoint, nsampling).Execute(sampling_points);
        }
    }

    /// Generate the background structure mesh
    static void GenerateStructuredMesh3D(std::vector<std::vector<std::vector<PointType> > >& sampling_points,
        const int& type,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::size_t>& nsampling)
    {
        if (type == 1)
        {
            GenerateStructuredMesh_Helper<3, 1>(StartPoint, EndPoint, nsampling).Execute(sampling_points);
        }
        else if (type == 2)
        {
            GenerateStructuredMesh_Helper<3, 2>(StartPoint, EndPoint, nsampling).Execute(sampling_points);
        }
        else if (type == 3)
        {
            GenerateStructuredMesh_Helper<3, 3>(StartPoint, EndPoint, nsampling).Execute(sampling_points);
        }
    }


    /// Create the quad elements based on given points list
    static std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> CreateQuadElements(ModelPart& r_model_part,
        const std::vector<std::vector<PointType> >& sampling_points,
        const std::string& sample_element_name,
        const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
        const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
        const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
        Properties::Pointer pProperties)
    {
        std::size_t last_node_id = FiniteCellAuxilliaryUtility::GetLastNodeId(r_model_part);
        std::size_t last_node_id_old = last_node_id;
        std::size_t num_division_1 = sampling_points.size() - 1;
        std::size_t num_division_2 = sampling_points[0].size() - 1;
        // KRATOS_WATCH(last_node_id)

        std::size_t num_1, num_2;
        if (close_dir == 1)
        {
            num_1 = num_division_1 + 1;
            num_2 = num_division_2;
        }
        else if (close_dir == 2)
        {
            num_1 = num_division_1;
            num_2 = num_division_2 + 1;
        }
        else
        {
            num_1 = num_division_1;
            num_2 = num_division_2;
        }

        Variable<int>& ACTIVATION_LEVEL_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("ACTIVATION_LEVEL"));

        // firstly create nodes and add to model_part
        ModelPart::NodesContainerType NewNodes;
        for (std::size_t i = 0; i < num_division_1+1; ++i)
        {
            for (std::size_t j = 0; j < num_division_2+1; ++j)
            {
                NodeType::Pointer pNewNode = r_model_part.CreateNewNode(++last_node_id,
                        sampling_points[i][j][0], sampling_points[i][j][1], sampling_points[i][j][2]);
                // std::cout << "node " << last_node_id << " is created at " << pNewNode->X0() << " " << pNewNode->Y0() << " " << pNewNode->Z0() << std::endl;
                NewNodes.push_back(pNewNode);
            }
        }

        // secondly create elements
        std::size_t last_element_id = FiniteCellAuxilliaryUtility::GetLastElementId(r_model_part);
        // KRATOS_WATCH(last_element_id)
        Element const& rCloneElement = KratosComponents<Element>::Get(sample_element_name);
        Element::NodesArrayType temp_element_nodes;
        ModelPart::ElementsContainerType NewElements;
        const std::string NodeKey("Node");
        std::vector<std::size_t> node;
        int activation_level;

        if (type == 1)
            node.resize(4);
        else if (type == 2)
            node.resize(8);
        else if (type == 2)
            node.resize(9);
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid type", type)

        if (activation_dir == 1) activation_level = -num_division_1;

        for (std::size_t i = 0; i < num_1; ++i)
        {
            if (activation_dir == 2) activation_level = -num_division_2;
            for (std::size_t j = 0; j < num_2; ++j)
            {
                if (type == 1)
                {
                    node[0] = last_node_id_old + i * (num_division_2 + 1) + j + 1;
                    node[2] = last_node_id_old + (i + 1) * (num_division_2 + 1) + j + 1;
                    if (j < num_division_2)
                    {
                        node[1] = last_node_id_old + i * (num_division_2 + 1) + j + 2;
                        node[3] = last_node_id_old + (i + 1) * (num_division_2 + 1) + j + 2;
                    }
                    else
                    {
                        node[1] = last_node_id_old + i * (num_division_2 + 1) + 1;
                        node[3] = last_node_id_old + (i + 1) * (num_division_2 + 1) + 1;
                    }
                    // std::cout << node[0] << " " << node[1] << " " << node[2] << " " << node[3] << std::endl;
                }
                else if (type == 2)
                {
                    // TODO
                }
                else if (type == 3)
                {
                    // TODO
                }

                temp_element_nodes.clear();
                temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node[0], NodeKey).base()));
                temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node[1], NodeKey).base()));
                temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node[3], NodeKey).base()));
                temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node[2], NodeKey).base()));

                Element::Pointer pNewElement = rCloneElement.Create(++last_element_id, temp_element_nodes, pProperties);
                // std::cout << "element " << pNewElement->Id() << " is created" << std::endl;
                pNewElement->Set(ACTIVE, true);
                pNewElement->SetValue(IS_INACTIVE, false);
                if (activation_dir != 0)
                    pNewElement->SetValue(ACTIVATION_LEVEL_var, activation_level);
                NewElements.push_back(pNewElement);

                if (activation_dir == 2) ++activation_level;
            }
            if (activation_dir == 1) ++activation_level;
        }

        for (ModelPart::ElementsContainerType::ptr_iterator it = NewElements.ptr_begin(); it != NewElements.ptr_end(); ++it)
        {
            r_model_part.Elements().push_back(*it);
        }

        r_model_part.Elements().Unique();

        std::cout << NewElements.size() << " " << sample_element_name << " elements are created and added to the model_part" << std::endl;

        return std::make_pair(NewNodes, NewElements);
    }


    /// Create the hex elements based on given points list
    static std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> CreateHexElements(ModelPart& r_model_part,
        const std::vector<std::vector<std::vector<PointType> > >& sampling_points,
        const std::string& sample_element_name,
        const int& type, // if 1: generate H8 elements; 2: H20 elements; 3: H27 elements
        const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir; 3: close on 3rd dir
        const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir; r: activation on 3rd dir
        Properties::Pointer pProperties)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling unimplemented", __FUNCTION__)
    }

    static ModelPart::NodesContainerType::iterator FindKey(ModelPart::NodesContainerType& ThisContainer,
            std::size_t& ThisKey, const std::string& ComponentName)
    {
        ModelPart::NodesContainerType::iterator i_result;
        if((i_result = ThisContainer.find(ThisKey)) == ThisContainer.end())
        {
            std::stringstream buffer;
            buffer << ComponentName << " #" << ThisKey << " is not found.";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }

        return i_result;
    }


    /// Create an element taking the same nodes as the original one, but taking the type of geometry from sample_element_name
    static Element::Pointer CreateParasiteElement(const std::string& sample_element_name,
        std::size_t& lastElementId,
        Element::Pointer pElement, Properties::Pointer pProperties )
    {
        Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);

        // REMARK: when creating the element here, the integration rule is not passed. Instead the default integration rule of this element_type is applied, which is not the same as the original element.
        Element::Pointer pNewElement = r_clone_element.Create(++lastElementId, pElement->pGetGeometry(), pProperties);

        std::cout << "1 element of type " << sample_element_name << " is created" << std::endl;

        return pNewElement;
    }


    /// Create the parasite element taking the same geometry of the original element, but the integration points are given.
    static Element::Pointer CreateParasiteElement(Element::Pointer pElement, // the parent element keep the geometry
        const std::string& sample_element_name,
        const int& RepresentativeIntegrationOrder,
        const GeometryType::IntegrationPointsArrayType& integration_points,
        std::size_t& lastElementId,
        Properties::Pointer pProperties)
    {
        Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);
        GeometryType& r_geom = *(pElement->pGetGeometry());

        GeometryData::IntegrationMethod RepresentativeIntegrationMethod
                = Function<double, double>::GetIntegrationMethod(RepresentativeIntegrationOrder);
        Variable<int>& INTEGRATION_ORDER_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("INTEGRATION_ORDER"));

        // create the new elements from sub-cell
        // here we make a clone of the geometry because we want to assign different geometry data later on
        // this also works with Bezier element, because Bezier geometry has implemented the Create method
        Element::Pointer pNewElement;
        pNewElement = r_clone_element.Create(++lastElementId, r_geom.Create(r_geom.Points()), pProperties);

        FiniteCellGeometryUtility::AssignGeometryData(pNewElement->GetGeometry(), RepresentativeIntegrationMethod, integration_points);
        pNewElement->SetValue(INTEGRATION_ORDER_var, RepresentativeIntegrationOrder);
        pNewElement->Initialize();

        return pNewElement;
    }


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
        return "Finite Cell Mesh Utility";
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
    FiniteCellMeshUtility& operator=(FiniteCellMeshUtility const& rOther);

    /// Copy constructor.
    FiniteCellMeshUtility(FiniteCellMeshUtility const& rOther);


    ///@}

}; // Class FiniteCellMeshUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, FiniteCellMeshUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const FiniteCellMeshUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.


#endif // KRATOS_FINITE_CELL_MESH_UTILITY_H_INCLUDED  defined
