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
/** class for auxilliary routines
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


    /// Create the quad-4 elements based on given points list
    static std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> CreateQ4ElementsClosedLoop(ModelPart& r_model_part,
        const std::vector<std::vector<PointType> >& sampling_points,
        const std::string& sample_element_name,
        Properties::Pointer pProperties)
    {
        std::size_t last_node_id = FiniteCellAuxilliaryUtility::GetLastNodeId(r_model_part);
        std::size_t last_node_id_old = last_node_id;
        std::size_t num_division_1 = sampling_points.size() - 1;
        std::size_t num_division_2 = sampling_points[0].size() - 1;
        // KRATOS_WATCH(last_node_id)

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
        std::size_t node_1, node_2, node_3, node_4;
        int activation_level = -num_division_1;
        for (std::size_t i = 0; i < num_division_1; ++i)
        {
            // KRATOS_WATCH(sampling_points[i].size())
            for (std::size_t j = 0; j < num_division_2+1; ++j)
            {
                node_1 = last_node_id_old + i * (num_division_2 + 1) + j + 1;
                node_3 = last_node_id_old + (i + 1) * (num_division_2 + 1) + j + 1;
                if (j < num_division_2)
                {
                    node_2 = last_node_id_old + i * (num_division_2 + 1) + j + 2;
                    node_4 = last_node_id_old + (i + 1) * (num_division_2 + 1) + j + 2;
                }
                else
                {
                    node_2 = last_node_id_old + i * (num_division_2 + 1) + 1;
                    node_4 = last_node_id_old + (i + 1) * (num_division_2 + 1) + 1;
                }
                // std::cout << node_1 << " " << node_2 << " " << node_3 << " " << node_4 << std::endl;

                temp_element_nodes.clear();
                temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node_1, NodeKey).base()));
                temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node_2, NodeKey).base()));
                temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node_4, NodeKey).base()));
                temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node_3, NodeKey).base()));

                Element::Pointer pNewElement = rCloneElement.Create(++last_element_id, temp_element_nodes, pProperties);
                // std::cout << "element " << pNewElement->Id() << " is created" << std::endl;
                pNewElement->Set(ACTIVE, true);
                pNewElement->SetValue(IS_INACTIVE, false);
                pNewElement->SetValue(ACTIVATION_LEVEL_var, activation_level);
                NewElements.push_back(pNewElement);
            }
            ++activation_level;
        }

        for (ModelPart::ElementsContainerType::ptr_iterator it = NewElements.ptr_begin(); it != NewElements.ptr_end(); ++it)
        {
            r_model_part.Elements().push_back(*it);
        }

        r_model_part.Elements().Unique();

        std::cout << NewElements.size() << " " << sample_element_name << " elements are created and added to the model_part" << std::endl;

        return std::make_pair(NewNodes, NewElements);
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
