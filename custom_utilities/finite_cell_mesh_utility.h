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

template<int TDim, int TType>
struct GenerateStructuredMesh_Helper
{
    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;
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


    /// Generate a uniform sampling points on 1D axis
    static void GenerateSampling(std::vector<double>& sampling,
        const double& s_min, const double& s_max, const std::size_t& nsampling);


    /// Generate a concentrated sampling points on 1D axis
    /// The points are sampled as following:
    ///            w1 N1                   w2 N2
    /// s =   --------------- s_min + --------------- s_max
    ///       (w1 N1 + w2 N2)         (w1 N1 + w2 N2)
    /// N1 and N2 are the first order basis functions:
    ///   N1 = 1-t
    ///   N2 = t
    static void GenerateSampling(std::vector<double>& sampling,
        const double& s_min, const double& s_max,
        const double& w1, const double& w2,
        const std::size_t& nsampling);


    /// Generate the points for background structure mesh
    static void GenerateStructuredMesh2D(std::vector<std::vector<PointType> >& sampling_points,
        const int& type,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::size_t>& nsampling);


    /// Generate the points for background structure mesh
    /// sampling vector contains the local coordinate s (0 <= s <= 1) of sampling points in each axis
    static void GenerateStructuredMesh2D(std::vector<std::vector<PointType> >& sampling_points,
        const int& type,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::vector<double> >& sampling);


    /// Generate the points for background structure mesh
    static void GenerateStructuredMesh3D(std::vector<std::vector<std::vector<PointType> > >& sampling_points,
        const int& type,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::size_t>& nsampling);


    /// Generate the points for background structure mesh
    /// sampling vector contains the local coordinate s (0 <= s <= 1) of sampling points in each axis
    static void GenerateStructuredMesh3D(std::vector<std::vector<std::vector<PointType> > >& sampling_points,
        const int& type,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::vector<double> >& sampling);


    /// Create the quad elements based on given points list
    static std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> CreateQuadElements(ModelPart& r_model_part,
        const std::vector<std::vector<PointType> >& sampling_points,
        const std::string& sample_element_name,
        const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
        const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
        const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
        Properties::Pointer pProperties);


    /// Create the hex elements based on given points list
    static std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> CreateHexElements(ModelPart& r_model_part,
        const std::vector<std::vector<std::vector<PointType> > >& sampling_points,
        const std::string& sample_element_name,
        const int& type, // if 1: generate H8 elements; 2: H20 elements; 3: H27 elements
        const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir; 3: close on 3rd dir
        const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir; r: activation on 3rd dir
        Properties::Pointer pProperties);


    /// Find the node by key
    static ModelPart::NodesContainerType::iterator FindKey(ModelPart::NodesContainerType& ThisContainer,
            std::size_t& ThisKey, const std::string& ComponentName);


    /// Create an element taking the same nodes as the original one, but taking the type of geometry from sample_element_name
    static Element::Pointer CreateParasiteElement(const std::string& sample_element_name,
        std::size_t& lastElementId,
        Element::Pointer pElement, Properties::Pointer pProperties );


    /// Create the parasite element taking the same geometry of the original element, but the integration points are given.
    static Element::Pointer CreateParasiteElement(Element::Pointer pElement, // the parent element keep the geometry
        const std::string& sample_element_name,
        const int& RepresentativeIntegrationOrder,
        const GeometryType::IntegrationPointsArrayType& integration_points,
        std::size_t& lastElementId,
        Properties::Pointer pProperties);


    /// Import the nodes from other model_part to this model_part
    static ModelPart::NodesContainerType ImportNodes(ModelPart& rThisModelPart, ModelPart& rOtherModelPart);


    /// Import the elements from list to the this model_part
    static ModelPart::ElementsContainerType ImportElements(ModelPart& rThisModelPart,
        ModelPart::ElementsContainerType& rOtherElements,
        const std::string& sample_element_name, Properties::Pointer pProperties);


    /// Import the conditions from list to the this model_part
    static ModelPart::ConditionsContainerType ImportConditions(ModelPart& rThisModelPart,
        ModelPart::ConditionsContainerType& rOtherConditions,
        const std::string& sample_cond_name, Properties::Pointer pProperties);


    /// Import the elements from list to the other list
    template<class TEntityType, class TEntityContainerType>
    static void ImportEntities(ModelPart& rThisModelPart,
        TEntityContainerType& rThisElements, // rThisModelPart.Elements() or rThisModelPart.Conditions()
        TEntityContainerType& rNewElements, // the added elements to rThisElements
        TEntityContainerType& rOtherElements,
        std::size_t& last_element_id,
        TEntityType const& r_clone_element,
        Properties::Pointer pProperties)
    {
        typename TEntityType::NodesArrayType temp_element_nodes;
        for (typename TEntityContainerType::ptr_iterator it = rOtherElements.ptr_begin(); it != rOtherElements.ptr_end(); ++it)
        {
            typename TEntityType::GeometryType r_geom = (*it)->GetGeometry();

            temp_element_nodes.clear();

            for (std::size_t i = 0; i < r_geom.size(); ++i)
            {
                std::size_t other_node_id = static_cast<std::size_t>(r_geom[i].GetValue(OTHER_NODE_ID));
                temp_element_nodes.push_back(*(FindKey(rThisModelPart.Nodes(), other_node_id, "Node").base()));
            }

            typename TEntityType::Pointer pNewElement;
            pNewElement = r_clone_element.Create(++last_element_id, temp_element_nodes, pProperties);
            rNewElements.push_back(pNewElement);
        }

        for (typename TEntityContainerType::ptr_iterator it = rNewElements.ptr_begin(); it != rNewElements.ptr_end(); ++it)
        {
            rThisElements.push_back(*it);
        }

        rThisElements.Unique();
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

#include "finite_cell_mesh_utility.hpp"

#endif // KRATOS_FINITE_CELL_MESH_UTILITY_H_INCLUDED  defined
