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


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "includes/element.h"
#include "includes/model_part.h"
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
            const double& max_weight) const
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
