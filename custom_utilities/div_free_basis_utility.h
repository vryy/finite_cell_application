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
//  Date:            11 Feb 2017
//


#if !defined(KRATOS_DIV_FREE_BASIS_UTILITY_H_INCLUDED )
#define  KRATOS_DIV_FREE_BASIS_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/quadrature_utility.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "geometries/line_2d_2.h"
#include "custom_geometries/finite_cell_geometry.h"
#include "custom_linear_solvers/least_square_lapack_solver.h"
#include "brep_application/custom_algebra/level_set/level_set.h"
#include "custom_utilities/finite_cell_geometry_utility.h"


#define ESTIMATE_RCOND


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
/** Detail class definition.
REF: M. Gehrke, Implementierung und Untersuchung numerischer Quadraturverfahren für diskontinuierliche Integranden in der Finiten-Zell-Methode, MSc thesis
*/
class DivFreeBasisUtility : public QuadratureUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DivFreeBasisUtility
    KRATOS_CLASS_POINTER_DEFINITION(DivFreeBasisUtility);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DivFreeBasisUtility()
    {}

    /// Destructor.
    virtual ~DivFreeBasisUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    template<std::size_t TDimension, std::size_t TDegree>
    static std::size_t BasisNumber()
    {
        if(TDimension == 2)
        {
            switch(TDegree)
            {
                case 0: return 2;
                case 1: return 5;
                case 2: return 9;
                case 3: return 14;
                default: return 0;
            }
        }
        else if(TDimension == 3)
        {
            // TODO
            return 0;
        }
    }


    template<std::size_t TDimension, std::size_t TDegree>
    static Matrix GetValues(const CoordinatesArrayType& rPoint)
    {
        Matrix Results;

        if(TDimension == 2)
        {
            std::size_t num_basis = BasisNumber<TDimension, TDegree>();
            if(num_basis == 0)
                return ZeroMatrix(0, TDimension);

            Results.resize(num_basis, TDimension, false);

            if(TDegree >= 0)
            {
                Results(0, 0) = 1.0;
                Results(0, 1) = 0.0;
                Results(1, 0) = 0.0;
                Results(1, 1) = 1.0;
            }

            if(TDegree >= 1)
            {
                Results(2, 0) = rPoint[1];
                Results(2, 1) = 0.0;
                Results(3, 0) = rPoint[0];
                Results(3, 1) = -rPoint[1];
                Results(4, 0) = 0.0;
                Results(4, 1) = rPoint[0];
            }

            if(TDegree >= 2)
            {
                Results(5, 0) = pow(rPoint[1], 2);
                Results(5, 1) = 0.0;
                Results(6, 0) = rPoint[0] * rPoint[1];
                Results(6, 1) = -0.5 * pow(rPoint[1], 2);
                Results(7, 0) = pow(rPoint[0], 2);
                Results(7, 1) = -2.0 * rPoint[0] * rPoint[1];
                Results(8, 0) = 0.0;
                Results(8, 1) = pow(rPoint[0], 2);
            }

            if(TDegree >= 3)
            {
                Results(9, 0) = pow(rPoint[1], 3);
                Results(9, 1) = 0.0;
                Results(10, 0) = rPoint[0]*pow(rPoint[1], 2);
                Results(10, 1) = -pow(rPoint[1], 3)/3;
                Results(11, 0) = pow(rPoint[0], 2)*rPoint[1];
                Results(11, 1) = -rPoint[0]*pow(rPoint[1], 2);
                Results(12, 0) = pow(rPoint[0], 3);
                Results(12, 1) = -3*pow(rPoint[0], 2)*rPoint[1];
                Results(13, 0) = 0.0;
                Results(13, 1) = pow(rPoint[0], 3);
            }
        }
        else if(TDimension == 3)
        {
            KRATOS_THROW_ERROR(std::logic_error, "Not implemented", "")
        }

        return Results;
    }

    /// Cut a geometry by a level set and construct the quadrature
    /// We assume that the geometry is a polytope, that said, in 2D the geometry is a triangle or quadrilateral. If it's second order, then the edges are still defined by the 3 (or 4) vertices.
    void AssignQuadrature2D(GeometryType& r_geom, const LevelSet& r_level_set,
            const unsigned int& integration_order, const unsigned int& div_free_order)
    {
        // pre-check
        if(r_geom.WorkingSpaceDimension() == 3)
            KRATOS_THROW_ERROR(std::logic_error, "This subroutine does not support 3D operation yet", "")

        if(    r_geom.GetGeometryType() == GeometryData::Kratos_Triangle2D6
            || r_geom.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8
            || r_geom.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9 )
        {
            KRATOS_THROW_ERROR(std::logic_error, "This subroutine does not support second order yet", "")
        }

        /* check if the cell is cut by the level set */

        std::size_t num_corners, num_edges;
        if(    r_geom.GetGeometryType() == GeometryData::Kratos_Triangle2D3
            || r_geom.GetGeometryType() == GeometryData::Kratos_Triangle2D6 )
        {
            num_corners = 3;
            num_edges = 3;
        }
        else if(r_geom.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4
             || r_geom.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8
             || r_geom.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9 )
        {
            num_corners = 4;
            num_edges = 4;
        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "This geometry is not supported", r_geom.GetGeometryType())
        }

        // check the cut status with the level set
        int configuration = 0; // check the cut status in the initial configuration
        int stat = r_level_set.CutStatus(r_geom, configuration);
//KRATOS_WATCH(stat)

        if(stat == -1) // the cell is cut
        {
            // obtain the quadrature rule, necessary to construct a cut-cell quadrature
            GeometryData::IntegrationMethod ThisIntegrationMethod
                    = LevelSet::GetIntegrationMethod(integration_order);

            const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints( ThisIntegrationMethod );

            /* compute matrix A */
            std::size_t num_basis;
            if(div_free_order == 0)      num_basis = BasisNumber<2, 0>();
            else if(div_free_order == 1) num_basis = BasisNumber<2, 1>();
            else if(div_free_order == 2) num_basis = BasisNumber<2, 2>();
            else if(div_free_order == 3) num_basis = BasisNumber<2, 3>();
            else if(div_free_order == 4) num_basis = BasisNumber<2, 4>();
            else num_basis = 0;

            Matrix fp(num_basis, 2);
            Matrix MA(num_basis, integration_points.size());

            CoordinatesArrayType global_coords;
            for(std::size_t i = 0; i < integration_points.size(); ++i)
            {
                // compute the global coordinates of the integration point
                r_geom.GlobalCoordinates(global_coords, integration_points[i]);

                // compute the normalized gradient of the level set at the integration point
                Vector grad_phi = r_level_set.GetGradient(integration_points[i]);
                Vector ni = grad_phi / norm_2(grad_phi);

                // compute the divergence free basis function
                if(div_free_order == 0)      noalias(fp) = GetValues<2, 0>(integration_points[i]);
                else if(div_free_order == 1) noalias(fp) = GetValues<2, 1>(integration_points[i]);
                else if(div_free_order == 2) noalias(fp) = GetValues<2, 2>(integration_points[i]);
                else if(div_free_order == 3) noalias(fp) = GetValues<2, 3>(integration_points[i]);
                else if(div_free_order == 4) noalias(fp) = GetValues<2, 4>(integration_points[i]);

                for(std::size_t b = 0; b < num_basis; ++b)
                {
                    MA(b, i) = fp(b, 0) * ni(0) + fp(b, 1) * ni(1);
                }
            }
//KRATOS_WATCH(MA)
            /* compute vector b */

            Vector Mb(num_basis);
            noalias(Mb) = ZeroVector(num_basis);

            // iterate through every edges of the cell, to compute the new quadrature weight
            for(std::size_t e = 0; e < num_edges; ++e)
            {
                std::vector<std::size_t> edge = this->GetEdge(e, r_geom.GetGeometryType());

                // compute the level set for two end of the edge
                double phi1 = r_level_set.GetValue(r_geom[edge[0]]);
                double phi2 = r_level_set.GetValue(r_geom[edge[1]]);

                /* here we identify 2 points of the edge using for integration */
                PointType P1, P2;
                bool need_integrate = false;

                // check if the edge is cut
                if(phi1 * phi2 < 0.0)
                {
                    // find the intersection point
                    r_level_set.Bisect(P1, r_geom[edge[0]], r_geom[edge[1]], r_level_set.Tolerance());

                    // the other point
                    P2 = (phi1 < 0.0 ? r_geom[edge[0]] : r_geom[edge[1]]);

                    need_integrate = true;
                }

                // check if the edge is inside the level set domain
                if(phi1 < 0.0 && phi2 < 0.0)
                {
                    P1 = r_geom[edge[0]];
                    P2 = r_geom[edge[1]];
                    need_integrate = true;
                }

                if(!need_integrate)
                    continue;

//KRATOS_WATCH(P1)
//KRATOS_WATCH(P2)

                // REMARK: we exclude the case that the level set cut the edge in 2 points. This is non trivial and cannot be handled by using this scheme
                // find the local coordinates of P1
                CoordinatesArrayType lP1, lP2;
                r_geom.PointLocalCoordinates(lP1, P1);
                r_geom.PointLocalCoordinates(lP2, P2);

                // create a quadrature rule on the line
                GeometryType::Pointer pLineGeom = GeometryType::Pointer( new Line2D2<NodeType>(Element::GeometryType::PointsArrayType( 2, Node<3>())) );

                const GeometryType::IntegrationPointsArrayType& line_integration_points = pLineGeom->IntegrationPoints( GeometryData::GI_GAUSS_4 ); // TODO parameterize this

                // integrate along P1-P2
                for(std::size_t i = 0; i < line_integration_points.size(); ++i)
                {
                    // compute integration point
                    double s1 = pLineGeom->ShapeFunctionValue(0, line_integration_points[i]);
                    double s2 = pLineGeom->ShapeFunctionValue(1, line_integration_points[i]);
                    CoordinatesArrayType lP = s1 * lP1 + s2 * lP2;
//KRATOS_WATCH(lP)
                    // compute fp
                    if(div_free_order == 0)      noalias(fp) = GetValues<2, 0>(lP);
                    else if(div_free_order == 1) noalias(fp) = GetValues<2, 1>(lP);
                    else if(div_free_order == 2) noalias(fp) = GetValues<2, 2>(lP);
                    else if(div_free_order == 3) noalias(fp) = GetValues<2, 3>(lP);
                    else if(div_free_order == 4) noalias(fp) = GetValues<2, 4>(lP);
//KRATOS_WATCH(fp)
                    // compute normal
                    double edge_length = norm_2(r_geom[edge[1]] - r_geom[edge[0]]);
                    double n0 = (r_geom[edge[1]][1] - r_geom[edge[0]][1]) / edge_length;
                    double n1 = -(r_geom[edge[1]][0] - r_geom[edge[0]][0]) / edge_length;

                    // compute dL
                    double dL = 0.5 * norm_2(P2 - P1);

                    // contribution to b
                    for(std::size_t b = 0; b < num_basis; ++b)
                    {
                        Mb(b) -= (fp(b, 0)*n0 + fp(b, 1)*n1) * dL * line_integration_points[i].Weight();
                    }
                }
            }
//KRATOS_WATCH(Mb)
            #ifdef ESTIMATE_RCOND
            double rcond = LeastSquareLAPACKSolver::EstimateRCond(MA);
            KRATOS_WATCH(rcond)
            #endif

            /* solve the non-square linear system by least square. NOTE: it can be very ill-conditioned */
            Vector Mw;
            LeastSquareLAPACKSolver::SolveDGELSY(MA, Mw, Mb);
//            LeastSquareLAPACKSolver::SolveDGELSS(MA, Mw, Mb);
            KRATOS_WATCH(Mw)
            KRATOS_WATCH(sum(Mw))


            /* create new quadrature and assign to the geometry */
            FiniteCellGeometryUtility::AssignGeometryData(r_geom, ThisIntegrationMethod, Mw);
        }
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
        return "Divergence-Free Basis";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {}


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


    std::vector<std::size_t> GetEdge(const std::size_t& e, const GeometryData::KratosGeometryType& Type)
    {
        std::vector<std::size_t> result;

        if(Type == GeometryData::Kratos_Triangle2D3 || Type == GeometryData::Kratos_Triangle3D3)
        {
            result.resize(2);
            if(e == 0)
            {
                result[0] = 0;
                result[1] = 1;
            }
            else if(e == 1)
            {
                result[0] = 1;
                result[1] = 2;
            }
            else if(e == 2)
            {
                result[0] = 2;
                result[1] = 0;
            }
        }
        else if(Type == GeometryData::Kratos_Triangle2D6 || Type == GeometryData::Kratos_Triangle3D6)
        {
            result.resize(3);
            if(e == 0)
            {
                result[0] = 0;
                result[1] = 1;
                result[2] = 3;
            }
            else if(e == 1)
            {
                result[0] = 1;
                result[1] = 2;
                result[2] = 4;
            }
            else if(e == 2)
            {
                result[0] = 2;
                result[1] = 0;
                result[2] = 5;
            }
        }
        else if(Type == GeometryData::Kratos_Quadrilateral2D4 || Type == GeometryData::Kratos_Quadrilateral3D4)
        {
            result.resize(2);
            if(e == 0)
            {
                result[0] = 0;
                result[1] = 1;
            }
            else if(e == 1)
            {
                result[0] = 1;
                result[1] = 2;
            }
            else if(e == 2)
            {
                result[0] = 2;
                result[1] = 3;
            }
            else if(e == 3)
            {
                result[0] = 3;
                result[1] = 0;
            }
        }
        else if(Type == GeometryData::Kratos_Quadrilateral2D8 || Type == GeometryData::Kratos_Quadrilateral3D8
             || Type == GeometryData::Kratos_Quadrilateral2D9 || Type == GeometryData::Kratos_Quadrilateral3D9)
        {
            result.resize(3);
            if(e == 0)
            {
                result[0] = 0;
                result[1] = 1;
                result[2] = 4;
            }
            else if(e == 1)
            {
                result[0] = 1;
                result[1] = 2;
                result[2] = 5;
            }
            else if(e == 2)
            {
                result[0] = 2;
                result[1] = 3;
                result[2] = 6;
            }
            else if(e == 3)
            {
                result[0] = 3;
                result[1] = 0;
                result[2] = 7;
            }
        }

        return result;
    }


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
    DivFreeBasisUtility& operator=(DivFreeBasisUtility const& rOther);

    /// Copy constructor.
    DivFreeBasisUtility(DivFreeBasisUtility const& rOther);


    ///@}

}; // Class DivFreeBasisUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                DivFreeBasisUtility& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const DivFreeBasisUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#undef ESTIMATE_RCOND

#endif // KRATOS_DIV_FREE_BASIS_UTILITY_H_INCLUDED  defined
