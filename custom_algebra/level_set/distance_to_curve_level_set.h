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
//  Date:            1 Feb 2018
//


#if !defined(KRATOS_DISTANCE_TO_CURVE_LEVEL_SET_H_INCLUDED )
#define  KRATOS_DISTANCE_TO_CURVE_LEVEL_SET_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_algebra/level_set/level_set.h"
#include "custom_utilities/finite_cell_mesh_utility.h"

#define PI 3.1415926535897932384626433832795028841971693

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
 * Level Set that compute the distance to a curve
 */
class DistanceToCurveLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DistanceToCurveLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(DistanceToCurveLevelSet);

    typedef LevelSet BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DistanceToCurveLevelSet(const FunctionR1R3::Pointer pAlignCurve, const double& R)
    : BaseType(), mpCurve(pAlignCurve), mR(R)
    {}

    /// Destructor.
    virtual ~DistanceToCurveLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual std::size_t WorkingSpaceDimension() const
    {
        return 3;
    }


    virtual double GetValue(const PointType& P) const
    {
        double t;
        // return DoNewtonRaphson(P, t);
        return DoBisection(P, t, -1.0, 2.0, 10);
    }


    virtual Vector GetGradient(const PointType& P) const
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "Not implemented")
    }


    /// Generate the sampling points on the level set surface
    std::vector<std::vector<PointType> > GeneratePoints(const std::size_t& nsampling_axial, const std::size_t& nsampling_radial) const
    {
        const double tmin = 0.0;
        const double tmax = 1.0;
        const double tol = 1.0e-10;
        // KRATOS_WATCH(nsampling_axial)
        // KRATOS_WATCH(nsampling_radial)

        std::vector<std::vector<PointType> > results;

        double t, d;
        PointType P, T, N, B, V, Up, Aux;
        Up[0] = 0.0; Up[1] = 0.0; Up[2] = 1.0;
        for (std::size_t i = 0; i < nsampling_axial; ++i)
        {
            // KRATOS_WATCH(t)
            t = tmin + i*(tmax-tmin)/(nsampling_axial-1);

            noalias(P) = mpCurve->GetValue(t);
            // KRATOS_WATCH(P)

            noalias(T) = mpCurve->GetDerivative(0, t);
            T *= 1.0/norm_2(T);
            // KRATOS_WATCH(T)

            /****************
            Using Frenet frame can cause twist
            // noalias(N) = mpCurve->GetSecondDerivative(0, 0, t);
            // N *= 1.0/norm_2(N);
            // KRATOS_WATCH(N)

            // noalias(B) = MathUtils<double>::CrossProduct(T, N);
            // KRATOS_WATCH(B)
            *****************/

            /****************/
            // this method is better, although the up vector must be defined, see nrbsweep.m
            if (i == 0)
            {
                noalias(Aux) = MathUtils<double>::CrossProduct(Up, T);
                noalias(B) = Aux / norm_2(Aux);
            }
            else
            {
                noalias(Aux) = B - inner_prod(B, T)*T;
                noalias(B) = Aux / norm_2(Aux);
            }
            // KRATOS_WATCH(B)

            noalias(N) = MathUtils<double>::CrossProduct(B, T);
            // KRATOS_WATCH(N)
            /****************/

            std::vector<PointType> radial_points(nsampling_radial);
            for (std::size_t j = 0; j < nsampling_radial; ++j)
            {
                d = j*2.0*PI/nsampling_radial;
                noalias(V) = std::cos(d)*N + std::sin(d)*B;

                noalias(radial_points[j]) = P + mR*V;
            }

            results.push_back(radial_points);
        }

        return results;
    }


    /// Create the elements based on sampling points on the surface
    std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> CreateQ4ElementsClosedLoop(ModelPart& r_model_part,
        const std::string& sample_element_name,
        Properties::Pointer pProperties,
        const std::size_t& nsampling_axial,
        const std::size_t& nsampling_radial) const
    {
        // firstly create the sampling points on surface
        std::vector<std::vector<PointType> > sampling_points = this->GeneratePoints(nsampling_axial, nsampling_radial);
        return FiniteCellMeshUtility::CreateQ4ElementsClosedLoop(r_model_part, sampling_points, sample_element_name, pProperties);
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
        return "Distance to Curve Level Set";
    }

    /// Print information about this object.
//    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "curve: " << *mpCurve;
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


    const FunctionR1R3::Pointer mpCurve;
    double mR;


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    // compute the distance based on Newton Raphson, quick but unstable
    double DoNewtonRaphson(const PointType& P, double& t) const
    {
        const double tol = 1.0e-10;
        const int max_iters = 300;

        // firstly compute the projection of point P to the curve
        int iter = 0;
        PointType Proj, dProj, ddProj;
        double rhs, drhs;

        do
        {
            noalias(Proj) = mpCurve->GetValue(t);
            noalias(dProj) = mpCurve->GetDerivative(0, t);
            noalias(ddProj) = mpCurve->GetSecondDerivative(0, 0, t);
            rhs = inner_prod(dProj, P - Proj);
            if (fabs(rhs) < tol) break;
            drhs = inner_prod(ddProj, P - Proj) - inner_prod(dProj, dProj);
            t -= rhs/drhs;
            ++iter;
            if (iter > max_iters) break;
        }
        while (fabs(rhs) > tol);

        if (iter > max_iters)
        {
            KRATOS_WATCH(P)
            KRATOS_WATCH(mpCurve->GetValue(0.0))
            KRATOS_WATCH(mpCurve->GetValue(1.0))
            KRATOS_WATCH(t)
            KRATOS_WATCH(rhs)
            KRATOS_THROW_ERROR(std::logic_error, "The local iteration does not converge", "")
        }
        // KRATOS_WATCH(t)
        // KRATOS_WATCH(P)
        // KRATOS_WATCH(Proj)

        // compute the distance
        return norm_2(P - Proj) - mR;
    }


    // compute the distance based on bisection algorithm, slow but stable
    double DoBisection(const PointType& P, double& t, const double& tmin, const double& tmax, const int& nsampling) const
    {
        const double tol = 1.0e-10;

        // firstly do the sampling
        std::vector<double> f(nsampling);
        PointType Proj, dProj;
        for (std::size_t i = 0; i < nsampling+1; ++i)
        {
            t = tmin + i*(tmax-tmin)/nsampling;
            noalias(Proj) = mpCurve->GetValue(t);
            noalias(dProj) = mpCurve->GetDerivative(0, t);
            f[i] = inner_prod(dProj, P - Proj);
        }

        bool found = false;
        for (std::size_t i = 0; i < nsampling; ++i)
        {
            if (fabs(f[i]) < tol)
            {
                found = true;
                t = tmin + i*(tmax-tmin)/nsampling;
                noalias(Proj) = mpCurve->GetValue(t);
                break;
            }

            if (f[i]*f[i+1] < 0.0)
            {
                // found the segment, do the bisection
                double left = tmin + i*(tmax-tmin)/nsampling;
                double right = tmin + (i+1)*(tmax-tmin)/nsampling;
                double mid;
                double fleft = f[i], fright = f[i+1], fmid;
                while ((right - left) > tol)
                {
                    mid = 0.5*(left + right);

                    noalias(Proj) = mpCurve->GetValue(mid);
                    noalias(dProj) = mpCurve->GetDerivative(0, mid);
                    fmid = inner_prod(dProj, P - Proj);

                    if (fabs(fmid) < tol)
                        break;

                    if (fmid * fleft > 0.0)
                    {
                        left = mid;
                        noalias(Proj) = mpCurve->GetValue(left);
                        noalias(dProj) = mpCurve->GetDerivative(0, left);
                        fleft = inner_prod(dProj, P - Proj);
                    }

                    if (fmid * fright > 0.0)
                    {
                        right = mid;
                        noalias(Proj) = mpCurve->GetValue(right);
                        noalias(dProj) = mpCurve->GetDerivative(0, right);
                        fright = inner_prod(dProj, P - Proj);
                    }
                }

                found = true;
                t = mid;
                break;
            }
        }

        if (!found)
        {
            KRATOS_WATCH(P)
            KRATOS_WATCH(mpCurve->GetValue(0.0))
            KRATOS_WATCH(mpCurve->GetValue(1.0))
            std::cout << "f:";
            for (std::size_t i = 0; i < f.size(); ++i)
                std::cout << " " << f[i];
            std::cout << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "Bisection error: there are no valid segment", "")
        }

        // compute the distance
        return norm_2(P - Proj) - mR;
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
    DistanceToCurveLevelSet& operator=(DistanceToCurveLevelSet const& rOther);

    /// Copy constructor.
    DistanceToCurveLevelSet(DistanceToCurveLevelSet const& rOther);


    ///@}

}; // Class DistanceToCurveLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                DistanceToCurveLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const DistanceToCurveLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SWEEP_LEVEL_SET_H_INCLUDED  defined
