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
//  Date:            10 Feb 2017
//


#if !defined(KRATOS_CIRCULAR_LEVEL_SET_H_INCLUDED )
#define  KRATOS_CIRCULAR_LEVEL_SET_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/level_set/level_set.h"
#include "custom_utilities/finite_cell_mesh_utility.h"


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
*/
class CircularLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CircularLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(CircularLevelSet);

    typedef LevelSet BaseType;

    static constexpr double PI = std::atan(1.0)*4;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CircularLevelSet(const double& cX, const double& cY, const double& R)
    : BaseType(), mcX(cX), mcY(cY), mR(R)
    {}

    /// Copy constructor.
    CircularLevelSet(CircularLevelSet const& rOther)
    : BaseType(rOther), mcX(rOther.mcX), mcY(rOther.mcY), mR(rOther.mR)
    {}

    /// Destructor.
    virtual ~CircularLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual LevelSet::Pointer CloneLevelSet() const
    {
        return LevelSet::Pointer(new CircularLevelSet(*this));
    }


    virtual std::size_t WorkingSpaceDimension() const
    {
        return 2;
    }


    virtual double GetValue(const PointType& P) const
    {
        return pow(P(0) - mcX, 2) + pow(P(1) - mcY, 2) - pow(mR, 2);
    }


    virtual Vector GetGradient(const PointType& P) const
    {
        Vector grad(2);
        grad(0) = 2.0 * (P(0) - mcX);
        grad(1) = 2.0 * (P(1) - mcY);
        return grad;
    }


    /// Generate the sampling points on the level set surface
    std::vector<PointType> GeneratePoints(const double& start_angle, const double& end_angle,
        const std::size_t& nsampling_radial) const
    {
        std::vector<PointType> radial_points(nsampling_radial);
        double small_angle = (end_angle - start_angle) / nsampling_radial;

        PointType V;
        V[2] = 0.0;
        double d;
        for (std::size_t j = 0; j < nsampling_radial; ++j)
        {
            d = start_angle + j*small_angle;
            V[0] = mcX + mR*std::cos(d);
            V[1] = mcY + mR*std::sin(d);
            noalias(radial_points[j]) = V;
        }

        return radial_points;
    }


    /// Generate the sampling points on the level set surface
    std::vector<PointType> GeneratePoints(const std::size_t& nsampling_radial) const
    {
        return GeneratePoints(0.0, 2*PI, nsampling_radial);
    }


    /// Create the elements based on sampling points on the line
    std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> CreateLineElements(ModelPart& r_model_part,
        const std::string& sample_element_name,
        Properties::Pointer pProperties,
        const double& start_angle,
        const double& end_angle,
        const std::size_t& nsampling_radial,
        const bool close = false) const
    {
        // firstly create the sampling points on surface
        std::vector<PointType> sampling_points = this->GeneratePoints(start_angle, end_angle, nsampling_radial);
        int order = 1;
        FiniteCellMeshUtility::MeshInfoType Info = FiniteCellMeshUtility::CreateLineElements(r_model_part, sampling_points, sample_element_name, order, close, pProperties);
        return std::make_pair(std::get<0>(Info), std::get<1>(Info));
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
        return "Circular Level Set";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << this->Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "cX: " << mcX << ", cY: " << mcY << ", R: " << mR;
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


    double mcX, mcY, mR;


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
    CircularLevelSet& operator=(CircularLevelSet const& rOther);

    ///@}

}; // Class CircularLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, CircularLevelSet& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const CircularLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " ";
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CIRCULAR_LEVEL_SET_H_INCLUDED  defined
