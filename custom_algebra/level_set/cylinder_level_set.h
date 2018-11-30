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
//  Date:            14 Feb 2017
//


#if !defined(KRATOS_CYLINDER_LEVEL_SET_H_INCLUDED )
#define  KRATOS_CYLINDER_LEVEL_SET_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/level_set/level_set.h"
#include "custom_utilities/finite_cell_mesh_utility.h"

/*#define PI 3.1415926535897932384626433832795028841971693*/

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
class CylinderLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CylinderLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(CylinderLevelSet);

    typedef LevelSet BaseType;

    static constexpr double PI = std::atan(1.0)*4;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CylinderLevelSet(const double& cX, const double& cY, const double& cZ, const double& dX, const double& dY, const double& dZ, const double& R)
    : BaseType(), mcX(cX), mcY(cY), mcZ(cZ), mR(R)
    {
        mLength = sqrt(pow(dX, 2) + pow(dY, 2) + pow(dZ, 2));

        if(mLength == 0.0)
            KRATOS_THROW_ERROR(std::logic_error, "The director vector can't be null", "")

        mdX = dX / mLength;
        mdY = dY / mLength;
        mdZ = dZ / mLength;
    }

    /// Copy constructor.
    CylinderLevelSet(CylinderLevelSet const& rOther)
    : BaseType(rOther), mcX(rOther.mcX), mcY(rOther.mcY), mcZ(rOther.mcZ)
    , mdX(rOther.mdX), mdY(rOther.mdY), mdZ(rOther.mdZ), mR(rOther.mR)
    {}

    /// Destructor.
    virtual ~CylinderLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual LevelSet::Pointer CloneLevelSet() const
    {
        return LevelSet::Pointer(new CylinderLevelSet(*this));
    }


    virtual std::size_t WorkingSpaceDimension() const
    {
        return 3;
    }


    virtual double GetValue(const PointType& P) const
    {
        double t = (P(0) - mcX) * mdX + (P(1) - mcY) * mdY + (P(2) - mcZ) * mdZ;
        double pX = mcX + t*mdX;
        double pY = mcY + t*mdY;
        double pZ = mcZ + t*mdZ;
//        double pX = (P(0) - mcX) * mdX;
//        double pY = (P(1) - mcY) * mdY;
//        double pZ = (P(2) - mcZ) * mdZ;
        return pow(P(0) - pX, 2) + pow(P(1) - pY, 2) + pow(P(2) - pZ, 2) - pow(mR, 2);
    }


    virtual Vector GetGradient(const PointType& P) const
    {
//        double pX = (P(0) - mcX) * mdX;
//        double pY = (P(1) - mcY) * mdY;
//        double pZ = (P(2) - mcZ) * mdZ;
//        Vector grad(3);
//        grad(0) = 2.0 * (P(0) - pX) * (1.0 - mdX);
//        grad(1) = 2.0 * (P(1) - pY) * (1.0 - mdY);
//        grad(2) = 2.0 * (P(2) - pZ) * (1.0 - mdZ);

        double t = (P(0) - mcX) * mdX + (P(1) - mcY) * mdY + (P(2) - mcZ) * mdZ;
        double pX = mcX + t*mdX;
        double pY = mcY + t*mdY;
        double pZ = mcZ + t*mdZ;
        Vector grad(3);
        grad(0) = 2.0 * (P(0) - pX) * (1.0 - mdX*mdX);
        grad(1) = 2.0 * (P(1) - pY) * (1.0 - mdY*mdY);
        grad(2) = 2.0 * (P(2) - pZ) * (1.0 - mdZ*mdZ);
        return grad;
    }


    /// Generate the sampling points on the level set surface
    std::vector<std::vector<PointType> > GeneratePoints(const std::size_t& nsampling_axial, const std::size_t& nsampling_radial,
        const double& start_angle, const double& end_angle) const
    {
        const double tmin = 0.0;
        const double tmax = 1.0;
        const double tol = 1.0e-10;
        // KRATOS_WATCH(nsampling_axial)
        // KRATOS_WATCH(nsampling_radial)

        std::vector<std::vector<PointType> > results;
        double small_angle = (end_angle - start_angle) / nsampling_radial;

        double t, d;
        PointType P, T, B, V, Up;
        Up[0] = 0.0; Up[1] = 0.0; Up[2] = 1.0;
        for (std::size_t i = 0; i < nsampling_axial; ++i)
        {
            t = tmin + i*(tmax-tmin)/(nsampling_axial-1);

            P[0] = mcX + t*mdX*mLength;
            P[1] = mcY + t*mdY*mLength;
            P[2] = mcZ + t*mdZ*mLength;
            // KRATOS_WATCH(P)

            T[0] = mdX;
            T[1] = mdY;
            T[2] = mdZ;
            // KRATOS_WATCH(T)

            noalias(B) = MathUtils<double>::CrossProduct(Up, T);
            // KRATOS_WATCH(B)

            std::vector<PointType> radial_points(nsampling_radial);
            for (std::size_t j = 0; j < nsampling_radial; ++j)
            {
                d = start_angle + j*small_angle;
                noalias(V) = std::cos(d)*Up + std::sin(d)*B;

                noalias(radial_points[j]) = P + mR*V;
            }

            results.push_back(radial_points);
        }

        return results;
    }


    /// Generate the sampling points on the level set surface
    std::vector<std::vector<PointType> > GeneratePoints(const std::size_t& nsampling_axial, const std::size_t& nsampling_radial) const
    {
        return GeneratePoints(nsampling_axial, nsampling_radial, 0.0, 2*PI);
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
        int order = 1;
        int close_dir = 2;
        int activation_dir = 1;
        FiniteCellMeshUtility::MeshInfoType Info = FiniteCellMeshUtility::CreateQuadElements(r_model_part, sampling_points, sample_element_name, order, close_dir, activation_dir, pProperties);
        return std::make_pair(std::get<0>(Info), std::get<1>(Info));
    }


    /// Create the elements based on sampling points on the surface
    std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> CreateQ4Elements(ModelPart& r_model_part,
        const std::string& sample_element_name,
        Properties::Pointer pProperties,
        const std::size_t& nsampling_axial,
        const std::size_t& nsampling_radial,
        const double& start_radial_angle,
        const double& end_radial_angle) const
    {
        // firstly create the sampling points on surface
        std::vector<std::vector<PointType> > sampling_points = this->GeneratePoints(nsampling_axial, nsampling_radial, start_radial_angle, end_radial_angle);
        int order = 1;
        int close_dir = 0;
        int activation_dir = 1;
        FiniteCellMeshUtility::MeshInfoType Info = FiniteCellMeshUtility::CreateQuadElements(r_model_part, sampling_points, sample_element_name, order, close_dir, activation_dir, pProperties);
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
        return "Cylinder Level Set";
    }

    /// Print information about this object.
//    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "cX: " << mcX << ", cY: " << mcY << ", cZ: " << mcZ
                 << "dX: " << mdX << ", dY: " << mdY << ", dZ: " << mdZ
                 << ", R: " << mR;
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


    double mcX, mcY, mcZ; // point on center line
    double mdX, mdY, mdZ; // director vector
    double mLength;
    double mR;


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
    CylinderLevelSet& operator=(CylinderLevelSet const& rOther);

    ///@}

}; // Class CylinderLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                CylinderLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const CylinderLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#undef PI

#endif // KRATOS_CYLINDER_LEVEL_SET_H_INCLUDED  defined
