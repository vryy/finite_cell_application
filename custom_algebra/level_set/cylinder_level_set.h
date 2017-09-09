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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CylinderLevelSet(const double& cX, const double& cY, const double& cZ, const double& dX, const double& dY, const double& dZ, const double& R)
    : BaseType(), mcX(cX), mcY(cY), mcZ(cZ), mR(R)
    {
        double norm_d = sqrt(pow(dX, 2) + pow(dY, 2) + pow(dZ, 2));

        if(norm_d == 0.0)
            KRATOS_THROW_ERROR(std::logic_error, "The director vector can't be null", "")

        mdX = dX / norm_d;
        mdY = dY / norm_d;
        mdZ = dZ / norm_d;
    }

    /// Destructor.
    virtual ~CylinderLevelSet() {}


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

    /// Copy constructor.
    CylinderLevelSet(CylinderLevelSet const& rOther);


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

#endif // KRATOS_CYLINDER_LEVEL_SET_H_INCLUDED  defined
