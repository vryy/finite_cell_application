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


#if !defined(KRATOS_PLANAR_LEVEL_SET_H_INCLUDED )
#define  KRATOS_PLANAR_LEVEL_SET_H_INCLUDED



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
class PlanarLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PlanarLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(PlanarLevelSet);

    typedef LevelSet BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PlanarLevelSet(const double& A, const double& B, const double& C, const double& D)
    : BaseType(), mA(A), mB(B), mC(C), mD(D)
    {}

    /// Copy constructor.
    PlanarLevelSet(PlanarLevelSet const& rOther)
    : BaseType(rOther), mA(rOther.mA), mB(rOther.mB), mC(rOther.mC), mD(rOther.mD)
    {}

    /// Destructor.
    virtual ~PlanarLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual LevelSet::Pointer CloneLevelSet() const
    {
        return LevelSet::Pointer(new PlanarLevelSet(*this));
    }


    virtual std::size_t WorkingSpaceDimension() const
    {
        return 3;
    }


    virtual double GetValue(const PointType& P) const
    {
        return mA*P(0) + mB*P(1) + mC*P(2) + mD;
    }


    virtual Vector GetGradient(const PointType& P) const
    {
        Vector grad(this->WorkingSpaceDimension());
        grad(0) = mA;
        grad(1) = mB;
        grad(2) = mC;
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
        return "Planar Level Set";
    }

    /// Print information about this object.
//    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "A: " << mA << ", B: " << mB << ", C: " << mC << ", D: " << mD;
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


    double mA, mB, mC, mD;


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
    PlanarLevelSet& operator=(PlanarLevelSet const& rOther);


    ///@}

}; // Class PlanarLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                PlanarLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const PlanarLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PLANAR_LEVEL_SET_H_INCLUDED  defined
