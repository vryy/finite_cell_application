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


#if !defined(KRATOS_LINEAR_LEVEL_SET_H_INCLUDED )
#define  KRATOS_LINEAR_LEVEL_SET_H_INCLUDED



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
class LinearLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(LinearLevelSet);

    typedef LevelSet BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LinearLevelSet(const double& A, const double& B, const double& C)
    : BaseType(), mA(A), mB(B), mC(C)
    {}

    /// Copy constructor.
    LinearLevelSet(LinearLevelSet const& rOther)
    : BaseType(rOther), mA(rOther.mA), mB(rOther.mB), mC(rOther.mC)
    {}

    /// Destructor.
    virtual ~LinearLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual LevelSet::Pointer CloneLevelSet() const
    {
        return LevelSet::Pointer(new LinearLevelSet(*this));
    }


    virtual std::size_t WorkingSpaceDimension() const
    {
        return 2;
    }


    virtual double GetValue(const PointType& P) const
    {
        return mA*P(0) + mB*P(1) + mC;
    }


    virtual Vector GetGradient(const PointType& P) const
    {
        Vector grad(2);
        grad(0) = mA;
        grad(1) = mB;
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
        return "Linear Level Set";
    }

    /// Print information about this object.
//    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "A: " << mA << ", B: " << mB << ", C: " << mC;
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


    double mA, mB, mC;


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
    LinearLevelSet& operator=(LinearLevelSet const& rOther);

    ///@}

}; // Class LinearLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                LinearLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const LinearLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LINEAR_LEVEL_SET_H_INCLUDED  defined
