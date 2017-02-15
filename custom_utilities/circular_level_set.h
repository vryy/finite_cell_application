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
#include "custom_utilities/level_set.h"


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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CircularLevelSet(const double& cX, const double& cY, const double& R)
    : BaseType(), mcX(cX), mcY(cY), mR(R)
    {}

    /// Destructor.
    virtual ~CircularLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


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
//    virtual void PrintInfo(std::ostream& rOStream) const;

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

    /// Copy constructor.
    CircularLevelSet(CircularLevelSet const& rOther);


    ///@}

}; // Class CircularLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                CircularLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const CircularLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CIRCULAR_LEVEL_SET_H_INCLUDED  defined
