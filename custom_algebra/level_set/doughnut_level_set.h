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
//  Date:            5 Jan 2018
//


#if !defined(KRATOS_DOUGHNUT_LEVEL_SET_H_INCLUDED )
#define  KRATOS_DOUGHNUT_LEVEL_SET_H_INCLUDED



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
 * REF: Massing et al, CutFEM: Discretizing geometry and partial differential equations
*/
class DoughnutLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DoughnutLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(DoughnutLevelSet);

    typedef LevelSet BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DoughnutLevelSet(const double& R, const double& r)
    : BaseType(), mR(R), mr(r)
    {
    }

    /// Destructor.
    virtual ~DoughnutLevelSet() {}


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
        return pow(mR - sqrt(pow(P(0), 2) + pow(P(1), 2)), 2) + pow(P(2), 2) - pow(mr, 2);
    }


    virtual Vector GetGradient(const PointType& P) const
    {
        Vector grad(3);
        double aux = pow(mR - sqrt(pow(P(0), 2) + pow(P(1), 2)), 2);
        grad(0) = 2.0 * aux * ( -P(0) / sqrt(pow(P(0), 2) + pow(P(1), 2)) );
        grad(1) = 2.0 * aux * ( -P(1) / sqrt(pow(P(0), 2) + pow(P(1), 2)) );
        grad(2) = 2.0 * P(2);
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
        return "Doughnut Level Set";
    }

    /// Print information about this object.
//    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << ", R: " << mR << ", r: " << mr;
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


    double mR, mr;


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
    DoughnutLevelSet& operator=(DoughnutLevelSet const& rOther);

    /// Copy constructor.
    DoughnutLevelSet(DoughnutLevelSet const& rOther);


    ///@}

}; // Class DoughnutLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                DoughnutLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const DoughnutLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DOUGHNUT_LEVEL_SET_H_INCLUDED  defined
