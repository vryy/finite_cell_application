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
//  Date:            20 Feb 2017
//


#if !defined(KRATOS_LOAD_FUNCTION_PLATE_WITH_THE_HOLE_H_INCLUDED )
#define  KRATOS_LOAD_FUNCTION_PLATE_WITH_THE_HOLE_H_INCLUDED



// System includes
#include <string>
#include <sstream>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/function/function.h"


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
///@name  LoadFunctionR3RnPlateWithTheHoles
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for the load apply on the sides of the standard plate with the hole problem
*/
template<std::size_t tload_side>
class LoadFunctionR3RnPlateWithTheHole : public FunctionR3Rn
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LoadFunctionR3RnPlateWithTheHole
    KRATOS_CLASS_POINTER_DEFINITION(LoadFunctionR3RnPlateWithTheHole);

    typedef FunctionR3Rn BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LoadFunctionR3RnPlateWithTheHole(const double P, const double r)
    : mP(P), mr(r)
    {}

    /// Copy constructor.
    LoadFunctionR3RnPlateWithTheHole(LoadFunctionR3RnPlateWithTheHole const& rOther)
    : BaseType(rOther), mP(rOther.mP), mr(rOther.mr)
    {}

    /// Destructor.
    virtual ~LoadFunctionR3RnPlateWithTheHole()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual BaseType::Pointer CloneFunction() const
    {
        return BaseType::Pointer(new LoadFunctionR3RnPlateWithTheHole(*this));
    }


    virtual OutputType GetValue(const InputType& P) const
    {
        double r = sqrt(pow(P[0], 2) + pow(P[1], 2));
        double theta = acos(P[0]/r);
//        KRATOS_WATCH(P)

        /// REF: Liu et al, Mesh Free Methods: Moving Beyond the Finite Element Method, example 6.11
        double o_xx = mP * (1.0 - pow(mr/r, 2) * (1.5*cos(2.0*theta) + cos(4.0*theta)) + 1.5*pow(mr/r, 4)*cos(4.0*theta));
        double o_yy = mP * (-pow(mr/r, 2) * (0.5*cos(2.0*theta) - cos(4.0*theta)) - 1.5*pow(mr/r, 4)*cos(4.0*theta));
        double o_xy = mP * (-pow(mr/r, 2) * (0.5*sin(2.0*theta) + sin(4.0*theta)) + 1.5*pow(mr/r, 4)*sin(4.0*theta)); // here the equation in the book is wrong, I have to manually fix this

        Vector Load(2);
        double nx, ny;
        if(tload_side == 0)
        {
            nx = 1.0;
            ny = 0.0;
        }
        else if(tload_side == 1)
        {
            nx = 0.0;
            ny = 1.0;
        }
        Load[0] = o_xx * nx + o_xy * ny;
        Load[1] = o_xy * nx + o_yy * ny;

        return Load;
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
        std::stringstream ss;
        ss << "Load Function for plate with the hole problem: P = " << mP << ", mr = " << mr;
        return ss.str();
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

    double mP, mr;

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
    LoadFunctionR3RnPlateWithTheHole& operator=(LoadFunctionR3RnPlateWithTheHole const& rOther);

    ///@}

}; // Class LoadFunctionR3RnPlateWithTheHole

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream LoadFunctionR3RnPlateWithTheHole
template<std::size_t tload_side>
inline std::istream& operator >> (std::istream& rIStream, LoadFunctionR3RnPlateWithTheHole<tload_side>& rThis)
{}

/// output stream LoadFunctionR3RnPlateWithTheHole
template<std::size_t tload_side>
inline std::ostream& operator << (std::ostream& rOStream, const LoadFunctionR3RnPlateWithTheHole<tload_side>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LOAD_FUNCTION_PLATE_WITH_THE_HOLE_H_INCLUDED  defined
