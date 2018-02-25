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
//  Date:            15 Feb 2017
//


#if !defined(KRATOS_PRODUCT_LEVEL_SET_H_INCLUDED )
#define  KRATOS_PRODUCT_LEVEL_SET_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
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
/** Class for product of two level sets
*/
class ProductLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ProductLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(ProductLevelSet);

    typedef LevelSet BaseType;

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ProductLevelSet(const BaseType::Pointer p_level_set_1, const BaseType::Pointer p_level_set_2)
    : BaseType(), mp_level_set_1(p_level_set_1), mp_level_set_2(p_level_set_2)
    {}

    /// Copy constructor.
    ProductLevelSet(ProductLevelSet const& rOther)
    : BaseType(rOther)
    , mp_level_set_1(rOther.mp_level_set_1->CloneLevelSet())
    , mp_level_set_2(rOther.mp_level_set_2->CloneLevelSet())
    {}

    /// Destructor.
    virtual ~ProductLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual LevelSet::Pointer CloneLevelSet() const
    {
        return LevelSet::Pointer(new ProductLevelSet(*this));
    }


    virtual std::size_t WorkingSpaceDimension() const
    {
        return mp_level_set_1->WorkingSpaceDimension();
    }


    virtual double GetValue(const PointType& P) const
    {
        return mp_level_set_1->GetValue(P) * mp_level_set_2->GetValue(P);
    }


    virtual Vector GetGradient(const PointType& P) const
    {
        double phi_1 = mp_level_set_1->GetValue(P);
        double phi_2 = mp_level_set_2->GetValue(P);

        Vector grad_1 = mp_level_set_1->GetGradient(P);
        Vector grad_2 = mp_level_set_2->GetGradient(P);

        const std::size_t dim = WorkingSpaceDimension();

        Vector result(dim);

        for(std::size_t i = 0; i < dim; ++i)
            result(i) = grad_1(i) * phi_2 + phi_1 * grad_2(i);

        return result;
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
        return "Product Level Set";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "LS1: " << *mp_level_set_1 << ", LS2: " << *mp_level_set_2;
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


    const BaseType::Pointer mp_level_set_1;
    const BaseType::Pointer mp_level_set_2;


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
    ProductLevelSet& operator=(ProductLevelSet const& rOther);

    ///@}

}; // Class ProductLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, ProductLevelSet& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const ProductLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " ";
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PRODUCT_LEVEL_SET_H_INCLUDED  defined
