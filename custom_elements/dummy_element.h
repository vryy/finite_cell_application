//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 28 Mar 17 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_DUMMY_ELEMENT_H_INCLUDED )
#define  KRATOS_DUMMY_ELEMENT_H_INCLUDED


// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{

/**
 * A zero contribution element to support for assigning condition
 */
class DummyElement : public Element
{
    public:
        // Counted pointer of DummyElement
        KRATOS_CLASS_POINTER_DEFINITION(DummyElement);
        
        /** 
         * Default constructor.
         */
        DummyElement();
        DummyElement( IndexType NewId, GeometryType::Pointer pGeometry);
        DummyElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        /**
         * Destructor.
         */
        virtual ~DummyElement();
  
        /**
         * Operations.
         */

        virtual Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const;

        virtual Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                                PropertiesType::Pointer pProperties) const;

        /**
         * Calculates the local system contributions for this contact element
         */
        void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                   VectorType& rRightHandSideVector, 
                                   ProcessInfo& rCurrentProcessInfo);

        void CalculateRightHandSide( VectorType& rRightHandSideVector, 
                                     ProcessInfo& rCurrentProcessInfo);

        void EquationIdVector( EquationIdVectorType& rResult, 
                               ProcessInfo& rCurrentProcessInfo);

        void GetDofList( DofsVectorType& ElementalDofList,
                         ProcessInfo& CurrentProcessInfo);

        void Initialize();
        /**
         * Turn back information as a string.
         * (DEACTIVATED)
         */
        //std::string Info();
  
        /**
         * Print information about this object.
         * (DEACTIVATED)
         */
        //virtual void PrintInfo(std::ostream& rOStream) const;

        /**
         * Print object's data.
         * (DEACTIVATED)
         */
        //virtual void PrintData(std::ostream& rOStream) const;
  
    protected:
    
    
    private:

        friend class Serializer;

        virtual void save ( Serializer& rSerializer ) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Element )
        }

        virtual void load ( Serializer& rSerializer )
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Element )
        }

        void CalculateAll( MatrixType& rLeftHandSideMatrix, 
                           VectorType& rRightHandSideVector,
                           ProcessInfo& rCurrentProcessInfo,
                           bool CalculateStiffnessMatrixFlag,
                           bool CalculateResidualVectorFlag);

}; // Class DummyElement 

}  // namespace Kratos.
  

#endif // KRATOS_DUMMY_ELEMENT_H_INCLUDED defined 

