//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 1 Mar 17 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_DUMMY_POINT_CONDITION_H_INCLUDED )
#define  KRATOS_DUMMY_POINT_CONDITION_H_INCLUDED


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
 */
class DummyPointCondition : public Condition
{
    public:
        // Counted pointer of DummyPointCondition
        KRATOS_CLASS_POINTER_DEFINITION(DummyPointCondition);
        
        /** 
         * Default constructor.
         */
        DummyPointCondition();
        DummyPointCondition( IndexType NewId, GeometryType::Pointer pGeometry);
        DummyPointCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        /**
         * Destructor.
         */
        virtual ~DummyPointCondition();
  
        /**
         * Operations.
         */

        virtual Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const;

        virtual Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
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

        void GetDofList( DofsVectorType& ConditionalDofList,
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
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Condition )
        }

        virtual void load ( Serializer& rSerializer )
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Condition )
        }

        void CalculateAll( MatrixType& rLeftHandSideMatrix, 
                           VectorType& rRightHandSideVector,
                           ProcessInfo& rCurrentProcessInfo,
                           bool CalculateStiffnessMatrixFlag,
                           bool CalculateResidualVectorFlag);

}; // Class DummyPointCondition 

}  // namespace Kratos.
  

#endif // KRATOS_DUMMY_POINT_CONDITION_H_INCLUDED defined 

