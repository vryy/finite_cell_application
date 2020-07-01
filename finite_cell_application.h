//  see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Feb 7, 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_FINITE_CELL_APPLICATION_H_INCLUDED)
#define KRATOS_FINITE_CELL_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_components.h"
#include "includes/kratos_application.h"
#include "custom_conditions/dummy_condition.h"
#include "custom_elements/dummy_element.h"


namespace Kratos
{

    ///@name Kratos Globals
    ///@{

    KRATOS_DEFINE_VARIABLE( Matrix, SUBCELL_WEIGHTS )
    KRATOS_DEFINE_VARIABLE( double, SUBCELL_DOMAIN_SIZE )
    KRATOS_DEFINE_VARIABLE( Vector, SUBCELL_DOMAIN_SIZES )
    KRATOS_DEFINE_VARIABLE( Vector, PHYSICAL_INTEGRATION_POINT_THREED_STRESSES )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(PHYSICAL_INTEGRATION_POINT_LOCAL)
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(PHYSICAL_INTEGRATION_POINT_GLOBAL)
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(PHYSICAL_INTEGRATION_POINT_DISPLACEMENT)
    KRATOS_DEFINE_VARIABLE( int, OTHER_NODE_ID )
    KRATOS_DEFINE_VARIABLE( int, OTHER_ID )
    KRATOS_DEFINE_VARIABLE( int, NUMBER_OF_PHYSICAL_POINTS )

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Enum's
    ///@{

    ///@}
    ///@name Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{

    /// Short class definition.
    /** Detail class definition.
    */
    class KratosFiniteCellApplication : public KratosApplication
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of KratosMultiphaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(KratosFiniteCellApplication);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        KratosFiniteCellApplication();

        /// Destructor.
        virtual ~KratosFiniteCellApplication(){}

        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        virtual void Register();

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
            return "Base application for Finite-Cell-based analysis";
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << Info();
            PrintData(rOStream);
        }

        ///// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const
        {
            rOStream << "in KratosFiniteCellApplication:";
            KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());
            rOStream << "Variables:" << std::endl;
            KratosComponents<VariableData>().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Elements:" << std::endl;
            KratosComponents<Element>().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Conditions:" << std::endl;
            KratosComponents<Condition>().PrintData(rOStream);
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

        const DummyCondition mDummySurfaceCondition2D3N;
        const DummyCondition mDummySurfaceCondition2D6N;
        const DummyCondition mDummySurfaceCondition2D4N;
        const DummyCondition mDummySurfaceCondition2D8N;
        const DummyCondition mDummySurfaceCondition2D9N;

        const DummyCondition mDummySurfaceCondition3D3N;
        const DummyCondition mDummySurfaceCondition3D6N;
        const DummyCondition mDummySurfaceCondition3D4N;
        const DummyCondition mDummySurfaceCondition3D8N;
        const DummyCondition mDummySurfaceCondition3D9N;

        const DummyCondition mDummyConditionPoint2D;
        const DummyCondition mDummyConditionPoint3D;
        const DummyCondition mDummyConditionLine2N;
        const DummyCondition mDummyConditionLine3N;
        const DummyCondition mDummyCondition2D3N;
        const DummyCondition mDummyCondition2D4N;
        const DummyCondition mDummyCondition2D6N;
        const DummyCondition mDummyCondition2D8N;
        const DummyCondition mDummyCondition2D9N;
        const DummyCondition mDummyCondition3D4N;
        const DummyCondition mDummyCondition3D10N;
        const DummyCondition mDummyCondition3D8N;
        const DummyCondition mDummyCondition3D20N;
        const DummyCondition mDummyCondition3D27N;
        const DummyCondition mDummyCondition3D6N;
        const DummyCondition mDummyCondition3D15N;

        const DummyElement mDummySurfaceElement2D3N;
        const DummyElement mDummySurfaceElement2D6N;
        const DummyElement mDummySurfaceElement2D4N;
        const DummyElement mDummySurfaceElement2D8N;
        const DummyElement mDummySurfaceElement2D9N;

        const DummyElement mDummySurfaceElement3D3N;
        const DummyElement mDummySurfaceElement3D6N;
        const DummyElement mDummySurfaceElement3D4N;
        const DummyElement mDummySurfaceElement3D8N;
        const DummyElement mDummySurfaceElement3D9N;

        const DummyElement mDummyVolumeElement3D4N;
        const DummyElement mDummyVolumeElement3D10N;
        const DummyElement mDummyVolumeElement3D8N;
        const DummyElement mDummyVolumeElement3D20N;
        const DummyElement mDummyVolumeElement3D27N;

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
        KratosFiniteCellApplication& operator=(KratosFiniteCellApplication const& rOther);

        /// Copy constructor.
        KratosFiniteCellApplication(KratosFiniteCellApplication const& rOther);


        ///@}

    }; // Class KratosFiniteCellApplication

    ///@}


    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}


} // namespace Kratos

#endif // KRATOS_FINITE_CELL_APPLICATION_H_INCLUDED defined

