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
//  Date:            29 Mar 2017
//


#if !defined(KRATOS_FINITE_CELL_AUXILLIARY_UTILITY_H_INCLUDED )
#define  KRATOS_FINITE_CELL_AUXILLIARY_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <fstream>


// External includes
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/foreach.hpp>
#include <boost/progress.hpp>


// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
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
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** class for auxilliary routines
*/
class FiniteCellAuxilliaryUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FiniteCellAuxilliaryUtility
    KRATOS_CLASS_POINTER_DEFINITION(FiniteCellAuxilliaryUtility);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FiniteCellAuxilliaryUtility() {}

    /// Destructor.
    virtual ~FiniteCellAuxilliaryUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Extract the element from the list of ids
    static ModelPart::ElementsContainerType GetElements(ModelPart& r_model_part, const std::set<std::size_t>& element_list)
    {
        ModelPart::ElementsContainerType pElements;

        GetElements(pElements, r_model_part, element_list);

        return pElements;
    }

    /// Extract the element from the list of ids
    static void GetElements(ModelPart::ElementsContainerType& rpElements,
        ModelPart& r_model_part, const std::set<std::size_t>& element_list)
    {
        rpElements.clear();

        boost::progress_display show_progress( element_list.size() );

        for(std::set<std::size_t>::iterator it = element_list.begin(); it != element_list.end(); ++it)
        {
            ModelPart::ElementsContainerType::iterator it_elem = r_model_part.Elements().find(*it);
            if(it_elem != r_model_part.Elements().end())
                rpElements.push_back(r_model_part.Elements()(*it));
            ++show_progress;
        }
    }

    /// Add an element to the container
    static void AddElement(ModelPart::ElementsContainerType& rpElements, Element::Pointer pElement)
    {
        rpElements.push_back(pElement);
    }

    /// Remove an element to the container
    static void RemoveElement(ModelPart::ElementsContainerType& rpElements, Element::Pointer pElement)
    {
        rpElements.erase(pElement->Id());
    }

    /// "Auxiliary method"
    /// Clean the set of conditions from the model_part (for the most used case, that is the linking conditions)
    void Clean(ModelPart& r_model_part, ModelPart::ConditionsContainerType& rConditions, const int& echo_level) const
    {
        // r_model_part.Conditions().erase(rConditions.begin(), rConditions.end());
        unsigned int cnt = 0, first_id = 0, last_id = 0;
        for(typename ModelPart::ConditionsContainerType::ptr_iterator it = rConditions.ptr_begin();
                it != rConditions.ptr_end(); ++it)
        {
            if(cnt == 0)
                first_id = (*it)->Id();
            last_id = (*it)->Id();

            r_model_part.RemoveCondition(*it);

            ++cnt;
        }

        if(echo_level > 0)
            std::cout << cnt << " conditions: " << first_id << "->" << last_id << " is removed from model_part" << std::endl;
    }


    static Condition::Pointer CreateCondition(ModelPart& r_model_part, const std::string& sample_cond_name,
        const std::size_t& Id, Properties::Pointer pProperties, const std::vector<std::size_t>& node_ids)
    {
        if(!KratosComponents<Condition>::Has(sample_cond_name))
            KRATOS_THROW_ERROR(std::logic_error, sample_cond_name, "is not registered to the KRATOS kernel")
        Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_cond_name);
//        KRATOS_WATCH(r_clone_condition)

        // create the points array
        GeometryType::PointsArrayType Points;
        for(std::size_t i = 0; i < node_ids.size(); ++i)
            Points.push_back(r_model_part.pGetNode(node_ids[i]));
//        GeometryType::Pointer pTempGeometry = r_clone_condition.GetGeometry().Create(Points);

        // create new condition
        Condition::Pointer pNewCond = r_clone_condition.Create(Id, Points, pProperties);
        pNewCond->SetValue(IS_INACTIVE, false);
        pNewCond->Set(ACTIVE, true);
        r_model_part.Conditions().push_back(pNewCond);

        return pNewCond;
    }


    /// Get the last node id of the model part
    static std::size_t GetLastNodeId(ModelPart& r_model_part)
    {
        std::size_t lastNodeId = 0;
        for(typename ModelPart::NodesContainerType::iterator it = r_model_part.Nodes().begin();
                it != r_model_part.Nodes().end(); ++it)
        {
            if(it->Id() > lastNodeId)
                lastNodeId = it->Id();
        }

        return lastNodeId;
    }


    /// Get the last element id of the model part
    static std::size_t GetLastElementId(ModelPart& r_model_part)
    {
        std::size_t lastElementId = 0;
        for(typename ModelPart::ElementsContainerType::ptr_iterator it = r_model_part.Elements().ptr_begin();
                it != r_model_part.Elements().ptr_end(); ++it)
        {
            if((*it)->Id() > lastElementId)
                lastElementId = (*it)->Id();
        }

        return lastElementId;
    }


    /// Get the last condition id of the model part
    static std::size_t GetLastConditionId(ModelPart& r_model_part)
    {
        std::size_t lastCondId = 0;
        for(typename ModelPart::ConditionsContainerType::ptr_iterator it = r_model_part.Conditions().ptr_begin();
                it != r_model_part.Conditions().ptr_end(); ++it)
        {
            if((*it)->Id() > lastCondId)
                lastCondId = (*it)->Id();
        }

        return lastCondId;
    }


    /// Get the last properties id of the model_part
    static std::size_t GetLastPropertiesId(ModelPart& r_model_part)
    {
        std::size_t lastPropId = 0;
        for(typename ModelPart::PropertiesContainerType::ptr_iterator it = r_model_part.rProperties().ptr_begin();
                it != r_model_part.rProperties().ptr_end(); ++it)
        {
            if((*it)->Id() > lastPropId)
                lastPropId = (*it)->Id();
        }

        return lastPropId;
    }


    /// Create a partitioning for multithreaded parallelization
    template<class TValueContainerType>
    static inline void CreatePartition(const unsigned int& number_of_threads,
            const std::size_t& number_of_elements, TValueContainerType& partitions)
    {
        partitions.resize(number_of_threads+1);
        int partition_size = number_of_elements / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_elements;
        for(unsigned int i = 1; i < number_of_threads; ++i)
            partitions[i] = partitions[i-1] + partition_size ;
    }


    /// Multithreaded refine the quadtree(-subcell) using a brep
    template<class TTreeType, class TBRepType>
    static void MultithreadedRefineBy(std::vector<typename TTreeType::Pointer>& r_trees, const TBRepType& r_brep)
    {
        unsigned int number_of_threads = 1;
#ifdef _OPENMP
        number_of_threads = omp_get_max_threads();
#endif
        std::cout << "number_of_threads for MultithreadedRefineBy: " << number_of_threads << std::endl;
        std::vector<unsigned int> tree_partition;
        CreatePartition(number_of_threads, r_trees.size(), tree_partition);

        boost::progress_display show_progress( r_trees.size() );

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for(int k = 0; k < number_of_threads; ++k)
        {
            typename std::vector<typename TTreeType::Pointer>::iterator it_first_tree = r_trees.begin() + tree_partition[k];
            typename std::vector<typename TTreeType::Pointer>::iterator it_last_tree = r_trees.begin() + tree_partition[k+1];

            for(typename std::vector<typename TTreeType::Pointer>::iterator it = it_first_tree; it != it_last_tree; ++it)
            {
                (*it)->RefineBy(r_brep);
                ++show_progress;
            }
        }
        std::cout << "\nMultithreadedRefineBy completed for " << r_trees.size() << " trees" << std::endl;
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

    void Print(GeometryType::Pointer pGeometry) const
    {
        for(std::size_t i = 0; i < pGeometry->size(); ++i)
            std::cout << " (" << (*pGeometry)[i].X()
                     << ", " << (*pGeometry)[i].Y()
                     << ", " << (*pGeometry)[i].Z() << "),";
        std::cout << std::endl;
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Finite Cell Auxiliary Utility";
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
    FiniteCellAuxilliaryUtility& operator=(FiniteCellAuxilliaryUtility const& rOther);

    /// Copy constructor.
    FiniteCellAuxilliaryUtility(FiniteCellAuxilliaryUtility const& rOther);


    ///@}

}; // Class FiniteCellAuxilliaryUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, FiniteCellAuxilliaryUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const FiniteCellAuxilliaryUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.


#endif // KRATOS_FINITE_CELL_AUXILLIARY_UTILITY_H_INCLUDED  defined
