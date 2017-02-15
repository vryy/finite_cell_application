// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 10 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_FINITE_CELL_GEOMETRY_H_INCLUDED )
#define  KRATOS_FINITE_CELL_GEOMETRY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes
#include <boost/array.hpp>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "containers/pointer_vector.h"
#include "integration/integration_point.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "utilities/math_utils.h"
#include "integration/quadrature.h"


namespace Kratos
{
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

///FiniteCellGeometry base class.
/** As a base class FiniteCellGeometry has all the common
 * interface of Kratos' geometries. Also it contains array of
 * pointers to its points, reference to shape functions values in
 * all integrations points and also local gradients of shape
 * functions evaluated in all integrations points.
 *
 * FiniteCellGeometry is a template class with just one template parameter:
 * - TPointType which reperesent the type of the point this geometry
 * type contain and build on.
 *
 * @see Point
 * @see Node
 * @see Formulation
 * @see GeometryAndFormulationElement
 */
template<class TBaseGeometryType>
class FiniteCellGeometry : public TBaseGeometryType
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// This FiniteCellGeometry type.
    typedef FiniteCellGeometry<TBaseGeometryType> GeometryType;

    /// Pointer definition of FiniteCellGeometry
    KRATOS_CLASS_POINTER_DEFINITION( FiniteCellGeometry );

    /** Base type for geometry.
    */
    typedef TBaseGeometryType BaseType;

    /** Redefinition of geometry point type.
     */
    typedef typename BaseType::PointType PointType;

    /** Array of counted pointers to point. This type used to hold
    geometry's points.
    */
    typedef PointerVector<PointType> PointsArrayType;

    /** Integration methods implemented in geometry.
    */
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /** A Vector of counted pointers to Geometries. Used for
    returning edges of the geometry.
     */
    typedef PointerVector<GeometryType> GeometriesArrayType;

    /** Type used for indexing in geometry class.std::size_t used for indexing
    point or integration point access methods and also all other
    methods which need point or integration point index.
    */
    typedef std::size_t IndexType;

    /** This typed used to return size or dimension in
    geometry. Dimension, WorkingDimension, PointsNumber and
    ... return this type as their results.
    */
    typedef std::size_t SizeType;

    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /** This type used for representing an integration point in
    geometry. This integration point is a point with an
    additional weight component.
    */
    typedef typename BaseType::IntegrationPointType IntegrationPointType;

    /** A Vector of IntegrationPointType which used to hold
    integration points related to an integration
    method. IntegrationPoints functions used this type to return
    their results.
    */
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /** A Vector of IntegrationPointsArrayType which used to hold
    integration points related to different integration method
    implemented in geometry.
    */
    typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

    /** A third order tensor used as shape functions' values
    continer.
    */
    typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;

    /** A fourth order tensor used as shape functions' local
    gradients container in geometry.
    */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    /** A third order tensor to hold jacobian matrices evaluated at
    integration points. Jacobian and InverseOfJacobian functions
    return this type as their result.
    */
    typedef typename BaseType::JacobiansType JacobiansType;

    /** A third order tensor to hold shape functions'  gradients.
    ShapefunctionsGradients function return this
    type as its result.
    */
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /** A third order tensor to hold shape functions' local second derivatives.
    ShapefunctionsLocalGradients function return this
    type as its result.
    */
    typedef typename BaseType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

    /** A fourth order tensor to hold shape functions' local third order derivatives
     */
    typedef typename BaseType::ShapeFunctionsThirdDerivativesType ShapeFunctionsThirdDerivativesType;

    /** Type of the normal vector used for normal to edges in geomety.
     */
    typedef typename BaseType::NormalType NormalType;
    
    typedef typename BaseType::NormalType ValuesContainerType;

    /**
     * Type of iterators
     */
    typedef typename BaseType::iterator                     iterator;
    typedef typename BaseType::const_iterator               const_iterator;
    typedef typename BaseType::reverse_iterator             reverse_iterator;
    typedef typename BaseType::const_reverse_iterator       const_reverse_iterator;
    typedef typename BaseType::ptr_iterator                 ptr_iterator;
    typedef typename BaseType::ptr_const_iterator           ptr_const_iterator;
    typedef typename BaseType::ptr_reverse_iterator         ptr_reverse_iterator;
    typedef typename BaseType::ptr_const_reverse_iterator   ptr_const_reverse_iterator;
    
    /**
     * Type of Matrix
     */
    typedef Matrix MatrixType;
    
    /**
     * Type of Vector
     */
    typedef Vector VectorType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    FiniteCellGeometry() : BaseType()
    {
    }

    /** Complete argument constructor. This constructor gives a
    complete set of arguments to pass all the initial value of
    all the member variables of geometry class. Also it has
    default value for integration variables to make it usefull
    in the case of constructing new geometry without mapping and
    integrating properties.

    @param ThisPoints Vector of pointers to points which this
    geometry constructing on them. Points must have dimension
    equal or greater than working space dimension though there
    is no control on it.

    TODO
    */
    FiniteCellGeometry( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints )
    {
    }

    /** Copy constructor.
    Construct this geometry as a copy of given geometry.

    @note This copy constructor don't copy the points and new
    geometry shares points with given source geometry. It's
    obvious that any change to this new geometry's point affect
    source geometry's points too.
    */
    FiniteCellGeometry( const FiniteCellGeometry& rOther )
        : BaseType( rOther )
    {
    }

    /** Copy constructor from a geometry with other point type.
    Construct this geometry as a copy of given geometry which
    has different type of points. The given goemetry's
    TOtherPointType* must be implicity convertible to this
    geometry PointType.

    @note This copy constructor don't copy the points and new
    geometry shares points with given source geometry. It's
    obvious that any change to this new geometry's point affect
    source geometry's points too.
    */
    template<class TOtherPointType> FiniteCellGeometry( FiniteCellGeometry<TOtherPointType> const & rOther )
        : BaseType( rOther.begin(), rOther.end() )
    {
    }

    /// Destructor. Do nothing!!!
    virtual ~FiniteCellGeometry() {}

    ///@}
    ///@name Operators
    ///@{

    /** Assignment operator.

    @note This operator don't copy the points and this
    geometry shares points with given source geometry. It's
    obvious that any change to this geometry's point affect
    source geometry's points too.

    @see Clone
    @see ClonePoints
    */
    FiniteCellGeometry& operator=( const FiniteCellGeometry& rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    /** Assignment operator for geometries with different point type.

    @note This operator don't copy the points and this
    geometry shares points with given source geometry. It's
    obvious that any change to this geometry's point affect
    source geometry's points too.

    @see Clone
    @see ClonePoints
    */
    template<class TOtherPointType>
    FiniteCellGeometry& operator=( FiniteCellGeometry<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /** Turn back information as a string.

    @return String contains information about this geometry.
    @see PrintData()
    @see PrintInfo()
    */
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << BaseType::Dimension()
               << " dimensional finite cell geometry in "
               << BaseType::WorkingSpaceDimension()
               << "D space";
        return buffer.str();
    }

    /** Print information about this object.

    @param rOStream Stream to print into it.
    @see PrintData()
    @see Info()
    */
    virtual void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << BaseType::Dimension()
                 << " dimensional finite cell geometry in "
                 << BaseType::WorkingSpaceDimension()
                 << "D space";
    }

    /** Print geometry's data into given stream. Prints it's points
    by the order they stored in the geometry and then center
    point of geometry.

    @param rOStream Stream to print into it.
    @see PrintInfo()
    @see Info()
    */
    virtual void PrintData( std::ostream& rOStream ) const
    {
        BaseType::PrintData(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:

    ///@}
    ///@name Protected static Member Variables
    ///@{
    

    ///@}
    ///@name Protected member Variables
    ///@{


    GeometryData::Pointer mpFiniteCellGeometryData;


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

    /** Protected Constructor.
    Avoids object to be created Except for derived classes
    */


    ///@}


private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
    }


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
    ///@name Private Friends
    ///@{

    template<class TOtherPointType> friend class FiniteCellGeometry;

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class FiniteCellGeometry

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TBaseGeometryType>
inline std::istream& operator >> ( std::istream& rIStream,
                                   FiniteCellGeometry<TBaseGeometryType>& rThis );

/// output stream function
template<class TBaseGeometryType>
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const FiniteCellGeometry<TBaseGeometryType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_FINITE_CELL_GEOMETRY_H_INCLUDED  defined 


