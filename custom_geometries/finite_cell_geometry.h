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
#ifdef SD_APP_FORWARD_COMPATIBILITY
#include "geometries/geometry_dimension.h"
#endif


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
 * - TBaseGeometryType which reperesent the type of the geometry this geometry
 * type build on.
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
    typedef FiniteCellGeometry<TBaseGeometryType> FiniteCellGeometryType;

    /// Pointer definition of FiniteCellGeometry
    KRATOS_CLASS_POINTER_DEFINITION( FiniteCellGeometry );

    /** Base type for base geometry.
    */
    typedef TBaseGeometryType BaseType;

    /** Redefinition of geometry point type.
     */
    typedef typename BaseType::PointType PointType;

    /** original geometry type.
    */
//    typedef typename BaseType::BaseType GeometryType;
    typedef Geometry<PointType> GeometryType;

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
    typedef PointerVector<FiniteCellGeometryType> GeometriesArrayType;

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
    typedef typename BaseType::ptr_iterator                 ptr_iterator;
    typedef typename BaseType::ptr_const_iterator           ptr_const_iterator;
#ifndef SD_APP_FORWARD_COMPATIBILITY
    typedef typename BaseType::reverse_iterator             reverse_iterator;
    typedef typename BaseType::const_reverse_iterator       const_reverse_iterator;
    typedef typename BaseType::ptr_reverse_iterator         ptr_reverse_iterator;
    typedef typename BaseType::ptr_const_reverse_iterator   ptr_const_reverse_iterator;
#endif

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
//        std::cout << "empty constructor FiniteCellGeometry is called" << std::endl;
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
//        std::cout << "constructor FiniteCellGeometry is called" << std::endl;
//        KRATOS_WATCH(typeid(*this).name())
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
//        std::cout << "copy constructor FiniteCellGeometry is called" << std::endl;
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

    ///@}
    ///@name Operations
    ///@{


    typename GeometryType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename GeometryType::Pointer( new FiniteCellGeometryType( ThisPoints ) );
    }


    /// Assign the new geometry_data to the geometry. Several considerations are made regarding the geometry for finite cell
    /// + the finite cell geometry uses the same quadrature location as the base geometry, however the weight are passed to the finite cell geometry
    /// + the finite cell geometry only supports single integration rule, which is passed in the arguments
    /// + the shape function and local gradients are the same as the base geometry
    void AssignGeometryData(const GeometryData::IntegrationMethod& ThisIntegrationMethod, const Vector& rWeights)
    {
        // copy the integration points data
        IntegrationPointsContainerType all_integration_points;
        all_integration_points[static_cast<int>(ThisIntegrationMethod)] = BaseType::IntegrationPoints(ThisIntegrationMethod);

        // modify the weight
        IntegrationPointsArrayType& integration_points = all_integration_points[static_cast<int>(ThisIntegrationMethod)];
        if (integration_points.size() != rWeights.size())
            KRATOS_THROW_ERROR(std::logic_error, "The number of integration points and the size of input array is not compatible", "")
            for (std::size_t i = 0; i < rWeights.size(); ++i)
            {
                integration_points[i].SetWeight(rWeights[i]);
            }

        // copy the shape function values data
        ShapeFunctionsValuesContainerType shape_functions_values;
        shape_functions_values[static_cast<int>(ThisIntegrationMethod)] = GeometryType::ShapeFunctionsValues(ThisIntegrationMethod);

        // copy the shape function local gradients data
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients;
//        shape_functions_local_gradients[static_cast<int>(ThisIntegrationMethod)] = BaseType::ShapeFunctionsLocalGradients(ThisIntegrationMethod); // I do not know why this does not compile. I skip this for now.
        shape_functions_local_gradients[static_cast<int>(ThisIntegrationMethod)] = GeometryType::ShapeFunctionsLocalGradients(ThisIntegrationMethod);

#ifdef SD_APP_FORWARD_COMPATIBILITY
        GeometryDimension::Pointer pGeometryDimension = GeometryDimension::Pointer(
                    new GeometryDimension(
                        BaseType::WorkingSpaceDimension(),
                        BaseType::LocalSpaceDimension())
                );

        // create new geometry data
        mpFiniteCellGeometryData = GeometryData::Pointer(
                                       new GeometryData(
                                           &(*pGeometryDimension),
                                           ThisIntegrationMethod,              //ThisDefaultMethod
                                           all_integration_points,             //ThisIntegrationPoints
                                           shape_functions_values,             //ThisShapeFunctionsValues
                                           shape_functions_local_gradients     //ThisShapeFunctionsLocalGradients
                                       )
                                   );

        // assign the geometry data back to the original geometry
        BaseType::SetGeometryData(&(*mpFiniteCellGeometryData));
#else
        // create new geometry data
        mpFiniteCellGeometryData = GeometryData::Pointer(
                                       new GeometryData(
                                           BaseType::Dimension(),
                                           BaseType::WorkingSpaceDimension(),
                                           BaseType::LocalSpaceDimension(),
                                           ThisIntegrationMethod,              //ThisDefaultMethod
                                           all_integration_points,             //ThisIntegrationPoints
                                           shape_functions_values,             //ThisShapeFunctionsValues
                                           shape_functions_local_gradients     //ThisShapeFunctionsLocalGradients
                                       )
                                   );

        // assign the geometry data back to the original geometry
        BaseType::mpGeometryData = &(*mpFiniteCellGeometryData);
#endif
    }


    /// Assign a list of integration points to the geometry. The ThisIntegrationMethod will
    /// ultimately becomes the default integration rule on the geometry.
    void AssignGeometryData(const GeometryData::IntegrationMethod& ThisIntegrationMethod,
                            const IntegrationPointsArrayType& integration_points)
    {
        // create the integration points data
        IntegrationPointsContainerType all_integration_points;
        all_integration_points[static_cast<int>(ThisIntegrationMethod)] = integration_points;

        // compute the shape_functions_values
        ShapeFunctionsValuesContainerType shape_functions_values;
        Matrix& rValues = shape_functions_values[static_cast<int>(ThisIntegrationMethod)];

        rValues.resize(integration_points.size(), this->size());
        Vector aux(this->size());
        for (std::size_t i = 0; i < integration_points.size(); ++i)
        {
            noalias( aux ) = BaseType::ShapeFunctionsValues( aux, integration_points[i] );
            noalias( row(rValues, i) ) = aux;
        }

        // compute the shape_functions_local_gradients
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients;
        ShapeFunctionsGradientsType& rLocalGradients = shape_functions_local_gradients[static_cast<int>(ThisIntegrationMethod)];

        rLocalGradients.resize(integration_points.size());
        for (std::size_t i = 0; i < integration_points.size(); ++i)
        {
            rLocalGradients[i].resize( this->size(), this->LocalSpaceDimension(), false );

            noalias( rLocalGradients[i] ) = BaseType::ShapeFunctionsLocalGradients( rLocalGradients[i], integration_points[i] );
        }

#ifdef SD_APP_FORWARD_COMPATIBILITY
        GeometryDimension::Pointer pGeometryDimension = GeometryDimension::Pointer(
                    new GeometryDimension(
                        BaseType::WorkingSpaceDimension(),
                        BaseType::LocalSpaceDimension())
                );

        // create new geometry data
        mpFiniteCellGeometryData = GeometryData::Pointer(
                                       new GeometryData(
                                           &(*pGeometryDimension),
                                           ThisIntegrationMethod,              //ThisDefaultMethod
                                           all_integration_points,             //ThisIntegrationPoints
                                           shape_functions_values,             //ThisShapeFunctionsValues
                                           shape_functions_local_gradients     //ThisShapeFunctionsLocalGradients
                                       )
                                   );

        // assign the geometry data back to the original geometry
        BaseType::SetGeometryData(&(*mpFiniteCellGeometryData));
#else
        // create new geometry data
        mpFiniteCellGeometryData = GeometryData::Pointer(
                                       new GeometryData(
                                           BaseType::Dimension(),
                                           BaseType::WorkingSpaceDimension(),
                                           BaseType::LocalSpaceDimension(),
                                           ThisIntegrationMethod,              //ThisDefaultMethod
                                           all_integration_points,             //ThisIntegrationPoints
                                           shape_functions_values,             //ThisShapeFunctionsValues
                                           shape_functions_local_gradients     //ThisShapeFunctionsLocalGradients
                                       )
                                   );

        // assign the geometry data back to the original geometry
        BaseType::mpGeometryData = &(*mpFiniteCellGeometryData);
#endif
    }


    ///@}
    ///@name Input and output
    ///@{

    /** Turn back information as a string.

    @return String contains information about this geometry.
    @see PrintData()
    @see PrintInfo()
    */
    std::string Info() const override
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
    void PrintInfo( std::ostream& rOStream ) const override
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
    void PrintData( std::ostream& rOStream ) const override
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

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
    }

    void load( Serializer& rSerializer ) override
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


