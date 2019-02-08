/*
see finite_cell_isogeometric_structural_application/LICENSE.txt
 */

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2014 Jan 29 $
//   Revision:            $Revision: 1.1 $
//
//
#if !defined(KRATOS_FINITE_CELL_GEO_2D_BEZIER_H_INCLUDED )
#define  KRATOS_FINITE_CELL_GEO_2D_BEZIER_H_INCLUDED

// System includes
#include <iostream>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "utilities/openmp_utils.h"
#include "integration/quadrature.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
//#include "finite_cell_application/custom_geometries/finite_cell_geometry.h"
#include "isogeometric_application/custom_geometries/geo_2d_bezier.h"
#include "isogeometric_application/custom_utilities/bspline_utils.h"
#include "isogeometric_application/custom_utilities/bezier_utils.h"

//#define DEBUG_LEVEL1
//#define DEBUG_LEVEL2
// #define DEBUG_LEVEL3
#define ENABLE_PROFILING

namespace Kratos
{

/**
 * Implementation of finite cell integration scheme for Bezier surface in 2D
 */

template<class TPointType>
class FiniteCellGeo2dBezier : public Geo2dBezier<TPointType>//, public FiniteCellGeometry<Geometry<TPointType> >
{
public:

    /**
     * Type Definitions
     */

    /**
     * IsogeometricGeometry as base class.
     */
    typedef Geo2dBezier<TPointType> BaseType;

    /**
     * The original geometry type
     */
    typedef typename BaseType::GeometryType GeometryType;

    /**
     * Pointer definition of FiniteCellGeo2dBezier
     */
    KRATOS_CLASS_POINTER_DEFINITION( FiniteCellGeo2dBezier );

    /**
     * Integration methods implemented in geometry.
     */
    typedef typename BaseType::IntegrationMethod IntegrationMethod;

    /**
     * A VectorType of counted pointers to Geometries. Used for
     * returning edges of the geometry.
     */
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    /**
     * Redefinition of template parameter TPointType.
     */
    typedef TPointType PointType;

    /**
     * Type used for indexing in geometry class.std::size_t used for indexing
     * point or integration point access methods and also all other
     * methods which need point or integration point index.
     */
    typedef typename BaseType::IndexType IndexType;

    /**
     * This typed used to return size or dimension in
     * geometry. Dimension, WorkingDimension, PointsNumber and
     * ... return this type as their results.
     */
    typedef typename BaseType::SizeType SizeType;

    /**
     * Array of counted pointers to point. This type used to hold
     * geometry's points.
     */
    typedef typename BaseType::PointsArrayType PointsArrayType;

    /**
     * This type used for representing an integration point in
     * geometry. This integration point is a point with an
     * additional weight component.
     */
    typedef typename BaseType::IntegrationPointType IntegrationPointType;

    /**
     * A VectorType of IntegrationPointType which used to hold
     * integration points related to an integration
     * method. IntegrationPoints functions used this type to return
     * their results.
     */
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /**
     * A VectorType of IntegrationPointsArrayType which used to hold
     * integration points related to different integration method
     * implemented in geometry.
     */
    typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

    /**
     * A third order tensor used as shape functions' values
     * container.
     */
    typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;

    /**
     * A fourth order tensor used as shape functions' local
     * gradients container in geometry.
     */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    /**
     * A third order tensor to hold jacobian matrices evaluated at
     * integration points. Jacobian and InverseOfJacobian functions
     * return this type as their result.
     */
    typedef typename BaseType::JacobiansType JacobiansType;

    /**
     * A third order tensor to hold shape functions' local
     * gradients. ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /**
     * A third order tensor to hold shape functions' local second derivatives.
     * ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

    /**
     * Type of the normal vector used for normal to edges in geomety.
     */
    typedef typename BaseType::NormalType NormalType;

    /**
     * Type of coordinates array
     */
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /**
     * Type of Matrix
     */
    typedef typename BaseType::MatrixType MatrixType;

    /**
     * Type of Vector
     */
    typedef typename BaseType::VectorType VectorType;

    /**
     * Type of values container
     */
    typedef typename BaseType::NormalType ValuesContainerType;

    /**
     * Life Cycle
     */

    FiniteCellGeo2dBezier()
    : BaseType( PointsArrayType() )
    {}

    FiniteCellGeo2dBezier( const PointsArrayType& ThisPoints )
    : BaseType( ThisPoints )
    {}

//    FiniteCellGeo2dBezier( const PointsArrayType& ThisPoints, const GeometryData* pGeometryData )
//    : BaseType( ThisPoints, pGeometryData )
//    {}

    /**
     * Copy constructor.
     * Construct this geometry as a copy of given geometry.
     *
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    FiniteCellGeo2dBezier( FiniteCellGeo2dBezier const& rOther )
    : BaseType( rOther )
    {}

    /**
     * Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicity convertible to this
     * geometry PointType.
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> FiniteCellGeo2dBezier( FiniteCellGeo2dBezier<TOtherPointType> const& rOther )
    : Geo2dBezier<TOtherPointType>( rOther )
    {}

    /**
     * Copy constructor.
     * Construct this geometry as a copy of given geometry.
     *
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    FiniteCellGeo2dBezier( BaseType const& rOther )
    : BaseType( rOther )
    {}

    /**
     * Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicity convertible to this
     * geometry PointType.
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> FiniteCellGeo2dBezier( Geo2dBezier<TOtherPointType> const& rOther )
    : Geo2dBezier<TOtherPointType>( rOther )
    {}

    /**
     * Destructor. Does nothing!!!
     */
    virtual ~FiniteCellGeo2dBezier()
    {}

    /**
     * Operators
     */

    /**
     * Assignment operator.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     * @see Clone
     * @see ClonePoints
     */
    FiniteCellGeo2dBezier& operator=( const FiniteCellGeo2dBezier& rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    /**
     * Assignment operator for geometries with different point type.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherPointType>
    FiniteCellGeo2dBezier& operator=( FiniteCellGeo2dBezier<TOtherPointType> const & rOther )
    {
        Geo2dBezier<TOtherPointType>::operator=( rOther );
        return *this;
    }

    /**
     * Assignment operator.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     * @see Clone
     * @see ClonePoints
     */
    FiniteCellGeo2dBezier& operator=( const BaseType& rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    /**
     * Assignment operator for geometries with different point type.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherPointType>
    FiniteCellGeo2dBezier& operator=( Geo2dBezier<TOtherPointType> const & rOther )
    {
        Geo2dBezier<TOtherPointType>::operator=( rOther );
        return *this;
    }

    /**
     * Operations
     */

    virtual typename GeometryType::Pointer Create( PointsArrayType const& ThisPoints ) const
    {
        typename BaseType::Pointer pBezierGeometry = boost::dynamic_pointer_cast<BaseType>(BaseType::Create( ThisPoints ));
        FiniteCellGeo2dBezier::Pointer pFiniteCellBezierGeometry = FiniteCellGeo2dBezier::Pointer(new FiniteCellGeo2dBezier(*pBezierGeometry));
        return pFiniteCellBezierGeometry;
    }

    virtual boost::shared_ptr< Geometry< Point<3> > > Clone() const
    {
//        Geometry< Point<3> >::PointsArrayType NewPoints;
//        //making a copy of the nodes TO POINTS (not Nodes!!!)

//        for ( IndexType i = 0; i < this->Points().size(); ++i )
//        NewPoints.push_back( this->Points()[i] );

//        //creating a geometry with the new points
//        boost::shared_ptr< Geometry< Point<3> > >
//        p_clone( new FiniteCellGeo2dBezier< Point<3> >( NewPoints ) );

//        p_clone->ClonePoints();

//        return p_clone;

        KRATOS_THROW_ERROR(std::logic_error, "FiniteCellGeo2dBezier geometry does not support for Clone", *this)
    }

    /**
     * Informations
     */

    virtual GeometryData::KratosGeometryFamily GetGeometryFamily()
    {
        return GeometryData::Kratos_NURBS;
    }

    virtual GeometryData::KratosGeometryType GetGeometryType()
    {
        return GeometryData::Kratos_Bezier2D;
    }

    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     *
     * @return String contains information about this geometry.
     * @see PrintData()
     * @see PrintInfo()
     */
    virtual std::string Info() const
    {
        return "2 dimensional Finite Cell Bezier decomposition surface in 2D space";
    }

    /**
     * Print information about this object.
     *
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    virtual void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << Info();
    }

    /**
     * Print geometry's data into given stream.
     * Prints it's points by the order they stored in the geometry
     * and then center point of geometry.
     *
     * @param rOStream Stream to print into it.
     * @see PrintInfo()
     * @see Info()
     */
    virtual void PrintData( std::ostream& rOStream ) const
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
//        MatrixType jacobian;
//        Jacobian( jacobian, PointType() );
//        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }

    /**
     * TO BE CALLED BY ELEMENT
     */
    /// Assign a list of integration points to the geometry. The ThisIntegrationMethod will
    /// ultimately becomes the default integration rule on the geometry.
    void AssignGeometryData(const GeometryData::IntegrationMethod ThisIntegrationMethod,
            const IntegrationPointsArrayType& integration_points)
    {
        BaseType::mpBezierGeometryData = BezierUtils::CreateIntegrationRule<2, 2, 2>(ThisIntegrationMethod, BaseType::mOrder1, BaseType::mOrder2, integration_points);
        GeometryType::mpGeometryData = &(*BaseType::mpBezierGeometryData);
    }

protected:

private:

    /**
     * Static Member Variables
     */

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointsArrayType );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointsArrayType );
    }

    /**
     * Private Operations
     */

    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class FiniteCellGeo2dBezier;

    /**
     * Un accessible methods
     */

};    // Class FiniteCellGeo2dBezier

/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >>(
        std::istream& rIStream, FiniteCellGeo2dBezier<TPointType>& rThis);

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator <<(
        std::ostream& rOStream, const FiniteCellGeo2dBezier<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}    // namespace Kratos.

#endif

