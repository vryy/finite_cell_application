//
//   Project Name:        Kratos
//   Last Modified by:    $Author:   hbui        $
//   Date:                $Date:   2018          $
//   Revision:            $Revision:         1.0 $
//
//

#if !defined(KRATOS_QUADRILATERAL_GAUSS_LOBATTO_ADDITIONAL_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_QUADRILATERAL_GAUSS_LOBATTO_ADDITIONAL_INTEGRATION_POINTS_H_INCLUDED


// System includes
#include <array>

// External includes

// Project includes
#include "integration/quadrature.h"

// add more Gauss-Lobatto quadrature for quadrilateral

namespace Kratos
{

class QuadrilateralGaussLobattoIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussLobattoIntegrationPoints3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 9> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 9;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
                IntegrationPointType( -1.00, -1.00, 1.00 / 9.00 ),
                IntegrationPointType( -1.00, 0.00, 4.00 / 9.00 ),
                IntegrationPointType( -1.00, 1.00, 1.00 / 9.00 ),
                IntegrationPointType( 0.00, -1.00, 4.00 / 9.00 ),
                IntegrationPointType( 0.00, 0.00, 16.00 / 9.00 ),
                IntegrationPointType( 0.00, 1.00, 4.00 / 9.00 ),
                IntegrationPointType( 1.00, -1.00, 1.00 / 9.00 ),
                IntegrationPointType( 1.00, 0.00, 4.00 / 9.00 ),
                IntegrationPointType( 1.00, 1.00, 1.00 / 9.00 )
            }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Gauss-Lobatto integration 3 ";
        return buffer.str();
    }
}; // Class QuadrilateralGaussLobattoIntegrationPoints3

class QuadrilateralGaussLobattoIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussLobattoIntegrationPoints4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 16> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 16;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
                IntegrationPointType(-1.00, -1.00, 1.00 / 36.00),
                IntegrationPointType(-std::sqrt(5.00) / 5.00, -1.00, 5.00 / 36.00),
                IntegrationPointType( std::sqrt(5.00) / 5.00, -1.00, 5.00 / 36.00),
                IntegrationPointType( 1.00, -1.00, 1.00 / 36.00),
                IntegrationPointType(-1.00, -std::sqrt(5.00) / 5.00, 5.00 / 36.00),
                IntegrationPointType(-std::sqrt(5.00) / 5.00, -std::sqrt(5.00) / 5.00, 25.00 / 36.00),
                IntegrationPointType( std::sqrt(5.00) / 5.00, -std::sqrt(5.00) / 5.00, 25.00 / 36.00),
                IntegrationPointType( 1.00, -std::sqrt(5.00) / 5.00, 5.00 / 36.00),
                IntegrationPointType(-1.00, std::sqrt(5.00) / 5.00, 5.00 / 36.00),
                IntegrationPointType(-std::sqrt(5.00) / 5.00, std::sqrt(5.00) / 5.00, 25.00 / 36.00),
                IntegrationPointType( std::sqrt(5.00) / 5.00, std::sqrt(5.00) / 5.00, 25.00 / 36.00),
                IntegrationPointType( 1.00, std::sqrt(5.00) / 5.00, 5.00 / 36.00),
                IntegrationPointType(-1.00, 1.00, 1.00 / 36.00),
                IntegrationPointType(-std::sqrt(5.00) / 5.00, 1.00, 5.00 / 36.00),
                IntegrationPointType( std::sqrt(5.00) / 5.00, 1.00, 5.00 / 36.00),
                IntegrationPointType( 1.00, 1.00, 1.00 / 36.00)
            }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Gauss-Lobatto integration 4 ";
        return buffer.str();
    }
}; // Class QuadrilateralGaussLobattoIntegrationPoints4

class QuadrilateralGaussLobattoIntegrationPoints5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussLobattoIntegrationPoints5);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 25> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 25;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
                IntegrationPointType(-1.00, -1.00, 1.00 / 10.00 * 1.00 / 10.00),
                IntegrationPointType(-std::sqrt(21.00) / 7.00, -1.00, 49.00 / 90.00 * 1.00 / 10.00),
                IntegrationPointType( 0.00, -1.00, 32.00 / 45.00 * 1.00 / 10.00),
                IntegrationPointType( std::sqrt(21.00) / 7.00, -1.00, 49.00 / 90.00 * 1.00 / 10.00),
                IntegrationPointType( 1.00, -1.00, 1.00 / 10.00 * 1.00 / 10.00),

                IntegrationPointType(-1.00, -std::sqrt(21.00) / 7.00, 1.00 / 10.00 * 49.00 / 90.00),
                IntegrationPointType(-std::sqrt(21.00) / 7.00, -std::sqrt(21.00) / 7.00, 49.00 / 90.00 * 49.00 / 90.00),
                IntegrationPointType( 0.00, -std::sqrt(21.00) / 7.00, 32.00 / 45.00 * 49.00 / 90.00),
                IntegrationPointType( std::sqrt(21.00) / 7.00, -std::sqrt(21.00) / 7.00, 49.00 / 90.00 * 49.00 / 90.00),
                IntegrationPointType( 1.00, -std::sqrt(21.00) / 7.00, 1.00 / 10.00 * 49.00 / 90.00),

                IntegrationPointType(-1.00, 0.00, 1.00 / 10.00 * 32.00 / 45.00),
                IntegrationPointType(-std::sqrt(21.00) / 7.00, 0.00, 49.00 / 90.00 * 32.00 / 45.00),
                IntegrationPointType( 0.00, 0.00, 32.00 / 45.00 * 32.00 / 45.00),
                IntegrationPointType( std::sqrt(21.00) / 7.00, 0.00, 49.00 / 90.00 * 32.00 / 45.00),
                IntegrationPointType( 1.00, 0.00, 1.00 / 10.00 * 32.00 / 45.00),

                IntegrationPointType(-1.00, std::sqrt(21.00) / 7.00, 1.00 / 10.00 * 49.00 / 90.00),
                IntegrationPointType(-std::sqrt(21.00) / 7.00, std::sqrt(21.00) / 7.00, 49.00 / 90.00 * 49.00 / 90.00),
                IntegrationPointType( 0.00, std::sqrt(21.00) / 7.00, 32.00 / 45.00 * 49.00 / 90.00),
                IntegrationPointType( std::sqrt(21.00) / 7.00, std::sqrt(21.00) / 7.00, 49.00 / 90.00 * 49.00 / 90.00),
                IntegrationPointType( 1.00, std::sqrt(21.00) / 7.00, 1.00 / 10.00 * 49.00 / 90.00),

                IntegrationPointType(-1.00, 1.00, 1.00 / 10.00 * 1.00 / 10.00),
                IntegrationPointType(-std::sqrt(21.00) / 7.00, 1.00, 49.00 / 90.00 * 1.00 / 10.00),
                IntegrationPointType( 0.00, 1.00, 32.00 / 45.00 * 1.00 / 10.00),
                IntegrationPointType( std::sqrt(21.00) / 7.00, 1.00, 49.00 / 90.00 * 1.00 / 10.00),
                IntegrationPointType( 1.00, 1.00, 1.00 / 10.00 * 1.00 / 10.00)
            }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Gauss-Lobatto integration 5 ";
        return buffer.str();
    }
}; // Class QuadrilateralGaussLobattoIntegrationPoints5

}

#endif // KRATOS_QUADRILATERAL_GAUSS_LOBATTO_ADDITIONAL_INTEGRATION_POINTS_H_INCLUDED defined


