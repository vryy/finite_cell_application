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
//  Date:            3 Feb 2018
//


#if !defined(KRATOS_FINITE_CELL_MESH_UTILITY_HPP_INCLUDED )
#define  KRATOS_FINITE_CELL_MESH_UTILITY_HPP_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"


namespace Kratos
{

template<>
struct GenerateStructuredMesh_Helper<2, 1> : public GenerateStructuredMesh_Helper<0, 0>
{
    static void Execute(std::vector<std::vector<PointType> >& rPoints,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::size_t>& nsampling)
    {
        double x, y;
        double dx = (EndPoint[0] - StartPoint[0]) / nsampling[0];
        double dy = (EndPoint[1] - StartPoint[1]) / nsampling[1];

        rPoints.resize(nsampling[1]+1);
        for (std::size_t j = 0; j < nsampling[1]+1; ++j)
        {
            rPoints[j].resize(nsampling[0]+1);
            y = StartPoint[1] + j*dy;

            for (std::size_t i = 0; i < nsampling[0]+1; ++i)
            {
                x = StartPoint[0] + i*dx;
                rPoints[j][i] = PointType(x, y, 0.0);
            }
        }
    }

    static void Execute(std::vector<std::vector<PointType> >& rPoints,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::vector<double> >& sampling)
    {
        std::stringstream ss;
        ss << __FUNCTION__ << " is not yet implemented for dim = 2 and type = 1";
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }
};

template<>
struct GenerateStructuredMesh_Helper<2, 2> : public GenerateStructuredMesh_Helper<0, 0>
{
    static void Execute(std::vector<std::vector<PointType> >& rPoints,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::size_t>& nsampling)
    {
        std::stringstream ss;
        ss << __FUNCTION__ << " is not yet implemented for dim = 2 and type = 2";
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    static void Execute(std::vector<std::vector<PointType> >& rPoints,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::vector<double> >& sampling)
    {
        std::stringstream ss;
        ss << __FUNCTION__ << " is not yet implemented for dim = 2 and type = 2";
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }
};

template<>
struct GenerateStructuredMesh_Helper<2, 3> : public GenerateStructuredMesh_Helper<0, 0>
{
    static void Execute(std::vector<std::vector<PointType> >& rPoints,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::size_t>& nsampling)
    {
        double x, y;
        double dx = 0.5 * (EndPoint[0] - StartPoint[0]) / nsampling[0];
        double dy = 0.5 * (EndPoint[1] - StartPoint[1]) / nsampling[1];

        rPoints.resize(2*nsampling[1]+1);
        for (std::size_t j = 0; j < 2*nsampling[1]+1; ++j)
        {
            rPoints[j].resize(2*nsampling[0]+1);
            y = StartPoint[1] + j*dy;

            for (std::size_t i = 0; i < 2*nsampling[0]+1; ++i)
            {
                x = StartPoint[0] + i*dx;
                rPoints[j][i] = PointType(x, y, 0.0);
            }
        }
    }

    static void Execute(std::vector<std::vector<PointType> >& rPoints,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::vector<double> >& sampling)
    {
        std::stringstream ss;
        ss << __FUNCTION__ << " is not yet implemented for dim = 2 and type = 3";
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }
};

template<>
struct GenerateStructuredMesh_Helper<3, 1> : public GenerateStructuredMesh_Helper<0, 0>
{
    static void Execute(std::vector<std::vector<std::vector<PointType> > >& rPoints,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::size_t>& nsampling)
    {
        double x, y, z;
        double dx = (EndPoint[0] - StartPoint[0]) / nsampling[0];
        double dy = (EndPoint[1] - StartPoint[1]) / nsampling[1];
        double dz = (EndPoint[2] - StartPoint[2]) / nsampling[2];

        rPoints.resize(nsampling[2]+1);
        for (std::size_t k = 0; k < nsampling[2]+1; ++k)
        {
            rPoints[k].resize(nsampling[1]+1);
            z = StartPoint[2] + k*dz;

            for (std::size_t j = 0; j < nsampling[1]+1; ++j)
            {
                rPoints[k][j].resize(nsampling[0]+1);
                y = StartPoint[1] + j*dy;

                for (std::size_t i = 0; i < nsampling[0]+1; ++i)
                {
                    x = StartPoint[0] + i*dx;
                    rPoints[k][j][i] = PointType(x, y, z);
                }
            }
        }
    }

    static void Execute(std::vector<std::vector<std::vector<PointType> > >& rPoints,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::vector<double> >& sampling)
    {
        double x, y, z;
        double dx = (EndPoint[0] - StartPoint[0]);
        double dy = (EndPoint[1] - StartPoint[1]);
        double dz = (EndPoint[2] - StartPoint[2]);

        rPoints.resize(sampling[2].size());
        for (std::size_t k = 0; k < sampling[2].size(); ++k)
        {
            z = StartPoint[2] + sampling[2][k]*dz;

            rPoints[k].resize(sampling[1].size());
            for (std::size_t j = 0; j < sampling[1].size(); ++j)
            {
                y = StartPoint[1] + sampling[1][j]*dy;

                rPoints[k][j].resize(sampling[0].size());
                for (std::size_t i = 0; i < sampling[0].size(); ++i)
                {
                    x = StartPoint[0] + sampling[0][i]*dx;
                    rPoints[k][j][i] = PointType(x, y, z);
                }
            }
        }
    }
};

template<>
struct GenerateStructuredMesh_Helper<3, 2> : public GenerateStructuredMesh_Helper<0, 0>
{
    static void Execute(std::vector<std::vector<std::vector<PointType> > >& rPoints,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::size_t>& nsampling)
    {
        std::stringstream ss;
        ss << __FUNCTION__ << " is not yet implemented for dim = 3 and type = 2";
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    static void Execute(std::vector<std::vector<std::vector<PointType> > >& rPoints,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::vector<double> >& sampling)
    {
        std::stringstream ss;
        ss << __FUNCTION__ << " is not yet implemented for dim = 3 and type = 2";
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }
};

template<>
struct GenerateStructuredMesh_Helper<3, 3> : public GenerateStructuredMesh_Helper<0, 0>
{
    static void Execute(std::vector<std::vector<std::vector<PointType> > >& rPoints,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::size_t>& nsampling)
    {
        double x, y, z;
        double dx = 0.5 * (EndPoint[0] - StartPoint[0]) / nsampling[0];
        double dy = 0.5 * (EndPoint[1] - StartPoint[1]) / nsampling[1];
        double dz = 0.5 * (EndPoint[2] - StartPoint[2]) / nsampling[2];

        rPoints.resize(2*nsampling[2]+1);
        for (std::size_t k = 0; k < 2*nsampling[2]+1; ++k)
        {
            rPoints[k].resize(2*nsampling[1]+1);
            z = StartPoint[2] + k*dz;

            for (std::size_t j = 0; j < 2*nsampling[1]+1; ++j)
            {
                rPoints[k][j].resize(2*nsampling[0]+1);
                y = StartPoint[1] + j*dy;

                for (std::size_t i = 0; i < 2*nsampling[0]+1; ++i)
                {
                    x = StartPoint[0] + i*dx;
                    rPoints[k][j][i] = PointType(x, y, z);
                }
            }
        }
    }

    static void Execute(std::vector<std::vector<std::vector<PointType> > >& rPoints,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const std::vector<std::vector<double> >& sampling)
    {
        double x, y, z;
        double dx = (EndPoint[0] - StartPoint[0]);
        double dy = (EndPoint[1] - StartPoint[1]);
        double dz = (EndPoint[2] - StartPoint[2]);

        rPoints.resize(2*sampling[2].size()-1);
        for (std::size_t k = 0; k < 2*sampling[2].size()-1; ++k)
        {
            if (k%2 == 0)
                z = StartPoint[2] + sampling[2][k/2]*dz;
            else
                z = StartPoint[2] + 0.5*(sampling[2][k/2] + sampling[2][k/2+1])*dz;

            rPoints[k].resize(2*sampling[1].size()-1);
            for (std::size_t j = 0; j < 2*sampling[1].size()-1; ++j)
            {
                if (j%2 == 0)
                    y = StartPoint[1] + sampling[1][j/2]*dy;
                else
                    y = StartPoint[1] + 0.5*(sampling[1][j/2] + sampling[1][j/2+1])*dy;

                rPoints[k][j].resize(2*sampling[0].size()-1);
                for (std::size_t i = 0; i < 2*sampling[0].size()-1; ++i)
                {
                    if (i%2 == 0)
                        x = StartPoint[0] + sampling[0][i/2]*dx;
                    else
                        x = StartPoint[0] + 0.5*(sampling[0][i/2]+sampling[0][i/2+1])*dx;
                    rPoints[k][j][i] = PointType(x, y, z);
                }
            }
        }
    }
};

}  // namespace Kratos.

#endif // KRATOS_FINITE_CELL_MESH_UTILITY_HPP_INCLUDED  defined

