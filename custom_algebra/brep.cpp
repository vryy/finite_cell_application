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
//  Date:            13 Mar 2017
//


// Project includes
#include "custom_algebra/brep.h"
#include "custom_utilities/finite_cell_mesh_utility.h"


namespace Kratos
{

BRep::BRep() : mTOL(1.0e-10) {}

BRep::BRep(BRep const& rOther) : mTOL(rOther.mTOL) {}

BRep::~BRep() {}

BRep::Pointer BRep::CloneBRep() const
{
    return BRep::Pointer(new BRep(*this));
}

void BRep::SetTolerance(const double& TOL) {mTOL = TOL;}

const double BRep::GetTolerance() const {return mTOL;}

std::size_t BRep::WorkingSpaceDimension() const
{
    KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
}

std::size_t BRep::LocalSpaceDimension() const
{
    KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
}

bool BRep::IsInside(const PointType& P) const
{
    KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
}

bool BRep::IsOnBoundary(const PointType& P, const double& tol) const
{
    KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
}

int BRep::CutStatus(Element::Pointer p_elem) const
{
    return this->CutStatus(p_elem->GetGeometry());
}

int BRep::CutStatus(GeometryType::Pointer p_geom) const
{
    return this->CutStatus(*p_geom);
}

int BRep::CutStatus(GeometryType& r_geom) const
{
    KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
}

int BRep::CutStatus(const std::vector<PointType>& r_points) const
{
    KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
}

int BRep::CutStatusBySampling(Element::Pointer p_elem, const std::size_t& nsampling) const
{
    return this->CutStatusBySampling(p_elem->GetGeometry(), nsampling);
}

int BRep::CutStatusBySampling(GeometryType::Pointer p_geom, const std::size_t& nsampling) const
{
    return this->CutStatusBySampling(*p_geom, nsampling);
}

int BRep::CutStatusBySampling(GeometryType& r_geom, const std::size_t& nsampling) const
{
    std::vector<PointType> SamplingPoints;
    FiniteCellMeshUtility::GenerateSamplingPoints(SamplingPoints, r_geom, nsampling);
    return this->CutStatus(SamplingPoints);
}

BRep::PointType BRep::Bisect(const BRep::PointType& P1, const BRep::PointType& P2, const double& tol) const
{
    KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
}

}  // namespace Kratos.

