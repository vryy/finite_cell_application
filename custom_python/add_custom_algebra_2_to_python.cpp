// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//



// Project includes
#include "includes/element.h"
#include "containers/array_1d.h"
#include "custom_python/add_custom_algebra_to_python.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/brep.h"
#include "custom_algebra/and_brep.h"
#include "custom_algebra/level_set/level_set.h"
#include "custom_algebra/level_set/circular_level_set.h"
#include "custom_algebra/level_set/doughnut_level_set.h"
#include "custom_algebra/level_set/spherical_level_set.h"
#include "custom_algebra/level_set/cylinder_level_set.h"
#include "custom_algebra/level_set/linear_level_set.h"
#include "custom_algebra/level_set/planar_level_set.h"
#include "custom_algebra/level_set/product_level_set.h"
#include "custom_algebra/level_set/inverse_level_set.h"
#include "custom_algebra/level_set/union_level_set.h"
#include "custom_algebra/level_set/intersection_level_set.h"
#include "custom_algebra/level_set/difference_level_set.h"
#include "custom_algebra/curve/parametric_curve.h"
#include "custom_algebra/surface/parametric_surface.h"
#include "custom_algebra/volume/parametric_volume.h"


namespace Kratos
{

const int BRep::_CUT;
const int BRep::_IN;
const int BRep::_OUT;

namespace Python
{

using namespace boost::python;

void FiniteCellApplication_AddBRepAndLevelSetToPython()
{
    /**************************************************************/
    /************* EXPORT INTERFACE FOR BREP **********************/
    /**************************************************************/

    int(BRep::*pointer_to_CutStatusElement)(Element::Pointer) const = &BRep::CutStatus;
    int(BRep::*pointer_to_CutStatusGeometry)(Element::GeometryType::Pointer) const = &BRep::CutStatus;
    int(BRep::*pointer_to_CutStatusBySamplingElement)(Element::Pointer, const std::size_t&) const = &BRep::CutStatusBySampling;
    int(BRep::*pointer_to_CutStatusBySamplingGeometry)(Element::GeometryType::Pointer, const std::size_t&) const = &BRep::CutStatusBySampling;

    class_<BRep, BRep::Pointer, boost::noncopyable>
    ( "BRep", init<>() )
    .def("SetTolerance", &BRep::SetTolerance)
    .def("GetTolerance", &BRep::GetTolerance)
    .def("CutStatus", pointer_to_CutStatusElement)
    .def("CutStatus", pointer_to_CutStatusGeometry)
    .def("CutStatusBySampling", pointer_to_CutStatusBySamplingElement)
    .def("CutStatusBySampling", pointer_to_CutStatusBySamplingGeometry)
    .def_readonly("_CUT", &BRep::_CUT)
    .def_readonly("_IN", &BRep::_IN)
    .def_readonly("_OUT", &BRep::_OUT)
    ;

    class_<ParametricCurve, ParametricCurve::Pointer, boost::noncopyable, bases<FunctionR1R3> >
    ("ParametricCurve", init<const FunctionR1R1::Pointer, const FunctionR1R1::Pointer, const FunctionR1R1::Pointer>())
    ;

    class_<ParametricSurface, ParametricSurface::Pointer, boost::noncopyable, bases<FunctionR2R3> >
    ("ParametricSurface", init<const FunctionR2R1::Pointer, const FunctionR2R1::Pointer, const FunctionR2R1::Pointer>())
    ;

    class_<ParametricVolume, ParametricVolume::Pointer, boost::noncopyable, bases<FunctionR3R3> >
    ("ParametricVolume", init<const FunctionR3R1::Pointer, const FunctionR3R1::Pointer, const FunctionR3R1::Pointer>())
    ;

    /**************************************************************/
    /************* EXPORT INTERFACE FOR LEVEL SET *****************/
    /**************************************************************/

    double(LevelSet::*LevelSet_pointer_to_GetValue)(const LevelSet::PointType&) const = &LevelSet::GetValue;

    class_<LevelSet, LevelSet::Pointer, boost::noncopyable, bases<FunctionR3R1, BRep> >
    ( "LevelSet", init<>() )
    .def("GetValue", LevelSet_pointer_to_GetValue)
    .def(self_ns::str(self))
    ;

    class_<CircularLevelSet, CircularLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "CircularLevelSet", init<const double&, const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<SphericalLevelSet, SphericalLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "SphericalLevelSet", init<const double&, const double&, const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<DoughnutLevelSet, DoughnutLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "DoughnutLevelSet", init<const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<CylinderLevelSet, CylinderLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "CylinderLevelSet", init<const double&, const double&, const double&, const double&, const double&, const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<LinearLevelSet, LinearLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "LinearLevelSet", init<const double&, const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<PlanarLevelSet, PlanarLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "PlanarLevelSet", init<const double&, const double&, const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<ProductLevelSet, ProductLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "ProductLevelSet", init<const LevelSet::Pointer, const LevelSet::Pointer>() )
    .def(self_ns::str(self))
    ;

    class_<InverseLevelSet, InverseLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "InverseLevelSet", init<const LevelSet::Pointer>() )
    .def(self_ns::str(self))
    ;

    class_<UnionLevelSet, UnionLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "UnionLevelSet", init<const LevelSet::Pointer, const LevelSet::Pointer>() )
    .def(self_ns::str(self))
    ;

    class_<IntersectionLevelSet, IntersectionLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "IntersectionLevelSet", init<const LevelSet::Pointer, const LevelSet::Pointer>() )
    .def(self_ns::str(self))
    ;

    class_<DifferenceLevelSet, DifferenceLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "DifferenceLevelSet", init<const LevelSet::Pointer, const LevelSet::Pointer>() )
    .def(self_ns::str(self))
    ;

    /**************************************************************/
    /************* EXPORT INTERFACE FOR PARTICULAR BREP ***********/
    /**************************************************************/

    class_<AndBRep, AndBRep::Pointer, bases<BRep>, boost::noncopyable>
    ( "AndBRep", init<BRep::Pointer, BRep::Pointer>() )
    ;

}
}  // namespace Python.
}  // namespace Kratos.

