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
//  Date:            23 Feb 2017
//


#if !defined(KRATOS_QUADRATURE_UTILITY_H_INCLUDED )
#define  KRATOS_QUADRATURE_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <fstream>


// External includes
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/numeric/ublas/lu.hpp> 


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "custom_geometries/finite_cell_geometry.h"
#include "custom_algebra/function.h"


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
/** Abstract class to compute/assign the quadrature weights for elements/conditions
*/
class QuadratureUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of QuadratureUtility
    KRATOS_CLASS_POINTER_DEFINITION(QuadratureUtility);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    QuadratureUtility() {}

    /// Destructor.
    virtual ~QuadratureUtility() {}


    void PyScaleQuadrature(Element::Pointer& p_elem,
            const int quadrature_order, const double ScaleFactor) const
    {
        GeometryData::IntegrationMethod ElementalIntegrationMethod
                = Function<double, double>::GetIntegrationMethod(quadrature_order);

        this->ScaleQuadrature(p_elem->GetGeometry(), ElementalIntegrationMethod, ScaleFactor);
    }


    void ScaleQuadrature(GeometryType& r_geom,
            const GeometryData::IntegrationMethod& ElementalIntegrationMethod,
            const double ScaleFactor) const
    {
        const GeometryType::IntegrationPointsArrayType& integration_points
                = r_geom.IntegrationPoints( ElementalIntegrationMethod );

        Vector Mw(integration_points.size());
        for(std::size_t i = 0; i < integration_points.size(); ++i)
        {
            Mw(i) = integration_points[i].Weight() * ScaleFactor;
        }

        /* create new quadrature and assign to the geometry */
        FiniteCellGeometry<GeometryType>::AssignGeometryData(r_geom, ElementalIntegrationMethod, Mw);
    }


    void PySaveQuadrature(boost::python::list& pyElemList, const std::string& fileName,
            const std::string& writeMode) const
    {
        std::ofstream myFile;

        int mode;
        if(writeMode == "binary")
        {
            myFile.open(fileName.c_str(), std::ios::out | std::ios::binary);
            mode = 0;
            myFile.write((char*)&mode, sizeof(int));
        }
        else if(writeMode == "binary-full")
        {
            myFile.open(fileName.c_str(), std::ios::out | std::ios::binary);
            mode = 1;
            myFile.write((char*)&mode, sizeof(int));
        }
        else if(writeMode == "ascii")
        {
            myFile.open(fileName.c_str(), std::ios::out);
            mode = 2;
            myFile << mode << std::endl;
        }
        else if(writeMode == "ascii-full")
        {
            myFile.open(fileName.c_str(), std::ios::out);
            mode = 3;
            myFile << mode << std::endl;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Unknown write mode", writeMode)

        myFile.precision(15);

        typedef boost::python::stl_input_iterator<Element::Pointer> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& p_elem,
                        std::make_pair(iterator_value_type(pyElemList), // begin
                        iterator_value_type() ) ) // end
        {
            switch(mode)
            {
                case 0:
                    SaveQuadrature<0>(myFile, p_elem, p_elem->GetGeometry().GetDefaultIntegrationMethod());
                    break;
                case 1:
                    SaveQuadrature<1>(myFile, p_elem, p_elem->GetGeometry().GetDefaultIntegrationMethod());
                    break;
                case 2:
                    SaveQuadrature<2>(myFile, p_elem, p_elem->GetGeometry().GetDefaultIntegrationMethod());
                    break;
                case 3:
                    SaveQuadrature<3>(myFile, p_elem, p_elem->GetGeometry().GetDefaultIntegrationMethod());
                    break;
            }
        }

        myFile.close();
    }


    template<int TMode>
    void SaveQuadrature(std::ofstream& rOStream, const Element::Pointer& p_elem,
            const GeometryData::IntegrationMethod& ElementalIntegrationMethod) const
    {
        const GeometryType::IntegrationPointsArrayType& integration_points
                = p_elem->GetGeometry().IntegrationPoints( ElementalIntegrationMethod );

        if(TMode == 0 || TMode == 1) // binary
        {
            std::size_t elem_id = p_elem->Id();
            std::size_t num_points = integration_points.size();
            rOStream.write((char*)&elem_id, sizeof(std::size_t));
            rOStream.write((char*)&ElementalIntegrationMethod, sizeof(int));
            rOStream.write((char*)&num_points, sizeof(std::size_t));

            double aux;
            if(TMode == 1)
            {
                for(std::size_t i = 0; i < integration_points.size(); ++i)
                {
                    aux = integration_points[i].X();
                    rOStream.write((char*)&aux, sizeof(double));

                    aux = integration_points[i].Y();
                    rOStream.write((char*)&aux, sizeof(double));

                    aux = integration_points[i].Z();
                    rOStream.write((char*)&aux, sizeof(double));
                }
            }

            for(std::size_t i = 0; i < integration_points.size(); ++i)
            {
                aux = integration_points[i].Weight();
                rOStream.write((char*)&aux, sizeof(double));
            }
        }
        else if(TMode == 2 || TMode == 3) // ascii
        {
            rOStream << p_elem->Id() << " " << ElementalIntegrationMethod
                     << " " << integration_points.size() << std::endl;

            for(std::size_t i = 0; i < integration_points.size(); ++i)
            {
                if(TMode == 3)
                {
                    rOStream << "  " << integration_points[i].X()
                             << " \t" << integration_points[i].Y()
                             << " \t" << integration_points[i].Z();
                }
                rOStream << " \t" << integration_points[i].Weight() << std::endl;
            }
        }
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


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
        return "Quadrature Utility";
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
    QuadratureUtility& operator=(QuadratureUtility const& rOther);

    /// Copy constructor.
    QuadratureUtility(QuadratureUtility const& rOther);


    ///@}

}; // Class QuadratureUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, QuadratureUtility& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const QuadratureUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_QUADRATURE_UTILITY_H_INCLUDED  defined
