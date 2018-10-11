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



// Project includes
#include "custom_utilities/finite_cell_mesh_utility.h"


namespace Kratos
{

void FiniteCellMeshUtility::GenerateSampling(std::vector<double>& sampling,
    const double& s_min, const double& s_max, const std::size_t& nsampling)
{
    if(sampling.size() != nsampling + 1)
        sampling.resize(nsampling + 1);

    double ds = (s_max - s_min) / nsampling;
    for (std::size_t i = 0; i < nsampling + 1; ++i)
        sampling[i] = s_min + i*ds;
}


void FiniteCellMeshUtility::GenerateSampling(std::vector<double>& sampling,
        const double& s_min, const double& s_max,
        const double& w1, const double& w2,
        const std::size_t& nsampling)
{
    if(sampling.size() != nsampling + 1)
        sampling.resize(nsampling + 1);

    double t, dt = 1.0 / nsampling, N1, N2;
    for (std::size_t i = 0; i < nsampling + 1; ++i)
    {
        t = i*dt;
        N1 = 1.0 - t;
        N2 = t;
        sampling[i] = w1*N1/(w1*N1+w2*N2)*s_min + w2*N2/(w1*N1+w2*N2)*s_max;
    }
}


void FiniteCellMeshUtility::GenerateStructuredPoints2D(std::vector<std::vector<PointType> >& sampling_points,
    const int& type,
    const PointType& StartPoint,
    const PointType& EndPoint,
    const std::vector<std::size_t>& nsampling)
{
    if (type == 1)
    {
        GenerateStructuredPoints_Helper<2, 1>::Execute(sampling_points, StartPoint, EndPoint, nsampling);
    }
    else if (type == 2)
    {
        GenerateStructuredPoints_Helper<2, 2>::Execute(sampling_points, StartPoint, EndPoint, nsampling);
    }
    else if (type == 3)
    {
        GenerateStructuredPoints_Helper<2, 3>::Execute(sampling_points, StartPoint, EndPoint, nsampling);
    }
}


/// Generate the points for background structure mesh
void FiniteCellMeshUtility::GenerateStructuredPoints2D(std::vector<std::vector<PointType> >& sampling_points,
    const int& type,
    const PointType& StartPoint,
    const PointType& EndPoint,
    const std::vector<std::vector<double> >& sampling)
{
    if (type == 1)
    {
        GenerateStructuredPoints_Helper<2, 1>::Execute(sampling_points, StartPoint, EndPoint, sampling);
    }
    else if (type == 2)
    {
        GenerateStructuredPoints_Helper<2, 2>::Execute(sampling_points, StartPoint, EndPoint, sampling);
    }
    else if (type == 3)
    {
        GenerateStructuredPoints_Helper<2, 3>::Execute(sampling_points, StartPoint, EndPoint, sampling);
    }
}


void FiniteCellMeshUtility::GenerateStructuredPoints3D(std::vector<std::vector<std::vector<PointType> > >& sampling_points,
    const int& type,
    const PointType& StartPoint,
    const PointType& EndPoint,
    const std::vector<std::size_t>& nsampling)
{
    if (type == 1)
    {
        GenerateStructuredPoints_Helper<3, 1>::Execute(sampling_points, StartPoint, EndPoint, nsampling);
    }
    else if (type == 2)
    {
        GenerateStructuredPoints_Helper<3, 2>::Execute(sampling_points, StartPoint, EndPoint, nsampling);
    }
    else if (type == 3)
    {
        GenerateStructuredPoints_Helper<3, 3>::Execute(sampling_points, StartPoint, EndPoint, nsampling);
    }
}


void FiniteCellMeshUtility::GenerateStructuredPoints3D(std::vector<std::vector<std::vector<PointType> > >& sampling_points,
    const int& type,
    const PointType& StartPoint,
    const PointType& EndPoint,
    const std::vector<std::vector<double> >& sampling)
{
    if (type == 1)
    {
        GenerateStructuredPoints_Helper<3, 1>::Execute(sampling_points, StartPoint, EndPoint, sampling);
    }
    else if (type == 2)
    {
        GenerateStructuredPoints_Helper<3, 2>::Execute(sampling_points, StartPoint, EndPoint, sampling);
    }
    else if (type == 3)
    {
        GenerateStructuredPoints_Helper<3, 3>::Execute(sampling_points, StartPoint, EndPoint, sampling);
    }
}


FiniteCellMeshUtility::MeshInfoType FiniteCellMeshUtility::CreateLineElements(ModelPart& r_model_part,
        const std::vector<PointType>& sampling_points,
        const std::string& sample_element_name,
        const int& type, // if 1: generate L2 elements; 2: L3 elements;
        const bool& close, // if false: open loop; true: close loop
        Properties::Pointer pProperties)
{
    std::size_t last_node_id = FiniteCellAuxilliaryUtility::GetLastNodeId(r_model_part);
    std::size_t last_node_id_old = last_node_id;
    std::size_t num_division_1 = sampling_points.size() - 1;
    // KRATOS_WATCH(last_node_id)

    std::size_t num_1, num_2;
    if (close)
    {
        num_1 = num_division_1 + 1;
    }
    else
    {
        num_1 = num_division_1;
    }

    // firstly create nodes and add to model_part
    ModelPart::NodesContainerType NewNodes;
    for (std::size_t i = 0; i < num_1; ++i)
    {
        NodeType::Pointer pNewNode = r_model_part.CreateNewNode(++last_node_id,
                sampling_points[i][0], sampling_points[i][1], sampling_points[i][2]);
        // std::cout << "node " << last_node_id << " is created at " << pNewNode->X0() << " " << pNewNode->Y0() << " " << pNewNode->Z0() << std::endl;
        NewNodes.push_back(pNewNode);
    }

    // secondly create elements
    std::size_t last_element_id = FiniteCellAuxilliaryUtility::GetLastElementId(r_model_part);
    // KRATOS_WATCH(last_element_id)
    Element const& rCloneElement = KratosComponents<Element>::Get(sample_element_name);
    Element::NodesArrayType temp_element_nodes;
    ModelPart::ElementsContainerType NewElements;
    const std::string NodeKey("Node");
    std::vector<std::size_t> node;
    int activation_level;

    if (type == 1)
        node.resize(2);
    else if (type == 2)
        node.resize(3);
    else
        KRATOS_THROW_ERROR(std::logic_error, "Invalid type", type)

    BoundaryLayerInfoType boundary_layers;
    BoundaryNodesInfoType boundary_nodes;

    for (std::size_t i = 0; i < num_1; ++i)
    {
        temp_element_nodes.clear();

        if (type == 1)
        {
            node[0] = last_node_id_old + i + 1;
            if (i < num_division_1)
            {
                node[1] = last_node_id_old + i + 2;
            }
            else
            {
                node[1] = last_node_id_old + 1;
            }

            temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node[0], NodeKey).base()));
            temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node[1], NodeKey).base()));
        }
        else if (type == 2)
        {
            // TODO
            KRATOS_THROW_ERROR(std::logic_error, "type == 2", "is not yet implemented")
        }

        Element::Pointer pNewElement = rCloneElement.Create(++last_element_id, temp_element_nodes, pProperties);
        // std::cout << "element " << pNewElement->Id() << " is created" << std::endl;
        pNewElement->Set(ACTIVE, true);
        pNewElement->SetValue(IS_INACTIVE, false);
        NewElements.push_back(pNewElement);
    }

    for (ModelPart::ElementsContainerType::ptr_iterator it = NewElements.ptr_begin(); it != NewElements.ptr_end(); ++it)
    {
        r_model_part.Elements().push_back(*it);
    }

    r_model_part.Elements().Unique();

    std::cout << NewElements.size() << " " << sample_element_name << " elements are created and added to the model_part" << std::endl;

    return std::make_tuple(NewNodes, NewElements, boundary_nodes, boundary_layers);
}


FiniteCellMeshUtility::MeshInfoType FiniteCellMeshUtility::CreateQuadElements(ModelPart& r_model_part,
    const std::vector<std::vector<PointType> >& sampling_points,
    const std::string& sample_element_name,
    const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
    const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
    const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
    Properties::Pointer pProperties)
{
    std::size_t last_node_id = FiniteCellAuxilliaryUtility::GetLastNodeId(r_model_part);
    std::size_t last_node_id_old = last_node_id;
    std::size_t num_division_1 = sampling_points.size() - 1;
    std::size_t num_division_2 = sampling_points[0].size() - 1;
    // KRATOS_WATCH(last_node_id)

    std::size_t num_1, num_2;
    if (close_dir == 1)
    {
        num_1 = num_division_1 + 1;
        num_2 = num_division_2;
    }
    else if (close_dir == 2)
    {
        num_1 = num_division_1;
        num_2 = num_division_2 + 1;
    }
    else
    {
        num_1 = num_division_1;
        num_2 = num_division_2;
    }

    Variable<int>& ACTIVATION_LEVEL_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("ACTIVATION_LEVEL"));

    // firstly create nodes and add to model_part
    ModelPart::NodesContainerType NewNodes;
    for (std::size_t i = 0; i < num_division_1+1; ++i)
    {
        for (std::size_t j = 0; j < num_division_2+1; ++j)
        {
            NodeType::Pointer pNewNode = r_model_part.CreateNewNode(++last_node_id,
                    sampling_points[i][j][0], sampling_points[i][j][1], sampling_points[i][j][2]);
            // std::cout << "node " << last_node_id << " is created at " << pNewNode->X0() << " " << pNewNode->Y0() << " " << pNewNode->Z0() << std::endl;
            NewNodes.push_back(pNewNode);
        }
    }

    // secondly create elements
    std::size_t last_element_id = FiniteCellAuxilliaryUtility::GetLastElementId(r_model_part);
    // KRATOS_WATCH(last_element_id)
    Element const& rCloneElement = KratosComponents<Element>::Get(sample_element_name);
    Element::NodesArrayType temp_element_nodes;
    ModelPart::ElementsContainerType NewElements;
    const std::string NodeKey("Node");
    std::vector<std::size_t> node;
    int activation_level;

    if (type == 1)
        node.resize(4);
    else if (type == 2)
        node.resize(8);
    else if (type == 3)
        node.resize(9);
    else
        KRATOS_THROW_ERROR(std::logic_error, "Invalid type", type)

    if (activation_dir == 1) activation_level = -num_division_1;

    BoundaryLayerInfoType boundary_layers;
    BoundaryNodesInfoType boundary_nodes;

    for (std::size_t i = 0; i < num_1; ++i)
    {
        if (activation_dir == 2) activation_level = -num_division_2;
        for (std::size_t j = 0; j < num_2; ++j)
        {
            temp_element_nodes.clear();

            if (type == 1)
            {
                node[0] = last_node_id_old + i * (num_division_2 + 1) + j + 1;
                node[2] = last_node_id_old + (i + 1) * (num_division_2 + 1) + j + 1;
                if (j < num_division_2)
                {
                    node[1] = last_node_id_old + i * (num_division_2 + 1) + j + 2;
                    node[3] = last_node_id_old + (i + 1) * (num_division_2 + 1) + j + 2;
                }
                else
                {
                    node[1] = last_node_id_old + i * (num_division_2 + 1) + 1;
                    node[3] = last_node_id_old + (i + 1) * (num_division_2 + 1) + 1;
                }
                // std::cout << node[0] << " " << node[1] << " " << node[2] << " " << node[3] << std::endl;

                temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node[0], NodeKey).base()));
                temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node[1], NodeKey).base()));
                temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node[3], NodeKey).base()));
                temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node[2], NodeKey).base()));
            }
            else if (type == 2)
            {
                // TODO
                KRATOS_THROW_ERROR(std::logic_error, "type == 2", "is not yet implemented")
            }
            else if (type == 3)
            {
                // TODO
                KRATOS_THROW_ERROR(std::logic_error, "type == 3", "is not yet implemented")
            }

            Element::Pointer pNewElement = rCloneElement.Create(++last_element_id, temp_element_nodes, pProperties);
            // std::cout << "element " << pNewElement->Id() << " is created" << std::endl;
            pNewElement->Set(ACTIVE, true);
            pNewElement->SetValue(IS_INACTIVE, false);
            if (activation_dir != 0)
                pNewElement->SetValue(ACTIVATION_LEVEL_var, activation_level);
            NewElements.push_back(pNewElement);

            if (activation_dir == 2) ++activation_level;
        }
        if (activation_dir == 1) ++activation_level;
    }

    for (ModelPart::ElementsContainerType::ptr_iterator it = NewElements.ptr_begin(); it != NewElements.ptr_end(); ++it)
    {
        r_model_part.Elements().push_back(*it);
    }

    r_model_part.Elements().Unique();

    std::cout << NewElements.size() << " " << sample_element_name << " elements are created and added to the model_part" << std::endl;

    return std::make_tuple(NewNodes, NewElements, boundary_nodes, boundary_layers);
}


FiniteCellMeshUtility::MeshInfoType FiniteCellMeshUtility::CreateHexElements(ModelPart& r_model_part,
    const std::vector<std::vector<std::vector<PointType> > >& sampling_points,
    const std::string& sample_element_name,
    const int& type, // if 1: generate H8 elements; 2: H20 elements; 3: H27 elements
    const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir; 3: close on 3rd dir
    const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir; r: activation on 3rd dir
    Properties::Pointer pProperties)
{
    std::size_t last_node_id = FiniteCellAuxilliaryUtility::GetLastNodeId(r_model_part);
    std::size_t last_node_id_old = last_node_id;
    std::size_t num_division_1, num_division_2, num_division_3;
    std::size_t num_1, num_2, num_3;

    if (type == 1)
    {
        num_division_1 = sampling_points.size() - 1;
        num_division_2 = sampling_points[0].size() - 1;
        num_division_3 = sampling_points[0][0].size() - 1;
        num_1 = num_division_1;
        num_2 = num_division_2;
        num_3 = num_division_3;
    }
    else if (type == 2)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
    }
    else if (type == 3)
    {
        num_division_1 = (sampling_points.size() - 1) / 2;
        num_division_2 = (sampling_points[0].size() - 1) / 2;
        num_division_3 = (sampling_points[0][0].size() - 1) / 2;
        num_1 = 2 * num_division_1;
        num_2 = 2 * num_division_2;
        num_3 = 2 * num_division_3;
    }
    // KRATOS_WATCH(last_node_id)

    // firstly create nodes and add to model_part
    ModelPart::NodesContainerType NewNodes;
    for (std::size_t i = 0; i < num_1+1; ++i)
    {
        for (std::size_t j = 0; j < num_2+1; ++j)
        {
            for (std::size_t k = 0; k < num_3+1; ++k)
            {
                NodeType::Pointer pNewNode = r_model_part.CreateNewNode(++last_node_id,
                        sampling_points[i][j][k][0], sampling_points[i][j][k][1], sampling_points[i][j][k][2]);
                // std::cout << "node " << last_node_id << " is created at " << pNewNode->X0() << " " << pNewNode->Y0() << " " << pNewNode->Z0() << std::endl;
                NewNodes.push_back(pNewNode);
            }
        }
    }

    // secondly create elements
    std::size_t last_element_id = FiniteCellAuxilliaryUtility::GetLastElementId(r_model_part);
    // KRATOS_WATCH(last_element_id)
    Element const& rCloneElement = KratosComponents<Element>::Get(sample_element_name);
    Element::NodesArrayType temp_element_nodes;
    ModelPart::ElementsContainerType NewElements;
    const std::string NodeKey("Node");
    std::vector<std::size_t> node;
    int activation_level;

    if (type == 1)
        node.resize(8);
    else if (type == 2)
        node.resize(20);
    else if (type == 3)
        node.resize(27);
    else
        KRATOS_THROW_ERROR(std::logic_error, "Invalid type", type)

//    KRATOS_WATCH(last_node_id_old)
    KRATOS_WATCH(num_division_1)
    KRATOS_WATCH(num_division_2)
    KRATOS_WATCH(num_division_3)

    BoundaryLayerInfoType boundary_layers;
    BoundaryNodesInfoType boundary_nodes;

    for (std::size_t i = 0; i < num_division_1; ++i)
    {
        for (std::size_t j = 0; j < num_division_2; ++j)
        {
            for (std::size_t k = 0; k < num_division_3; ++k)
            {
                temp_element_nodes.clear();

                if (type == 1)
                {
                    node[0] = last_node_id_old + (i * (num_2 + 1) + j) * (num_3 + 1) + k + 1;
                    node[1] = last_node_id_old + (i * (num_2 + 1) + j + 1) * (num_3 + 1) + k + 1;
                    node[2] = last_node_id_old + ((i + 1) * (num_2 + 1) + j + 1) * (num_3 + 1) + k + 1;
                    node[3] = last_node_id_old + ((i + 1) * (num_2 + 1) + j) * (num_3 + 1) + k + 1;
                    node[4] = node[0] + 1;
                    node[5] = node[1] + 1;
                    node[6] = node[2] + 1;
                    node[7] = node[3] + 1;

                    for (int n = 0; n < 8; ++n)
                    {
//                        std::cout << "node " << n << ": " << node[n] << std::endl;
                        temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node[n], NodeKey).base()));
                    }

                    // extract the layer information
                    if (k == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[3], node[2], node[1]};
                        boundary_layers["xmin"].push_back(layer_cond);
                        boundary_nodes["xmin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (k == num_division_3-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[4], node[5], node[6], node[7]};
                        boundary_layers["xmax"].push_back(layer_cond);
                        boundary_nodes["xmax"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (j == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[4], node[7], node[3]};
                        boundary_layers["ymin"].push_back(layer_cond);
                        boundary_nodes["ymin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (j == num_division_2-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[1], node[2], node[6], node[5]};
                        boundary_layers["ymax"].push_back(layer_cond);
                        boundary_nodes["ymax"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (i == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[1], node[5], node[4]};
                        boundary_layers["zmin"].push_back(layer_cond);
                        boundary_nodes["zmin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (i == num_division_1-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[2], node[3], node[7], node[6]};
                        boundary_layers["zmax"].push_back(layer_cond);
                        boundary_nodes["zmax"].insert(layer_cond.begin(), layer_cond.end());
                    }
                }
                else if (type == 2)
                {
                    node[0] = last_node_id_old + 2*((i * (num_2 + 1) + j) * (num_3 + 1) + k) + 1;
                    node[1] = last_node_id_old + 2*((i * (num_2 + 1) + j + 1) * (num_3 + 1) + k) + 1;
                    node[2] = last_node_id_old + 2*(((i + 1) * (num_2 + 1) + j + 1) * (num_3 + 1) + k) + 1;
                    node[3] = last_node_id_old + 2*(((i + 1) * (num_2 + 1) + j) * (num_3 + 1) + k) + 1;
                    node[4] = node[0] + 2;
                    node[5] = node[1] + 2;
                    node[6] = node[2] + 2;
                    node[7] = node[3] + 2;

                    node[8] = (node[0] + node[1]) / 2;
                    node[9] = (node[1] + node[2]) / 2;
                    node[10] = (node[2] + node[3]) / 2;
                    node[11] = (node[0] + node[3]) / 2;

                    node[12] = (node[0] + node[4]) / 2;
                    node[13] = (node[1] + node[5]) / 2;
                    node[14] = (node[2] + node[6]) / 2;
                    node[15] = (node[3] + node[7]) / 2;

                    node[16] = (node[4] + node[5]) / 2;
                    node[17] = (node[5] + node[6]) / 2;
                    node[18] = (node[6] + node[7]) / 2;
                    node[19] = (node[4] + node[7]) / 2;

                    for (int n = 0; n < 20; ++n)
                    {
//                        std::cout << "node " << n << ": " << node[n] << std::endl;
                        temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node[n], NodeKey).base()));
                    }

//                    std::cout << "element " << i << " " << j << " " << k << ":" << std::endl;
//                    for (int n = 0; n < 27; ++n)
//                        std::cout << " " << node[n];
//                    std::cout << std::endl;

                    // extract the layer information
                    if (k == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[3], node[2], node[1], node[11], node[10], node[9], node[8]};
                        boundary_layers["xmin"].push_back(layer_cond);
                        boundary_nodes["xmin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (k == num_division_3-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[4], node[5], node[6], node[7], node[16], node[17], node[18], node[19]};
                        boundary_layers["xmax"].push_back(layer_cond);
                        boundary_nodes["xmax"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (j == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[4], node[7], node[3], node[12], node[19], node[15], node[11]};
                        boundary_layers["ymin"].push_back(layer_cond);
                        boundary_nodes["ymin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (j == num_division_2-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[1], node[2], node[6], node[5], node[9], node[14], node[17], node[13]};
                        boundary_layers["ymax"].push_back(layer_cond);
                        boundary_nodes["ymax"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (i == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[1], node[5], node[4], node[8], node[13], node[16], node[12]};
                        boundary_layers["zmin"].push_back(layer_cond);
                        boundary_nodes["zmin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (i == num_division_1-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[2], node[3], node[7], node[6], node[10], node[15], node[18], node[14]};
                        boundary_layers["zmax"].push_back(layer_cond);
                        boundary_nodes["zmax"].insert(layer_cond.begin(), layer_cond.end());
                    }
                }
                else if (type == 3)
                {
                    node[0] = last_node_id_old + 2*((i * (num_2 + 1) + j) * (num_3 + 1) + k) + 1;
                    node[1] = last_node_id_old + 2*((i * (num_2 + 1) + j + 1) * (num_3 + 1) + k) + 1;
                    node[2] = last_node_id_old + 2*(((i + 1) * (num_2 + 1) + j + 1) * (num_3 + 1) + k) + 1;
                    node[3] = last_node_id_old + 2*(((i + 1) * (num_2 + 1) + j) * (num_3 + 1) + k) + 1;
                    node[4] = node[0] + 2;
                    node[5] = node[1] + 2;
                    node[6] = node[2] + 2;
                    node[7] = node[3] + 2;

                    node[8] = (node[0] + node[1]) / 2;
                    node[9] = (node[1] + node[2]) / 2;
                    node[10] = (node[2] + node[3]) / 2;
                    node[11] = (node[0] + node[3]) / 2;

                    node[12] = (node[0] + node[4]) / 2;
                    node[13] = (node[1] + node[5]) / 2;
                    node[14] = (node[2] + node[6]) / 2;
                    node[15] = (node[3] + node[7]) / 2;

                    node[16] = (node[4] + node[5]) / 2;
                    node[17] = (node[5] + node[6]) / 2;
                    node[18] = (node[6] + node[7]) / 2;
                    node[19] = (node[4] + node[7]) / 2;

                    node[20] = (node[0] + node[1] + node[2] + node[3]) / 4;
                    node[21] = (node[0] + node[1] + node[4] + node[5]) / 4;
                    node[22] = (node[1] + node[2] + node[5] + node[6]) / 4;
                    node[23] = (node[2] + node[3] + node[6] + node[7]) / 4;
                    node[24] = (node[0] + node[3] + node[4] + node[7]) / 4;
                    node[25] = (node[4] + node[5] + node[6] + node[7]) / 4;

                    node[26] = (node[0] + node[1] + node[2] + node[3] + node[4] + node[5] + node[6] + node[7]) / 8;

                    for (int n = 0; n < 27; ++n)
                    {
//                        std::cout << "node " << n << ": " << node[n] << std::endl;
                        temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), node[n], NodeKey).base()));
                    }

//                    std::cout << "element " << i << " " << j << " " << k << ":" << std::endl;
//                    for (int n = 0; n < 27; ++n)
//                        std::cout << " " << node[n];
//                    std::cout << std::endl;

                    // extract the layer information
                    if (k == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[3], node[2], node[1], node[11], node[10], node[9], node[8], node[20]};
                        boundary_layers["xmin"].push_back(layer_cond);
                        boundary_nodes["xmin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (k == num_division_3-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[4], node[5], node[6], node[7], node[16], node[17], node[18], node[19], node[25]};
                        boundary_layers["xmax"].push_back(layer_cond);
                        boundary_nodes["xmax"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (j == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[4], node[7], node[3], node[12], node[19], node[15], node[11], node[24]};
                        boundary_layers["ymin"].push_back(layer_cond);
                        boundary_nodes["ymin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (j == num_division_2-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[1], node[2], node[6], node[5], node[9], node[14], node[17], node[13], node[22]};
                        boundary_layers["ymax"].push_back(layer_cond);
                        boundary_nodes["ymax"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (i == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[1], node[5], node[4], node[8], node[13], node[16], node[12], node[21]};
                        boundary_layers["zmin"].push_back(layer_cond);
                        boundary_nodes["zmin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (i == num_division_1-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[2], node[3], node[7], node[6], node[10], node[15], node[18], node[14], node[23]};
                        boundary_layers["zmax"].push_back(layer_cond);
                        boundary_nodes["zmax"].insert(layer_cond.begin(), layer_cond.end());
                    }
                }

                Element::Pointer pNewElement = rCloneElement.Create(++last_element_id, temp_element_nodes, pProperties);
                // std::cout << "element " << pNewElement->Id() << " is created" << std::endl;
                pNewElement->Set(ACTIVE, true);
                pNewElement->SetValue(IS_INACTIVE, false);
                NewElements.push_back(pNewElement);
            }
        }
    }

    for (ModelPart::ElementsContainerType::ptr_iterator it = NewElements.ptr_begin(); it != NewElements.ptr_end(); ++it)
    {
        r_model_part.Elements().push_back(*it);
    }

    r_model_part.Elements().Unique();

    std::cout << NewElements.size() << " " << sample_element_name << " elements are created and added to the model_part" << std::endl;

    return std::make_tuple(NewNodes, NewElements, boundary_nodes, boundary_layers);
}


ModelPart::NodesContainerType::iterator FiniteCellMeshUtility::FindKey(ModelPart::NodesContainerType& ThisContainer,
        std::size_t& ThisKey, const std::string& ComponentName)
{
    ModelPart::NodesContainerType::iterator i_result;
    if((i_result = ThisContainer.find(ThisKey)) == ThisContainer.end())
    {
        std::stringstream buffer;
        buffer << ComponentName << " #" << ThisKey << " is not found.";
        KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
    }

    return i_result;
}


Element::Pointer FiniteCellMeshUtility::CreateParasiteElement(const std::string& sample_element_name,
    std::size_t& lastElementId,
    Element::Pointer pElement, Properties::Pointer pProperties )
{
    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);

    // REMARK: when creating the element here, the integration rule is not passed. Instead the default integration rule of this element_type is applied, which is not the same as the original element.
    Element::Pointer pNewElement = r_clone_element.Create(++lastElementId, pElement->pGetGeometry(), pProperties);

    std::cout << "1 element of type " << sample_element_name << " is created" << std::endl;

    return pNewElement;
}


Element::Pointer FiniteCellMeshUtility::CreateParasiteElement(Element::Pointer pElement, // the parent element keep the geometry
    const std::string& sample_element_name,
    const int& RepresentativeIntegrationOrder,
    const GeometryType::IntegrationPointsArrayType& integration_points,
    std::size_t& lastElementId,
    Properties::Pointer pProperties)
{
    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);
    GeometryType& r_geom = *(pElement->pGetGeometry());

    GeometryData::IntegrationMethod RepresentativeIntegrationMethod
            = Function<double, double>::GetIntegrationMethod(RepresentativeIntegrationOrder);
    Variable<int>& INTEGRATION_ORDER_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("INTEGRATION_ORDER"));

    // create the new elements from sub-cell
    // here we make a clone of the geometry because we want to assign different geometry data later on
    // this also works with Bezier element, because Bezier geometry has implemented the Create method
    Element::Pointer pNewElement;
    pNewElement = r_clone_element.Create(++lastElementId, r_geom.Create(r_geom.Points()), pProperties);

    FiniteCellGeometryUtility::AssignGeometryData(pNewElement->GetGeometry(), RepresentativeIntegrationMethod, integration_points);
    pNewElement->SetValue(INTEGRATION_ORDER_var, RepresentativeIntegrationOrder);
    pNewElement->Initialize();

    return pNewElement;
}


ModelPart::NodesContainerType FiniteCellMeshUtility::ImportNodes(ModelPart& rThisModelPart, ModelPart& rOtherModelPart)
{
    std::size_t last_node_id = FiniteCellAuxilliaryUtility::GetLastNodeId(rThisModelPart);

    // create nodes and add to model_part
    ModelPart::NodesContainerType NewNodes;
    for(ModelPart::NodeIterator it = rOtherModelPart.NodesBegin(); it != rOtherModelPart.NodesEnd(); ++it)
    {
        std::size_t new_node_id = ++last_node_id;
        NodeType::Pointer pNewNode = rThisModelPart.CreateNewNode(new_node_id, it->X(), it->Y(), it->Z());
        it->SetValue(OTHER_NODE_ID, new_node_id);
        NewNodes.push_back(pNewNode);
    }

    std::cout << NewNodes.size() << " nodes from " << rOtherModelPart.Name()
              << " are added to the model_part " << rThisModelPart.Name() << std::endl;
    return NewNodes;
}


ModelPart::ElementsContainerType FiniteCellMeshUtility::ImportElements(ModelPart& rThisModelPart,
    ModelPart::ElementsContainerType& rOtherElements,
    const std::string& sample_element_name, Properties::Pointer pProperties)
{
    std::size_t last_element_id = FiniteCellAuxilliaryUtility::GetLastElementId(rThisModelPart);
    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);

    ModelPart::ElementsContainerType NewElements;
    ImportEntities(rThisModelPart, rThisModelPart.Elements(), NewElements, rOtherElements, last_element_id, r_clone_element, pProperties);

    std::cout << NewElements.size() << " " << sample_element_name
              << " elements are created and added to the model_part " << rThisModelPart.Name() << std::endl;

    return NewElements;
}


ModelPart::ConditionsContainerType FiniteCellMeshUtility::ImportConditions(ModelPart& rThisModelPart,
    ModelPart::ConditionsContainerType& rOtherConditions,
    const std::string& sample_cond_name, Properties::Pointer pProperties)
{
    std::size_t last_cond_id = FiniteCellAuxilliaryUtility::GetLastConditionId(rThisModelPart);
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_cond_name);

    ModelPart::ConditionsContainerType NewConditions;
    ImportEntities(rThisModelPart, rThisModelPart.Conditions(), NewConditions, rOtherConditions, last_cond_id, r_clone_condition, pProperties);

    std::cout << NewConditions.size() << " " << sample_cond_name
              << " conditions are created and added to the model_part " << rThisModelPart.Name() << std::endl;

    return NewConditions;
}

}  // namespace Kratos.

