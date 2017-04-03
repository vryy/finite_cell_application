import math
import pprint
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.FiniteCellApplication import *

class FiniteCellSimulator:

    def __init__(self, params):
        self.params = params

    ###COMPUTE GLOBAL DISPLACEMENT (L2) ERROR###
    def compute_L2_error(self, elements, process_info, solution, P):
        nom = 0.0
        denom = 0.0
        for element in elements:
            if element.GetValue(IS_INACTIVE) == False:
                u = element.GetValuesOnIntegrationPoints(DISPLACEMENT, process_info)
                J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
                W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                for i in range(0, len(u)):
                    ana_u = solution.get_displacement(P, Q[i][0], Q[i][1], Q[i][2])
                    nom = nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2) + pow(u[i][2] - ana_u[2], 2)) * W[i][0] * J0[i][0]
                    denom = denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2) + pow(ana_u[2], 2)) * W[i][0] * J0[i][0]
        print("nom:", nom)
        print("denom:", denom)
        if denom == 0.0:
            if nom == 0.0:
                return 0.0
            else:
                return float('nan');
        else:
            return math.sqrt(abs(nom / denom))

    ###COMPUTE GLOBAL DISPLACEMENT (H1) ERROR###
    def compute_H1_error(self, elements, process_info, solution, P):
        nom = 0.0
        denom = 0.0
        for element in elements:
            if element.GetValue(IS_INACTIVE) == False:
                o = element.GetValuesOnIntegrationPoints(THREED_STRESSES, process_info)
                J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
                W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                for i in range(0, len(o)):
                    ana_o = solution.get_stress_3d(P, Q[i][0], Q[i][1], Q[i][2])
                    nom = nom + (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2))) * W[i][0] * J0[i][0]
                    denom = denom + (pow(ana_o[0], 2) + pow(ana_o[1], 2) + pow(ana_o[2], 2) + 2.0*(pow(ana_o[3], 2) + pow(ana_o[4], 2) + pow(ana_o[5], 2))) * W[i][0] * J0[i][0]
        print("nom:", nom)
        print("denom:", denom)
        if denom == 0.0:
            if nom == 0.0:
                return 0.0
            else:
                return float('nan');
        else:
            return math.sqrt(abs(nom / denom))

    ###TREES & FOREST CREATION#############
    def CreateForest(self, elements):
        ###############################################################
        ######### CREATE TREES AND FOREST
        ###############################################################
        self.forest = {}
        for elem in elements:
            self.forest[elem.Id] = QuadTreeLocal(elem)

        #######QUAD TREE POST PROCESSING###########
#        qt_depth = self.params["qt_depth"]
#        for qi, qt in self.forest.iteritems():
#            for i in range(0, qt_depth):
#                qt.RefineBy(self.brep)
    #    sample_element = "KinematicLinear3D8N"
    #    tmodel_part = ModelPart("QuadTree")
    #    last_node_id = 0x0000000000000000 + len(model.model_part.Nodes)
    #    last_element_id = 0x0000000000000000
    #    for qi, qt in self.forest.iteritems():
    #        last_ids = qt.AddToModelPart(tmodel_part, sample_element, last_node_id, last_element_id)
    #        last_node_id = last_ids[0]
    #        last_element_id = last_ids[1]
    #    print(tmodel_part)
    #    model.gid_io.InitializeMesh( 1234.0 )
    #    mesh = tmodel_part.GetMesh()
    #    model.gid_io.WriteMesh( mesh )
    #    model.gid_io.FinalizeMesh()
    #    sys.exit(0)
        ###########################################

    def CreateForestSubCell(self, elements):
        self.forest = {}
        subcell_fit_mode = self.params["subcell_fit_mode"]
        cut_cell_quadrature_order = self.params["cut_cell_quadrature_order"]
        for elem in elements:
            if   subcell_fit_mode == "subcell fit gauss":
                self.forest[elem.Id] = MomentFittedQuadTreeSubCell(elem, "gauss", cut_cell_quadrature_order)
#            elif subcell_fit_mode == "subcell fit equal-dist":
#                self.forest[elem.Id] = MomentFittedQuadTreeSubCell(elem, "equal-dist", 3, 3)
            elif subcell_fit_mode == "subcell fit two-layer":
                self.forest[elem.Id] = MomentFittedQuadTreeSubCell(elem, "gauss", cut_cell_quadrature_order)
            elif subcell_fit_mode == "subcell nonfit":
                self.forest[elem.Id] = MomentFittedQuadTreeSubCell(elem)
            else:
                print("unknown subcell_fit_mode", subcell_fit_mode)
                sys.exit(0)

    ###FITTING FUNCTIONS#############
    def CreateFittingFunctions(self):
        fit_funcs = []
        fit_space_dim = self.params["fitting_space_dimension"]
        fit_degree = self.params["fitting_function_degree"]

        if fit_space_dim == 2:
            f_1 = ScalarFunctionR3R1(1.0)
            f_x = MonomialFunctionR3R1X()
            f_y = MonomialFunctionR3R1Y()
            f_xy = MonomialFunctionR3R1XY()

            f_x2 = MonomialFunctionR3R1X2()
            f_y2 = MonomialFunctionR3R1Y2()
            f_x2y = MonomialFunctionR3R1X2Y()
            f_xy2 = MonomialFunctionR3R1XY2()
            f_x2y2 = MonomialFunctionR3R1X2Y2()

            f_x3 = MonomialFunctionR3R1X3()
            f_y3 = MonomialFunctionR3R1Y3()
            f_x3y = MonomialFunctionR3R1X3Y()
            f_xy3 = MonomialFunctionR3R1XY3()
            f_x3y2 = MonomialFunctionR3R1X3Y2()
            f_x2y3 = MonomialFunctionR3R1X2Y3()
            f_x3y3 = MonomialFunctionR3R1X3Y3()

            if fit_degree >= 1:
                fit_funcs.append(f_1)
                fit_funcs.append(f_x)
                fit_funcs.append(f_y)
                fit_funcs.append(f_xy)

            if fit_degree >= 2:
                fit_funcs.append(f_x2)
                fit_funcs.append(f_y2)
                fit_funcs.append(f_x2y)
                fit_funcs.append(f_xy2)
                fit_funcs.append(f_x2y2)

            if fit_degree >= 3:
                fit_funcs.append(f_x3)
                fit_funcs.append(f_y3)
                fit_funcs.append(f_x3y)
                fit_funcs.append(f_xy3)
                fit_funcs.append(f_x3y2)
                fit_funcs.append(f_x2y3)
                fit_funcs.append(f_x3y3)

            if fit_degree >= 4:
                print("Unsupported fit_degree " + str(fit_degree))

        elif fit_space_dim == 3:
            f_1 = ScalarFunctionR3R1(1.0)
            f_x = MonomialFunctionR3R1X()
            f_y = MonomialFunctionR3R1Y()
            f_z = MonomialFunctionR3R1Z()
            f_xy = MonomialFunctionR3R1XY()
            f_yz = MonomialFunctionR3R1YZ()
            f_xz = MonomialFunctionR3R1XZ()
            f_xyz = MonomialFunctionR3R1XYZ()

            if fit_degree >= 1:
                fit_funcs.append(f_1)
                fit_funcs.append(f_x)
                fit_funcs.append(f_y)
                fit_funcs.append(f_z)
                fit_funcs.append(f_xy)
                fit_funcs.append(f_yz)
                fit_funcs.append(f_xz)
                fit_funcs.append(f_xyz)
            
            if fit_degree >= 2:
                print("Unsupported fit_degree " + str(fit_degree))

        return fit_funcs

    ###FITTING DRIVER#############
    def MomentFit(self, bulk_elements):
        print("fit parameters:")
        pprint.pprint(self.params)
        qt_depth = self.params["qt_depth"]
        ###########################################
        ##QUADTREE QUADRATURE
        ###########################################
        self.CreateForest(bulk_elements)
        self.elements = {}
        for elem in bulk_elements:
            self.elements[elem.Id] = elem
        cut_cell_quadrature_method = self.params["cut_cell_quadrature_method"]
        ###########################################
        ##MOMENT FIT
        ###########################################
        fit_util = MomentFittingUtility()
        fit_funcs = self.CreateFittingFunctions()

        cut_elems = []
        exclude_elems = []

        fit_mode = self.params["fit_mode"]
        integrator_quadrature_method = self.params["integrator_quadrature_method"]
        solver_type = self.params["fit_solver_type"]
        echo_level = self.params["fit_echo_level"]

        if fit_mode == "serial":
            if echo_level == -2:
                num_elems = len(self.forest)
                cnt = 0
                cnt2 = 0
                print("moment-fitting begins")
            for qi, qt in self.forest.iteritems():
                elem = self.elements[qi]
                stat = self.brep.CutStatus(elem)
                if stat == -1:
                    for i in range(0, qt_depth):
                        qt.RefineBy(self.brep)
                    fit_util.FitQuadrature(elem, fit_funcs, self.brep, qt, cut_cell_quadrature_method, integrator_quadrature_method, solver_type, echo_level)
                    cut_elems.append(elem)
                elif stat == 1:
                    exclude_elems.append(elem)
                if echo_level == -2:
                    if cnt > 0.01*num_elems:
                        cnt = 0
                        cnt2 = cnt2 + 1
                        sys.stdout.write(str(cnt2) + " ")
                        sys.stdout.flush()
                    else:
                        cnt = cnt + 1
            if echo_level == -2:
                sys.stdout.write("\n")
        elif fit_mode == "multithread":
            integrators = []
            for qi, qt in self.forest.iteritems():
                elem = self.elements[qi]
                stat = self.brep.CutStatus(elem)
                if stat == -1:
                    cut_elems.append(elem)
                    integrators.append(qt)
                elif stat == 1:
                    exclude_elems.append(elem)
            aux_util = FiniteCellAuxilliaryUtility()
            for i in range(0, qt_depth):
                aux_util.MultithreadedQuadTreeRefineBy(integrators, self.brep)
            fit_util.MultithreadedFitQuadrature(cut_elems, fit_funcs, self.brep, integrators, cut_cell_quadrature_method, integrator_quadrature_method, solver_type, echo_level)
        else:
            print("Unknown fit_mode", fit_mode)
            sys.exit(0)

        print("construct quadrature using moment fitting successfully")
        quad_filename = self.params["quad_filename"]
        quad_filetype = self.params["quad_filetype"]
        accuracy = self.params["quad_accuracy"]
        fit_util.SaveQuadrature(quad_filename, quad_filetype, cut_elems, exclude_elems, accuracy)

    def MomentFitSubCell(self, bulk_elements):
        print("fit parameters:")
        pprint.pprint(self.params)
        min_qt_depth = self.params["min_qt_depth"]
        max_qt_depth = self.params["max_qt_depth"]
        ###########################################
        ##QUADTREE QUADRATURE
        ###########################################
        self.CreateForestSubCell(bulk_elements)
        ###########################################
        ##MOMENT FIT
        ###########################################
        fit_util = MomentFittingUtility()
        fit_funcs = self.CreateFittingFunctions()

        cut_elems = []
        exclude_elems = []

#        fit_mode = self.params["fit_mode"]
        integrator_quadrature_method = self.params["integrator_quadrature_method"]
        solver_type = self.params["fit_solver_type"]
        echo_level = self.params["fit_echo_level"]

        for qi, qs in self.forest.iteritems():
            elem = qs.GetElement()
            stat = self.brep.CutStatus(elem)
            elem.SetValue(CUT_STATUS, stat)
            if stat == -1:
                cut_elems.append(qs)
            elif stat == 1:
                exclude_elems.append(qs)
        aux_util = FiniteCellAuxilliaryUtility()
        for i in range(0, min_qt_depth):
            aux_util.MultithreadedQuadTreeSubCellRefineBy(cut_elems, self.brep)
        # recursive refine the subcell that are cut but not integrable by the quadtree # this is required if Gauss-Legendre quadrature is used
#        for qs in cut_elems:
#            for i in range(0, qs.NumberOfSubCells()):
#                qt = qs.CreateQuadTree(i)
#                stat = self.brep.CutStatus(qt.pGetGeometry())
#                if stat == -1:
#                    cnt = min_qt_depth
#                    found = False
#                    while cnt < max_qt_depth:
#                        if qt.DomainSize(self.brep) == 0.0:
#                            qt.RefineBy(self.brep)
#                            cnt = cnt + 1
#                        else:
#                            found = True
#                            break
#                    if (not found) and (cnt >= max_qt_depth):
#                        print("can't find the physical point within subcell " + str(i) + " at " + str(max_qt_depth) + " quadtree level")
#                        sys.exit(0)
        fit_util.MultithreadedFitQuadratureSubCell(cut_elems, fit_funcs, self.brep, integrator_quadrature_method, solver_type, echo_level)

        print("construct quadrature using moment fitting for subcell successfully")
        quad_filename = self.params["quad_filename"]
        quad_filetype = self.params["quad_filetype"]
        accuracy = self.params["quad_accuracy"]
        fit_util.SaveQuadratureSubCell(quad_filename, quad_filetype, cut_elems, exclude_elems, accuracy)

    ###SIMULATION DRIVER#############
    def Initialize(self, model, bulk_elements):
        self.quadrature_method = self.params["quadrature_method"]
        self.mpu = self.params["material_properties_utility"]
        self.mat_type = self.params["mat_type"]
        print("simulator parameters:")
        pprint.pprint(self.params)
        quad_util = QuadratureUtility()

        if self.quadrature_method == "quadtree":
            qt_depth = self.params["qt_depth"]
            small_weight = self.params["small_weight"]
            ###########################################
            ##QUADTREE
            ###########################################
            self.CreateForest(bulk_elements)
            cut_cell_quadrature_method = self.params["cut_cell_quadrature_method"]
            cut_cell_quadrature_order = quad_util.GetQuadratureOrder(cut_cell_quadrature_method)
            total_add_quadrature_pnts = 0
            for qi, qt in self.forest.iteritems():
                elem = model.model_part.Elements[qi]
                stat = self.brep.CutStatus(elem)
                if stat == -1:
                    for i in range(0, qt_depth):
                        qt.RefineBy(self.brep)
                    np = qt.ConstructQuadrature(self.brep, cut_cell_quadrature_method, small_weight)
                    total_add_quadrature_pnts = total_add_quadrature_pnts + np
                    elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                    elem.Initialize()
                elif stat == 1:
                    elem.SetValue(ACTIVATION_LEVEL, -1)
                    elem.SetValue(IS_INACTIVE, True)
            print("construct quadrature using quad-tree successfully, " + str(total_add_quadrature_pnts) + " quadrature points are added")
            for elem in bulk_elements:
                self.mpu.SetMaterialProperties(model.model_part, elem, self.mat_type)
                physical_constitutive_laws = elem.GetValuesOnIntegrationPoints(CONSTITUTIVE_LAW, model.model_part.ProcessInfo)
                for i in range(0, len(physical_constitutive_laws)):
                    physical_constitutive_laws[i].SetValue(PARENT_ELEMENT_ID, elem.Id, model.model_part.ProcessInfo)
                    physical_constitutive_laws[i].SetValue(INTEGRATION_POINT_INDEX, i, model.model_part.ProcessInfo)
            print("reading " + self.mat_type + " material completed")

        elif self.quadrature_method == "moment-fit quadtree":
            cut_cell_quadrature_order = self.params["cut_cell_quadrature_order"]
            quadrature_data = self.params["quadrature_data"]
            ###########################################
            ##OBTAIN THE PRE-CALCULATED QUADRATURE FROM MOMENT FIT
            ###########################################
            cut_elems = quadrature_data.GetCutElements()
            exclude_elems = quadrature_data.GetExcludeElements()
            cut_cell_quadrature = quadrature_data.GetCutCellQuadrature()
            for elem_id in cut_elems:
                elem = model.model_part.Elements[elem_id]
                quad_util.SetQuadrature(elem, cut_cell_quadrature_order, cut_cell_quadrature[elem_id])
                elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                elem.Initialize()
            for elem_id in exclude_elems:
                elem = model.model_part.Elements[elem_id]
                elem.SetValue(ACTIVATION_LEVEL, -1)
                elem.SetValue(IS_INACTIVE, True)
            print("obtain moment-fit quadrature successfully")
            for elem in bulk_elements:
                self.mpu.SetMaterialProperties(model.model_part, elem, self.mat_type)
                physical_constitutive_laws = elem.GetValuesOnIntegrationPoints(CONSTITUTIVE_LAW, model.model_part.ProcessInfo)
                for i in range(0, len(physical_constitutive_laws)):
                    physical_constitutive_laws[i].SetValue(PARENT_ELEMENT_ID, elem.Id, model.model_part.ProcessInfo)
                    physical_constitutive_laws[i].SetValue(INTEGRATION_POINT_INDEX, i, model.model_part.ProcessInfo)
            print("reading " + self.mat_type + " material completed")

        elif self.quadrature_method == "moment-fit subcell":
            cut_cell_quadrature_order = self.params["cut_cell_quadrature_order"]
            quadrature_data = self.params["quadrature_data"]
            subcell_element_type = self.params["subcell_element_type"]
            quad_util = QuadratureUtility()
            aux_util = FiniteCellAuxilliaryUtility()
            ###########################################
            ##QUADTREE
            ###########################################
            self.CreateForestSubCell(bulk_elements)
            cut_elems = quadrature_data.GetCutElements()
            exclude_elems = quadrature_data.GetExcludeElements()
            cut_cell_quadrature = quadrature_data.GetCutCellQuadrature()
            cut_cell_full_quadrature = quadrature_data.GetCutCellFullQuadrature()
            subcell_weights = quadrature_data.GetCutCellSubCellWeights()
            lastElementId = aux_util.GetLastElementId(model.model_part)
            lastCondId = aux_util.GetLastConditionId(model.model_part)
            for elem in bulk_elements:
                if elem.Id in cut_elems:
#                    print("at cut element " + str(elem.Id))
                    if not cut_cell_quadrature[elem.Id]: # if cut cell quadrature is empty then skip
                        elem.SetValue(ACTIVATION_LEVEL, -1)
                        elem.SetValue(IS_INACTIVE, True)
                        elem.Set(ACTIVE, False)
                        continue

                    quad_util.SetQuadrature(elem, cut_cell_quadrature_order, cut_cell_quadrature[elem.Id])
                    elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                    elem.Initialize()

                    output = self.forest[elem.Id].CreateSubCellElements(model.model_part, subcell_element_type, cut_cell_quadrature_order, cut_cell_full_quadrature[elem.Id], subcell_weights[elem.Id], lastElementId, lastCondId)
                    lastElementId = output[0]
                    lastCondId = output[1]
                    new_sub_elems = output[2]

                    self.mpu.SetMaterialProperties(model.model_part, elem, self.mat_type)
                    elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                    elem.Initialize()

                    physical_constitutive_laws = elem.GetValuesOnIntegrationPoints(CONSTITUTIVE_LAW, model.model_part.ProcessInfo)
                    for i in range(0, len(physical_constitutive_laws)):
                        physical_constitutive_laws[i].SetValue(PARENT_ELEMENT_ID, elem.Id, model.model_part.ProcessInfo)
                        physical_constitutive_laws[i].SetValue(INTEGRATION_POINT_INDEX, i, model.model_part.ProcessInfo)
                    physical_integration_points = elem.GetValuesOnIntegrationPoints(INTEGRATION_POINT_LOCAL, model.model_part.ProcessInfo)
                    cnt = 0
                    for sub_elem in new_sub_elems:
                        sub_elem.SetValue(CONSTITUTIVE_LAW, physical_constitutive_laws[cnt])
                        sub_elem.SetValue(INTEGRATION_POINT_LOCAL, physical_integration_points[cnt])
                        cnt = cnt + 1
#                    print("cut element " + str(elem.Id) + " completed")
                elif elem.Id in exclude_elems:
                    elem.SetValue(ACTIVATION_LEVEL, -1)
                    elem.SetValue(IS_INACTIVE, True)
                    elem.Set(ACTIVE, False)
                else:
                    self.mpu.SetMaterialProperties(model.model_part, elem, self.mat_type)
            print("obtain moment-fit subcell quadrature successfully")
        else:
            print("Unknown quadrature_method", self.quadrature_method)
            sys.exit(0)

    ###CHECKING FUNCTIONS#############
    def InspectQuadTreeSubCell(self, bulk_elements, qt_depth):
        self.CreateForestSubCell(bulk_elements)

        cut_elems = []

        for qi, qs in self.forest.iteritems():
            elem = qs.GetElement()
            stat = self.brep.CutStatus(elem)
            elem.SetValue(CUT_STATUS, stat)
            if stat == -1:
                cut_elems.append(qs)

        aux_util = FiniteCellAuxilliaryUtility()
        for i in range(0, qt_depth):
            aux_util.MultithreadedQuadTreeSubCellRefineBy(cut_elems, self.brep)

        for qs in cut_elems:
            print("checking quadtree subcell of element " + str(qs.GetElement().Id) + ":")
            print("    number of subcells: " + str(qs.NumberOfSubCells()))
            for i in range(0, qs.NumberOfSubCells()):
                qt = qs.CreateQuadTree(i)
                stat = self.brep.CutStatus(qt.pGetGeometry())
                print("    subcell " + str(i) + " info:")
                if stat == -1:
                    print("        cut: yes")
                elif stat == 1:
                    print("        cut: out")
                else:
                    print("        cut: in")
                print("        domain size: " + str(qt.DomainSize(self.brep)))
                print("        center of gravity: " + str(qt.CenterOfGravity(self.brep)))
                print("        occupied geometry:")
                aux_util.Print(qt.pCreateGeometry())

    def ExportQuadTreeSubCell(self, model, bulk_elements, qt_depth, sample_element_name, sample_cond_name, selected_elems = "all"):
        self.CreateForestSubCell(bulk_elements)

        cut_elems = []

        if selected_elems == "all":
            for qi, qs in self.forest.iteritems():
                elem = qs.GetElement()
                stat = self.brep.CutStatus(elem)
                elem.SetValue(CUT_STATUS, stat)
                if stat == -1:
                    cut_elems.append(qs)
        else:
            for qi, qs in self.forest.iteritems():
                elem = qs.GetElement()
                if elem.Id in selected_elems:
                    stat = self.brep.CutStatus(elem)
                    elem.SetValue(CUT_STATUS, stat)
                    if stat == -1:
                        cut_elems.append(qs)

        aux_util = FiniteCellAuxilliaryUtility()
        for i in range(0, qt_depth):
            aux_util.MultithreadedQuadTreeSubCellRefineBy(cut_elems, self.brep)

        # add to the model_part
        cog_points = []
        for qs in cut_elems:
            lastNodeId = aux_util.GetLastNodeId(model.model_part)
            lastElementId = aux_util.GetLastElementId(model.model_part)
            qs.DeepAddToModelPart(model.model_part, sample_element_name, lastNodeId, lastElementId)
            for i in range(0, qs.NumberOfSubCells()):
                qt = qs.CreateQuadTree(i)
                cog_points.append(qt.CenterOfGravity(self.brep, 0x22))

        print("len(cog_points):", len(cog_points))
        quad_util = QuadratureUtility()
        quad_util.CreateConditionFromPoint(model.model_part, cog_points, sample_cond_name)

        print(model.model_part)
#        model.gid_io.InitializeMesh( 1234.0 )
#        mesh = model.model_part.GetMesh()
#        model.gid_io.WriteMesh( mesh )
#        model.gid_io.FinalizeMesh()
        model.WriteOutput( 1234.0 )


    ## Check if any quadtree subcell is integrable by refinement to some depth
    def CheckQuadTreeSubCell(self, bulk_elements, qt_depth):
        self.CreateForestSubCell(bulk_elements)

        cut_elems = []

        for qi, qs in self.forest.iteritems():
            elem = qs.GetElement()
            stat = self.brep.CutStatus(elem)
            elem.SetValue(CUT_STATUS, stat)
            if stat == -1:
                cut_elems.append(qs)

        aux_util = FiniteCellAuxilliaryUtility()
        for i in range(0, qt_depth):
            aux_util.MultithreadedQuadTreeSubCellRefineBy(cut_elems, self.brep)

        Hf = HeavisideFunctionR3R1(self.brep)

        print("checking subcell")
        for qs in cut_elems:
            for i in range(0, qs.NumberOfSubCells()):
                qt = qs.CreateQuadTree(i)
                stat = self.brep.CutStatus(qt.pCreateGeometry())
                if stat == -1:
#                    domain_size = qt.DomainSize(self.brep)
#                    domain_size = qt.Integrate(Hf, 0x12) # Gauss Legendre order 2
#                    domain_size = qt.Integrate(Hf, 0x13) # Gauss Legendre order 3
                    domain_size = qt.Integrate(Hf, 0x22) # Gauss Lobatto order 2
#                    domain_size = qt.Integrate(Hf, 0x23) # Gauss Lobatto order 3
#                    print("element " + str(qs.GetElement().Id) + ", subcell " + str(i) + " domain_size: " + str(domain_size))
                    if domain_size == 0.0:
                        print("at quadtree subcell of element " + str(qs.GetElement().Id) + ", subcell " + str(i) + ": the domain size is zero although cut")
                        # that means the quadtree can't integrate the domain size of the geometry occupied by the quadtree
        print("checking subcell done")

