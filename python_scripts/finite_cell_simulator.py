import math
import pprint
import time as time_module
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.FiniteCellApplication import *
#from KratosMultiphysics.FiniteCellStructuralApplication import *
CheckForPreviousImport()

###FITTING FUNCTIONS#############

## Create the list of function used for fitting
## if rank is 2 and degree 2 this subroutine will generate [1, x, x^2] x [1, y, y^2]
## if rank is 3 and degree 3 this subroutine will generate [1, x, x^2, x^3] x [1, y, y^2, y^3] x [1, z, z^2, z^3]
def CreateTensorProductFittingFunctions(rank, degree):
    funcs = []
    if rank == 2:
        for i in range(0, degree+1):
            if i > 0:
                expr1 = "pow(y, " + str(i) + ")"
            else:
                expr1 = '1'
            for j in range(0, degree+1):
                if j > 0:
                    expr2 = "pow(x, " + str(j) + ")"
                else:
                    expr2 = '1'
                funcs.append(MathPressoFunctionR3R1(expr1 + '*' + expr2))
    elif rank == 3:
        for i in range(0, degree+1):
            if i > 0:
                expr1 = "pow(z, " + str(i) + ")"
            else:
                expr1 = '1'
            for j in range(0, degree+1):
                if j > 0:
                    expr2 = "pow(y, " + str(j) + ")"
                else:
                    expr2 = '1'
                for k in range(0, degree+1):
                    if k > 0:
                        expr3 = "pow(x, " + str(k) + ")"
                    else:
                        expr3 = '1'
                    funcs.append(MathPressoFunctionR3R1(expr1 + '*' + expr2 + '*' + expr3))
    else:
        print("rank " + str(rank) + " is not supported")
        sys.exit(1)

    return funcs

def CreateFittingFunctions(fit_space_dim, fit_degree):
    fit_funcs = []

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

#            if fit_degree >= 3:
#                fit_funcs.append(f_x3)
#                fit_funcs.append(f_y3)
#                fit_funcs.append(f_x3y)
#                fit_funcs.append(f_xy3)
#                fit_funcs.append(f_x3y2)
#                fit_funcs.append(f_x2y3)
#                fit_funcs.append(f_x3y3)

        if fit_degree >= 3:
            fit_funcs = CreateTensorProductFittingFunctions(2, fit_degree)
#                print("Unsupported fit_degree " + str(fit_degree))

    elif fit_space_dim == 3:
        f_1 = ScalarFunctionR3R1(1.0)
        f_x = MonomialFunctionR3R1X()
        f_y = MonomialFunctionR3R1Y()
        f_z = MonomialFunctionR3R1Z()
        f_yz = MonomialFunctionR3R1YZ()
        f_xz = MonomialFunctionR3R1XZ()
        f_xy = MonomialFunctionR3R1XY()
        f_xyz = MonomialFunctionR3R1XYZ()

        f_x2 = MonomialFunctionR3R1X2()
        f_y2 = MonomialFunctionR3R1Y2()
        f_z2 = MonomialFunctionR3R1Z2()
        f_y2z = MonomialFunctionR3R1Y2Z()
        f_yz2 = MonomialFunctionR3R1YZ2()
        f_x2z = MonomialFunctionR3R1X2Z()
        f_xz2 = MonomialFunctionR3R1XZ2()
        f_x2y = MonomialFunctionR3R1X2Y()
        f_xy2 = MonomialFunctionR3R1XY2()
        f_y2z2 = MonomialFunctionR3R1Y2Z2()
        f_x2z2 = MonomialFunctionR3R1X2Z2()
        f_x2y2 = MonomialFunctionR3R1X2Y2()
        f_x2yz = MonomialFunctionR3R1X2YZ()
        f_xy2z = MonomialFunctionR3R1XY2Z()
        f_xyz2 = MonomialFunctionR3R1XYZ2()
        f_xy2z2 = MonomialFunctionR3R1XY2Z2()
        f_x2yz2 = MonomialFunctionR3R1X2YZ2()
        f_x2y2z = MonomialFunctionR3R1X2Y2Z()
        f_x2y2z2 = MonomialFunctionR3R1X2Y2Z2()

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
            fit_funcs.append(f_x2)
            fit_funcs.append(f_y2)
            fit_funcs.append(f_z2)
            fit_funcs.append(f_y2z)
            fit_funcs.append(f_yz2)
            fit_funcs.append(f_x2z)
            fit_funcs.append(f_xz2)
            fit_funcs.append(f_x2y)
            fit_funcs.append(f_xy2)
            fit_funcs.append(f_y2z2)
            fit_funcs.append(f_x2z2)
            fit_funcs.append(f_x2y2)
            fit_funcs.append(f_x2yz)
            fit_funcs.append(f_xy2z)
            fit_funcs.append(f_xyz2)
            fit_funcs.append(f_xy2z2)
            fit_funcs.append(f_x2yz2)
            fit_funcs.append(f_x2y2z)
            fit_funcs.append(f_x2y2z2)

        if fit_degree >= 3:
            fit_funcs = CreateTensorProductFittingFunctions(3, fit_degree)
#                print("Unsupported fit_degree " + str(fit_degree))

    return fit_funcs

def CreateQuadTree(element, nsampling):
    if nsampling == 1:
        return QuadTreeLocalReference(element)
    elif nsampling == 2:
        return QuadTreeLocalReference2(element)
    elif nsampling == 3:
        return QuadTreeLocalReference3(element)
    elif nsampling == 4:
        return QuadTreeLocalReference4(element)
    elif nsampling == 5:
        return QuadTreeLocalReference5(element)
    elif nsampling == 6:
        return QuadTreeLocalReference6(element)
    elif nsampling == 7:
        return QuadTreeLocalReference7(element)
    elif nsampling == 8:
        return QuadTreeLocalReference8(element)
    elif nsampling == 9:
        return QuadTreeLocalReference9(element)
    elif nsampling == 10:
        return QuadTreeLocalReference10(element)

def CreateQuadTreeSubCell(element, nsampling):
    if nsampling == 1:
        return MomentFittedQuadTreeSubCellReference(element)
    elif nsampling == 2:
        return MomentFittedQuadTreeSubCellReference2(element)
    elif nsampling == 3:
        return MomentFittedQuadTreeSubCellReference3(element)
    elif nsampling == 4:
        return MomentFittedQuadTreeSubCellReference4(element)
    elif nsampling == 5:
        return MomentFittedQuadTreeSubCellReference5(element)
    elif nsampling == 6:
        return MomentFittedQuadTreeSubCellReference6(element)
    elif nsampling == 7:
        return MomentFittedQuadTreeSubCellReference7(element)
    elif nsampling == 8:
        return MomentFittedQuadTreeSubCellReference8(element)
    elif nsampling == 9:
        return MomentFittedQuadTreeSubCellReference9(element)
    elif nsampling == 10:
        return MomentFittedQuadTreeSubCellReference10(element)

class FiniteCellSimulator:

    def __init__(self, params = None):
        self.params = params

        if 'enable_ghost_penalty' not in self.params: # turn off ghost penalty by default
            self.params['enable_ghost_penalty'] = False

        if 'enable_skeleton_penalty' not in self.params: # turn off skeleton penalty by default
            self.params['enable_skeleton_penalty'] = False

        if 'export_quadtree_cell' not in self.params:
            self.params['export_quadtree_cell'] = False

        if self.params['export_quadtree_cell'] == True:
            if 'sample_quad_element_name' not in self.params:
                self.params["sample_quad_element_name"] = "DummySurfaceElement2D4N"

    ###TREES & FOREST CREATION#############
    def CreateForest(self, elements, nsampling = 1):
        ###############################################################
        ######### CREATE TREES AND FOREST
        ###############################################################
        self.forest = {}
        for elem in elements:
            self.forest[elem.Id] = CreateQuadTree(elem, nsampling)

    def CreateForestSubCell(self, elements, nsampling = 1):
        self.forest = {}
        subcell_fit_mode = self.params["subcell_fit_mode"]
        cut_cell_quadrature_method = self.params["cut_cell_quadrature_method"]
        quad_util = QuadratureUtility()
        for elem in elements:
            self.forest[elem.Id] = CreateQuadTreeSubCell(elem, nsampling)
            if   subcell_fit_mode == "subcell fit gauss":
                self.forest[elem.Id].ConstructSubCellsBasedOnGaussQuadrature(cut_cell_quadrature_method)
            elif subcell_fit_mode == "subcell fit equal-dist":
                cut_cell_quadrature_order = quad_util.GetQuadratureOrder(cut_cell_quadrature_method)
                self.forest[elem.Id].ConstructSubCellsBasedOnEqualDistribution(cut_cell_quadrature_order)
            elif subcell_fit_mode == "subcell fit two-layer": # for backward compatibility
                self.forest[elem.Id].ConstructSubCellsBasedOnGaussQuadrature(cut_cell_quadrature_method)
            elif subcell_fit_mode == "subcell nonfit": # for backward compatibility
                pass
            else:
                print("unknown subcell_fit_mode", subcell_fit_mode)
                sys.exit(0)

    ###FITTING DRIVER#############
    def MomentFit(self, bulk_elements):
        print("fit parameters:")
        pprint.pprint(self.params)
        qt_depth = self.params["qt_depth"]
        nsampling = self.params["number_of_samplings"] if ("number_of_samplings" in self.params) else 1
        ###########################################
        ##QUADTREE QUADRATURE
        ###########################################
        self.CreateForest(bulk_elements, nsampling)
        self.elements = {}
        for elem in bulk_elements:
            self.elements[elem.Id] = elem
        cut_cell_quadrature_method = self.params["cut_cell_quadrature_method"]
        ###########################################
        ##MOMENT FIT
        ###########################################
        fit_util = MomentFittingUtility()
        fit_space_dim = self.params["fitting_space_dimension"]
        fit_degree = self.params["fitting_function_degree"]
        fit_funcs = CreateFittingFunctions(fit_space_dim, fit_degree)

        cut_elems = []
        exclude_elems = []

        fit_mode = self.params["fit_mode"]
        integrator_quadrature_method = self.params["integrator_quadrature_method"]
        solver_type = self.params["fit_solver_type"]
        echo_level = self.params["fit_echo_level"]
        small_weight = self.params["fit_small_weight"]

        if fit_mode == "serial":
            if echo_level == -2:
                num_elems = len(self.forest)
                cnt = 0
                cnt2 = 0
                print("moment-fitting begins")
            for qi, qt in self.forest.iteritems():
                elem = self.elements[qi]
                stat = self.brep.CutStatusBySampling(elem, nsampling)
                if stat == self.brep._CUT:
                    for i in range(0, qt_depth):
                        qt.RefineBy(self.brep)
                    fit_util.FitQuadrature(elem, fit_funcs, self.brep, qt, cut_cell_quadrature_method, integrator_quadrature_method, solver_type, echo_level, small_weight)
                    cut_elems.append(elem)
                elif stat == self.brep._OUT:
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
                stat = self.brep.CutStatusBySampling(elem, nsampling)
                if stat == self.brep._CUT:
                    cut_elems.append(elem)
                    integrators.append(qt)
                elif stat == self.brep._OUT:
                    exclude_elems.append(elem)
            aux_util = FiniteCellAuxilliaryUtility()
            for i in range(0, qt_depth):
                print("Refine level " + str(i+1) + " ...")
#                aux_util.MultithreadedQuadTreeRefineBy(integrators, self.brep)
                aux_util.MultithreadedRefineBy(integrators, self.brep)
                print("###################################################")
            fit_util.MultithreadedFitQuadrature(cut_elems, fit_funcs, self.brep, integrators, cut_cell_quadrature_method, integrator_quadrature_method, solver_type, echo_level, small_weight)
        else:
            print("Unknown fit_mode", fit_mode)
            sys.exit(0)

        print("construct quadrature using moment fitting successfully")
        quad_filename = self.params["quad_filename"]
        quad_filetype = self.params["quad_filetype"]
        accuracy = self.params["quad_accuracy"]
        fit_util.SaveQuadrature(quad_filename, quad_filetype, cut_elems, exclude_elems, accuracy)
        fid = open(quad_filename, "a")
        fid.write("\ndef GetParameter():\n")
        fid.write("    params = " + str(self.params) + "\n")
        fid.write("    return params\n")
        fid.close()

    def MomentFitSubCell(self, model_part, bulk_elements):
        print("fit parameters:")
        pprint.pprint(self.params)
        qt_depth = self.params["qt_depth"]
        max_qt_depth = self.params["max_qt_depth"]
        extra_qt_depth = self.params["extra_qt_depth"]
        small_domain_size = self.params["small_domain_size"]
        number_of_minimum_physical_points = self.params["number_of_minimum_physical_points"]
        nsampling = self.params["number_of_samplings"] if ("number_of_samplings" in self.params) else 1
        quad_util = QuadratureUtility()
        aux_util = FiniteCellAuxilliaryUtility()
        ###########################################
        ##QUADTREE QUADRATURE
        ###########################################
        self.CreateForestSubCell(bulk_elements, nsampling)
        ###########################################
        ##MOMENT FIT
        ###########################################
        fit_util = MomentFittingUtility()
        fit_space_dim = self.params["fitting_space_dimension"]
        fit_degree = self.params["fitting_function_degree"]
        fit_funcs = CreateFittingFunctions(fit_space_dim, fit_degree)
        print("list of fitting functions for fit_degree = " + str(fit_degree) + ":")
        for func in fit_funcs:
            print(func)

        cut_qs_elems = []
        exclude_qs_elems = []
        small_qs_elems = []

        integrator_quadrature_method = self.params["integrator_quadrature_method"]
        solver_type = self.params["fit_solver_type"]
        echo_level = self.params["fit_echo_level"]
        fit_small_weight = self.params["fit_small_weight"]
        action_on_small_cut_cell = self.params["action_on_small_cut_cell"]
        fict_weight = self.params["fictitious_weight"] if ("fictitious_weight" in self.params) else 0.0

        subcell_fit_func_callback = getattr(fit_util, "MultithreadedFitQuadratureSubCell")
        small_subcell_fit_func_callback = getattr(fit_util, "MultithreadedFitQuadratureSubCellUnique")
        save_quadrature_subcell_callback = getattr(fit_util, "SaveQuadratureSubCell")
        generate_physical_integration_points_callback = getattr(aux_util, "MultithreadedGeneratePhysicalIntegrationPoints")

        for qi, qs in self.forest.iteritems():
            elem = qs.GetElement()
            stat = self.brep.CutStatusBySampling(elem, nsampling)
            elem.SetValue(CUT_STATUS, stat)
            if stat == self.brep._CUT:
                cut_qs_elems.append(qs)
            elif stat == self.brep._OUT:
                exclude_qs_elems.append(qs)
        aux_util = FiniteCellAuxilliaryUtility()
        for i in range(0, qt_depth):
            print("Refine level " + str(i+1) + " ...")
#            aux_util.MultithreadedQuadTreeSubCellRefineBy(cut_qs_elems, self.brep)
            aux_util.MultithreadedRefineBy(cut_qs_elems, self.brep)
            print("###################################################")
        # recursive refine the subcell that are cut but not integrable by the quadtree # this is required if Gauss-Legendre quadrature is used
        if max_qt_depth > qt_depth:
            for qs in cut_qs_elems:
                for i in range(0, qs.NumberOfSubCells()):
                    qt = qs.CreateQuadTree(i)
                    stat = self.brep.CutStatusBySampling(qt.pCreateGeometry(), nsampling)
                    if stat == self.brep._CUT:
                        cnt = qt_depth
                        found = False
                        while cnt < max_qt_depth:
                            if qt.DomainSize(self.brep, integrator_quadrature_method) == 0.0:
                                qt.RefineBy(self.brep)
                                cnt = cnt + 1
                                print("subcell " + str(i) + " of element " + str(qs.GetElement().Id) + " is cut but domain size is zero, refine to level " + str(cnt))
                            else:
                                found = True
                                if cnt > qt_depth:
                                    for j in range(0, extra_qt_depth):
                                        qt.RefineBy(self.brep)
                                break
                        if (not found) and (cnt >= max_qt_depth):
                            print("can't find the physical point within subcell " + str(i) + " at " + str(max_qt_depth) + " quadtree level")
                            sys.exit(0)
        print("Refine completed")
        sys.stdout.flush()

        # generate the physical integration points in each cut-cell
        #### serial version
#        print("Generating physical integration points for sub-cells...")
#        cnt = 0
#        for qs in cut_qs_elems:
#            elem_id = qs.GetElement().Id

#            # for cut element, we generate and store the physical integration points
#            qs.GeneratePhysicalIntegrationPoints(self.brep, integrator_quadrature_method)

#            # if we want to add the fictitious element later on, the weight shall be set
#            if fict_weight != 0.0:
#                qs.SetFictitiousWeight(fict_weight)

#            cnt = cnt+1
#            print(" element " + str(elem_id) + " is done with physical integration points generation..." + str(cnt*100/len(cut_qs_elems)) + "% done")
#            sys.stdout.flush()
#        print("Generate physical integration points completed for " + str(len(cut_qs_elems)) + " cut-cells")
#        sys.stdout.flush()
        ###### multithreaded version
        print("Generating physical integration points for sub-cells...")
        generate_physical_integration_points_callback(cut_qs_elems, self.brep, integrator_quadrature_method)

        # if we want to add the fictitious element later on, the weight shall be set
        if fict_weight != 0.0:
            for qs in cut_qs_elems:
                qs.SetFictitiousWeight(fict_weight)

        print("Generate physical integration points completed for " + str(len(cut_qs_elems)) + " cut-cells")
        sys.stdout.flush()
        #########

        # here we check the small subcell; the small subcell is determined by either domain_size less than a threshold, or number of physical points less than specific number
        if small_domain_size > 0.0:
            criteria_check_small_domain_size = True
        else:
            criteria_check_small_domain_size = False

        if number_of_minimum_physical_points != 0:
            criteria_check_number_of_physical_points = True
        else:
            criteria_check_number_of_physical_points = False

        proper_cut_qs_elems = []
        cnt = 0
        cnt2 = 0
        cnt3 = 0
        num_cut_qs_elems = len(cut_qs_elems)

        if not(criteria_check_small_domain_size or criteria_check_number_of_physical_points):
            proper_cut_qs_elems = cut_qs_elems
            print("\nNo criteria is selected. No quadtree subcell is removed")
        else:
            for qs in cut_qs_elems:

                # here we check if the cut element shall be excluded or not. The criteria for exclusion is either small size or small number of physical integration points
                include_element = True

                if include_element:
                    if criteria_check_small_domain_size:
                        domain_size = qs.DomainSize(self.brep, integrator_quadrature_method)
                        if domain_size < small_domain_size:
                            include_element = False
                            print("\nthe quadtree subcell of element " + str(elem_id) + " is too small, domain size = " + str(domain_size) + ". It will be skipped.")

                if include_element:
                    if criteria_check_number_of_physical_points:
                        print("checking the number of physical point for element " + str(elem_id))
                        num_physical_points = len(qs.GetPhysicalIntegrationPoints())
                        print("the quadtree subcell of element " + str(elem_id) + " has number of physical point(s) = " + str(num_physical_points))
                        if num_physical_points < number_of_minimum_physical_points:
                            include_element = False
                            print("the quadtree subcell of element " + str(elem_id) + " has number of physical point(s) = " + str(num_physical_points) + " < " + str(number_of_minimum_physical_points) + ". It will be considered as small subcell. The action on it is: " + str(action_on_small_cut_cell) + "\n")

                if include_element:
                    proper_cut_qs_elems.append(qs)
                else:
                    if action_on_small_cut_cell == 'eliminate':
                        exclude_qs_elems.append(qs)
                        cnt = cnt + 1
                    elif action_on_small_cut_cell == 'replace by subcell quadtree':
                        small_qs_elems.append(qs)
                        cnt = cnt + 1
                    elif action_on_small_cut_cell == 'replace by quadtree':
                        small_qs_elems.append(qs)
                        cnt = cnt + 1
                    elif action_on_small_cut_cell == 'replace by moment-fit quadtree':
                        small_qs_elems.append(qs)
                        cnt = cnt + 1
                    else:
                        small_qs_elems.append(qs)
                        cnt = cnt + 1

                if cnt2 > 0.01*num_cut_qs_elems:
                    cnt2 = 0
                    sys.stdout.write(str(cnt3) + " ")
                    sys.stdout.flush()
                    cnt3 = cnt3 + 1
                else:
                    cnt2 = cnt2 + 1

            if action_on_small_cut_cell == 'eliminate':
                print("\nRemoval completed, " + str(cnt) + " quadtree subcell is removed")
            elif (action_on_small_cut_cell == 'replace by subcell quadtree') or (action_on_small_cut_cell == 'replace by quadtree'):
                print("\nReplace completed, " + str(cnt) + " quadtree subcell will be filled with quadtree quadrature")

        start_fit_t = time_module.time()
        subcell_fit_func_callback(proper_cut_qs_elems, fit_funcs, self.brep, integrator_quadrature_method, solver_type, echo_level, fit_small_weight)
        end_fit_t = time_module.time()
        print("construct quadrature using moment fitting for subcell successfully in " + str(end_fit_t - start_fit_t) + " s")
#        sys.exit(0)
        quad_filename = self.params["quad_filename"]
        quad_filetype = self.params["quad_filetype"]
        accuracy = self.params["quad_accuracy"]

        ## export the physical integration points for debugging if needed
        if self.params["export_physical_integration_point"]:
            cog_physical_points = []
            for qs in cut_qs_elems:
                elem = qs.GetElement()
                points = qs.GetPhysicalIntegrationPoints()
                for point in points:
                    cog_physical_points.append(quad_util.CreatePoint(elem, point[0], point[1], point[2]))
            print("len(cog_physical_points):", len(cog_physical_points))
            prop_id = self.params["physical_integration_point_prop_id"]
            quad_util.CreateConditionFromPoint(model_part, cog_physical_points, "DummyConditionPoint3D", model_part.Properties[prop_id])

        ## export the physical integration points for debugging if needed
        if self.params["export_fictitious_integration_point"]:
            cog_fictitious_points = []
            for qs in cut_qs_elems:
                elem = qs.GetElement()
                points = qs.GetFictitiousIntegrationPoints()
                for point in points:
                    cog_fictitious_points.append(quad_util.CreatePoint(elem, point[0], point[1], point[2]))
            print("len(cog_fictitious_points):", len(cog_fictitious_points))
            prop_id = self.params["fictitious_integration_point_prop_id"]
            quad_util.CreateConditionFromPoint(model_part, cog_fictitious_points, "DummyConditionPoint3D", model_part.Properties[prop_id])

        ## export the quadtree for debugging if needed
        if self.params["export_quadtree_cell"]:
            lastNodeId = aux_util.GetLastNodeId(model_part)
            lastElementId = aux_util.GetLastElementId(model_part)
            sample_quad_element_name = self.params["sample_quad_element_name"]
            sample_quad_element_prop_id = self.params["sample_quad_element_prop_id"]
            for qs in cut_qs_elems:
                last_ids = qs.DeepAddToModelPart(model_part, sample_quad_element_name, lastNodeId, lastElementId, sample_quad_element_prop_id)
                lastNodeId = last_ids[0]
                lastElementId = last_ids[1]

        ## handle the small cells
        if action_on_small_cut_cell == 'eliminate':
            save_quadrature_subcell_callback(quad_filename, quad_filetype, proper_cut_qs_elems, exclude_qs_elems, small_qs_elems, accuracy)
        elif action_on_small_cut_cell == 'replace by quadtree':
            # perform quadtree refinement for "eliminated" cell
            total_add_quadrature_pnts = 0
            small_cut_cell_quadrature_method = self.params["small_cut_cell_quadrature_method"]
            quadtree_small_weight = self.params["quadtree_small_weight"]
            small_cut_cell_qt_depth = self.params["small_cut_cell_qt_depth"]
            valid_quadtree_qs_elems = []
            for qs in small_qs_elems:
                qs.Clear()
                qs.ConstructSubCellsBasedOnGaussQuadrature(1)
                for i in range(0, small_cut_cell_qt_depth):
                    qs.RefineBy(self.brep)
                np = qs.ConstructQuadrature(self.brep, small_cut_cell_quadrature_method, quadtree_small_weight)
                total_add_quadrature_pnts = total_add_quadrature_pnts + np
                if np == 0:
                    exclude_qs_elems.append(qs)
                else:
                    valid_quadtree_qs_elems.append(qs)
            save_quadrature_subcell_callback(quad_filename, quad_filetype, proper_cut_qs_elems, exclude_qs_elems, valid_quadtree_qs_elems, accuracy)
            print("generate quadtree quadrature for small cell completed, " + str(total_add_quadrature_pnts) + " quadrature points are added")
            ## export the quadtree for debugging if needed
            if self.params["export_small_quadtree_cell"]:
                lastNodeId = aux_util.GetLastNodeId(model_part)
                lastElementId = aux_util.GetLastElementId(model_part)
                sample_quad_element_name = self.params["sample_quad_element_name"]
                for qs in small_qs_elems:
                    last_ids = qs.DeepAddToModelPart(model_part, sample_quad_element_name, lastNodeId, lastElementId, 51)
                    lastNodeId = last_ids[0]
                    lastElementId = last_ids[1]
        elif action_on_small_cut_cell == 'replace by subcell quadtree':
            # perform quadtree refinement for "eliminated" cell
            total_add_quadrature_pnts = 0
            small_cut_cell_quadrature_method = self.params["small_cut_cell_quadrature_method"]
            quadtree_small_weight = self.params["quadtree_small_weight"]
            small_cut_cell_qt_depth = self.params["small_cut_cell_qt_depth"]
            valid_quadtree_qs_elems = []
            for qs in small_qs_elems:
                qs.Clear()
                for i in range(0, small_cut_cell_qt_depth):
                    qs.RefineBy(self.brep)
                np = qs.ConstructQuadrature(self.brep, small_cut_cell_quadrature_method, quadtree_small_weight)
                total_add_quadrature_pnts = total_add_quadrature_pnts + np
                if np == 0:
                    exclude_qs_elems.append(qs)
                else:
                    valid_quadtree_qs_elems.append(qs)
            save_quadrature_subcell_callback(quad_filename, quad_filetype, proper_cut_qs_elems, exclude_qs_elems, valid_quadtree_qs_elems, accuracy)
            print("generate subcell quadtree quadrature for small subcell completed, " + str(total_add_quadrature_pnts) + " quadrature points are added")
        elif action_on_small_cut_cell == 'replace by moment-fit quadtree':
            small_cut_cell_quadrature_method = self.params["small_cut_cell_quadrature_method"]
            small_cut_cell_qt_depth = self.params["small_cut_cell_qt_depth"]
            small_elems = []
            small_qt_elems = []
            print("len(small_qs_elems):", len(small_qs_elems))
            for qs in small_qs_elems:
                qs.Clear()
                qt = CreateQuadTree(qs.GetElement(), nsampling)
                for i in range(0, small_cut_cell_qt_depth):
                    qt.RefineBy(self.brep)
                small_qt_elems.append(qt)
                small_elems.append(qs.GetElement())
            small_fit_solver_type = self.params["small_fit_solver_type"]
            small_fit_echo_level = self.params["small_fit_echo_level"]
            small_fit_small_weight = self.params["small_fit_small_weight"]
            small_subcell_fit_func_callback(small_elems, fit_funcs, self.brep, small_qs_elems, small_cut_cell_quadrature_method, integrator_quadrature_method, small_fit_solver_type, small_fit_echo_level, small_fit_small_weight)
            save_quadrature_subcell_callback(quad_filename, quad_filetype, proper_cut_qs_elems, exclude_qs_elems, small_qs_elems, accuracy)
            print("generate moment-fit quadrature for small cell completed")

        fid = open(quad_filename, "a")
        fid.write("\ndef GetParameter():\n")
        fid.write("    params = " + str(self.params) + "\n")
        fid.write("    return params\n")
        fid.close()

    ###SIMULATION DRIVER#############
    def Initialize(self, model, bulk_elements):
        self.quadrature_method = self.params["quadrature_method"]
        self.mpu = self.params["material_properties_utility"]
        self.mat_type = self.params["mat_type"]
        nsampling = self.params["number_of_samplings"] if ("number_of_samplings" in self.params) else 1
        print("simulator parameters:")
        pprint.pprint(self.params)
        quad_util = QuadratureUtility()
        aux_util = FiniteCellAuxilliaryUtility()
        mesh_util = FiniteCellMeshUtility()

        if self.quadrature_method == "quadtree":
            qt_depth = self.params["qt_depth"]
            small_weight = self.params["small_weight"]
            ###########################################
            ##QUADTREE
            ###########################################
            self.CreateForest(bulk_elements, nsampling)
            cut_cell_quadrature_method = self.params["cut_cell_quadrature_method"]
            cut_cell_quadrature_order = quad_util.GetQuadratureOrder(cut_cell_quadrature_method)
            total_add_quadrature_pnts = 0
            self.proper_cut_elems = ElementsArray()
            cut_elems = []
            exclude_elems = []
            for qi, qt in self.forest.iteritems():
                elem = model.model_part.Elements[qi]
                stat = self.brep.CutStatusBySampling(elem, nsampling)
                if stat == self.brep._CUT:
                    for i in range(0, qt_depth):
                        qt.RefineBy(self.brep)
                    np = qt.ConstructQuadrature(self.brep, cut_cell_quadrature_method, small_weight)
                    total_add_quadrature_pnts = total_add_quadrature_pnts + np
                    elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                    elem.Initialize()
                    aux_util.AddElement(self.proper_cut_elems, elem)
                    cut_elems.append(elem)
                elif stat == self.brep._OUT:
                    elem.SetValue(ACTIVATION_LEVEL, -1)
                    elem.SetValue(IS_INACTIVE, True)
                    elem.Set(ACTIVE, False)
                    exclude_elems.append(elem)
                    print("element " + str(elem.Id) + " is outside the physical domain. It will be deactivated.")
            print("construct quadrature using quad-tree successfully, " + str(total_add_quadrature_pnts) + " quadrature points are added")
            for elem in bulk_elements:
                self.mpu.SetMaterialProperties(model.model_part, elem, self.mat_type)
                physical_constitutive_laws = elem.GetValuesOnIntegrationPoints(CONSTITUTIVE_LAW, model.model_part.ProcessInfo)
                for i in range(0, len(physical_constitutive_laws)):
                    physical_constitutive_laws[i].SetValue(PARENT_ELEMENT_ID, elem.Id, model.model_part.ProcessInfo)
                    physical_constitutive_laws[i].SetValue(INTEGRATION_POINT_INDEX, i, model.model_part.ProcessInfo)
            print("reading " + self.mat_type + " material completed")

            ## export the physical integration points for debugging if needed
            if self.params["export_physical_integration_point"]:
                cog_points = []
                for elem in self.proper_cut_elems:
                    points = elem.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
                    for point in points:
                        cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))
                print("len(cog_points):", len(cog_points))
                prop_id = self.params["physical_integration_point_prop_id"]
                quad_util.CreateConditionFromPoint(model.model_part, cog_points, "DummyConditionPoint3D", model.model_part.Properties[prop_id])

            ## export the quadtree for debugging if needed
            if self.params["export_quadtree_cell"]:
                lastNodeId = aux_util.GetLastNodeId(model.model_part)
                lastElementId = aux_util.GetLastElementId(model.model_part)
                sample_quad_element_name = self.params["sample_quad_element_name"]
                for qi, qt in self.forest.iteritems():
                    last_ids = qt.AddToModelPart(model.model_part, sample_quad_element_name, lastNodeId, lastElementId)
                    lastNodeId = last_ids[0]
                    lastElementId = last_ids[1]

            ## export the quadrature if user needs
            if self.params["write_quadrature_to_file"] == True:
                quad_filename = self.params["quad_filename"]
                quad_filetype = self.params["quad_filetype"]
                accuracy = self.params["quad_accuracy"]
                self.params["material_properties_utility"] = ""
                quad_util.SaveQuadrature(quad_filename, quad_filetype, cut_elems, exclude_elems, accuracy)
                fid = open(quad_filename, "a")
                fid.write("\ndef GetParameter():\n")
                fid.write("    params = " + str(self.params) + "\n")
                fid.write("    return params\n")
                fid.close()
            # end if self.quadrature_method == "quadtree":
        elif self.quadrature_method == "quadtree preload":
            cut_cell_quadrature_method = self.params["cut_cell_quadrature_method"]
            cut_cell_quadrature_order = quad_util.GetQuadratureOrder(cut_cell_quadrature_method)
            quadrature_data = self.params["quadrature_data"]
            cut_elems = quadrature_data.GetCutElements()
            exclude_elems = quadrature_data.GetExcludeElements()
            cut_cell_quadrature = quadrature_data.GetCutCellQuadrature()
            ###########################################
            ##QUADTREE PRELOAD FROM FILE
            ###########################################
            total_add_quadrature_pnts = 0
            self.proper_cut_elems = ElementsArray()
            for elem_id in cut_elems:
                elem = model.model_part.Elements[elem_id]
                quad_util.SetQuadrature(elem, cut_cell_quadrature_order, cut_cell_quadrature[elem_id])
                total_add_quadrature_pnts = total_add_quadrature_pnts + len(cut_cell_quadrature[elem_id])
                elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                elem.Initialize()
                aux_util.AddElement(self.proper_cut_elems, elem)
            for elem_id in exclude_elems:
                elem = model.model_part.Elements[elem_id]
                elem.SetValue(ACTIVATION_LEVEL, -1)
                elem.SetValue(IS_INACTIVE, True)
                elem.Set(ACTIVE, False)
            print("obtain quadtree quadrature successfully, " + str(total_add_quadrature_pnts) + " quadrature points are added")
            for elem_id in cut_elems:
                elem = model.model_part.Elements[elem_id]
                self.mpu.SetMaterialProperties(model.model_part, elem, self.mat_type)
                physical_constitutive_laws = elem.GetValuesOnIntegrationPoints(CONSTITUTIVE_LAW, model.model_part.ProcessInfo)
                for i in range(0, len(physical_constitutive_laws)):
                    physical_constitutive_laws[i].SetValue(PARENT_ELEMENT_ID, elem.Id, model.model_part.ProcessInfo)
                    physical_constitutive_laws[i].SetValue(INTEGRATION_POINT_INDEX, i, model.model_part.ProcessInfo)
            print("reading " + self.mat_type + " material completed")
            if self.params["export_physical_integration_point"]:
                cog_points = []
                for elem in self.proper_cut_elems:
                    points = elem.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
                    for point in points:
                        cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))
                print("len(cog_points):", len(cog_points))
                prop_id = self.params["physical_integration_point_prop_id"]
                quad_util.CreateConditionFromPoint(model.model_part, cog_points, "DummyConditionPoint3D", model.model_part.Properties[prop_id])
            #end elif self.quadrature_method == "quadtree":
        elif self.quadrature_method == "moment-fit quadtree":
            cut_cell_quadrature_order = self.params["cut_cell_quadrature_order"]
            quadrature_data = self.params["quadrature_data"]
            ###########################################
            ##OBTAIN THE PRE-CALCULATED QUADRATURE FROM MOMENT FIT
            ###########################################
            cut_elems = quadrature_data.GetCutElements()
            exclude_elems = quadrature_data.GetExcludeElements()
            cut_cell_quadrature = quadrature_data.GetCutCellQuadrature()
            self.proper_cut_elems = ElementsArray()
            for elem_id in cut_elems:
                elem = model.model_part.Elements[elem_id]
                quad_util.SetQuadrature(elem, cut_cell_quadrature_order, cut_cell_quadrature[elem_id])
                elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                elem.Initialize()
                aux_util.AddElement(self.proper_cut_elems, elem)
            for elem_id in exclude_elems:
                elem = model.model_part.Elements[elem_id]
                elem.SetValue(ACTIVATION_LEVEL, -1)
                elem.SetValue(IS_INACTIVE, True)
                elem.Set(ACTIVE, False)
            print("obtain moment-fit quadrature successfully")
            for elem in bulk_elements:
                self.mpu.SetMaterialProperties(model.model_part, elem, self.mat_type)
                physical_constitutive_laws = elem.GetValuesOnIntegrationPoints(CONSTITUTIVE_LAW, model.model_part.ProcessInfo)
                for i in range(0, len(physical_constitutive_laws)):
                    physical_constitutive_laws[i].SetValue(PARENT_ELEMENT_ID, elem.Id, model.model_part.ProcessInfo)
                    physical_constitutive_laws[i].SetValue(INTEGRATION_POINT_INDEX, i, model.model_part.ProcessInfo)
            print("reading " + self.mat_type + " material completed")
            if self.params["export_physical_integration_point"]:
                cog_points = []
                for elem in self.proper_cut_elems:
                    points = elem.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
                    for point in points:
                        cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))
                print("len(cog_points):", len(cog_points))
                prop_id = self.params["physical_integration_point_prop_id"]
                quad_util.CreateConditionFromPoint(model.model_part, cog_points, "DummyConditionPoint3D", model.model_part.Properties[prop_id])
            # end elif self.quadrature_method == "moment-fit quadtree":
        elif self.quadrature_method == "moment-fit subcell":
            cut_cell_quadrature_method = self.params["cut_cell_quadrature_method"]
            cut_cell_quadrature_order = quad_util.GetQuadratureOrder(cut_cell_quadrature_method)
            quadrature_data = self.params["quadrature_data"]
            extrapolation_mode = self.params["extrapolation_mode"]
            if extrapolation_mode == "extrapolated constant stress element":
                subcell_element_type = "ExtrapolatedConstantStress" + self.params["subcell_element_basetype"]
            elif extrapolation_mode == "extrapolated element":
                subcell_element_type = "Extrapolated" + self.params["subcell_element_basetype"]
            else:
                subcell_element_type = self.params["subcell_element_basetype"]
            if "extrapolated material" in extrapolation_mode:
                extrapolated_material_mode = extrapolation_mode.split()[2] # must be either "implicit", "explicit" or "implex"
            ###########################################
            ##QUADTREE
            ###########################################
            self.CreateForestSubCell(bulk_elements, nsampling)
            cut_elems = quadrature_data.GetCutElements()
            exclude_elems = quadrature_data.GetExcludeElements()
            quadtree_elems = quadrature_data.GetQuadTreeElements()
            self.fict_elems = ElementsArray()
            cut_cell_quadrature = quadrature_data.GetCutCellQuadrature()
            cut_cell_fictitious_quadrature = quadrature_data.GetCutCellFictitiousQuadrature()
            cut_cell_full_quadrature = quadrature_data.GetCutCellFullQuadrature()
            subcell_weights = quadrature_data.GetCutCellSubCellWeights()
            subcell_domain_sizes = quadrature_data.GetCutCellSubCellDomainSizes()
            quadtree_quadrature = quadrature_data.GetQuadTreeQuadrature()
            self.proper_cut_elems = ElementsArray()
            lastElementId = aux_util.GetLastElementId(model.model_part)
            lastCondId = aux_util.GetLastConditionId(model.model_part)
            self.all_subcell_elems = []
            extrapolated_prop = model.model_part.Properties[self.params["extrapolated_prop_id"]]
            if self.params["export_physical_integration_point"]:
                cog_points = []
            if self.params["export_fictitious_integration_point"]:
                fict_cog_points = []
            if self.params["add_fictitious_element"]:
                sample_fictitious_element_name = self.params["fictitious_element_name"]
                if "fictitious_prop" not in self.params:
                    fict_prop = "same as parent element"
                else:
                    fict_prop = self.params["fictitious_prop"]
            for elem in bulk_elements:
                if elem.Id in cut_elems:
#                    print("at cut element " + str(elem.Id))
                    if len(cut_cell_quadrature[elem.Id]) == 0: # if cut cell quadrature is empty then skip
                        elem.SetValue(ACTIVATION_LEVEL, -1)
                        elem.SetValue(IS_INACTIVE, True)
                        elem.Set(ACTIVE, False)
                        print("element " + str(elem.Id) + " has no physical integration points. It will be excluded.")
                    else:
                        aux_util.AddElement(self.proper_cut_elems, elem)

                        quad_util.SetQuadrature(elem, cut_cell_quadrature_order, cut_cell_quadrature[elem.Id])
                        elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                        elem.Initialize()

                        output = self.forest[elem.Id].CreateSubCellElements(model.model_part, subcell_element_type, cut_cell_quadrature_order, cut_cell_full_quadrature[elem.Id], subcell_weights[elem.Id], lastElementId, lastCondId, extrapolated_prop)
                        lastElementId = output[0]
                        lastCondId = output[1]
                        new_sub_elems = output[2] # !!!IMPORTANT!!! here it's assumed that the sub-elements are Extrapolated ones

                        self.mpu.SetMaterialProperties(model.model_part, elem, self.mat_type)
                        elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                        elem.Initialize()

                        physical_constitutive_laws = elem.GetValuesOnIntegrationPoints(CONSTITUTIVE_LAW, model.model_part.ProcessInfo)
                        for i in range(0, len(physical_constitutive_laws)):
                            physical_constitutive_laws[i].SetValue(PARENT_ELEMENT_ID, elem.Id, model.model_part.ProcessInfo)
                            physical_constitutive_laws[i].SetValue(INTEGRATION_POINT_INDEX, i, model.model_part.ProcessInfo)
    #                        physical_constitutive_laws[i].SetValue(REPRESENTATIVE_WEIGHT, subcell_domain_sizes[elem.Id][i] / len(cut_cell_full_quadrature[elem.Id]), model.model_part.ProcessInfo) # here we divide to number of integration point of the full quadrature. Because all the subcell element will be added up when calculating the error

                        if extrapolation_mode == "extrapolated element" or extrapolation_mode == "extrapolated constant stress element":
                            physical_integration_points = elem.GetValuesOnIntegrationPoints(INTEGRATION_POINT_LOCAL, model.model_part.ProcessInfo)
                            cnt = 0
                            for sub_elem in new_sub_elems:
                                self.all_subcell_elems.append(sub_elem)
                                sub_elem.SetValue(CONSTITUTIVE_LAW, physical_constitutive_laws[cnt])
                                sub_elem.SetValue(PHYSICAL_INTEGRATION_POINT_LOCAL, physical_integration_points[cnt])
                                sub_elem.SetValue(SUBCELL_DOMAIN_SIZE, subcell_domain_sizes[elem.Id][cnt])
                                cnt = cnt + 1
                        elif "extrapolated material" in extrapolation_mode:
                            representative_weight = subcell_domain_sizes[elem.Id][i] / len(cut_cell_full_quadrature[elem.Id]) # here we divide to number of integration point of the full quadrature. Because all the subcell element will be added up when calculating the error
                            cnt = 0
                            for sub_elem in new_sub_elems:
                                self.all_subcell_elems.append(sub_elem)
                                self.mpu.SetExtrapolatedMaterial(model.model_part, sub_elem, physical_constitutive_laws[cnt], representative_weight, extrapolated_material_mode)
                                cnt = cnt + 1

                        if self.params["export_physical_integration_point"]:
                            points = elem.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
                            for point in points:
                                cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))

                    if self.params["add_fictitious_element"]:
                        if len(cut_cell_fictitious_quadrature[elem.Id]) > 0:
                            if fict_prop == "same as parent element": # this uses the material proper from the element, even it's nonlinear
                                fict_elem = mesh_util.CreateParasiteElement(sample_fictitious_element_name, lastElementId, elem, cut_cell_quadrature_order, cut_cell_fictitious_quadrature[elem.Id], elem.Properties)
                                self.mpu.SetMaterialProperties(model.model_part, fict_elem, self.mat_type)
                                fict_elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                                fict_elem.Initialize()
                            else: # this allows to set the material property of fictitious element DIFFERENT to the parent element
                                fict_elem = mesh_util.CreateParasiteElement(sample_fictitious_element_name, lastElementId, elem, cut_cell_quadrature_order, cut_cell_fictitious_quadrature[elem.Id], fict_prop)
                            print("fictitious element " + str(fict_elem.Id) + " is created out from element " + str(elem.Id))
                            aux_util.AddElement(self.fict_elems, fict_elem)
                            lastElementId = fict_elem.Id

                            if self.params["export_fictitious_integration_point"]:
                                points = fict_elem.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
                                for point in points:
                                    fict_cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))

                elif elem.Id in exclude_elems:
                    elem.SetValue(ACTIVATION_LEVEL, -1)
                    elem.SetValue(IS_INACTIVE, True)
                    elem.Set(ACTIVE, False)
                    print("element " + str(elem.Id) + " is outside. It will be excluded")

                elif elem.Id in quadtree_elems:
                    quad_util.SetQuadrature(elem, cut_cell_quadrature_order, quadtree_quadrature[elem.Id])
                    elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                    elem.Initialize()
                    self.mpu.SetMaterialProperties(model.model_part, elem, self.mat_type)
#                    aux_util.AddElement(self.proper_cut_elems, elem) # a quadtree element is NOT considered as proper_cut_elems
                    print("element " + str(elem.Id) + " is assigned " + str(len(quadtree_quadrature[elem.Id])) + " quadrature points")
                    if self.params["export_physical_integration_point"]:
                        points = elem.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
                        for point in points:
                            cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))
                else:
                    self.mpu.SetMaterialProperties(model.model_part, elem, self.mat_type)
            print("obtain moment-fit subcell quadrature successfully")
            if self.params["export_physical_integration_point"]:
                prop_id = self.params["physical_integration_point_prop_id"]
                quad_util.CreateConditionFromPoint(model.model_part, cog_points, "DummyConditionPoint3D", model.model_part.Properties[prop_id])
            if self.params["export_fictitious_integration_point"]:
                prop_id = self.params["fictitious_integration_point_prop_id"]
                quad_util.CreateConditionFromPoint(model.model_part, fict_cog_points, "DummyConditionPoint3D", model.model_part.Properties[prop_id])

            if self.params["add_fictitious_element"] == True:
                for fict_elem in self.fict_elems:
                    model.model_part.AddElement(fict_elem)
                    print("fictitious element " + str(fict_elem.Id) + " is added to the model_part")

            # end elif self.quadrature_method == "moment-fit subcell":
        elif self.quadrature_method == "moment-fit subcell no-preload":
            ###########################################
            ##QUADTREE
            ###########################################
            self.CreateForestSubCell(bulk_elements, nsampling)
            self.proper_cut_elems = ElementsArray()
        else:
            print("Unknown quadrature_method", self.quadrature_method)
            sys.exit(0)

        if self.params["enable_ghost_penalty"] == True:
            ghost_penalty_util = GhostPenaltyUtility()

        if self.params["enable_skeleton_penalty"] == True:
            ghost_penalty_util = SkeletonPenaltyUtility()

        if self.params["enable_ghost_penalty"] == True and self.params["enable_skeleton_penalty"] == True:
            print("Both ghost_penalty and skeleton_penalty are activated. You should choose only one.")
            sys.exit(0)

        if self.params["enable_ghost_penalty"] == True or self.params["enable_skeleton_penalty"] == True:
            space_dim = self.params["space_dim"]
            estimated_neighbours = self.params["estimated_number_of_neighbours"]

            FindElementalNeighboursProcess(model.model_part, space_dim, estimated_neighbours).Execute()

            sample_cond = self.params["sample_ghost_penalty_condition"]
            ghost_prop = self.params["ghost_penalty_properties"]
            lastCondId = aux_util.GetLastConditionId(model.model_part)
            ghost_penalty_conds = ghost_penalty_util.SetUpSurfacePenaltyConditions(model.model_part, bulk_elements, sample_cond, self.brep, lastCondId, ghost_prop)

    # end Initialize

    ###COMPUTE GLOBAL DISPLACEMENT (L2) ERROR###
    def compute_L2_error(self, elements, process_info, solution, P):
        print("!!!compute_L2_error:Please turn off the MoveMeshFlag in order to have correct results")
        nom = 0.0
        denom = 0.0
        for element in elements:
            if element.Is(ACTIVE):
                u = element.GetValuesOnIntegrationPoints(DISPLACEMENT, process_info)
                J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
                W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
#                print("number of integration point on element " + str(element.Id) + ":" + str(len(Q)))
#                if element.Id == 6:
#                    print("at element ", element.Id)
#                    print("u:", u)
#                    print("Q:", Q)
#                    print("W:", W)
#                    print("J0:", J0)
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

    def compute_L2_error_mfsc(self, elements, all_subcell_elems, process_info, solution, P):
        print("!!!compute_L2_error_mfsc:Please turn off the MoveMeshFlag in order to have correct results")
        nom = 0.0
        denom = 0.0
        for element in elements:
            if element.Is(ACTIVE):
                u = element.GetValuesOnIntegrationPoints(DISPLACEMENT, process_info)
                J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
                W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                for i in range(0, len(u)):
                    ana_u = solution.get_displacement(P, Q[i][0], Q[i][1], Q[i][2])
                    nom = nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2) + pow(u[i][2] - ana_u[2], 2)) * W[i][0] * J0[i][0]
                    denom = denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2) + pow(ana_u[2], 2)) * W[i][0] * J0[i][0]
        for element in all_subcell_elems:
            if element.Is(ACTIVE):
                u = element.GetValuesOnIntegrationPoints(PHYSICAL_INTEGRATION_POINT_DISPLACEMENT, process_info)
#                J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.GetValuesOnIntegrationPoints(PHYSICAL_INTEGRATION_POINT_GLOBAL, process_info)
                W = element.GetValuesOnIntegrationPoints(SUBCELL_DOMAIN_SIZE, process_info)
                for i in range(0, len(u)):
                    ana_u = solution.get_displacement(P, Q[i][0], Q[i][1], Q[i][2])
                    nom = nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2) + pow(u[i][2] - ana_u[2], 2)) * W[i][0]# * J0[i][0]
                    denom = denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2) + pow(ana_u[2], 2)) * W[i][0]# * J0[i][0]
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
        print("!!!compute_H1_error:Please turn off the MoveMeshFlag in order to have correct results")
        nom = 0.0
        denom = 0.0
        for element in elements:
            if element.Is(ACTIVE):
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

    def compute_H1_error_mfsc(self, elements, all_subcell_elems, process_info, solution, P):
        print("!!!compute_H1_error_mfsc:Please turn off the MoveMeshFlag in order to have correct results")
        nom = 0.0
        denom = 0.0
        for element in elements:
            if element.Is(ACTIVE):
                o = element.GetValuesOnIntegrationPoints(THREED_STRESSES, process_info)
                J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
                W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                for i in range(0, len(o)):
                    ana_o = solution.get_stress_3d(P, Q[i][0], Q[i][1], Q[i][2])
                    nom = nom + (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2))) * W[i][0] * J0[i][0]
                    denom = denom + (pow(ana_o[0], 2) + pow(ana_o[1], 2) + pow(ana_o[2], 2) + 2.0*(pow(ana_o[3], 2) + pow(ana_o[4], 2) + pow(ana_o[5], 2))) * W[i][0] * J0[i][0]
        for element in all_subcell_elems:
            if element.Is(ACTIVE):
                o = element.GetValuesOnIntegrationPoints(PHYSICAL_INTEGRATION_POINT_THREED_STRESSES, process_info)
#                J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.GetValuesOnIntegrationPoints(PHYSICAL_INTEGRATION_POINT_GLOBAL, process_info)
                W = element.GetValuesOnIntegrationPoints(SUBCELL_DOMAIN_SIZE, process_info)
                for i in range(0, len(o)):
                    ana_o = solution.get_stress_3d(P, Q[i][0], Q[i][1], Q[i][2])
                    nom = nom + (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2))) * W[i][0]# * J0[i][0]
                    denom = denom + (pow(ana_o[0], 2) + pow(ana_o[1], 2) + pow(ana_o[2], 2) + 2.0*(pow(ana_o[3], 2) + pow(ana_o[4], 2) + pow(ana_o[5], 2))) * W[i][0]# * J0[i][0]
        print("nom:", nom)
        print("denom:", denom)
        if denom == 0.0:
            if nom == 0.0:
                return 0.0
            else:
                return float('nan');
        else:
            return math.sqrt(abs(nom / denom))

    ###COMPUTE DOMAIN SIZE##############
    def compute_domain_size(self, elements, process_info):
        domain_size = 0.0
        for element in elements:
            if element.Is(ACTIVE):
                J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
                W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                for i in range(0, len(W)):
                    domain_size = domain_size + W[i][0] * J0[i][0]
        return domain_size

    def compute_domain_size_subcell(self, all_subcell_elems, process_info):
        domain_size = 0.0
        for element in all_subcell_elems:
            if element.Is(ACTIVE):
                W = element.GetValuesOnIntegrationPoints(SUBCELL_DOMAIN_SIZE, process_info)
                for i in range(0, len(W)):
                    domain_size = domain_size + W[i][0]
        return domain_size

    ###CHECKING FUNCTIONS#############
    def ExportQuadTree(self, model, sample_element, group = None):
        nsampling = self.params["number_of_samplings"] if ("number_of_samplings" in self.params) else 1
        self.CreateForest(model.model_part.Elements, nsampling)
        #######QUAD TREE POST PROCESSING###########
        qt_depth = self.params["qt_depth"]
        for qi, qt in self.forest.iteritems():
            if qi in group:
                for i in range(0, qt_depth):
                    qt.RefineBy(self.brep)
        tmodel_part = ModelPart("QuadTree")
        last_node_id = 0x0000000000000000 + len(model.model_part.Nodes)
        last_element_id = 0x0000000000000000
        for qi, qt in self.forest.iteritems():
            last_ids = qt.AddToModelPart(tmodel_part, sample_element, last_node_id, last_element_id)
            last_node_id = last_ids[0]
            last_element_id = last_ids[1]
        print(tmodel_part)
        model.gid_io.InitializeMesh( 1234.0 )
        mesh = tmodel_part.GetMesh()
        model.gid_io.WriteMesh( mesh )
        model.gid_io.FinalizeMesh()
        ###########################################

    ## check if the moment fitting can reproduce the constant function
    def IntegrateConstantFunction(self, elements, process_info):
        domain_size = 0.0
        for element in elements:
            if element.Is(ACTIVE):
                J0 = element.GetValuesOnIntegrationPoints(JACOBIAN_0, process_info)
                W = element.GetValuesOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                for i in range(0, len(W)):
                    domain_size = domain_size + W[i][0] * J0[i][0]
        return domain_size

    def InspectQuadTreeSubCell(self, bulk_elements, qt_depth):
        nsampling = self.params["number_of_samplings"] if ("number_of_samplings" in self.params) else 1
        self.CreateForestSubCell(bulk_elements, nsampling)

        cut_elems = []

        for qi, qs in self.forest.iteritems():
            elem = qs.GetElement()
            stat = self.brep.CutStatusBySampling(elem, nsampling)
            elem.SetValue(CUT_STATUS, stat)
            if stat == self.brep._CUT:
                cut_elems.append(qs)

        aux_util = FiniteCellAuxilliaryUtility()
        for i in range(0, qt_depth):
#            aux_util.MultithreadedQuadTreeSubCellRefineBy(cut_elems, self.brep)
            aux_util.MultithreadedRefineBy(cut_elems, self.brep)

        for qs in cut_elems:
            print("checking quadtree subcell of element " + str(qs.GetElement().Id) + ":")
            print("    number of subcells: " + str(qs.NumberOfSubCells()))
            for i in range(0, qs.NumberOfSubCells()):
                qt = qs.CreateQuadTree(i)
                stat = self.brep.CutStatusBySampling(qt.pGetGeometry(), nsampling)
                print("    subcell " + str(i) + " info:")
                if stat == self.brep._CUT:
                    print("        cut: yes")
                elif stat == self.brep._OUT:
                    print("        cut: out")
                else:
                    print("        cut: in")
                print("        domain size: " + str(qt.DomainSize(self.brep)))
                print("        center of gravity: " + str(qt.CenterOfGravity(self.brep)))
                print("        occupied geometry:")
                aux_util.Print(qt.pCreateGeometry())

    def ExportQuadTreeSubCell(self, model, bulk_elements, qt_depth, sample_element_name, sample_cond_name, integrator_quadrature_method, selected_elems = "all"):
        nsampling = self.params["number_of_samplings"] if ("number_of_samplings" in self.params) else 1
        self.CreateForestSubCell(bulk_elements, nsampling)

        cut_elems = []

        if selected_elems == "all":
            for qi, qs in self.forest.iteritems():
                elem = qs.GetElement()
                stat = self.brep.CutStatusBySampling(elem, nsampling)
                elem.SetValue(CUT_STATUS, stat)
                if stat == self.brep._CUT:
                    cut_elems.append(qs)
        else:
            for qi, qs in self.forest.iteritems():
                elem = qs.GetElement()
                if elem.Id in selected_elems:
                    stat = self.brep.CutStatusBySampling(elem, nsampling)
                    elem.SetValue(CUT_STATUS, stat)
                    if stat == self.brep._CUT:
                        cut_elems.append(qs)

        aux_util = FiniteCellAuxilliaryUtility()
        for i in range(0, qt_depth):
#            aux_util.MultithreadedQuadTreeSubCellRefineBy(cut_elems, self.brep)
            aux_util.MultithreadedRefineBy(cut_elems, self.brep)

        # add to the model_part
        cog_points = []
        for qs in cut_elems:
            lastNodeId = aux_util.GetLastNodeId(model.model_part)
            lastElementId = aux_util.GetLastElementId(model.model_part)
            # print("qs.NumberOfSubCells():", qs.NumberOfSubCells())
            qs.DeepAddToModelPart(model.model_part, sample_element_name, lastNodeId, lastElementId)
            for i in range(0, qs.NumberOfSubCells()):
                qt = qs.CreateQuadTree(i)
                cog_points.append(qt.CenterOfGravity(self.brep, integrator_quadrature_method)) # remark: the zero cog still be exported, that happen when the integrator can't integrate deep enough to find the point

        # print("len(cog_points):", len(cog_points))
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
        nsampling = self.params["number_of_samplings"] if ("number_of_samplings" in self.params) else 1
        self.CreateForestSubCell(bulk_elements, nsampling)

        cut_elems = []

        for qi, qs in self.forest.iteritems():
            elem = qs.GetElement()
            stat = self.brep.CutStatusBySampling(elem, nsampling)
            elem.SetValue(CUT_STATUS, stat)
            if stat == self.brep._CUT:
                cut_elems.append(qs)

        aux_util = FiniteCellAuxilliaryUtility()
        for i in range(0, qt_depth):
#            aux_util.MultithreadedQuadTreeSubCellRefineBy(cut_elems, self.brep)
            aux_util.MultithreadedRefineBy(cut_elems, self.brep)

        Hf = HeavisideFunctionR3R1(self.brep)

        print("checking subcell")
        for qs in cut_elems:
            for i in range(0, qs.NumberOfSubCells()):
                qt = qs.CreateQuadTree(i)
                stat = self.brep.CutStatusBySampling(qt.pCreateGeometry(), nsampling)
                if stat == self.brep._CUT:
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

