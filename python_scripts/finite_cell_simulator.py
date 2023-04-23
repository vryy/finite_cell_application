import math
import pprint
import time as time_module
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.BRepApplication import *
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

def CreateQuadTree(element, nsampling, quadtree_frame = 0, quadtree_configuration = 0):
    if quadtree_frame == 0:
        if nsampling == 1:
            return QuadTreeLocalReference(element, quadtree_configuration)
        elif nsampling == 2:
            return QuadTreeLocalReference2(element, quadtree_configuration)
        elif nsampling == 3:
            return QuadTreeLocalReference3(element, quadtree_configuration)
        elif nsampling == 4:
            return QuadTreeLocalReference4(element, quadtree_configuration)
        elif nsampling == 5:
            return QuadTreeLocalReference5(element, quadtree_configuration)
        elif nsampling == 6:
            return QuadTreeLocalReference6(element, quadtree_configuration)
        elif nsampling == 7:
            return QuadTreeLocalReference7(element, quadtree_configuration)
        elif nsampling == 8:
            return QuadTreeLocalReference8(element, quadtree_configuration)
        elif nsampling == 9:
            return QuadTreeLocalReference9(element, quadtree_configuration)
        elif nsampling == 10:
            return QuadTreeLocalReference10(element, quadtree_configuration)
    elif quadtree_frame == 1:
        if nsampling == 1:
            return QuadTreeLocalCurrent(element, quadtree_configuration)
        elif nsampling == 2:
            return QuadTreeLocalCurrent2(element, quadtree_configuration)
        elif nsampling == 3:
            return QuadTreeLocalCurrent3(element, quadtree_configuration)
        elif nsampling == 4:
            return QuadTreeLocalCurrent4(element, quadtree_configuration)
        elif nsampling == 5:
            return QuadTreeLocalCurrent5(element, quadtree_configuration)
        elif nsampling == 6:
            return QuadTreeLocalCurrent6(element, quadtree_configuration)
        elif nsampling == 7:
            return QuadTreeLocalCurrent7(element, quadtree_configuration)
        elif nsampling == 8:
            return QuadTreeLocalCurrent8(element, quadtree_configuration)
        elif nsampling == 9:
            return QuadTreeLocalCurrent9(element, quadtree_configuration)
        elif nsampling == 10:
            return QuadTreeLocalCurrent10(element, quadtree_configuration)

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
        if params == None:
            self.params = {}
        else:
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

        if 'write_quadrature_to_file' not in self.params:
            self.params['write_quadrature_to_file'] = False
            self.params['export_quadrature_in_reference_frame'] = False
            self.params['export_quadrature_in_current_frame'] = False
        else:
            if 'export_quadrature_in_reference_frame' not in self.params:
                self.params['export_quadrature_in_reference_frame'] = False
            if 'export_quadrature_in_current_frame' not in self.params:
                self.params['export_quadrature_in_current_frame'] = False

        if 'number_of_samplings' not in self.params:
            self.params['number_of_samplings'] = 1

        if 'quadtree_configuration' not in self.params:
            self.params['quadtree_configuration'] = 0

        if 'quadtree_frame' not in self.params:
            self.params['quadtree_frame'] = 0

        if 'export_physical_integration_point' not in self.params:
            self.params['export_physical_integration_point'] = False

        if 'export_fictitious_integration_point' not in self.params:
            self.params['export_fictitious_integration_point'] = False

        self.forest = {}

        self.aux_util = FiniteCellAuxiliaryUtility()

    ###TREES & FOREST CREATION#############
    def CleanForest(self):
        self.forest = {}

    def CreateForest(self, elements, nsampling, quadtree_frame, quadtree_configuration):
        print("###############################################################")
        print("######### CREATE TREES AND FOREST FOR ADAPTIVE QUADRATURE")
        print("Frame of the quadtree: " + str(quadtree_frame))
        print("###############################################################")
        for elem in elements:
            self.forest[elem.Id] = CreateQuadTree(elem, nsampling, quadtree_frame, quadtree_configuration)

    def CreateForestSubCell(self, elements, nsampling = 1):
        subcell_fit_mode = self.params["subcell_fit_mode"]
        cut_cell_quadrature_method = self.params["cut_cell_quadrature_method"]
        quad_util = QuadratureUtility()
        for elem in elements:
            self.forest[elem.Id] = CreateQuadTreeSubCell(elem, nsampling)
            if   subcell_fit_mode == "subcell fit gauss":
                self.forest[elem.Id].ConstructSubCellsBasedOnGaussQuadrature(cut_cell_quadrature_method)
            elif subcell_fit_mode == "subcell fit equal-dist":
                self.forest[elem.Id].ConstructSubCellsBasedOnEqualDistribution(cut_cell_quadrature_method)
            elif subcell_fit_mode == "subcell fit two-layer": # for backward compatibility
                self.forest[elem.Id].ConstructSubCellsBasedOnGaussQuadrature(cut_cell_quadrature_method)
            elif subcell_fit_mode == "subcell nonfit": # for backward compatibility
                pass
            else:
                print("unknown subcell_fit_mode", subcell_fit_mode)
                sys.exit(0)

    ###FITTING DRIVER#############
    def MomentFit(self, bulk_elements, process_info = ProcessInfo()):
        print("fit parameters:")
        pprint.pprint(self.params)
        qt_depth = self.params["qt_depth"]
        nsampling = self.params["number_of_samplings"]
        quadtree_frame = self.params["quadtree_frame"]
        quadtree_configuration = self.params["quadtree_configuration"]
        ###########################################
        ##QUADTREE QUADRATURE
        ###########################################
        self.CreateForest(bulk_elements, nsampling, quadtree_frame, quadtree_configuration)
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
        configuration = 0 # reference configuration

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
                stat = self.brep.CutStatusBySampling(elem, nsampling, configuration)
                if stat == BRep._CUT:
                    for i in range(0, qt_depth):
                        qt.RefineBy(self.brep)
                    fit_util.FitQuadrature(elem, fit_funcs, self.brep, qt, cut_cell_quadrature_method, integrator_quadrature_method, solver_type, echo_level, small_weight)
                    cut_elems.append(elem)
                elif stat == BRep._OUT:
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
                stat = self.brep.CutStatusBySampling(elem, nsampling, configuration)
                if stat == BRep._CUT:
                    cut_elems.append(elem)
                    integrators.append(qt)
                elif stat == BRep._OUT:
                    exclude_elems.append(elem)
            for i in range(0, qt_depth):
                print("Refine level " + str(i+1) + " ...")
#                self.aux_util.MultithreadedQuadTreeRefineBy(integrators, self.brep)
                self.aux_util.MultithreadedRefineBy(integrators, self.brep)
                print("###################################################")
            fit_util.MultithreadedFitQuadrature(cut_elems, fit_funcs, self.brep, integrators, cut_cell_quadrature_method, integrator_quadrature_method, solver_type, echo_level, small_weight)
        else:
            print("Unknown fit_mode", fit_mode)
            sys.exit(0)

        print("construct quadrature using moment fitting successfully")
        if self.params["write_quadrature_to_file"] == True:
            quad_filename = self.params["quad_filename"]
            quad_filetype = self.params["quad_filetype"]
            accuracy = self.params["quad_accuracy"]
            fit_util.SaveQuadrature(quad_filename, quad_filetype, cut_elems, exclude_elems, accuracy)
            fid = open(quad_filename, "a")
            fid.write("\ndef GetParameter():\n")
            fid.write("    params = " + str(self.params) + "\n")
            fid.write("    return params\n")
            fid.close()
        else:
            quad_data = {}
            for elem in cut_elems:
                elem.Initialize(process_info)
                ipoints = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_LOCAL, process_info)
                iweights = elem.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                # print(ipoints)
                # print(iweights)
                data = []
                for i in range(0, len(iweights)):
                    data.append([ipoints[i][0], ipoints[i][1], ipoints[i][2], iweights[i][0]])
                quad_data[elem.Id] = data
            return quad_data

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
        start_fit_t0 = time_module.time()
        model_part.ProcessInfo[RESET_CONFIGURATION] = 0 # this to ensure that the default initialization behaviour is invoked
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
        fit_mode = self.params["fit_mode"]
        fit_funcs = CreateFittingFunctions(fit_space_dim, fit_degree)
        print("list of fitting functions for fit_degree = " + str(fit_degree) + ":")
        for func in fit_funcs:
            print(func)
        configuration = 0 # reference configuration

        cut_qs_elems = []
        exclude_qs_elems = []
        small_qs_elems = []

        integrator_quadrature_method = self.params["integrator_quadrature_method"]
        solver_type = self.params["fit_solver_type"]
        echo_level = self.params["fit_echo_level"]
        fit_small_weight = self.params["fit_small_weight"]
        action_on_small_cut_cell = self.params["action_on_small_cut_cell"]
        fict_weight = self.params["fictitious_weight"] if ("fictitious_weight" in self.params) else 0.0

        if fit_mode == 'serial':
            subcell_fit_func_callback = getattr(fit_util, "FitQuadratureSubCell")
            small_subcell_fit_func_callback = getattr(fit_util, "FitQuadratureSubCellUnique")
        elif fit_mode == "multithread":
            subcell_fit_func_callback = getattr(fit_util, "MultithreadedFitQuadratureSubCell")
            small_subcell_fit_func_callback = getattr(fit_util, "MultithreadedFitQuadratureSubCellUnique")
        generate_physical_integration_points_callback = getattr(self.aux_util, "MultithreadedGeneratePhysicalIntegrationPoints")
        save_quadrature_subcell_callback = getattr(fit_util, "SaveQuadratureSubCell")

        for qi, qs in self.forest.iteritems():
            elem = qs.GetElement()
            stat = self.brep.CutStatusBySampling(elem, nsampling, configuration)
            elem.SetValue(CUT_STATUS, stat)
            if stat == BRep._CUT:
                cut_qs_elems.append(qs)
            elif stat == BRep._OUT:
                exclude_qs_elems.append(qs)
        for i in range(0, qt_depth):
            print("Refine level " + str(i+1) + " ...")
#            self.aux_util.MultithreadedQuadTreeSubCellRefineBy(cut_qs_elems, self.brep)
            self.aux_util.MultithreadedRefineBy(cut_qs_elems, self.brep)
            print("###################################################")
        # recursive refine the subcell that are cut but not integrable by the quadtree # this is required if Gauss-Legendre quadrature is used
        if max_qt_depth > qt_depth:
            for qs in cut_qs_elems:
                for i in range(0, qs.NumberOfSubCells()):
                    qt = qs.CreateQuadTree(i)
                    stat = self.brep.CutStatusBySampling(qt.pCreateGeometry(), nsampling, configuration)
                    if stat == BRep._CUT:
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

        start_fit_t1 = time_module.time()
        if fit_mode == 'serial':
            for qs in proper_cut_qs_elems:
                subcell_fit_func_callback(qs, fit_funcs, self.brep, integrator_quadrature_method, solver_type, echo_level, fit_small_weight, model_part.ProcessInfo)
        elif fit_mode == 'multithread':
            subcell_fit_func_callback(proper_cut_qs_elems, fit_funcs, self.brep, integrator_quadrature_method, solver_type, echo_level, fit_small_weight, model_part.ProcessInfo)
        end_fit_t1 = time_module.time()
        print("construct quadrature using moment fitting for subcell successfully in " + str(end_fit_t1 - start_fit_t1) + " s")
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
            lastNodeId = self.aux_util.GetLastNodeId(model_part)
            lastElementId = self.aux_util.GetLastElementId(model_part)
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
                lastNodeId = self.aux_util.GetLastNodeId(model_part)
                lastElementId = self.aux_util.GetLastElementId(model_part)
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
            small_subcell_fit_func_callback(small_elems, fit_funcs, self.brep, small_qs_elems, small_cut_cell_quadrature_method, integrator_quadrature_method, small_fit_solver_type, small_fit_echo_level, small_fit_small_weight, model_part.ProcessInfo)
            save_quadrature_subcell_callback(quad_filename, quad_filetype, proper_cut_qs_elems, exclude_qs_elems, small_qs_elems, accuracy)
            print("generate moment-fit quadrature for small cell completed")

        end_fit_t0 = time_module.time()
        print("Total fitting time: " + str(end_fit_t0 - start_fit_t0))
        fid = open(quad_filename, "a")
        fid.write("\ndef GetParameter():\n")
        fid.write("    params = " + str(self.params) + "\n")
        fid.write("    return params\n")
        fid.close()

    ###SIMULATION DRIVER#############
    def Initialize(self, model, bulk_elements):
        """
        Initialize the finite cell model
        """
        self.quadrature_method = self.params["quadrature_method"]
        self.mpu = self.params["material_properties_utility"]
        nsampling = self.params["number_of_samplings"]
        print("simulator parameters:")
        pprint.pprint(self.params)
        quad_util = QuadratureUtility()
        mesh_util = FiniteCellMeshUtility()

        # this to ensure that the default initialization behaviour is set
        model.model_part.ProcessInfo[RESET_CONFIGURATION] = 0

        if self.quadrature_method == "quadtree":
            qt_depth = self.params["qt_depth"]
            small_weight = self.params["small_weight"]
            if "set_small_weight" in self.params:
                set_small_weight = self.params["set_small_weight"]
            else:
                set_small_weight = True # set the small weight by default
            configuration = 0 # reference configuration
            quadtree_frame = self.params['quadtree_frame']
            quadtree_configuration = self.params['quadtree_configuration']
            ## parameters for small cut cell handling
            if 'action_on_small_cut_cell' in self.params:
                action_on_small_cut_cell = self.params['action_on_small_cut_cell']
            else:
                action_on_small_cut_cell = None
            ##
            if 'small_domain_size' in self.params:
                small_domain_size = self.params["small_domain_size"]
            else:
                small_domain_size = 0.0
            ##
            if 'small_domain_size_ratio' in self.params:
                small_domain_size_ratio = self.params["small_domain_size_ratio"]
            else:
                small_domain_size_ratio = 0.0
            ##
            if 'number_of_minimum_physical_points' in self.params:
                number_of_minimum_physical_points = self.params["number_of_minimum_physical_points"]
            else:
                number_of_minimum_physical_points = 0
            ##
            if small_domain_size > 0.0:
                criteria_check_small_domain_size = True
            else:
                criteria_check_small_domain_size = False
            #
            if small_domain_size_ratio > 0.0:
                criteria_check_small_domain_size_ratio = True
            else:
                criteria_check_small_domain_size_ratio = False
            #
            if number_of_minimum_physical_points != 0:
                criteria_check_number_of_physical_points = True
            else:
                criteria_check_number_of_physical_points = False
            ##
            ###########################################
            ##QUADTREE
            ###########################################
            print("criteria_check_small_domain_size: " + str(criteria_check_small_domain_size))
            print("criteria_check_small_domain_size_ratio: " + str(criteria_check_small_domain_size_ratio))
            print("criteria_check_number_of_physical_points: " + str(criteria_check_number_of_physical_points))
            self.CreateForest(bulk_elements, nsampling, quadtree_frame, quadtree_configuration)
            cut_cell_quadrature_method = self.params["cut_cell_quadrature_method"]
            cut_cell_quadrature_order = quad_util.GetQuadratureOrder(cut_cell_quadrature_method)
            total_add_quadrature_pnts = 0
            self.proper_cut_elems = ElementsArray()
            cut_elems = []
            exclude_elems = []
            for elem in bulk_elements:
                qt = self.forest[elem.Id]
                if criteria_check_small_domain_size_ratio:
                    cell_size = qt.DomainSize()
                stat = self.brep.CutStatusBySampling(elem, nsampling, configuration)
                elem.SetValue(CUT_STATUS, stat)
                if stat == BRep._CUT:
                    for i in range(0, qt_depth):
                        qt.RefineBy(self.brep)
                    if set_small_weight:
                        np = qt.ConstructQuadrature(self.brep, cut_cell_quadrature_method, small_weight)
                    else:
                        np = qt.ConstructQuadrature(cut_cell_quadrature_method)
                    print("cut element " + str(elem.Id) + " has " + str(np) + " quadrature points")
                    print("cut element " + str(elem.Id) + " has " + str(qt.NumberOfCells()) + " leaf cells")
                    total_add_quadrature_pnts = total_add_quadrature_pnts + np[0]
                    include_element = True
                    if action_on_small_cut_cell != None:
                        # here we check the small cutcell; the small subcell is determined by either domain_size less than a threshold, or number of physical points less than specific number

                        if criteria_check_small_domain_size:
                            domain_size = qt.DomainSize(self.brep, cut_cell_quadrature_method)
                            if domain_size < small_domain_size:
                                include_element = False

                        if criteria_check_small_domain_size_ratio:
                            domain_size = qt.DomainSize(self.brep, cut_cell_quadrature_method)
                            print("element " + str(elem.Id) + " cut ratio: " + str(domain_size/cell_size))
                            if domain_size/cell_size < small_domain_size_ratio:
                                include_element = False

                        if criteria_check_number_of_physical_points:
                            num_physical_points = np[0]
                            if num_physical_points < number_of_minimum_physical_points:
                                include_element = False

                    if include_element:
                        elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                        elem.Initialize(model.model_part.ProcessInfo)
                        self.aux_util.AddElement(self.proper_cut_elems, elem)
                        cut_elems.append(elem)
                    else:
                        if action_on_small_cut_cell == 'eliminate':
                            elem.SetValue(ACTIVATION_LEVEL, -1)
                            elem.SetValue(IS_INACTIVE, True)
                            elem.Set(ACTIVE, False)
                            exclude_elems.append(elem)
                            print("element " + str(elem.Id) + " is detected as small cut cell. It will be deactivated.")
                        else:
                            raise Exception("Invalid action_on_small_cut_cell " + str(action_on_small_cut_cell))

                elif stat == BRep._OUT:
                    elem.SetValue(ACTIVATION_LEVEL, -1)
                    elem.SetValue(IS_INACTIVE, True)
                    elem.Set(ACTIVE, False)
                    exclude_elems.append(elem)
                    print("element " + str(elem.Id) + " is outside the physical domain. It will be deactivated.")
            print("construct quadrature using quad-tree successfully, " + str(total_add_quadrature_pnts) + " quadrature points are added")
            for elem in bulk_elements:
                self.mpu.SetMaterialProperties(model.model_part, elem)
                physical_constitutive_laws = elem.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, model.model_part.ProcessInfo)
                for i in range(0, len(physical_constitutive_laws)):
                    physical_constitutive_laws[i].SetValue(PARENT_ELEMENT_ID, elem.Id, model.model_part.ProcessInfo)
                    physical_constitutive_laws[i].SetValue(INTEGRATION_POINT_INDEX, i, model.model_part.ProcessInfo)
            print("reading material completed")

            ## export the physical integration points for debugging if needed
#             if self.params["export_physical_integration_point"]:
# #                cog_points = []
# #                for elem in self.proper_cut_elems:
# #                    points = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
# #                    for point in points:
# #                        cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))
# #                print("len(cog_points):", len(cog_points))
# #                prop_id = self.params["physical_integration_point_prop_id"]
# #                quad_util.CreateConditionFromPoint(model.model_part, cog_points, "DummyConditionPoint3D", model.model_part.Properties[prop_id])
#                 self.ExportPhysicalIntegrationPoints(model.model_part, self.proper_cut_elems, self.params["physical_integration_point_prop_id"])

            if self.params["export_physical_integration_point"] or self.params['export_fictitious_integration_point']:
                cog_points = []
                fict_cog_points = []
                for elem in self.proper_cut_elems:
                    points = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL_IN_REFERENCE_CONFIGURATION, model.model_part.ProcessInfo)
                    weights = elem.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, model.model_part.ProcessInfo)
                    num_points = len(points)
                    for i in range(0, num_points):
                        point = points[i]
                        weight = weights[i][0]
                        if weight > 1.1*small_weight: # coefficient 1.1 is used to improve the stability of the numerical comparison
                            cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))
                        else:
                            fict_cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))
                print("len(cog_points):", len(cog_points))
                print("len(fict_cog_points):", len(fict_cog_points))
                if self.params["export_physical_integration_point"]:
                    prop_id = self.params["physical_integration_point_prop_id"]
                    cog_conds = quad_util.CreateConditionFromPoint(model.model_part, cog_points, "DummyConditionPoint3D", model.model_part.Properties[prop_id])
                    # model.model_part.Properties[prop_id].SetValue(LAYER_NAME, "physical quadrature points")
                if self.params["export_fictitious_integration_point"]:
                    prop_id = self.params["fictitious_integration_point_prop_id"]
                    fict_cog_conds = quad_util.CreateConditionFromPoint(model.model_part, fict_cog_points, "DummyConditionPoint3D", model.model_part.Properties[prop_id])
                    # model.model_part.Properties[prop_id].SetValue(LAYER_NAME, "fictitious quadrature points")

            ## export the quadtree for debugging if needed
            if self.params["export_quadtree_cell"]:
                lastNodeId = self.aux_util.GetLastNodeId(model.model_part)
                lastElementId = self.aux_util.GetLastElementId(model.model_part)
                sample_quad_element_name = self.params["sample_quad_element_name"]
                if not ("starting_level" in self.params):
                    starting_level = 1
                else:
                    starting_level = self.params["starting_level"]
                for elem in bulk_elements:
                    qt = self.forest[elem.Id]
                    if starting_level == 1:
                        last_ids = qt.AddToModelPart(model.model_part, sample_quad_element_name, lastNodeId, lastElementId)
                    else:
                        last_ids = qt.AddToModelPartWithLevel(model.model_part, sample_quad_element_name, lastNodeId, lastElementId, starting_level)
                    lastNodeId = last_ids[0]
                    lastElementId = last_ids[1]

            ## export the quadrature if user needs
            if self.params["write_quadrature_to_file"] == True:
                quad_filename = self.params["quad_filename"]
                quad_filetype = self.params["quad_filetype"]
                accuracy = self.params["quad_accuracy"]
                quad_util.SaveQuadrature(quad_filename, quad_filetype, cut_elems, exclude_elems, accuracy)
                fid = open(quad_filename, "a")
                fid.write("\ndef GetParameter():\n")
                fid.write("    params = " + str(self.params) + "\n")
                fid.write("    return params\n")
                fid.close()
                if self.params["export_quadrature_in_reference_frame"] == True:
                    quad_filename = self.params["quad_filename_reference"]
                    quad_filetype = self.params["quad_filetype_reference"]
                    accuracy = self.params["quad_accuracy"]
                    quad_util.ExportQuadratureInReferenceFrame(quad_filename, quad_filetype, cut_elems, exclude_elems, accuracy)
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
            for elem in model.model_part.Elements:
                elem.SetValue(CUT_STATUS, BRep._IN)
            for elem_id in cut_elems:
                elem = model.model_part.Elements[elem_id]
                quad_util.SetQuadrature(elem, cut_cell_quadrature_order, cut_cell_quadrature[elem_id])
                total_add_quadrature_pnts = total_add_quadrature_pnts + len(cut_cell_quadrature[elem_id])
                elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                elem.SetValue(CUT_STATUS, BRep._CUT)
                elem.Initialize(model.model_part.ProcessInfo)
                self.aux_util.AddElement(self.proper_cut_elems, elem)
            for elem_id in exclude_elems:
                elem = model.model_part.Elements[elem_id]
                elem.SetValue(ACTIVATION_LEVEL, -1)
                elem.SetValue(IS_INACTIVE, True)
                elem.SetValue(CUT_STATUS, BRep._OUT)
                elem.Set(ACTIVE, False)
            print("obtain quadtree quadrature successfully, " + str(total_add_quadrature_pnts) + " quadrature points are added")
            for elem_id in cut_elems:
                elem = model.model_part.Elements[elem_id]
                self.mpu.SetMaterialProperties(model.model_part, elem)
                physical_constitutive_laws = elem.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, model.model_part.ProcessInfo)
                for i in range(0, len(physical_constitutive_laws)):
                    physical_constitutive_laws[i].SetValue(PARENT_ELEMENT_ID, elem.Id, model.model_part.ProcessInfo)
                    physical_constitutive_laws[i].SetValue(INTEGRATION_POINT_INDEX, i, model.model_part.ProcessInfo)
            print("reading material completed")
            if self.params["export_physical_integration_point"]:
                cog_points = []
                for elem in self.proper_cut_elems:
                    points = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
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
            for elem in model.model_part.Elements:
                elem.SetValue(CUT_STATUS, BRep._IN)
            for elem_id in cut_elems:
                elem = model.model_part.Elements[elem_id]
                quad_util.SetQuadrature(elem, cut_cell_quadrature_order, cut_cell_quadrature[elem_id])
                elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                elem.SetValue(CUT_STATUS, BRep._CUT)
                elem.Initialize(model.model_part.ProcessInfo)
                self.aux_util.AddElement(self.proper_cut_elems, elem)
            for elem_id in exclude_elems:
                elem = model.model_part.Elements[elem_id]
                elem.SetValue(ACTIVATION_LEVEL, -1)
                elem.SetValue(IS_INACTIVE, True)
                elem.SetValue(CUT_STATUS, BRep._OUT)
                elem.Set(ACTIVE, False)
            print("obtain moment-fit quadrature successfully")
            for elem in bulk_elements:
                self.mpu.SetMaterialProperties(model.model_part, elem)
                physical_constitutive_laws = elem.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, model.model_part.ProcessInfo)
                for i in range(0, len(physical_constitutive_laws)):
                    physical_constitutive_laws[i].SetValue(PARENT_ELEMENT_ID, elem.Id, model.model_part.ProcessInfo)
                    physical_constitutive_laws[i].SetValue(INTEGRATION_POINT_INDEX, i, model.model_part.ProcessInfo)
            print("reading material completed")
            if self.params["export_physical_integration_point"]:
                cog_points = []
                for elem in self.proper_cut_elems:
                    points = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
                    for point in points:
                        cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))
                print("len(cog_points):", len(cog_points))
                prop_id = self.params["physical_integration_point_prop_id"]
                quad_util.CreateConditionFromPoint(model.model_part, cog_points, "DummyConditionPoint3D", model.model_part.Properties[prop_id])
            # end elif self.quadrature_method == "moment-fit quadtree":
        elif self.quadrature_method == "moment-fit subcell":
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
            lastElementId = self.aux_util.GetLastElementId(model.model_part)
            lastCondId = self.aux_util.GetLastConditionId(model.model_part)
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
            number_of_physical_integration_points = 0
            for elem in bulk_elements:
                if elem.Id in cut_elems:
                    print("at cut element " + str(elem.Id))
                    cut_cell_quadrature_order = self.forest[elem.Id].GetRepresentativeIntegrationOrder()
                    elem.SetValue(CUT_STATUS, BRep._CUT)
                    if len(cut_cell_quadrature[elem.Id]) == 0: # if cut cell quadrature is empty then skip
                        elem.SetValue(ACTIVATION_LEVEL, -1)
                        elem.SetValue(IS_INACTIVE, True)
                        elem.Set(ACTIVE, False)
                        print("element " + str(elem.Id) + " has no physical integration points. It will be excluded.")
                    else:
                        self.aux_util.AddElement(self.proper_cut_elems, elem)

                        quad_util.SetQuadrature(elem, cut_cell_quadrature_order, cut_cell_quadrature[elem.Id])
                        elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                        elem.Initialize(model.model_part.ProcessInfo)

                        output = self.forest[elem.Id].CreateSubCellElements(model.model_part, subcell_element_type, cut_cell_quadrature_order, cut_cell_full_quadrature[elem.Id], subcell_weights[elem.Id], lastElementId, lastCondId, extrapolated_prop)
                        lastElementId = output[0]
                        lastCondId = output[1]
                        new_sub_elems = output[2] # !!!IMPORTANT!!! here it's assumed that the sub-elements are Extrapolated ones

                        for sub_elem in new_sub_elems:
                            sub_elem.SetValue(IS_INACTIVE, False)
                            sub_elem.Set(ACTIVE, True)
                            sub_elem.SetValue(PARENT_ELEMENT_ID, elem.Id)

                        self.mpu.SetMaterialProperties(model.model_part, elem)
                        elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                        elem.Initialize(model.model_part.ProcessInfo)

                        physical_constitutive_laws = elem.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, model.model_part.ProcessInfo)
                        for i in range(0, len(physical_constitutive_laws)):
                            physical_constitutive_laws[i].SetValue(PARENT_ELEMENT_ID, elem.Id, model.model_part.ProcessInfo)
                            physical_constitutive_laws[i].SetValue(INTEGRATION_POINT_INDEX, i, model.model_part.ProcessInfo)
                            # physical_constitutive_laws[i].SetValue(REPRESENTATIVE_WEIGHT, subcell_domain_sizes[elem.Id][i] / len(cut_cell_full_quadrature[elem.Id]), model.model_part.ProcessInfo) # here we divide to number of integration point of the full quadrature. Because all the subcell element will be added up when calculating the error
                        print("len(physical_constitutive_laws)", len(physical_constitutive_laws))
                        print("len(new_sub_elems)", len(new_sub_elems))
                        number_of_physical_integration_points = number_of_physical_integration_points + len(physical_constitutive_laws)

                        if extrapolation_mode == "extrapolated element" or extrapolation_mode == "extrapolated constant stress element":
                            physical_integration_points = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_LOCAL, model.model_part.ProcessInfo)
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
                            points = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
                            for point in points:
                                cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))

                    if self.params["add_fictitious_element"]:
                        if len(cut_cell_fictitious_quadrature[elem.Id]) > 0:
                            cut_cell_quadrature_order = 1
                            if fict_prop == "same as parent element": # this uses the material proper from the element, even it's nonlinear
                                fict_elem = mesh_util.CreateParasiteElement(sample_fictitious_element_name, lastElementId, elem, cut_cell_quadrature_order, cut_cell_fictitious_quadrature[elem.Id], elem.Properties)
                                fict_elem.Initialize(model.model_part.ProcessInfo)
                                self.mpu.SetMaterialProperties(model.model_part, fict_elem)
                                fict_elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                                fict_elem.SetValue(IS_INACTIVE, False)
                                fict_elem.Set(ACTIVE, True)
                                fict_elem.SetValue(PARENT_ELEMENT_ID, elem.Id)
                                fict_elem.Initialize(model.model_part.ProcessInfo)
                            else: # this allows to set the material property of fictitious element DIFFERENT to the parent element
                                fict_elem = mesh_util.CreateParasiteElement(sample_fictitious_element_name, lastElementId, elem, cut_cell_quadrature_order, cut_cell_fictitious_quadrature[elem.Id], fict_prop)
                            print("fictitious element " + str(fict_elem.Id) + " is created out from element " + str(elem.Id))
                            self.aux_util.AddElement(self.fict_elems, fict_elem)
                            lastElementId = fict_elem.Id

                            if self.params["export_fictitious_integration_point"]:
                                points = fict_elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
                                for point in points:
                                    fict_cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))

                elif elem.Id in exclude_elems:
                    elem.SetValue(ACTIVATION_LEVEL, -1)
                    elem.SetValue(IS_INACTIVE, True)
                    elem.SetValue(CUT_STATUS, BRep._OUT)
                    elem.Set(ACTIVE, False)
                    print("element " + str(elem.Id) + " is outside. It will be excluded")

                elif elem.Id in quadtree_elems:
                    cut_cell_quadrature_order = 1
                    quad_util.SetQuadrature(elem, cut_cell_quadrature_order, quadtree_quadrature[elem.Id])
                    elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                    elem.SetValue(CUT_STATUS, BRep._CUT)
                    elem.Initialize(model.model_part.ProcessInfo)
                    self.mpu.SetMaterialProperties(model.model_part, elem)
#                    self.aux_util.AddElement(self.proper_cut_elems, elem) # a quadtree element is NOT considered as proper_cut_elems
                    print("element " + str(elem.Id) + " is assigned " + str(len(quadtree_quadrature[elem.Id])) + " quadrature points")
                    if self.params["export_physical_integration_point"]:
                        points = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
                        for point in points:
                            cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))
                else:
                    elem.SetValue(CUT_STATUS, BRep._IN)
                    self.mpu.SetMaterialProperties(model.model_part, elem)
            print("obtain moment-fit subcell quadrature successfully")
            print("number of physical integration points: " + str(number_of_physical_integration_points))
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

        if self.params["enable_ghost_penalty"] and self.params["enable_skeleton_penalty"]:
            print("Both ghost_penalty and skeleton_penalty are activated. You should choose only one.")
            sys.exit(0)

        if self.params["enable_ghost_penalty"] or self.params["enable_skeleton_penalty"]:
            space_dim = self.params["space_dim"]
            estimated_neighbours = self.params["estimated_number_of_neighbours"]

            FindElementalNeighboursProcess(model.model_part, space_dim, estimated_neighbours).Execute()

            sample_cond = self.params["sample_ghost_penalty_condition"]
            # ghost_prop = self.params["ghost_penalty_properties"]
            ghost_prop_id = self.params["ghost_penalty_prop_id"]
            ghost_prop = model.model_part.Properties[ghost_prop_id]
            self.params['ghost_penalty_setting_callback'](ghost_prop)
            lastCondId = self.aux_util.GetLastConditionId(model.model_part)
            if "ghost_echo_level" in self.params:
                ghost_echo_level = self.params["ghost_echo_level"]
            else:
                ghost_echo_level = 1
            ghost_penalty_conds = ghost_penalty_util.SetUpSurfacePenaltyConditions(model.model_part, bulk_elements, sample_cond, self.brep, lastCondId, ghost_prop, ghost_echo_level)
            for cond in ghost_penalty_conds:
                cond.Initialize(model.model_part.ProcessInfo)

    # end Initialize

    ### Utility function to export the integration points to the mesh for post-processing ###
    def ExportPhysicalIntegrationPoints(self, model_part, elements, prop_id, current_frame = False):
        quad_util = QuadratureUtility()
        cog_points = []
        for elem in elements:
            points = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model_part.ProcessInfo)
            if current_frame:
                disps = elem.CalculateOnIntegrationPoints(DISPLACEMENT, model_part.ProcessInfo)
            for i in range(0, len(points)):
                point = points[i]
                if current_frame:
                    disp = disps[i]
                    point[0] = point[0] + disp[0]
                    point[1] = point[1] + disp[1]
                    point[2] = point[2] + disp[2]
            for point in points:
                cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))
        print("len(cog_points):", len(cog_points))
        quad_util.CreateConditionFromPoint(model_part, cog_points, "DummyConditionPoint3D", model_part.Properties[prop_id])

    ###COMPUTE GLOBAL DISPLACEMENT (L2) ERROR###
    def compute_L2_error(self, elements, process_info, solution, P):
        print("!!!compute_L2_error:Please turn off the MoveMeshFlag in order to have correct results")
        nom = 0.0
        denom = 0.0
        for element in elements:
            if element.Is(ACTIVE):
                u = element.CalculateOnIntegrationPoints(DISPLACEMENT, process_info)
                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
                W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
#                print("number of integration point on element " + str(element.Id) + ":" + str(len(Q)))
                # if element.Id == 6:
                #     print("at element ", element.Id)
                #     print("u:", u)
                #     print("Q:", Q)
                #     print("W:", W)
                #     print("J0:", J0)
                elem_nom = 0.0
                elem_denom = 0.0
                for i in range(0, len(u)):
                    ana_u = solution.get_displacement(P, Q[i][0], Q[i][1], Q[i][2])
                    elem_nom = elem_nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2) + pow(u[i][2] - ana_u[2], 2)) * W[i][0] * J0[i][0]
                    elem_denom = elem_denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2) + pow(ana_u[2], 2)) * W[i][0] * J0[i][0]
                # if (element.GetValue(CUT_STATUS) == BRep._CUT):
                #     print("l2 error on cut element " + str(element.Id) + ": " + str(elem_nom))
                nom = nom + elem_nom
                denom = denom + elem_denom
        print("nom:", nom)
        print("denom:", denom)
        if denom == 0.0:
            if nom == 0.0:
                return 0.0
            else:
                return float('nan');
        else:
            return math.sqrt(abs(nom / denom))

    ###COMPUTE GLOBAL DISPLACEMENT (L2) ERROR FOR MFSC SCHEME###
    def compute_L2_error_mfsc(self, elements, all_subcell_elems, fict_elems, process_info, solution, P):
        print("!!!compute_L2_error_mfsc:Please turn off the MoveMeshFlag in order to have correct results")
        nom = 0.0
        denom = 0.0
        for element in elements:
            if element.Is(ACTIVE):
                u = element.CalculateOnIntegrationPoints(DISPLACEMENT, process_info)
                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
                W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                for i in range(0, len(u)):
                    ana_u = solution.get_displacement(P, Q[i][0], Q[i][1], Q[i][2])
                    nom = nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2) + pow(u[i][2] - ana_u[2], 2)) * W[i][0] * J0[i][0]
                    denom = denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2) + pow(ana_u[2], 2)) * W[i][0] * J0[i][0]
        for element in all_subcell_elems:
            if element.Is(ACTIVE):
                elem_nom = 0.0
                elem_denom = 0.0
                #####method 1
#                 u = element.CalculateOnIntegrationPoints(PHYSICAL_INTEGRATION_POINT_DISPLACEMENT, process_info)
# #                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
#                 Q = element.CalculateOnIntegrationPoints(PHYSICAL_INTEGRATION_POINT_GLOBAL, process_info)
#                 W = element.CalculateOnIntegrationPoints(SUBCELL_DOMAIN_SIZE, process_info)
#                 for i in range(0, len(u)):
#                     ana_u = solution.get_displacement(P, Q[i][0], Q[i][1], Q[i][2])
#                     elem_nom = elem_nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2) + pow(u[i][2] - ana_u[2], 2)) * W[i][0]# * J0[i][0]
#                     elem_denom = elem_denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2) + pow(ana_u[2], 2)) * W[i][0]# * J0[i][0]
                #####method 2
                u = element.CalculateOnIntegrationPoints(DISPLACEMENT, process_info)
                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
                W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                # print("Q:", Q)
                # print("W:", W)
                for i in range(0, len(u)):
                    ana_u = solution.get_displacement(P, Q[i][0], Q[i][1], Q[i][2])
                    elem_nom = elem_nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2) + pow(u[i][2] - ana_u[2], 2)) * W[i][0] * J0[i][0]
                    elem_denom = elem_denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2) + pow(ana_u[2], 2)) * W[i][0] * J0[i][0]
                #############
                nom = nom + elem_nom
                denom = denom + elem_denom
                # print("l2 error on mfsc cut element " + str(element.Id) + ", parent = " + str(element.GetValue(PARENT_ELEMENT_ID)) + ": " + str(elem_nom))
        for element in fict_elems:
            if element.Is(ACTIVE):
                u = element.CalculateOnIntegrationPoints(DISPLACEMENT, process_info)
                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
                W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                elem_nom = 0.0
                elem_denom = 0.0
                for i in range(0, len(u)):
                    ana_u = solution.get_displacement(P, Q[i][0], Q[i][1], Q[i][2])
                    elem_nom = elem_nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2) + pow(u[i][2] - ana_u[2], 2)) * W[i][0] * J0[i][0]
                    elem_denom = elem_denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2) + pow(ana_u[2], 2)) * W[i][0] * J0[i][0]
                nom = nom + elem_nom
                denom = denom + elem_denom
                # print("l2 error on mfsc fictitious element " + str(element.Id) + ", parent = " + str(element.GetValue(PARENT_ELEMENT_ID)) + ": " + str(elem_nom))
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
                o = element.CalculateOnIntegrationPoints(THREED_STRESSES, process_info)
                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
                W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
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

    def compute_H1_error_mfsc(self, elements, all_subcell_elems, fict_elems, process_info, solution, P):
        print("!!!compute_H1_error_mfsc:Please turn off the MoveMeshFlag in order to have correct results")
        nom = 0.0
        denom = 0.0
        for element in elements:
            if element.Is(ACTIVE):
                o = element.CalculateOnIntegrationPoints(THREED_STRESSES, process_info)
                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
                W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                for i in range(0, len(o)):
                    ana_o = solution.get_stress_3d(P, Q[i][0], Q[i][1], Q[i][2])
                    nom = nom + (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2))) * W[i][0] * J0[i][0]
                    denom = denom + (pow(ana_o[0], 2) + pow(ana_o[1], 2) + pow(ana_o[2], 2) + 2.0*(pow(ana_o[3], 2) + pow(ana_o[4], 2) + pow(ana_o[5], 2))) * W[i][0] * J0[i][0]
        for element in all_subcell_elems:
            if element.Is(ACTIVE):
                o = element.CalculateOnIntegrationPoints(PHYSICAL_INTEGRATION_POINT_THREED_STRESSES, process_info)
#                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.CalculateOnIntegrationPoints(PHYSICAL_INTEGRATION_POINT_GLOBAL, process_info)
                W = element.CalculateOnIntegrationPoints(SUBCELL_DOMAIN_SIZE, process_info)
                for i in range(0, len(o)):
                    ana_o = solution.get_stress_3d(P, Q[i][0], Q[i][1], Q[i][2])
                    nom = nom + (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2))) * W[i][0]# * J0[i][0]
                    denom = denom + (pow(ana_o[0], 2) + pow(ana_o[1], 2) + pow(ana_o[2], 2) + 2.0*(pow(ana_o[3], 2) + pow(ana_o[4], 2) + pow(ana_o[5], 2))) * W[i][0]# * J0[i][0]
        for element in fict_elems:
            if element.Is(ACTIVE):
                o = element.CalculateOnIntegrationPoints(THREED_STRESSES, process_info)
                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
                W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
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

    ###COMPUTE DOMAIN SIZE###
    def compute_domain_size(self, elements, process_info):
        domain_size = 0.0
        for element in elements:
            if element.Is(ACTIVE):
                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                elem_size = 0.0
                for i in range(0, len(W)):
                    elem_size = elem_size + W[i][0] * J0[i][0]
                if (element.GetValue(CUT_STATUS) == BRep._CUT):
                    print("domain size of cut element " + str(element.Id) + ": " + str(elem_size))
                else:
                    print("domain size of uncut element " + str(element.Id) + ": " + str(elem_size))
                domain_size = domain_size + elem_size
        return domain_size

    ###COMPUTE DOMAIN SIZE###
    def compute_domain_size_subcell(self, all_subcell_elems, process_info):
        domain_size = 0.0
        for element in all_subcell_elems:
            if element.Is(ACTIVE):
                W = element.CalculateOnIntegrationPoints(SUBCELL_DOMAIN_SIZE, process_info)
                for i in range(0, len(W)):
                    domain_size = domain_size + W[i][0]
        return domain_size

    ###COMPUTE DOMAIN SIZE###
    def compute_domain_size_mfsc(self, elements, all_subcell_elems, fict_elems, process_info):
        domain_size = 0.0
        for element in elements:
            if element.Is(ACTIVE):
                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                elem_size = 0.0
                for i in range(0, len(W)):
                    elem_size = elem_size + W[i][0] * J0[i][0]
                domain_size = domain_size + elem_size
                print("domain size of uncut element " + str(element.Id) + ": " + str(elem_size))
        mfsc_elems = {}
        for element in all_subcell_elems:
            if element.Is(ACTIVE):
                parent_element_id = element.GetValue(PARENT_ELEMENT_ID)
                if parent_element_id not in mfsc_elems:
                    mfsc_elems[parent_element_id] = 0.0
                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                for i in range(0, len(W)):
                    mfsc_elems[parent_element_id] = mfsc_elems[parent_element_id] + W[i][0] * J0[i][0]
        for element in fict_elems:
            if element.Is(ACTIVE):
                parent_element_id = element.GetValue(PARENT_ELEMENT_ID)
                if parent_element_id not in mfsc_elems:
                    mfsc_elems[parent_element_id] = 0.0
                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                for i in range(0, len(W)):
                    mfsc_elems[parent_element_id] = mfsc_elems[parent_element_id] + W[i][0] * J0[i][0]
        for elem_id, elem_size in mfsc_elems.iteritems():
            print("domain size of cut element " + str(elem_id) + ": " + str(elem_size))
            domain_size = domain_size + elem_size
        return domain_size

    ###CHECKING FUNCTIONS#############
    def ExportQuadTree(self, model, sample_element, group = None):
        nsampling = self.params["number_of_samplings"]
        quadtree_frame = self.params["quadtree_frame"]
        quadtree_configuration = self.params["quadtree_configuration"]
        self.CreateForest(model.model_part.Elements, nsampling, quadtree_frame, quadtree_configuration)
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
                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                for i in range(0, len(W)):
                    domain_size = domain_size + W[i][0] * J0[i][0]
        return domain_size

    def InspectQuadTreeSubCell(self, bulk_elements, qt_depth):
        nsampling = self.params["number_of_samplings"] if ("number_of_samplings" in self.params) else 1
        configuration = 1 # we check the cut status in the current configuration
        self.CreateForestSubCell(bulk_elements, nsampling)

        cut_elems = []

        for qi, qs in self.forest.iteritems():
            elem = qs.GetElement()
            stat = self.brep.CutStatusBySampling(elem, nsampling, configuration)
            elem.SetValue(CUT_STATUS, stat)
            if stat == BRep._CUT:
                cut_elems.append(qs)

        for i in range(0, qt_depth):
#            self.aux_util.MultithreadedQuadTreeSubCellRefineBy(cut_elems, self.brep)
            self.aux_util.MultithreadedRefineBy(cut_elems, self.brep)

        for qs in cut_elems:
            print("checking quadtree subcell of element " + str(qs.GetElement().Id) + ":")
            print("    number of subcells: " + str(qs.NumberOfSubCells()))
            for i in range(0, qs.NumberOfSubCells()):
                qt = qs.CreateQuadTree(i)
                stat = self.brep.CutStatusBySampling(qt.pGetGeometry(), nsampling, configuration)
                print("    subcell " + str(i) + " info:")
                if stat == BRep._CUT:
                    print("        cut: yes")
                elif stat == BRep._OUT:
                    print("        cut: out")
                else:
                    print("        cut: in")
                print("        domain size: " + str(qt.DomainSize(self.brep)))
                print("        center of gravity: " + str(qt.CenterOfGravity(self.brep)))
                print("        occupied geometry:")
                self.aux_util.Print(qt.pCreateGeometry())

    def ExportQuadTreeSubCell(self, model, bulk_elements, qt_depth, sample_element_name, sample_cond_name, integrator_quadrature_method, selected_elems = "all"):
        nsampling = self.params["number_of_samplings"] if ("number_of_samplings" in self.params) else 1
        configuration = 1 # we check the cut status in the current configuration
        self.CreateForestSubCell(bulk_elements, nsampling)

        cut_elems = []

        if selected_elems == "all":
            for qi, qs in self.forest.iteritems():
                elem = qs.GetElement()
                stat = self.brep.CutStatusBySampling(elem, nsampling, configuration)
                elem.SetValue(CUT_STATUS, stat)
                if stat == BRep._CUT:
                    cut_elems.append(qs)
        else:
            for qi, qs in self.forest.iteritems():
                elem = qs.GetElement()
                if elem.Id in selected_elems:
                    stat = self.brep.CutStatusBySampling(elem, nsampling, configuration)
                    elem.SetValue(CUT_STATUS, stat)
                    if stat == BRep._CUT:
                        cut_elems.append(qs)

        for i in range(0, qt_depth):
#            self.aux_util.MultithreadedQuadTreeSubCellRefineBy(cut_elems, self.brep)
            self.aux_util.MultithreadedRefineBy(cut_elems, self.brep)

        # add to the model_part
        cog_points = []
        for qs in cut_elems:
            lastNodeId = self.aux_util.GetLastNodeId(model.model_part)
            lastElementId = self.aux_util.GetLastElementId(model.model_part)
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
        configuration = 1 # we check the cut status in the current configuration
        self.CreateForestSubCell(bulk_elements, nsampling)

        cut_elems = []

        for qi, qs in self.forest.iteritems():
            elem = qs.GetElement()
            stat = self.brep.CutStatusBySampling(elem, nsampling, configuration)
            elem.SetValue(CUT_STATUS, stat)
            if stat == BRep._CUT:
                cut_elems.append(qs)

        for i in range(0, qt_depth):
#            self.aux_util.MultithreadedQuadTreeSubCellRefineBy(cut_elems, self.brep)
            self.aux_util.MultithreadedRefineBy(cut_elems, self.brep)

        Hf = HeavisideFunctionR3R1(self.brep)

        print("checking subcell")
        for qs in cut_elems:
            for i in range(0, qs.NumberOfSubCells()):
                qt = qs.CreateQuadTree(i)
                stat = self.brep.CutStatusBySampling(qt.pCreateGeometry(), nsampling, configuration)
                if stat == BRep._CUT:
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

