import math
import pprint
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.FiniteCellApplication import *
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

class FiniteCellSimulator:

    def __init__(self, params = None):
        self.params = params

    ###TREES & FOREST CREATION#############
    def CreateForest(self, elements, nsampling = 1):
        ###############################################################
        ######### CREATE TREES AND FOREST
        ###############################################################
        self.forest = {}
        for elem in elements:
            if nsampling == 1:
                self.forest[elem.Id] = QuadTreeLocal(elem)
            elif nsampling == 2:
                self.forest[elem.Id] = QuadTreeLocal2(elem)
            elif nsampling == 3:
                self.forest[elem.Id] = QuadTreeLocal3(elem)
            elif nsampling == 4:
                self.forest[elem.Id] = QuadTreeLocal4(elem)
            elif nsampling == 5:
                self.forest[elem.Id] = QuadTreeLocal5(elem)

    def CreateForestSubCell(self, elements, nsampling = 1):
        self.forest = {}
        subcell_fit_mode = self.params["subcell_fit_mode"]
        cut_cell_quadrature_method = self.params["cut_cell_quadrature_method"]
        quad_util = QuadratureUtility()
        if nsampling > 1:
            subcell_class = getattr(FiniteCellApplication, "MomentFittedQuadTreeSubCell" + str(nsampling))
        else:
            subcell_class = getattr(FiniteCellApplication, "MomentFittedQuadTreeSubCell")
        for elem in elements:
            self.forest[elem.Id] = subcell_class(elem)
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
                stat = self.brep.CutStatus(elem)
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
                stat = self.brep.CutStatus(elem)
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
            if nsampling == 1:
                fit_util.MultithreadedFitQuadrature(cut_elems, fit_funcs, self.brep, integrators, cut_cell_quadrature_method, integrator_quadrature_method, solver_type, echo_level, small_weight)
            elif nsampling == 2:
                fit_util.MultithreadedFitQuadrature2(cut_elems, fit_funcs, self.brep, integrators, cut_cell_quadrature_method, integrator_quadrature_method, solver_type, echo_level, small_weight)
            elif nsampling == 3:
                fit_util.MultithreadedFitQuadrature3(cut_elems, fit_funcs, self.brep, integrators, cut_cell_quadrature_method, integrator_quadrature_method, solver_type, echo_level, small_weight)
            elif nsampling == 4:
                fit_util.MultithreadedFitQuadrature4(cut_elems, fit_funcs, self.brep, integrators, cut_cell_quadrature_method, integrator_quadrature_method, solver_type, echo_level, small_weight)
            elif nsampling == 5:
                fit_util.MultithreadedFitQuadrature5(cut_elems, fit_funcs, self.brep, integrators, cut_cell_quadrature_method, integrator_quadrature_method, solver_type, echo_level, small_weight)
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

    def MomentFitSubCell(self, bulk_elements):
        print("fit parameters:")
        pprint.pprint(self.params)
        qt_depth = self.params["qt_depth"]
        max_qt_depth = self.params["max_qt_depth"]
        extra_qt_depth = self.params["extra_qt_depth"]
        small_domain_size = self.params["small_domain_size"]
        number_of_minimum_physical_points = self.params["number_of_minimum_physical_points"]
        nsampling = self.params["number_of_samplings"] if ("number_of_samplings" in self.params) else 1
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

        cut_elems = []
        exclude_elems = []
        quadtree_elems = []

#        fit_mode = self.params["fit_mode"]
        integrator_quadrature_method = self.params["integrator_quadrature_method"]
        solver_type = self.params["fit_solver_type"]
        echo_level = self.params["fit_echo_level"]
        fit_small_weight = self.params["fit_small_weight"]
        action_on_small_subcell = self.params["action_on_small_subcell"]

        if nsampling > 1:
            fit_func_callback = getattr(fit_util, "MultithreadedFitQuadratureSubCell" + str(nsampling))
            save_quadrature_subcell_callback = getattr(fit_util, "SaveQuadratureSubCell" + str(nsampling))
        else:
            fit_func_callback = getattr(fit_util, "MultithreadedFitQuadratureSubCell")
            save_quadrature_subcell_callback = getattr(fit_util, "SaveQuadratureSubCell")

        for qi, qs in self.forest.iteritems():
            elem = qs.GetElement()
            stat = self.brep.CutStatus(elem)
            elem.SetValue(CUT_STATUS, stat)
            if stat == self.brep._CUT:
                cut_elems.append(qs)
            elif stat == self.brep._OUT:
                exclude_elems.append(qs)
        aux_util = FiniteCellAuxilliaryUtility()
        for i in range(0, qt_depth):
            print("Refine level " + str(i+1) + " ...")
#            aux_util.MultithreadedQuadTreeSubCellRefineBy(cut_elems, self.brep)
            aux_util.MultithreadedRefineBy(cut_elems, self.brep)
            print("###################################################")
        # recursive refine the subcell that are cut but not integrable by the quadtree # this is required if Gauss-Legendre quadrature is used
        if max_qt_depth > qt_depth:
            for qs in cut_elems:
                for i in range(0, qs.NumberOfSubCells()):
                    qt = qs.CreateQuadTree(i)
                    stat = self.brep.CutStatus(qt.pCreateGeometry())
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
        # check again, if any quadtree subcell is too small, we shall eliminate them
        print("Refine completed")
        sys.stdout.flush()

        # here we check the small subcell; the small subcell is determined by either domain_size less than a threshold, or number of physical points less than specific number
        if small_domain_size > 0.0:
            criteria_check_small_domain_size = True
        else:
            criteria_check_small_domain_size = False

        if number_of_minimum_physical_points != 0:
            criteria_check_number_of_physical_points = True
        else:
            criteria_check_number_of_physical_points = False

        proper_cut_elems = []
        cnt = 0
        cnt2 = 0
        cnt3 = 0
        num_cut_elems = len(cut_elems)

        if not(criteria_check_small_domain_size or criteria_check_number_of_physical_points):
            proper_cut_elems = cut_elems
            print("\nNo criteria is selected. No quadtree subcell is removed")
        else:
            for qs in cut_elems:
                elem_id = qs.GetElement().Id

                include_element = True

                if include_element:
                    if criteria_check_small_domain_size:
                        domain_size = qs.DomainSize(self.brep, integrator_quadrature_method)
                        if domain_size < small_domain_size:
                            include_element = False
                            print("\nthe quadtree subcell of element " + str(elem_id) + " is too small, domain size = " + str(domain_size) + ". It will be skipped.")

                if include_element:
                    if criteria_check_number_of_physical_points:
                        num_physical_points = qs.GetNumberOfPhysicalIntegrationPoint(self.brep, integrator_quadrature_method)
                        if num_physical_points < number_of_minimum_physical_points:
                            include_element = False
                            print("\nthe quadtree subcell of element " + str(elem_id) + " has number of physical point(s) = " + str(num_physical_points) + " < " + str(number_of_minimum_physical_points) + ". It will be skipped.")

#                print("include_element:", include_element)
                if include_element:
                    proper_cut_elems.append(qs)
                else:
                    if action_on_small_subcell == 'eliminate':
                        exclude_elems.append(qs)
                        cnt = cnt + 1
                    elif action_on_small_subcell == 'replace by quadtree':
                        quadtree_elems.append(qs)
                        cnt = cnt + 1

                if cnt2 > 0.01*num_cut_elems:
                    cnt2 = 0
                    sys.stdout.write(str(cnt3) + " ")
                    sys.stdout.flush()
                    cnt3 = cnt3 + 1
                else:
                    cnt2 = cnt2 + 1

            if action_on_small_subcell == 'eliminate':
                print("\nRemoval completed, " + str(cnt) + " quadtree subcell is removed")
            elif action_on_small_subcell == 'replace by quadtree':
                print("\nReplace completed, " + str(cnt) + " quadtree subcell will be filled with quadtree quadrature")


#        if small_domain_size > 0.0:
#            proper_cut_elems = []
#            cnt = 0
#            cnt2 = 0
#            cnt3 = 0
#            num_cut_elems = len(cut_elems)
#            for qs in cut_elems:
#                elem_id = qs.GetElement().Id
#                domain_size = qs.DomainSize(self.brep, integrator_quadrature_method)
#                if domain_size > small_domain_size:
#                    if number_of_minimum_physical_points != 0:
#                        num_physical_points = qs.GetNumberOfPhysicalIntegrationPoint(self.brep, integrator_quadrature_method)
#                        if num_physical_points >= number_of_minimum_physical_points:
#                            proper_cut_elems.append(qs)
#                        else:
#                            cnt = cnt + 1
#                            print("\nthe quadtree subcell of element " + str(elem_id) + " has number of physical point(s) = " + str(num_physical_points) + " < " + str(number_of_minimum_physical_points) + ". It will be skipped.")
#                else:
#                    cnt = cnt + 1
#                    exclude_elems.append(qs)
#                    print("\nthe quadtree subcell of element " + str(elem_id) + " is too small, domain size = " + str(domain_size) + ". It will be skipped.")
#                if cnt2 > 0.01*num_cut_elems:
#                    cnt2 = 0
#                    sys.stdout.write(str(cnt3) + " ")
#                    sys.stdout.flush()
#                    cnt3 = cnt3 + 1
#                else:
#                    cnt2 = cnt2 + 1
#            print("\nRemoval completed, " + str(cnt) + " quadtree subcell is removed")
#        else:
#            if number_of_minimum_physical_points != 0:
#                proper_cut_elems = []
#                cnt = 0
#                cnt2 = 0
#                cnt3 = 0
#                num_cut_elems = len(cut_elems)
#                for qs in cut_elems:
#                    elem_id = qs.GetElement().Id
#                    num_physical_points = qs.GetNumberOfPhysicalIntegrationPoint(self.brep, integrator_quadrature_method)
#                    if num_physical_points >= number_of_minimum_physical_points:
#                        proper_cut_elems.append(qs)
#                    else:
#                        cnt = cnt + 1
#                        print("\nthe quadtree subcell of element " + str(elem_id) + " has number of physical point(s) = " + str(num_physical_points) + " < " + str(number_of_minimum_physical_points) + ". It will be skipped.")
#                    if cnt2 > 0.01*num_cut_elems:
#                        cnt2 = 0
#                        sys.stdout.write(str(cnt3) + " ")
#                        sys.stdout.flush()
#                        cnt3 = cnt3 + 1
#                    else:
#                        cnt2 = cnt2 + 1
#                print("\nRemoval completed, " + str(cnt) + " quadtree subcell is removed")
#            else:
#                proper_cut_elems = cut_elems

        fit_func_callback(proper_cut_elems, fit_funcs, self.brep, integrator_quadrature_method, solver_type, echo_level, fit_small_weight)
        print("construct quadrature using moment fitting for subcell successfully")
        quad_filename = self.params["quad_filename"]
        quad_filetype = self.params["quad_filetype"]
        accuracy = self.params["quad_accuracy"]

        if action_on_small_subcell == 'eliminate':
            save_quadrature_subcell_callback(quad_filename, quad_filetype, proper_cut_elems, exclude_elems, quadtree_elems, accuracy)
        elif action_on_small_subcell == 'replace by quadtree':
            # perform quadtree refinement for "eliminated" cell
            total_add_quadrature_pnts = 0
            small_cut_cell_quadrature_method = self.params["small_cut_cell_quadrature_method"]
            quadtree_small_weight = self.params["quadtree_small_weight"]
            valid_quadtree_elems = []
            for qs in quadtree_elems:
                np = qs.ConstructQuadrature(self.brep, small_cut_cell_quadrature_method, quadtree_small_weight)
                total_add_quadrature_pnts = total_add_quadrature_pnts + np
                if np == 0:
                    exclude_elems.append(qs)
                else:
                    valid_quadtree_elems.append(qs)
            save_quadrature_subcell_callback(quad_filename, quad_filetype, proper_cut_elems, exclude_elems, valid_quadtree_elems, accuracy)
            print("generate quadtree quadrature for small subcell completed, " + str(total_add_quadrature_pnts) + " quadrature points are added")

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
                stat = self.brep.CutStatus(elem)
                if stat == self.brep._CUT:
                    for i in range(0, qt_depth):
                        qt.RefineBy(self.brep)
                    np = qt.ConstructQuadrature(self.brep, cut_cell_quadrature_method, small_weight)
                    total_add_quadrature_pnts = total_add_quadrature_pnts + np
                    elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                    elem.Initialize()
                    aux_util.AddElement(self.proper_cut_elems, elem)
                    cut_elems.append(elem)
                elif stat == 1:
                    elem.SetValue(ACTIVATION_LEVEL, -1)
                    elem.SetValue(IS_INACTIVE, True)
                    elem.Set(ACTIVE, False)
                    exclude_elems.append(elem)
            print("construct quadrature using quad-tree successfully, " + str(total_add_quadrature_pnts) + " quadrature points are added")
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
                quad_util.CreateConditionFromPoint(model.model_part, cog_points, "DummyConditionPoint3D")

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
                quad_util.CreateConditionFromPoint(model.model_part, cog_points, "DummyConditionPoint3D")

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
            cut_cell_quadrature = quadrature_data.GetCutCellQuadrature()
            cut_cell_full_quadrature = quadrature_data.GetCutCellFullQuadrature()
            subcell_weights = quadrature_data.GetCutCellSubCellWeights()
            subcell_domain_sizes = quadrature_data.GetCutCellSubCellDomainSizes()
            quadtree_quadrature = quadrature_data.GetQuadTreeQuadrature()
            self.proper_cut_elems = ElementsArray()
            lastElementId = aux_util.GetLastElementId(model.model_part)
            lastCondId = aux_util.GetLastConditionId(model.model_part)
            self.all_subcell_elems = []
            if self.params["export_physical_integration_point"]:
                cog_points = []
            for elem in bulk_elements:
                if elem.Id in cut_elems:
#                    print("at cut element " + str(elem.Id))
                    if not cut_cell_quadrature[elem.Id]: # if cut cell quadrature is empty then skip
                        elem.SetValue(ACTIVATION_LEVEL, -1)
                        elem.SetValue(IS_INACTIVE, True)
                        elem.Set(ACTIVE, False)
                        continue

                    aux_util.AddElement(self.proper_cut_elems, elem)

                    quad_util.SetQuadrature(elem, cut_cell_quadrature_order, cut_cell_quadrature[elem.Id])
                    elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                    elem.Initialize()

                    output = self.forest[elem.Id].CreateSubCellElements(model.model_part, subcell_element_type, cut_cell_quadrature_order, cut_cell_full_quadrature[elem.Id], subcell_weights[elem.Id], lastElementId, lastCondId)
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
                elif elem.Id in exclude_elems:
                    elem.SetValue(ACTIVATION_LEVEL, -1)
                    elem.SetValue(IS_INACTIVE, True)
                    elem.Set(ACTIVE, False)
                elif elem.Id in quadtree_elems:
                    quad_util.SetQuadrature(elem, cut_cell_quadrature_order, quadtree_quadrature[elem.Id])
                    elem.SetValue(INTEGRATION_ORDER, cut_cell_quadrature_order)
                    elem.Initialize()
                    aux_util.AddElement(self.proper_cut_elems, elem) # a quadtree element is also considered as proper_cut_elems
                    print("element " + str(elem.Id) + " is assigned " + str(len(quadtree_quadrature[elem.Id])) + " quadrature points")
                    if self.params["export_physical_integration_point"]:
                        points = elem.GetValuesOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model.model_part.ProcessInfo)
                        for point in points:
                            cog_points.append(quad_util.CreatePoint(point[0], point[1], point[2]))
                else:
                    self.mpu.SetMaterialProperties(model.model_part, elem, self.mat_type)
            print("obtain moment-fit subcell quadrature successfully")
            if self.params["export_physical_integration_point"]:
                quad_util.CreateConditionFromPoint(model.model_part, cog_points, "DummyConditionPoint3D")
        else:
            print("Unknown quadrature_method", self.quadrature_method)
            sys.exit(0)

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
            stat = self.brep.CutStatus(elem)
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
                stat = self.brep.CutStatus(qt.pGetGeometry())
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
                stat = self.brep.CutStatus(elem)
                elem.SetValue(CUT_STATUS, stat)
                if stat == self.brep._CUT:
                    cut_elems.append(qs)
        else:
            for qi, qs in self.forest.iteritems():
                elem = qs.GetElement()
                if elem.Id in selected_elems:
                    stat = self.brep.CutStatus(elem)
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
            stat = self.brep.CutStatus(elem)
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
                stat = self.brep.CutStatus(qt.pCreateGeometry())
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

