import math
import pprint

import finite_cell_simulator
from finite_cell_simulator import *

from KratosMultiphysics import *
from KratosMultiphysics.mpi import *

class FiniteCellMpiSimulator(FiniteCellSimulator, object):

    def __init__(self, params):
        super(FiniteCellMpiSimulator, self).__init__(params)

    ###COMPUTE GLOBAL DISPLACEMENT (L2) ERROR###
    def compute_L2_error(self, elements, process_info, solution, P):
        if mpi.rank == 0:
            print("computing L2 error at load " + str(P))
        nom = 0.0
        denom = 0.0
        for element in elements:
            if element.GetValue(IS_INACTIVE) == False:
                u = element.CalculateOnIntegrationPoints(DISPLACEMENT, process_info)
                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
                W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                for i in range(0, len(u)):
                    ana_u = solution.get_displacement(P, Q[i][0], Q[i][1], Q[i][2])
                    nom = nom + (pow(u[i][0] - ana_u[0], 2) + pow(u[i][1] - ana_u[1], 2) + pow(u[i][2] - ana_u[2], 2)) * W[i][0] * J0[i][0]
                    denom = denom + (pow(ana_u[0], 2) + pow(ana_u[1], 2) + pow(ana_u[2], 2)) * W[i][0] * J0[i][0]
        list_nom = [nom]
        list_denom = [denom]
        if mpi.rank == 0:
            print("computing L2 error: integration completed, now start to gather the results")

        mpi.world.barrier()

        all_list_nom = mpi.allgather_list_double(mpi.world, list_nom)
        all_nom = 0.0
        for nom in all_list_nom:
            all_nom = all_nom + nom
        print("all_list_nom:", all_list_nom)
        print("all_nom:", all_nom)

        all_list_denom = mpi.allgather_list_double(mpi.world, list_denom)
        all_denom = 0.0
        for denom in all_list_denom:
            all_denom = all_denom + denom
        print("all_list_denom:", all_list_denom)
        print("all_denom:", all_denom)

        if all_denom == 0.0:
            if all_nom == 0.0:
                return 0.0
            else:
                return float('nan');
        else:
            return math.sqrt(abs(all_nom / all_denom))

    ###COMPUTE GLOBAL DISPLACEMENT (H1) ERROR###
    def compute_H1_error(self, elements, process_info, solution, P):
        if mpi.rank == 0:
            print("computing H1 error at load " + str(P))
        nom = 0.0
        denom = 0.0
        for element in elements:
            if element.GetValue(IS_INACTIVE) == False:
                o = element.CalculateOnIntegrationPoints(THREED_STRESSES, process_info)
                J0 = element.CalculateOnIntegrationPoints(JACOBIAN_0, process_info)
                Q = element.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, process_info)
                W = element.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, process_info)
                for i in range(0, len(o)):
                    ana_o = solution.get_stress_3d(P, Q[i][0], Q[i][1], Q[i][2])
                    nom = nom + (pow(o[i][0] - ana_o[0], 2) + pow(o[i][1] - ana_o[1], 2) + pow(o[i][2] - ana_o[2], 2) + 2.0*(pow(o[i][3] - ana_o[3], 2) + pow(o[i][4] - ana_o[4], 2) + pow(o[i][5] - ana_o[5], 2))) * W[i][0] * J0[i][0]
                    denom = denom + (pow(ana_o[0], 2) + pow(ana_o[1], 2) + pow(ana_o[2], 2) + 2.0*(pow(ana_o[3], 2) + pow(ana_o[4], 2) + pow(ana_o[5], 2))) * W[i][0] * J0[i][0]
        list_nom = [nom]
        list_denom = [denom]
        if mpi.rank == 0:
            print("Computing H1 error: integration completed, now gather the results")

        mpi.world.barrier()

        all_list_nom = mpi.allgather_list_double(mpi.world, list_nom)
        all_nom = 0.0
        for nom in all_list_nom:
            all_nom = all_nom + nom

        all_list_denom = mpi.allgather_list_double(mpi.world, list_denom)
        all_denom = 0.0
        for denom in all_list_denom:
            all_denom = all_denom + denom

        if all_denom == 0.0:
            if all_nom == 0.0:
                return 0.0
            else:
                return float('nan');
        else:
            return math.sqrt(abs(all_nom / all_denom))

