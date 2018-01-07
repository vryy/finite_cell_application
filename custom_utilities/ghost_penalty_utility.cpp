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
//  Date:            6 Jan 2018
//


#include "custom_utilities/ghost_penalty_utility.h"


namespace Kratos
{

const int GhostPenaltyUtility::msEdgesT3[][2] = { {0, 1}, {1, 2}, {2, 0} };
const int GhostPenaltyUtility::msEdgesT6[][3] = { {0, 1, 3}, {1, 2, 4}, {2, 0, 5} };
const int GhostPenaltyUtility::msEdgesQ4[][2] = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
const int GhostPenaltyUtility::msEdgesQ8[][3] = { {0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7} };
const int GhostPenaltyUtility::msEdgesQ9[][3] = { {0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7} };
const int GhostPenaltyUtility::msFacesT4[][3] = { {0, 2, 1}, {0, 3, 2}, {0, 1, 3}, {2, 3, 1} };
const int GhostPenaltyUtility::msFacesT10[][6] = { {0, 2, 1, 6, 5, 4}, {0, 3, 2, 7, 9, 6}, {0, 1, 3, 4, 8, 7}, {2, 3, 1, 9, 8, 5} };
const int GhostPenaltyUtility::msFacesH8[][4] = { {3, 2, 1, 0}, {0, 1, 5, 4}, {2, 6, 5, 1}, {7, 6, 2, 3}, {7, 3, 0, 4}, {4, 5, 6, 7} };
const int GhostPenaltyUtility::msFacesH20[][8] = { {3, 2, 1, 0, 10, 9, 8, 11}, {0, 1, 5, 4, 8, 13, 16, 12}, {2, 6, 5, 1, 14, 17, 13, 9}, {7, 6, 2, 3, 14, 18, 10, 15}, {7, 3, 0, 4, 15, 11, 12, 19}, {4, 5, 6, 7, 16, 17, 18, 19} };
const int GhostPenaltyUtility::msFacesH27[][9] = { {3, 2, 1, 0, 10, 9, 8, 11, 20}, {0, 1, 5, 4, 8, 13, 16, 12, 21}, {2, 6, 5, 1, 14, 17, 13, 9, 22}, {7, 6, 2, 3, 14, 18, 10, 15, 23}, {7, 3, 0, 4, 15, 11, 12, 19, 24}, {4, 5, 6, 7, 16, 17, 18, 19, 25} };

}  // namespace Kratos.

