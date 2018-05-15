// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//


// Project includes
#include "includes/define.h"
#include "custom_utilities/quad_tree_node.h"


namespace Kratos
{

template class QuadTreeNode<GLOBAL_REFERENCE>;
template class QuadTreeNode<GLOBAL_CURRENT>;
template class QuadTreeNodeQ4<GLOBAL_REFERENCE>;
template class QuadTreeNodeQ4<GLOBAL_CURRENT>;
#ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
template class QuadTreeNodeBezier2D<GLOBAL_REFERENCE>;
template class QuadTreeNodeBezier2D<GLOBAL_CURRENT>;
template class QuadTreeNodeBezier3D<GLOBAL_REFERENCE>;
template class QuadTreeNodeBezier3D<GLOBAL_CURRENT>;
#endif
template class QuadTreeNodeH8<GLOBAL_REFERENCE>;
template class QuadTreeNodeH8<GLOBAL_CURRENT>;
template class QuadTreeNodeT3<GLOBAL_REFERENCE>;
template class QuadTreeNodeT3<GLOBAL_CURRENT>;
template class QuadTreeNodeT4<GLOBAL_REFERENCE>;
template class QuadTreeNodeT4<GLOBAL_CURRENT>;

}

