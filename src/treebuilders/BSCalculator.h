/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "TreeCalculator.h"

namespace mrcpp {

class BSCalculator final : public TreeCalculator<2> {
public:
    explicit BSCalculator(const ScalingBasis &basis, int n);

private:
    const int diff_order;
    Eigen::MatrixXd S_m1;
    Eigen::MatrixXd S_0;
    Eigen::MatrixXd S_p1;

    void calcNode(MWNode<2> &node);
    void readSMatrix(const ScalingBasis &basis, char n);
};

} // namespace mrcpp
