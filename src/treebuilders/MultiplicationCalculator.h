/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "TreeCalculator.h"
#include "trees/FunctionTreeVector.h"

namespace mrcpp {

template <int D, typename T> class MultiplicationCalculator final : public TreeCalculator<D, T> {
public:
    MultiplicationCalculator(const FunctionTreeVector<D, T> &inp, bool conjugate = false)
            : prod_vec(inp)
            , conj(conjugate) {}

private:
    FunctionTreeVector<D, T> prod_vec;
    bool conj;

    void calcNode(MWNode<D, double> &node_o) {
        const NodeIndex<D> &idx = node_o.getNodeIndex();
        double *coefs_o = node_o.getCoefs();
        for (int j = 0; j < node_o.getNCoefs(); j++) { coefs_o[j] = 1.0; }
        for (int i = 0; i < this->prod_vec.size(); i++) {
            double c_i = get_coef(this->prod_vec, i);
            FunctionTree<D, double> &func_i = get_func(this->prod_vec, i);
            // This generates missing nodes
            MWNode<D, double> node_i = func_i.getNode(idx); // Copy node
            node_i.mwTransform(Reconstruction);
            node_i.cvTransform(Forward);
            const double *coefs_i = node_i.getCoefs();
            int n_coefs = node_i.getNCoefs();
            for (int j = 0; j < n_coefs; j++) { coefs_o[j] *= c_i * coefs_i[j]; }
        }
        node_o.cvTransform(Backward);
        node_o.mwTransform(Compression);
        node_o.setHasCoefs();
        node_o.calcNorms();
    }
    void calcNode(MWNode<D, ComplexDouble> &node_o) {
        const NodeIndex<D> &idx = node_o.getNodeIndex();
        ComplexDouble *coefs_o = node_o.getCoefs();
        for (int j = 0; j < node_o.getNCoefs(); j++) { coefs_o[j] = 1.0; }
        for (int i = 0; i < this->prod_vec.size(); i++) {
            ComplexDouble c_i = get_coef(this->prod_vec, i);
            FunctionTree<D, ComplexDouble> &func_i = get_func(this->prod_vec, i);
            // ComplexDoublehis generates missing nodes
            MWNode<D, ComplexDouble> node_i = func_i.getNode(idx); // Copy node
            node_i.mwTransform(Reconstruction);
            node_i.cvTransform(Forward);
            const ComplexDouble *coefs_i = node_i.getCoefs();
            int n_coefs = node_i.getNCoefs();
            if (func_i.conjugate() xor (conj and i == 0)) {
                for (int j = 0; j < n_coefs; j++) { coefs_o[j] *= c_i * std::conj(coefs_i[j]); }
            } else {
                for (int j = 0; j < n_coefs; j++) { coefs_o[j] *= c_i * coefs_i[j]; }
            }
        }
        node_o.cvTransform(Backward);
        node_o.mwTransform(Compression);
        node_o.setHasCoefs();
        node_o.calcNorms();
    }
};

} // namespace mrcpp
