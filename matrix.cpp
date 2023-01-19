#include"matrix.h"

state matrix::eigenvector1 () const {
    if (at(0,1) == complex<double>(0,0))
        return state(1,0);
    return state((eigenvalue2()-at(0,0))/at(0,1), 1, true);
}
state matrix::eigenvector2 () const {
    if (at(0,1) == complex<double>(0,0))
        return state(0,1);
    return state((eigenvalue1()-at(0,0))/at(0,1), 1, true);
}

void matrix::dagger () {
    elem[0][0] = conj(elem[0][0]);
    elem[1][1] = conj(elem[1][1]);
    complex<double> tmp = elem[0][1];
    elem[0][1] = conj(elem[1][0]);
    elem[1][0] = conj(tmp);
}

matrix matrix::operator* (matrix B) const {
    matrix C;
    for (int i=0;i<2;i++) {
        for (int j=0;j<2;j++) {
            C.set_elem(i,j, elem[i][0]*B.at(0,j) + elem[i][1]*B.at(1,j));
        }
    }
    return C;
}