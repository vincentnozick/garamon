#include <cc3ga/Mvec.hpp>


// make a point in g2ga
template<typename T>
cc3ga::Mvec<T> point(const T &x, const T &y){

    cc3ga::Mvec<T> mv;
    mv[cc3ga::E1] = x;
    mv[cc3ga::E2] = y;
    mv[cc3ga::Ei1] = 0.5 * x*x;
    mv[cc3ga::Ei2] = 0.5 * y*y;
    mv[cc3ga::Ei3] = x*y;
    mv[cc3ga::Ei4] = x*x*x;
    mv[cc3ga::Ei5] = x*x*y;
    mv[cc3ga::Ei6] = x*y*y;
    mv[cc3ga::Ei7] = y*y*y;
    mv[cc3ga::Eo1] = 1.0;
    mv[cc3ga::Eo2] = 1.0;

    return mv;
}
