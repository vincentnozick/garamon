#include <cc2ga/Mvec.hpp>


// make a point in cc2ga
template<typename T>
cc2ga::Mvec<T> point(const T &x, const T &y){

    cc2ga::Mvec<T> mv;
    mv[cc2ga::E1] = x;
    mv[cc2ga::E2] = y;
    mv[cc2ga::Ei1] = 0.5 * x*x;
    mv[cc2ga::Ei2] = 0.5 * y*y;
    mv[cc2ga::Ei3] = x*y;
    mv[cc2ga::Ei4] = x*x*x;
    mv[cc2ga::Ei5] = x*x*y;
    mv[cc2ga::Ei6] = x*y*y;
    mv[cc2ga::Ei7] = y*y*y;
    mv[cc2ga::Eo1] = 1.0;
    mv[cc2ga::Eo2] = 1.0;

    return mv;
}
