#include <iostream>
#include <cstdlib>
#include <c3ga/Mvec.hpp>


int constructorDefaultTest(){
    c3ga::Mvec<double> mv;
    if(!(mv.isEmpty())) return EXIT_FAILURE;
    return EXIT_SUCCESS;
}

int constructorCopyTest(){

    // empty vector
    {
        c3ga::Mvec<double> mv;
        c3ga::Mvec<double> mv2(mv);
        if(!(mv2.isEmpty())) return EXIT_FAILURE;
    }

    // any multivector
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 1;
        mv[c3ga::E12i] = 3;
        mv[c3ga::E0123] = 4;
        c3ga::Mvec<double> mv2(mv);
        if(mv2[c3ga::E1] != 1)    return EXIT_FAILURE;
        if(mv2[c3ga::E12i] != 3)  return EXIT_FAILURE;
        if(mv2[c3ga::E0123] != 4) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int constructorTemplateConverterTest(){

    c3ga::Mvec<double> mv;
    mv[c3ga::E1] = 1;
    mv[c3ga::E12i] = 3;
    mv[c3ga::E0123] = 4;
    c3ga::Mvec<double> mv2(mv);
    if(mv2[c3ga::E1] != 1)    return EXIT_FAILURE;
    if(mv2[c3ga::E12i] != 3)  return EXIT_FAILURE;
    if(mv2[c3ga::E0123] != 4) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int constructorScalarTest(){
    c3ga::Mvec<double> mv(5);

    if(mv[c3ga::scalar] != 5) return EXIT_FAILURE;
    if(mv.grade() != 0)      return EXIT_FAILURE;

    double a = 5;
    c3ga::Mvec<double> mv2(a);
    if(mv2[c3ga::scalar] != 5) return EXIT_FAILURE;
    if(mv2.grade() != 0)      return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int castTest(){

    // empty vector
    {
        c3ga::Mvec<double> mv;
        if( int(mv) != 0) return EXIT_FAILURE;
    }

    // scalar
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::scalar] = 42;
        if( int(mv) != 42) return EXIT_FAILURE;
    }

    // vector
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 42;
        if( int(mv) != 0) return EXIT_FAILURE;
    }

    // tri-vector
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E12i] = 42;
        if( double(mv) != 0.0f) return EXIT_FAILURE;
    }

    // multivrctor with scalar
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::scalar] = 42;
        mv[c3ga::E12i] = 42;
        if( int(mv) != 42) return EXIT_FAILURE;
    }

    // multivector without scalar
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E12]  = 42;
        mv[c3ga::E12i] = 42;
        if( double(mv) != 0.0) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int operatorEqualTest(){

    // operator =
    c3ga::Mvec<double> a;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] = -6.0;
    c3ga::Mvec<double> b = a;
    if(b[c3ga::E12] != 5.0)  return EXIT_FAILURE;
    if(b[c3ga::E13] != -6.0) return EXIT_FAILURE;

    // operator = (const Mvec)
    const c3ga::Mvec<double> c(a);
    c3ga::Mvec<double> d = c;
    if(d[c3ga::E12] != 5.0)  return EXIT_FAILURE;
    if(d[c3ga::E13] != -6.0) return EXIT_FAILURE;

    // operator = (non ref Mvec)
    c3ga::Mvec<double> e = c3ga::Mvec<double>(c);
    if(e[c3ga::E12] != 5.0)  return EXIT_FAILURE;
    if(e[c3ga::E13] != -6.0) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int operatorEqualEqualTest(){

    c3ga::Mvec<double> a,b;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] = -6.0;
    b[c3ga::E12] =  a[c3ga::E12];
    b[c3ga::E13] =  a[c3ga::E13];
    if(!(a==b)) return EXIT_FAILURE;

    b[c3ga::E23] = 2.0;
    if(a==b) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int operatorNotEqualTest(){

    c3ga::Mvec<double> a,b;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] = -6.0;
    b[c3ga::E12] =  a[c3ga::E12];
    b[c3ga::E13] =  a[c3ga::E13];
    if(a!=b) return EXIT_FAILURE;

    b[c3ga::E23] = 2.0;
    if(!(a!=b)) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int operatorPlusEqualTest(){
    c3ga::Mvec<double> a;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] = -6.0;
    c3ga::Mvec<double> b;
    b[c3ga::E13] = 1.0;
    b[c3ga::E1] =  2.0;
    b[c3ga::E3] =  1.0;
    a += b;
    if(a[c3ga::E12] !=  5.0) return EXIT_FAILURE;
    if(a[c3ga::E13] != -5.0) return EXIT_FAILURE;
    if(a[c3ga::E1]  !=  2.0) return EXIT_FAILURE;
    if(a[c3ga::E3]  !=  1.0) return EXIT_FAILURE;
    std::vector<unsigned int> g = a.grades();
    if(g[0] != 1) return EXIT_FAILURE;
    if(g[1] != 2) return EXIT_FAILURE;
    if(a.grade() != 2) return EXIT_FAILURE;


    a += 8.0; // from constructor with scalar

    if(a[c3ga::scalar] !=  8.0) return EXIT_FAILURE;
    g = a.grades();
    if(g[0] != 0) return EXIT_FAILURE;
    if(g[1] != 1) return EXIT_FAILURE;
    if(g[2] != 2) return EXIT_FAILURE;
    if(a.grade() != 2) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}


int operatorPlusTest(){

    c3ga::Mvec<double> a;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] = -6.0;

    c3ga::Mvec<double> b;
    b[c3ga::E13] = 1.0;
    b[c3ga::E1] =  2.0;
    b[c3ga::E3] =  1.0;

    c3ga::Mvec<double> c;
    c = a + b;

    if(c[c3ga::E12] !=  5.0) return EXIT_FAILURE;
    if(c[c3ga::E13] != -5.0) return EXIT_FAILURE;
    if(c[c3ga::E1]  !=  2.0) return EXIT_FAILURE;
    if(c[c3ga::E3]  !=  1.0) return EXIT_FAILURE;

    std::vector<unsigned int> g = c.grades();
    if(g[0] != 1) return EXIT_FAILURE;
    if(g[1] != 2) return EXIT_FAILURE;
    if(c.grade() != 2) return EXIT_FAILURE;

    c3ga::Mvec<double> d;
    d = a + 42; // int
    if(d[c3ga::E12] !=  5.0)    return EXIT_FAILURE;
    if(d[c3ga::E13] != -6.0)    return EXIT_FAILURE;
    if(d[c3ga::scalar] != 42.0) return EXIT_FAILURE;
    g = d.grades();
    if(g[0] != 0) return EXIT_FAILURE;
    if(g[1] != 2) return EXIT_FAILURE;
    if(d.grade() != 2) return EXIT_FAILURE;

    c3ga::Mvec<double> e;
    e = a + 42.0f; // double
    if(e[c3ga::E12] !=  5.0)    return EXIT_FAILURE;
    if(e[c3ga::E13] != -6.0)    return EXIT_FAILURE;
    if(e[c3ga::scalar] != 42.0) return EXIT_FAILURE;
    g = e.grades();
    if(g[0] != 0) return EXIT_FAILURE;
    if(g[1] != 2) return EXIT_FAILURE;
    if(e.grade() != 2) return EXIT_FAILURE;

    c3ga::Mvec<double> f;
    f = 42 + a; // int
    if(f[c3ga::E12] !=  5.0)    return EXIT_FAILURE;
    if(f[c3ga::E13] != -6.0)    return EXIT_FAILURE;
    if(f[c3ga::scalar] != 42.0) return EXIT_FAILURE;
    g = f.grades();
    if(g[0] != 0) return EXIT_FAILURE;
    if(g[1] != 2) return EXIT_FAILURE;
    if(f.grade() != 2) return EXIT_FAILURE;

    c3ga::Mvec<double> h;
    h = 42.0f + a;  // double
    if(h[c3ga::E12] !=  5.0)    return EXIT_FAILURE;
    if(h[c3ga::E13] != -6.0)    return EXIT_FAILURE;
    if(h[c3ga::scalar] != 42.0) return EXIT_FAILURE;
    g = h.grades();
    if(g[0] != 0) return EXIT_FAILURE;
    if(g[1] != 2) return EXIT_FAILURE;
    if(h.grade() != 2) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int operatorMinusTest(){

    c3ga::Mvec<double> a;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] = -6.0;

    c3ga::Mvec<double> b;
    b[c3ga::E13] = 1.0;
    b[c3ga::E1] =  2.0;
    b[c3ga::E3] =  1.0;

    c3ga::Mvec<double> c;
    c = a - b;

    if(c[c3ga::E12] !=  5.0) return EXIT_FAILURE;
    if(c[c3ga::E13] != -7.0) return EXIT_FAILURE;
    if(c[c3ga::E1]  != -2.0) return EXIT_FAILURE;
    if(c[c3ga::E3]  != -1.0) return EXIT_FAILURE;

    std::vector<unsigned int> g = c.grades();
    if(g[0] != 1) return EXIT_FAILURE;
    if(g[1] != 2) return EXIT_FAILURE;
    if(c.grade() != 2) return EXIT_FAILURE;

    c3ga::Mvec<double> d;
    d = a - 42; // int
    if(d[c3ga::E12] !=  5.0)    return EXIT_FAILURE;
    if(d[c3ga::E13] != -6.0)    return EXIT_FAILURE;
    if(d[c3ga::scalar] != -42.0) return EXIT_FAILURE;
    g = d.grades();
    if(g[0] != 0) return EXIT_FAILURE;
    if(g[1] != 2) return EXIT_FAILURE;
    if(d.grade() != 2) return EXIT_FAILURE;

    c3ga::Mvec<double> e;
    e = a - 42.0f; // double
    if(e[c3ga::E12] !=  5.0)    return EXIT_FAILURE;
    if(e[c3ga::E13] != -6.0)    return EXIT_FAILURE;
    if(e[c3ga::scalar] != -42.0) return EXIT_FAILURE;
    g = e.grades();
    if(g[0] != 0) return EXIT_FAILURE;
    if(g[1] != 2) return EXIT_FAILURE;
    if(e.grade() != 2) return EXIT_FAILURE;

    c3ga::Mvec<double> f;
    f = 42 - a; // int
    if(f[c3ga::E12] != -5.0)    return EXIT_FAILURE;
    if(f[c3ga::E13] !=  6.0)    return EXIT_FAILURE;
    if(f[c3ga::scalar] != 42.0) return EXIT_FAILURE;
    g = f.grades();
    if(g[0] != 0) return EXIT_FAILURE;
    if(g[1] != 2) return EXIT_FAILURE;
    if(f.grade() != 2) return EXIT_FAILURE;

    c3ga::Mvec<double> h;
    h = 42.0f - a;  // double
    if(h[c3ga::E12] != -5.0)    return EXIT_FAILURE;
    if(h[c3ga::E13] !=  6.0)    return EXIT_FAILURE;
    if(h[c3ga::scalar] != 42.0) return EXIT_FAILURE;
    g = h.grades();
    if(g[0] != 0) return EXIT_FAILURE;
    if(g[1] != 2) return EXIT_FAILURE;
    if(h.grade() != 2) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int operatorMinusEqualTest(){
    c3ga::Mvec<double> a;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] = -6.0;
    c3ga::Mvec<double> b;
    b[c3ga::E1]  =  2.0;
    b[c3ga::E3]  =  1.0;
    b[c3ga::E13] =  1.0;

    a -= b;
    if(a[c3ga::E12] !=  5.0) return EXIT_FAILURE;
    if(a[c3ga::E13] != -7.0) return EXIT_FAILURE;
    if(a[c3ga::E1]  != -2.0) return EXIT_FAILURE;
    if(a[c3ga::E3]  != -1.0) return EXIT_FAILURE;
    std::vector<unsigned int> g = a.grades();
    if(g[0] != 1) return EXIT_FAILURE;
    if(g[1] != 2) return EXIT_FAILURE;
    if(a.grade() != 2) return EXIT_FAILURE;

    a-= 8.0;
    if(a[c3ga::scalar] != -8.0) return EXIT_FAILURE;
    g = a.grades();
    if(g[0] != 0) return EXIT_FAILURE;
    if(g[1] != 1) return EXIT_FAILURE;
    if(g[2] != 2) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int operatorUnaryMinusTest(){

    c3ga::Mvec<double> a;
    a[c3ga::E1] =  2.0;
    a[c3ga::E3] =  1.0;
    a[c3ga::E12] =  5.0;
    a[c3ga::E03] = -6.0;
    a[c3ga::E23i] = -7.0;
    c3ga::Mvec<double> b;
    b = - a;

    if(b[c3ga::E1]   != -2.0) return EXIT_FAILURE;
    if(b[c3ga::E3]   != -1.0) return EXIT_FAILURE;
    if(b[c3ga::E12]  != -5.0) return EXIT_FAILURE;
    if(b[c3ga::E03]  !=  6.0) return EXIT_FAILURE;
    if(b[c3ga::E23i] !=  7.0) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int gradeTest(){

    c3ga::Mvec<double> a;
    if(a.grade() != 0) return EXIT_FAILURE;

    a[c3ga::scalar] =  5.0;
    if(a.grade() != 0) return EXIT_FAILURE;

    a[c3ga::E1] =  5.0;
    a[c3ga::E3] = -6.0;
    if(a.grade() != 1) return EXIT_FAILURE;

    a[c3ga::E12] =  5.0;
    a[c3ga::E03] = -6.0;
    if(a.grade() != 2) return EXIT_FAILURE;

    a[c3ga::E12i] =  5.0;
    a[c3ga::E013] = -6.0;
    if(a.grade() != 3) return EXIT_FAILURE;

    a[c3ga::E012i] =  5.0;
    a[c3ga::E0123] = -6.0;
    if(a.grade() != 4) return EXIT_FAILURE;

    a[c3ga::E0123i] = 5.0; // ne fonctionne pas
    if(a.grade() != 5) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int isGradeTest(){

    c3ga::Mvec<double> a;
    if(a.isGrade(0) == false) return EXIT_FAILURE;
    if(a.isGrade(1) == true)  return EXIT_FAILURE;

    a[c3ga::scalar] =  5.0;
    if(a.isGrade(0) == false) return EXIT_FAILURE;
    if(a.isGrade(1) == true)  return EXIT_FAILURE;

    a[c3ga::E1] =  5.0;
    a[c3ga::E3] = -6.0;
    if(a.isGrade(1) == false) return EXIT_FAILURE;
    if(a.isGrade(2) == true)  return EXIT_FAILURE;

    a[c3ga::E12] =  5.0;
    a[c3ga::E03] = -6.0;
    if(a.isGrade(2) == false) return EXIT_FAILURE;
    if(a.isGrade(3) == true)  return EXIT_FAILURE;

    a[c3ga::E12i] =  5.0;
    a[c3ga::E013] = -6.0;
    if(a.isGrade(3) == false) return EXIT_FAILURE;
    if(a.isGrade(4) == true)  return EXIT_FAILURE;

    a[c3ga::E012i] =  5.0;
    a[c3ga::E0123] = -6.0;
    if(a.isGrade(4) == false) return EXIT_FAILURE;
    if(a.isGrade(5) == true)  return EXIT_FAILURE;

    a[c3ga::E0123i] = 5.0;
    if(a.isGrade(5) == false) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int accessOperatorTest(){

    c3ga::Mvec<double> a;

    a[c3ga::scalar] =  5.0;

    a[c3ga::E0] =  5.0;
    a[c3ga::E1] =  5.0;
    a[c3ga::E2] =  5.0;
    a[c3ga::E3] =  5.0;
    a[c3ga::Ei] =  5.0;

    a[c3ga::E01] =  5.0;
    a[c3ga::E02] =  5.0;
    a[c3ga::E03] =  5.0;
    a[c3ga::E0i] =  5.0;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] =  5.0;
    a[c3ga::E1i] =  5.0;
    a[c3ga::E23] =  5.0;
    a[c3ga::E2i] =  5.0;
    a[c3ga::E3i] =  5.0;

    a[c3ga::E012] =  5.0;
    a[c3ga::E013] =  5.0;
    a[c3ga::E01i] =  5.0;
    a[c3ga::E023] =  5.0;
    a[c3ga::E02i] =  5.0;
    a[c3ga::E03i] =  5.0;
    a[c3ga::E123] =  5.0;
    a[c3ga::E12i] =  5.0;
    a[c3ga::E13i] =  5.0;
    a[c3ga::E23i] =  5.0;

    a[c3ga::E0123] =  5.0;
    a[c3ga::E012i] =  5.0;
    a[c3ga::E013i] =  5.0;
    a[c3ga::E023i] =  5.0;
    a[c3ga::E123i] =  5.0;

    a[c3ga::E0123i] =  5.0;

    return EXIT_SUCCESS;
}

int staticGetTest(){

    // bivector
    {
        c3ga::Mvec<double> mv = c3ga::e12<double>();
        if(mv.grade()   != 2) return EXIT_FAILURE;
        if(mv[c3ga::E12] != 1) return EXIT_FAILURE;
    }

    // trivector
    {
        c3ga::Mvec<double> mv = 2 * c3ga::e12i<double>();
        if(mv.grade()    != 3) return EXIT_FAILURE;
        if(mv[c3ga::E12i] != 2) return EXIT_FAILURE;
    }

    // quadvector
    {
        c3ga::Mvec<double> mv =  c3ga::e1<double>() ^ c3ga::e2<double>() ^ c3ga::e3i<double>();
        if(mv.grade()     != 4) return EXIT_FAILURE;
        if(mv[c3ga::E123i] != 1) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int methodGetTest(){

    // bivector
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E12] = 2;
        a[c3ga::E23] = 3;
        b = a.e12();
        if(b.grade()   != 2) return EXIT_FAILURE;
        if(b[c3ga::E12] != 2) return EXIT_FAILURE;
        if(b[c3ga::E23] != 0) return EXIT_FAILURE;
    }

    // empty
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E12i] = 2;
        a[c3ga::E23i] = 3;
        b = a.e12();
        if(b.grade() != 0) return EXIT_FAILURE;
    }

    // trivector
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E012] = 2;
        a[c3ga::E023] = 3;
        b = 2 * a.e012() ^ c3ga::ei<double>();
        if(b.grade()     != 4) return EXIT_FAILURE;
        if(b[c3ga::E012i] != 4) return EXIT_FAILURE;
        if(b[c3ga::E023i] != 0) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


int gradesTest(){

    c3ga::Mvec<double> a;

    // empty multivector
    std::vector<unsigned int> g = a.grades();
    if(g.size() != 0) return EXIT_FAILURE;

    // grade 0
    a[c3ga::scalar] =  5.0;
    g = a.grades();
    if(g.size() != 1) return EXIT_FAILURE;
    if(g[0] != 0)     return EXIT_FAILURE;

    // grade 1
    a[c3ga::E2] =  5.0;
    g = a.grades();
    if(g.size() != 2) return EXIT_FAILURE;
    if(g[0] != 0)     return EXIT_FAILURE;
    if(g[1] != 1)     return EXIT_FAILURE;

    // grade 2
    a[c3ga::E13] =  5.0;
    g = a.grades();
    if(g.size() != 3) return EXIT_FAILURE;
    if(g[0] != 0)     return EXIT_FAILURE;
    if(g[1] != 1)     return EXIT_FAILURE;
    if(g[2] != 2)     return EXIT_FAILURE;

    // grade 3
    a[c3ga::E13i] =  5.0;
    g = a.grades();
    if(g.size() != 4) return EXIT_FAILURE;
    if(g[0] != 0)     return EXIT_FAILURE;
    if(g[1] != 1)     return EXIT_FAILURE;
    if(g[2] != 2)     return EXIT_FAILURE;
    if(g[3] != 3)     return EXIT_FAILURE;

    // grade 4
    a[c3ga::E123i] =  5.0;
    g = a.grades();
    if(g.size() != 5) return EXIT_FAILURE;
    if(g[0] != 0)     return EXIT_FAILURE;
    if(g[1] != 1)     return EXIT_FAILURE;
    if(g[2] != 2)     return EXIT_FAILURE;
    if(g[3] != 3)     return EXIT_FAILURE;
    if(g[4] != 4)     return EXIT_FAILURE;

    // grade 4
    a[c3ga::E0123i] =  5.0;
    g = a.grades();
    if(g.size() != 6) return EXIT_FAILURE;
    if(g[0] != 0)     return EXIT_FAILURE;
    if(g[1] != 1)     return EXIT_FAILURE;
    if(g[2] != 2)     return EXIT_FAILURE;
    if(g[3] != 3)     return EXIT_FAILURE;
    if(g[4] != 4)     return EXIT_FAILURE;
    if(g[5] != 5)     return EXIT_FAILURE;

    // non consecutives grades
    c3ga::Mvec<double> b;
    b[c3ga::E1]     =  5.0;
    b[c3ga::E123]   =  5.0;
    b[c3ga::E0123i] =  5.0;
    g = b.grades();
    if(g.size() != 3) return EXIT_FAILURE;
    if(g[0] != 1)     return EXIT_FAILURE;
    if(g[1] != 3)     return EXIT_FAILURE;
    if(g[2] != 5)     return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int sameGradesTest(){

    // empty vectors
    {
        c3ga::Mvec<double> a,b;
        std::vector<unsigned int> g = a.grades();
        if(g.size() != 0) return EXIT_FAILURE;
    }

    // scalars
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::scalar] = 42;
        b[c3ga::scalar] = 2;
        if(a.sameGrade(b) == false) return EXIT_FAILURE;
    }

    // scalars
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::scalar] = 42;
        b[c3ga::E1] = 2;
        if(a.sameGrade(b) == true) return EXIT_FAILURE;
    }

    // vectors
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E2] = 42;
        b[c3ga::E1] = 2;
        if(a.sameGrade(b) == false) return EXIT_FAILURE;
    }

    // vectors
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E2] = 42;
        b[c3ga::E12] = 2;
        if(a.sameGrade(b) == true) return EXIT_FAILURE;
    }

    // any k-vector
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E2] = 42;
        a[c3ga::E23] = 42;
        a[c3ga::E23i] = 42;
        b[c3ga::E123] = 2;
        if(a.sameGrade(b) == false) return EXIT_FAILURE;
    }

    // any k-vector
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E2] = 42;
        a[c3ga::E23] = 42;
        a[c3ga::E23i] = 42;
        b[c3ga::E123i] = 2;
        if(a.sameGrade(b) == true) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int isHomogeneousTest()
{
    // grade 0
    {
        c3ga::Mvec<double> a;
        if(a.isHomogeneous() == false) return EXIT_FAILURE;
    }

    // grade 0
    {
        c3ga::Mvec<double> a;
        a[c3ga::scalar] = 42;
        if(a.isHomogeneous() == false) return EXIT_FAILURE;
    }

    // grade 1
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 42;
        if(a.isHomogeneous() == false) return EXIT_FAILURE;
    }

    // grade 2
    {
        c3ga::Mvec<double> a;
        a[c3ga::E13] = 42;
        if(a.isHomogeneous() == false) return EXIT_FAILURE;
    }

    // grade 3
    {
        c3ga::Mvec<double> a;
        a[c3ga::E13i] = 42;
        if(a.isHomogeneous() == false) return EXIT_FAILURE;
    }

    // grade 4
    {
        c3ga::Mvec<double> a;
        a[c3ga::E013i] = 42;
        if(a.isHomogeneous() == false) return EXIT_FAILURE;
    }

    // grade 5
    {
        c3ga::Mvec<double> a;
        a[c3ga::E0123i] = 42;
        if(a.isHomogeneous() == false) return EXIT_FAILURE;
    }

    // non homogeneous 1
    {
        c3ga::Mvec<double> a;
        a[c3ga::E012] = 42;
        a[c3ga::E1]   = 42;
        if(a.isHomogeneous() == true) return EXIT_FAILURE;
    }

    // non homogeneous 2
    {
        c3ga::Mvec<double> a;
        a[c3ga::E012] = 42;
        a[c3ga::scalar] = 42;
        a[c3ga::E12]  = 42;
        if(a.isHomogeneous() == true) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int roundZeroTest(){

    c3ga::Mvec<double> mv;
    mv[c3ga::E12]  = 1.0e-10;
    mv[c3ga::E123] =  5.0;
    mv[c3ga::E12i] = -5.0;
    mv[c3ga::E13i] = -1.0e-10;

    mv.roundZero(1.0e-6);

    std::vector<unsigned int> g = mv.grades();
    if(g.size() != 1){return EXIT_FAILURE;}
    if(g[0] != 3)          return EXIT_FAILURE;


    if(mv[c3ga::E12]  != 0.0f) return EXIT_FAILURE; // when doing this, create a zero array for bivectors
    if(mv[c3ga::E123] == 0.0f) return EXIT_FAILURE;
    if(mv[c3ga::E12i] == 0.0f) return EXIT_FAILURE;
    if(mv[c3ga::E13i] != 0.0f) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int isEmptyTest(){
    c3ga::Mvec<double> mv;
    if(!mv.isEmpty()) return EXIT_FAILURE;

    mv[c3ga::E3] = 42.0;
    if(mv.isEmpty()) return EXIT_FAILURE;

    mv[c3ga::E123] = 42.0;
    if(mv.isEmpty()) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int clearTest(){

    c3ga::Mvec<double> mv;
    mv[c3ga::E3]   = 42.0;
    mv[c3ga::E123] = 42.0;


    // erase specific grade
    mv.clear(3);
    std::vector<unsigned  int> g = mv.grades();
    if(g.size() != 1) return EXIT_FAILURE;
    if(g[0] != 1)     return EXIT_FAILURE;

    // erase non existing grade
    mv.clear(2);
    g = mv.grades();
    if(g.size() != 1) return EXIT_FAILURE;
    if(g[0] != 1)     return EXIT_FAILURE;

    // erase all
    mv[c3ga::E123]  = 42.0;
    mv[c3ga::E0123] = 42.0;
    mv.clear();
    g = mv.grades();
    if(g.size() != 0) return EXIT_FAILURE;
    if(!mv.isEmpty()) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int wedgeTest(){

    double epsilon = 1.0e-7;

     // test1
     {
         c3ga::Mvec<double> a;
         a[c3ga::E1] = 1.0;
         a[c3ga::E2] = 1.0;
         c3ga::Mvec<double> b;
         b[c3ga::E1] = 1.0;
         b[c3ga::E3] = 1.0;
         c3ga::Mvec<double> c = a ^b;
         if (fabs(c[c3ga::E13] - 1.0) > epsilon) return EXIT_FAILURE;
         if (fabs(c[c3ga::E12] + 1.0) > epsilon) return EXIT_FAILURE;
         if (fabs(c[c3ga::E23] - 1.0) > epsilon) return EXIT_FAILURE;
     }

     // test2
     {
         c3ga::Mvec<double> a;
         a[c3ga::E1] = 1.0;
         a[c3ga::E2] = 1.0;
         c3ga::Mvec<double> b;
         b[c3ga::E12] = 1.0;
         b[c3ga::E23] = 1.0;
         c3ga::Mvec<double> c = a ^b;
         if (fabs(c[c3ga::E123] - 1.0) > epsilon) return EXIT_FAILURE;
     }

    // test3 (homogeneous multivectors)
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        c3ga::Mvec<double> b;
        b[c3ga::E12] = 4.0;
        b[c3ga::E23] = 5.0;
        c3ga::Mvec<double> c = a^b;
        if(fabs(c[c3ga::E123] - 10.0) > epsilon) return EXIT_FAILURE;
    }

    // test4
    {
        c3ga::Mvec<double> a;
        a[c3ga::scalar] = 2.0;
        a[c3ga::E3] = 3.0;
        c3ga::Mvec<double> b;
        b[c3ga::E1]  = 1.0;
        b[c3ga::E23] = 1.0;
        c3ga::Mvec<double> c = a^b;
        if(fabs(c[c3ga::E1]  - 2.0) > epsilon) return EXIT_FAILURE;
        if(fabs(c[c3ga::E23] - 2.0) > epsilon) return EXIT_FAILURE;
        if(fabs(c[c3ga::E13] + 3.0) > epsilon) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int wedgeEqualTest(){

    double epsilon = 1.0e-7;

    // test1
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 1.0;
        a[c3ga::E2] = 1.0;
        c3ga::Mvec<double> b;
        b[c3ga::E1] = 1.0;
        b[c3ga::E3] = 1.0;
        a ^= b;
        if (fabs(a[c3ga::E13] - 1.0) > epsilon) return EXIT_FAILURE;
        if (fabs(a[c3ga::E12] + 1.0) > epsilon) return EXIT_FAILURE;
        if (fabs(a[c3ga::E23] - 1.0) > epsilon) return EXIT_FAILURE;
    }

    // test2
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 1.0;
        a[c3ga::E2] = 1.0;
        c3ga::Mvec<double> b;
        b[c3ga::E12] = 1.0;
        b[c3ga::E23] = 1.0;
        a ^= b;
        if (fabs(a[c3ga::E123] - 1.0) > epsilon) return EXIT_FAILURE;
    }

    // test3
    {
        c3ga::Mvec<double> a;
        a[c3ga::scalar] = 2.0;
        a[c3ga::E3] = 3.0;
        c3ga::Mvec<double> b;
        b[c3ga::E1]  = 1.0;
        b[c3ga::E23] = 1.0;
        a ^= b;
        if(fabs(a[c3ga::E1]  - 2.0) > epsilon) return EXIT_FAILURE;
        if(fabs(a[c3ga::E23] - 2.0) > epsilon) return EXIT_FAILURE;
        if(fabs(a[c3ga::E13] + 3.0) > epsilon) return EXIT_FAILURE;
    }

    // test4
    {
        c3ga::Mvec<double> a;
        a[c3ga::scalar] = 2.0;
        a[c3ga::E3] = 3.0;
        a ^= 42.0f; // double
        if(fabs(a[c3ga::scalar] - 84.0) > epsilon) return EXIT_FAILURE;
        if(fabs(a[c3ga::E3] - 126.0) > epsilon) return EXIT_FAILURE;
    }

    // test5
    {
        c3ga::Mvec<double> a;
        a[c3ga::scalar] = 2.0;
        a[c3ga::E3] = 3.0;
        a ^= 42; // int
        if(fabs(a[c3ga::scalar] - 84.0) > epsilon) return EXIT_FAILURE;
        if(fabs(a[c3ga::E3] - 126.0) > epsilon) return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


/// test the outer product between the primal form of a multivector and the dual of another
int wedgePrimalDualTest(){
    double epsilon = 1.0e-7;

    // test1
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 4.0;
        a[c3ga::E2] = 3.0;
        c3ga::Mvec<double> b; // we consider this as b and will compute a ^ dual(b)
        b[c3ga::E1] = 2.0;
        b[c3ga::E3] = 5.0;
        b[c3ga::E0] = 6.0;
        c3ga::Mvec<double> c = a.outerPrimalDual(b);
        if (fabs(c[c3ga::E0123i] - 8.0) > epsilon) {std::cout << "wedge outer primal dual: test 1"  << " failed "<< "mv3 = "<< c <<std::endl; return EXIT_FAILURE;} // first, all verifications are now wrong (dual instead of primal) I miss one verification here, I have to put it on the todo list
    }

    // test2
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 4.0;
        a[c3ga::E2] = 3.0;
        c3ga::Mvec<double> b;
        b[c3ga::E12] = 2.0;
        b[c3ga::E23] = 5.0;
        b[c3ga::E2i] = 6.0;
        c3ga::Mvec<double> c = a.outerPrimalDual(b);
        if (fabs(c[c3ga::E0123] - 18.0) > epsilon) {std::cout << "wedge outer primal dual: test 2"  << " failed "<< "mv3 = "<< c <<std::endl; return EXIT_FAILURE;}
        if (fabs(c[c3ga::E012i] + 15.0) > epsilon) {std::cout << "wedge outer primal dual: test 2"  << " failed "<< "mv3 = "<< c <<std::endl; return EXIT_FAILURE;}
        if (fabs(c[c3ga::E013i] + 8.0) > epsilon) {std::cout << "wedge outer primal dual: test 2"  << " failed "<< "mv3 = "<< c <<std::endl; return EXIT_FAILURE;}
        if (fabs(c[c3ga::E023i] + 6.0) > epsilon) {std::cout << "wedge outer primal dual: test 2"  << " failed "<< "mv3 = "<< c <<std::endl; return EXIT_FAILURE;}
    }


    return EXIT_SUCCESS;
}



int innerProductTest(){

    const double epsilon = 1.0e-7;
    int numTest = 0;
    // grade 1
    {
        numTest++;
        c3ga::Mvec<double> mv1;
        mv1[c3ga::E1] = 1.2;
        mv1[c3ga::E2] = 2.2;
        c3ga::Mvec<double> mv2;
        mv2[c3ga::E1] = 1;
        mv2[c3ga::E3] = 2;
        c3ga::Mvec<double> mv3 = mv1 | mv2;
        //std::cout << std::endl;
        //std::cout << "InnerProduct: test " << numTest << "mv1 = "<< mv1 <<std::endl;
        //std::cout << "InnerProduct: test " << numTest << "mv2 = "<< mv2 <<std::endl;
        //std::cout << "InnerProduct: test " << numTest << "mv3 = "<< mv3 <<std::endl;
        if (mv3.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::scalar] - 1.2) > epsilon){
            std::cout << "InnerProduct: test " << numTest << " failed "<< "mv3 = "<< mv3 <<std::endl;
            return EXIT_FAILURE;
        }
    }

    // grade 1 c3ga with eo and ei
    {
        //std::cout << std::endl;
        numTest++;
        c3ga::Mvec<double> mv1,mv2,mv3;
        mv1[c3ga::E0] = 2.0;
        mv1[c3ga::E2] = 3.0;
        mv1[c3ga::Ei] = 4.0;
        mv2[c3ga::E0] = 1;
        mv3 = mv1 | mv2;
        //std::cout << "InnerProduct: test " << numTest << "mv1 = "<< mv1 <<std::endl;
        //std::cout << "InnerProduct: test " << numTest << "mv2 = "<< mv2 <<std::endl;
        //std::cout << "InnerProduct: test " << numTest << "mv3 = "<< mv3 << std::endl;
        if (mv3.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::scalar] + 4.0) > epsilon){
            std::cout << "InnerProduct: test " << numTest << " failed "<< " mv3"<< mv3 <<std::endl;
            return EXIT_FAILURE;
        }
    }

    // grade 3 c3ga
    {
        numTest++;
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E012] = 2.0;
        mv1[c3ga::E023] = 3.0;
        mv1[c3ga::E23i] = 4.0;
        mv2[c3ga::E023] = 1;
        mv3 = mv1 | mv2;
        if (mv3.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::scalar] - 4.0) > epsilon){
            std::cout << "InnerProduct: test " << numTest << " failed "<< "mv3 = "<< mv3 <<std::endl;
            return EXIT_FAILURE;
        }
    }

    // different grades
    {
        numTest++;
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E02] = 1;
        mv3 = mv1 | mv2;
        if (mv3.grade() != 0){
            std::cout << "InnerProduct: test " << numTest << " failed "<< "grade is "<< mv3.grade()<<"; result is: " << mv3 <<std::endl;
            return EXIT_FAILURE;
        }
        if (fabs(mv3[c3ga::scalar]) > epsilon){
            return EXIT_FAILURE;
        }

    }

    // scalar
    {
        numTest++;
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::scalar] = 42;
        mv3 = mv1 | mv2;
        if (mv3.grade() != 3){
            std::cout << "InnerProduct: test " << numTest << " failed "<< "grade is "<< mv3.grade() <<std::endl;
            return EXIT_FAILURE;
        }
        if (fabs(mv3[c3ga::scalar]) > epsilon){
            return EXIT_FAILURE;
        }
    }

    // scalar
    {
        numTest++;
        c3ga::Mvec<double> mv1, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv3 = mv1 | 42;
        if (mv3.grade() != 3) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::scalar]) > epsilon){
            std::cout << "InnerProduct: test " << numTest << " failed "<< ", mv3 is "<< mv3 <<std::endl;
            return EXIT_FAILURE;
        }

    }

    // scalar
    {
        numTest++;
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::scalar] = 42;
        mv3 = mv2 | mv1;
        if (mv3.grade() != 3){
            std::cout << "InnerProduct: test " << numTest << " failed "<< "grade is "<< mv3.grade() <<std::endl;
            return EXIT_FAILURE;
        }
//        std::cout << "InnerProduct: test " << numTest << " success"<<std::endl;
    }

    // scalar
    {
        numTest++;
        c3ga::Mvec<double> mv1, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv3 = 42 | mv1;
        if (mv3.grade() != 3){
            std::cout << "InnerProduct: test " << numTest << " failed "<< "grade is "<< mv3.grade() <<std::endl;
            return EXIT_FAILURE;
        }
//        std::cout << "InnerProduct: test " << numTest << " success"<<std::endl;
    }

    // scalar, vector and bivector
    {
        numTest++;
        c3ga::Mvec<double> mv1, mv2;
        mv1[c3ga::E12]  = 2.0;
        mv1[c3ga::E2]   = 3.0;
        mv2[c3ga::E12]  = 4.0;
        mv2[c3ga::E23]  = 6.0;
        mv2[c3ga::E123] = 5.0;
        c3ga::Mvec<double> mv3 = mv1 | mv2; // -8.00 - 12.00*e1 + 8.00*e3 - 15.00*e1^e3
//        std::cout << "mv1 : "  << mv1 << std::endl;
//        std::cout << "mv2 : "  << mv2 << std::endl;
//        std::cout << "mv3 : "  << mv3 << std::endl;
        if(mv3.grade() != 2) {std::cout << "InnerProduct: test " << numTest << " failed "<< "grade is "<< mv3.grade() <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv3[c3ga::scalar]  + 8.0) > epsilon) {std::cout << "InnerProduct: test " << numTest << " failed "<< ", mv3 is "<< mv3 <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv3[c3ga::E1]  + 12.0) > epsilon) {std::cout << "InnerProduct: test " << numTest << " failed "<< ", mv3 is "<< mv3 <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv3[c3ga::E3]  -  8.0) > epsilon) {std::cout << "InnerProduct: test " << numTest << " failed "<< ", mv3 is "<< mv3 <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv3[c3ga::E13] + 15.0) > epsilon) {std::cout << "InnerProduct: test " << numTest << " failed "<< ", mv3 is "<< mv3 <<std::endl; return EXIT_FAILURE;}
//        std::cout << "test " << numTest << " OK" << std::endl;
    }

    {
        numTest++;
        c3ga::Mvec<double> mv1,mv2;
        mv1[c3ga::scalar]  = 5.0;
        mv1[c3ga::E1]  = 2.0;
        mv1[c3ga::Ei]  = 3.0;
        mv2[c3ga::E12] = 4.0;
        mv2[c3ga::E1i] = 5.0;
        c3ga::Mvec<double> mv3 = mv1 | mv2; // 8.00*e2 + 10.00*ni + 20.00*e1^e2 + 25.00*e1^ni

        if(mv3.grade() != 2) {std::cout << "InnerProduct: test " << numTest << " failed "<< "grade is "<< mv3.grade() <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv3[c3ga::E2] -8.0)  > epsilon) {std::cout << "InnerProduct: test " << numTest << " failed "<< ", mv3 is "<< mv3 <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv3[c3ga::Ei] -10.0) > epsilon) {std::cout << "InnerProduct: test " << numTest << " failed "<< ", mv3 is "<< mv3 <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv3[c3ga::E12]-20.0) > epsilon) {std::cout << "InnerProduct: test " << numTest << " failed "<< ", mv3 is "<< mv3 <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv3[c3ga::E1i]-25.0) > epsilon) {std::cout << "InnerProduct: test " << numTest << " failed "<< ", mv3 is "<< mv3 <<std::endl; return EXIT_FAILURE;}
//        std::cout << "test " << numTest << " OK" << std::endl;
    }

    return EXIT_SUCCESS;
}

int innerProductEqualTest(){

    const double epsilon = 1.0e-7;

    // grade 1
    {
        c3ga::Mvec<double> mv1;
        mv1[c3ga::E1] = 1.2;
        mv1[c3ga::E2] = 2.2;
        c3ga::Mvec<double> mv2;
        mv2[c3ga::E1] = 1;
        mv2[c3ga::E3] = 2;
        mv1 |= mv2;
        if (mv1.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv1[c3ga::scalar] - 1.2) > epsilon) return EXIT_FAILURE;
    }

    // grade 1 c3ga
    {
        c3ga::Mvec<double> mv1,mv2,mv3;
        mv1[c3ga::E0] = 2.0;
        mv1[c3ga::E2] = 3.0;
        mv1[c3ga::Ei] = 4.0;
        mv2[c3ga::E0] = 1;
        mv1 |= mv2;
        if (mv1.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv1[c3ga::scalar] + 4.0) > epsilon) return EXIT_FAILURE;
    }


    // grade 3 c3ga
    {
        c3ga::Mvec<double> mv1, mv2;
        mv1[c3ga::E012] = 2.0;
        mv1[c3ga::E023] = 3.0;
        mv1[c3ga::E23i] = 4.0;
        mv2[c3ga::E023] = 1;
        mv1 |= mv2;
        if (mv1.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv1[c3ga::scalar] - 4.0) > epsilon) return EXIT_FAILURE;
    }


    // different grades
    {
        c3ga::Mvec<double> mv1, mv2;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E02] = 1;
        mv1 |= mv2;
        if (mv1.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv1[c3ga::scalar]) > epsilon) return EXIT_FAILURE;
    }


    // scalar
    {
        c3ga::Mvec<double> mv1, mv2;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::scalar] = 42;
        mv1 |= mv2;
        if (mv1.grade() != 3) return EXIT_FAILURE;
    }

    // scalar
    {
        c3ga::Mvec<double> mv1;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv1 |= 42;
        if (mv1.grade() != 3) return EXIT_FAILURE;
    }

    // scalar
    {
        c3ga::Mvec<double> mv1, mv2;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::scalar] = 42;
        mv1 |= mv2;
        if (mv1.grade() != 3) return EXIT_FAILURE;

    }

    return EXIT_SUCCESS;
}

int leftContractionScalarTest() {

    const double epsilon = 1.0e-7;

    // mv < scalar = mv.scalar x scalar
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv1[c3ga::scalar] = 2.0;
        mv2[c3ga::scalar] = 42;
        mv3 = mv1 < mv2;

        if(mv3.grade() != 0) return EXIT_FAILURE;
        if(fabs(mv3[c3ga::scalar] - (2.0*42.0)) > epsilon) return EXIT_FAILURE;
    }

    // mv < scalar = mv.scalar x scalar
    {
        c3ga::Mvec<double> mv1, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv1[c3ga::scalar] = 5.0;
        mv3 = mv1 < 42;

        if(mv3.grade() != 0) return EXIT_FAILURE;
        if(fabs(mv3[c3ga::scalar] - (5.0*42.0)) > epsilon) return EXIT_FAILURE;
    }

    // scalar < mv = mv x scalar
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::scalar] = 2;
        mv2[c3ga::E123] = 2.0;
        mv2[c3ga::E013] = 3.0;
        mv2[c3ga::E023] = 4.0;
        mv3 = mv1 < mv2;

        if (mv3.grade() != 3) return EXIT_FAILURE;
        if (mv3 != (2 ^ mv2)) return EXIT_FAILURE;
    }

    // scalar < mv
    {
        c3ga::Mvec<double> mv1, mv3;
        mv1[c3ga::scalar] = 2.0;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv3 = 2 < mv1;

        if (mv3.grade() != 3) return EXIT_FAILURE;
        if (mv3 != (2 ^ mv1)) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int leftContractionTest() {

    const double epsilon = 1.0e-7;

    // grade 1
    {
        c3ga::Mvec<double> mv1;
        mv1[c3ga::E1] = 1.2;
        mv1[c3ga::E2] = 2.2;
        c3ga::Mvec<double> mv2;
        mv2[c3ga::E1] = 1;
        mv2[c3ga::E3] = 2;
        c3ga::Mvec<double> mv3 = (mv1 < mv2);

        if (mv3.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::scalar] - 1.2) > epsilon) return EXIT_FAILURE;
    }

    // grade 1 c3ga
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E0] = 2.0;
        mv1[c3ga::E2] = 3.0;
        mv1[c3ga::Ei] = 4.0;
        mv2[c3ga::E0] = 1;
        mv3 = (mv2 < mv1);
        if (mv3.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::scalar] + 4.0) > epsilon) return EXIT_FAILURE;
    }

    // grade 3 c3ga
//    {
//        c3ga::Mvec<double> mv1, mv2, mv3;
//        mv1[c3ga::E012] = 2.0;
//        mv1[c3ga::E023] = 3.0;
//        mv1[c3ga::E23i] = 4.0;
//        mv2[c3ga::E023] = 1;
//        mv3 = mv2 < mv1;
//        std::cout << mv3 << std::endl;
//        if (mv3.grade() != 0) return EXIT_FAILURE;
//        if (fabs(mv3[c3ga::scalar] - 4.0) > epsilon) return EXIT_FAILURE;
//    }

    // different grades
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E02] = 1;
        mv3 = mv1 < mv2;



        if (mv3.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::scalar]) > epsilon) return EXIT_FAILURE;
    }

    // different grades
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E23] = 1;
        mv3 = mv2 < mv1;


        if (mv3.grade() != 1) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::E0] + 4.0) > epsilon) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::E1] + 2.0) > epsilon) return EXIT_FAILURE;
    }

    // different grades with "leftContraction"
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E23] = 1;
        mv3 = c3ga::leftContraction(mv2,mv1);


        if (mv3.grade() != 1) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::E0] + 4.0) > epsilon) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::E1] + 2.0) > epsilon) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


int rightContractionScalarTest() {

    const double epsilon = 1.0e-7;

    // scalar > mv = mv.scalar x scalar
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::scalar] = 42;
        mv2[c3ga::scalar] = 3.0;
        mv2[c3ga::E123] = 2.0;
        mv2[c3ga::E013] = 3.0;
        mv2[c3ga::E023] = 4.0;
        mv3 = mv1 > mv2;
        if (mv3.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::scalar] - (42*3.0)) > epsilon) return EXIT_FAILURE; // grade(mv1) > grade(mv2)
    }

    // scalar > mv = mv.scalar x scalar
    {
        c3ga::Mvec<double> mv1, mv3;
        mv1[c3ga::scalar] = 5.0;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv3 = 42 > mv1;
        if (mv3.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::scalar] - (5.0*42)) > epsilon) return EXIT_FAILURE; // grade(mv1) > grade(mv2)
    }

    // mv > scalar = mv x scalar
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::scalar] = 1.0;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::scalar] = 2;
        mv3 = mv1 > mv2;
        if (mv3.grade() != 3) return EXIT_FAILURE;
        if (mv3 != (2^mv1))   return EXIT_FAILURE;
    }

    // mv > scalar = mv x scalar
    {
        c3ga::Mvec<double> mv1, mv3;
        mv1[c3ga::scalar] = 1.0;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv3 = mv1 > 2;
        if (mv3.grade() != 3) return EXIT_FAILURE;
        if (mv3 != (2^mv1))   return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


int rightContractionTest(){

    const double epsilon = 1.0e-7;

    // grade 1
    {
        c3ga::Mvec<double> mv1;
        mv1[c3ga::E1] = 1.2;
        mv1[c3ga::E2] = 2.2;
        c3ga::Mvec<double> mv2;
        mv2[c3ga::E1] = 1;
        mv2[c3ga::E3] = 2;
        c3ga::Mvec<double> mv3 = (mv1 > mv2);


        if (mv3.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::scalar] - 1.2) > epsilon) return EXIT_FAILURE;
    }

    // grade 1 c3ga
    {
        c3ga::Mvec<double> mv1,mv2,mv3;
        mv1[c3ga::E0] = 2.0;
        mv1[c3ga::E2] = 3.0;
        mv1[c3ga::Ei] = 4.0;
        mv2[c3ga::E0] = 1;
        mv3 = (mv2 > mv1);
        if (mv3.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::scalar] + 4.0) > epsilon) return EXIT_FAILURE;
    }

    // grade 3 c3ga
//    {
//        c3ga::Mvec<double> mv1, mv2, mv3;
//        mv1[c3ga::E012] = 2.0;
//        mv1[c3ga::E023] = 3.0;
//        mv1[c3ga::E23i] = 4.0;
//        mv2[c3ga::E023] = 1;
//        mv3 = mv2 < mv1;
//        std::cout << mv3 << std::endl;
//        if (mv3.grade() != 0) return EXIT_FAILURE;
//        if (fabs(mv3[c3ga::scalar] - 4.0) > epsilon) return EXIT_FAILURE;
//    }

    // different grades
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E02] = 1;
        mv3 = mv2 > mv1;
        if (mv3.grade() != 0) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::scalar]) > epsilon) return EXIT_FAILURE; // grade(mv1) > grade(mv2)
    }

    // different grades
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E23] = 1;
        mv3 = mv1 > mv2;


        if (mv3.grade() != 1) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::E0] + 4.0) > epsilon) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::E1] + 2.0) > epsilon) return EXIT_FAILURE;
    }

    // different grades with "rightContraction"
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E23] = 1;
        mv3 = c3ga::rightContraction(mv1, mv2);


        if (mv3.grade() != 1) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::E0] + 4.0) > epsilon) return EXIT_FAILURE;
        if (fabs(mv3[c3ga::E1] + 2.0) > epsilon) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int geometricProductTest(){

    const double epsilon = 1.0e-7;

    // test 1
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        b[c3ga::E2] = 4.0;
        b[c3ga::E3] = 5.0;
        c = a * b; // 12.00 + 15.00*e2^e3 + -10.00*e3^e1 + 8.00*e1^e2
        if(c.grade() != 2) return EXIT_FAILURE;
        if(fabs(c[c3ga::scalar] - 12.0) > epsilon) return EXIT_FAILURE;
        if(fabs(c[c3ga::E12] -  8.0) > epsilon)    return EXIT_FAILURE;
        if(fabs(c[c3ga::E13] - 10.0) > epsilon)    return EXIT_FAILURE;
        if(fabs(c[c3ga::E23] - 15.0) > epsilon)    return EXIT_FAILURE;
    }

    // test 2
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        b[c3ga::E12] = 4.0;
        b[c3ga::E13] = 5.0;
        c = a * b;  // -12.00*e1 + 8.00*e2 + 10.00*e3 + -15.00*e1^e2^e3
        if(c.grade() != 3) return EXIT_FAILURE;
        if(fabs(c[c3ga::E1] + 12.0) > epsilon)    return EXIT_FAILURE;
        if(fabs(c[c3ga::E2] -  8.0) > epsilon)    return EXIT_FAILURE;
        if(fabs(c[c3ga::E3] - 10.0) > epsilon)    return EXIT_FAILURE;
        if(fabs(c[c3ga::E123] + 15.0) > epsilon)  return EXIT_FAILURE;
    }

    // test 3
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        c = a * 2.0f;  // double
        if(c.grade() != 1) return EXIT_FAILURE;
        if(fabs(c[c3ga::E1] - 4.0) > epsilon)    return EXIT_FAILURE;
        if(fabs(c[c3ga::E2] - 6.0) > epsilon)    return EXIT_FAILURE;
    }

    // test 4
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        c = a * 2;  // int
        if(c.grade() != 1) return EXIT_FAILURE;
        if(fabs(c[c3ga::E1] - 4.0) > epsilon) return EXIT_FAILURE;
        if(fabs(c[c3ga::E2] - 6.0) > epsilon) return EXIT_FAILURE;
    }

    // test 5
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        c = 2 * a;  // int
        if(c.grade() != 1) return EXIT_FAILURE;
        if(fabs(c[c3ga::E1] - 4.0) > epsilon)    return EXIT_FAILURE;
        if(fabs(c[c3ga::E2] - 6.0) > epsilon)    return EXIT_FAILURE;
    }

    // test 6
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        c = 2 * a;  // int
        if(c.grade() != 1) return EXIT_FAILURE;
        if(fabs(c[c3ga::E1] - 4.0) > epsilon)    return EXIT_FAILURE;
        if(fabs(c[c3ga::E2] - 6.0) > epsilon)    return EXIT_FAILURE;
    }

    // test 7
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        a[c3ga::E23] = 4.0;
        a[c3ga::E23i] = 5.0;
        a[c3ga::E023i] = 6.0;
        a[c3ga::E0123i] = 7.0;
        c = (a * a);  // -88.00 + 84.00*e1 + -40.00*ni + 70.00*e1^ni + 30.00*e3^ni + -48.00*no^ni + 16.00*e1^e2^e3 + 56.00*e1^no^ni + 42.00*e1^e3^no^ni + -28.00*e2^e3^no^ni + 24.00*e1^e2^e3^no^ni
//        std::cout << "a * a = " << c;
        if(c.grade() != 5) return EXIT_FAILURE;
        if(fabs(c[c3ga::scalar] + 88.0) > epsilon) return EXIT_FAILURE;
    }

    // test 8  => very important test (do not restrict to outer and inner)
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1i] = 2.0;
        b[c3ga::E02] = 3.0;
        c = (a * b);  // - 6.00*e1^e2 - 6.00*eo^e1^e2^ei
        if(c.grade() != 4) return EXIT_FAILURE;
        if(fabs(c[c3ga::E12] + 6.0) > epsilon)   return EXIT_FAILURE;
        if(fabs(c[c3ga::E012i] + 6.0) > epsilon) return EXIT_FAILURE;
    }

    // test 9 => very important test (do not restrict to outer and inner)
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1i] = 2.0;
        b[c3ga::E023] = 3.0;
        b[c3ga::E012] = 3.0;
        c = (a * b);  // - 6.00*e2 - 6.00*e1^e2^e3 + 6.00*e0^e2^ei + 6.00*eo^e1^e2^e3^ei
        if(c.grade() != 5) return EXIT_FAILURE;
        if(fabs(c[c3ga::E2] + 6.0) > epsilon)     return EXIT_FAILURE;
        if(fabs(c[c3ga::E123] + 6.0) > epsilon)   return EXIT_FAILURE;
        if(fabs(c[c3ga::E02i] - 6.0) > epsilon)   return EXIT_FAILURE;
        if(fabs(c[c3ga::E0123i] - 6.0) > epsilon) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int geometricProductEqualTest(){

    const double epsilon = 1.0e-7;

    // test 1
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        b[c3ga::E2] = 4.0;
        b[c3ga::E3] = 5.0;
        a *= b; // 12.00 + 15.00*e2^e3 + -10.00*e3^e1 + 8.00*e1^e2
        if(a.grade() != 2) return EXIT_FAILURE;
        if(fabs(a[c3ga::scalar] - 12.0) > epsilon) return EXIT_FAILURE;
        if(fabs(a[c3ga::E12] -  8.0) > epsilon)    return EXIT_FAILURE;
        if(fabs(a[c3ga::E13] - 10.0) > epsilon)    return EXIT_FAILURE;
        if(fabs(a[c3ga::E23] - 15.0) > epsilon)    return EXIT_FAILURE;
    }

    // test 2
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        b[c3ga::E12] = 4.0;
        b[c3ga::E13] = 5.0;
        a *= b;  // -12.00*e1 + 8.00*e2 + 10.00*e3 + -15.00*e1^e2^e3
        if(a.grade() != 3) return EXIT_FAILURE;
        if(fabs(a[c3ga::E1] + 12.0) > epsilon)    return EXIT_FAILURE;
        if(fabs(a[c3ga::E2] -  8.0) > epsilon)    return EXIT_FAILURE;
        if(fabs(a[c3ga::E3] - 10.0) > epsilon)    return EXIT_FAILURE;
        if(fabs(a[c3ga::E123] + 15.0) > epsilon)  return EXIT_FAILURE;
    }

    // test 3
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        a *= 2.0f;  // double
        if(a.grade() != 1) return EXIT_FAILURE;
        if(fabs(a[c3ga::E1] - 4.0) > epsilon)    return EXIT_FAILURE;
        if(fabs(a[c3ga::E2] - 6.0) > epsilon)    return EXIT_FAILURE;
    }

    // test 4
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        a *= 2;  // int
        if(a.grade() != 1) return EXIT_FAILURE;
        if(fabs(a[c3ga::E1] - 4.0) > epsilon) return EXIT_FAILURE;
        if(fabs(a[c3ga::E2] - 6.0) > epsilon) return EXIT_FAILURE;
    }

    // test 5
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        a[c3ga::E23] = 4.0;
        a[c3ga::E23i] = 5.0;
        a[c3ga::E023i] = 6.0;
        a[c3ga::E0123i] = 7.0;
        a *= a;  // -88.00 + 84.00*e1 + -40.00*ni + 70.00*e1^ni + 30.00*e3^ni + -48.00*no^ni + 16.00*e1^e2^e3 + 56.00*e1^no^ni + 42.00*e1^e3^no^ni + -28.00*e2^e3^no^ni + 24.00*e1^e2^e3^no^ni
//        std::cout << "a *= a " << a;
        if(a.grade() != 5) return EXIT_FAILURE;
        if(fabs(a[c3ga::scalar] + 88.0) > epsilon) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


int scalarProductTest(){

    const double epsilon = 1.0e-7;

    // test 1
    {
        c3ga::Mvec<double> mv1, mv2;
        mv1[c3ga::E12]  = 2.0;
        mv1[c3ga::E2]   = 3.0;
        mv2[c3ga::E12]  = 4.0;
        mv2[c3ga::E23]  = 6.0;
        mv2[c3ga::E123] = 5.0;
        c3ga::Mvec<double> mv3 = mv1.scalarProduct(mv2); //
        if(mv3.grade() != 0) return EXIT_FAILURE;
        if(fabs(mv3[c3ga::scalar] + 2.0*4.0) > epsilon) return EXIT_FAILURE;
    }

    // test 2
    {
        c3ga::Mvec<double> mv1,mv2;
        mv1[c3ga::E1]  = 2.0;
        mv1[c3ga::Ei]  = 3.0;
        mv2[c3ga::E12] = 4.0;
        mv2[c3ga::E1i] = 5.0;
        c3ga::Mvec<double> mv3 = mv1.scalarProduct(mv2); //
        if(mv3.grade() != 0) return EXIT_FAILURE;
        if(fabs(mv3[c3ga::scalar]) > epsilon) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


int hestenesProductTest(){

    const double epsilon = 1.0e-7;

    // test 1
    {
        c3ga::Mvec<double> mv1, mv2;
        mv1[c3ga::E12]  = 2.0;
        mv1[c3ga::E2]   = 3.0;
        mv2[c3ga::E12]  = 4.0;
        mv2[c3ga::E23]  = 6.0;
        mv2[c3ga::E123] = 5.0;
        c3ga::Mvec<double> mv3 = mv1.hestenesProduct(mv2); // -12.00*e1 + 8.00*e3 - 15.00*e1^e3
        if(mv3.grade() != 2) return EXIT_FAILURE;
        if(fabs(mv3[c3ga::E1]  + 12.0) > epsilon) return EXIT_FAILURE;
        if(fabs(mv3[c3ga::E3]  -  8.0) > epsilon) return EXIT_FAILURE;
        if(fabs(mv3[c3ga::E13] + 15.0) > epsilon) return EXIT_FAILURE;
    }

    // test 2
    {
        c3ga::Mvec<double> mv1,mv2;
        mv1[c3ga::scalar]  = 5.0;
        mv1[c3ga::E1]  = 2.0;
        mv1[c3ga::Ei]  = 3.0;

        mv2[c3ga::E12] = 4.0;
        mv2[c3ga::E1i] = 5.0;
        c3ga::Mvec<double> mv3 = mv1.hestenesProduct(mv2); // 8e2 + 10ei

        if(mv3.grade() != 1) return EXIT_FAILURE;
        if(fabs(mv3[c3ga::E2] -8.0)  > epsilon) return EXIT_FAILURE;
        if(fabs(mv3[c3ga::Ei] -10.0) > epsilon) return EXIT_FAILURE;
        if(fabs(mv3[c3ga::E12]) > epsilon) return EXIT_FAILURE;
        if(fabs(mv3[c3ga::E1i]) > epsilon) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


int invertTest(){

    const double epsilon = 1.0e-6;

    int numTest =0;

    // test 1
    {
        numTest++;
        c3ga::Mvec<double> a,b;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        b = 1.0 / a;


        if(b.grade() != 1) return EXIT_FAILURE;
        if(fabs(b[c3ga::E1] - 0.153846) > epsilon) return EXIT_FAILURE;
        if(fabs(b[c3ga::E2] - 0.230769) > epsilon) return EXIT_FAILURE;
    }

    // test 2
    {
        numTest++;
        c3ga::Mvec<double> a,b;
        a[c3ga::E13] = 2.0;
        a[c3ga::E23] = 3.0;
        b = 1.0 / a;

        if(b.grade() != 2) return EXIT_FAILURE;
        if(fabs(b[c3ga::E13] + 0.153846) > epsilon) return EXIT_FAILURE;
        if(fabs(b[c3ga::E23] + 0.230769) > epsilon) return EXIT_FAILURE;
    }

    // test 2 with "inv(mv)"
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E13] = 2.0;
        a[c3ga::E23] = 3.0;
        b = a.inv();

        if(b.grade() != 2) return EXIT_FAILURE;
        if(fabs(b[c3ga::E13] + 0.153846) > epsilon) return EXIT_FAILURE;
        if(fabs(b[c3ga::E23] + 0.230769) > epsilon) return EXIT_FAILURE;
    }

    // test 3
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E13] = 6.0;
        a[c3ga::E23] = 4.0;
        b = a / 2;
        if(b.grade() != 2) return EXIT_FAILURE;
        if(fabs(b[c3ga::E13] - 3) > epsilon) return EXIT_FAILURE;
        if(fabs(b[c3ga::E23] - 2) > epsilon) return EXIT_FAILURE;
    }

    // test 4
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        a[c3ga::E13] = 2.0;
        a[c3ga::E23] = 3.0;
        b = 1.0 / a;
//        std::cout << "inv: " << b << "    with norm : " << b.norm()<< "    with sqrt quadratic norm : "<< sqrt(b.quadraticNorm()) << "   "<<std::endl;  // norm should be 5.1 (error with the inner product)
        if(b.grade() != 2) return EXIT_FAILURE;
        if(fabs(b[c3ga::E13] + 0.08) > 1e-2) return EXIT_FAILURE; // not sure about the result: gaviewer : inverse() or general_inverse() does not provide the same answer
        if(fabs(b[c3ga::E23] + 0.12) > 1e-2) return EXIT_FAILURE;
        if(fabs(b[c3ga::E1]  - 0.08) > 1e-2) return EXIT_FAILURE;
        if(fabs(b[c3ga::E2]  - 0.12) > 1e-2) return EXIT_FAILURE;
    }

    // test 5
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        b[c3ga::E13] = 2.0;
        b[c3ga::E23] = 3.0;
        c = a / b;
        c.roundZero(); // else get epsilon component in e123
        if(c.grade() != 1) return EXIT_FAILURE;
        if(fabs(c[c3ga::E3] + 1.0) > epsilon) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int invertEqualTest(){

    const double epsilon = 1.0e-6;

    // test 1
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 6.0;
        a[c3ga::E2] = 4.0;
        a /= 2.0;
        if(a.grade() != 1) return EXIT_FAILURE;
        if(fabs(a[c3ga::E1] - 3) > epsilon) return EXIT_FAILURE;
        if(fabs(a[c3ga::E2] - 2) > epsilon) return EXIT_FAILURE;
    }

    // test 2
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E1] = 6.0;
        a[c3ga::E2] = 4.0;
        b[c3ga::scalar] = 2;
        a /= b;
        if(a.grade() != 1) return EXIT_FAILURE;
        if(fabs(a[c3ga::E1] - 3) > epsilon) return EXIT_FAILURE;
        if(fabs(a[c3ga::E2] - 2) > epsilon) return EXIT_FAILURE;
    }

    // test 3
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        b[c3ga::E13] = 2.0;
        b[c3ga::E23] = 3.0;
        a /= b;
        a.roundZero(); // else get epsilon component in e123
        if(a.grade() != 1) return EXIT_FAILURE;
        if(fabs(a[c3ga::E3] + 1.0) > epsilon) return EXIT_FAILURE;
    }

    // test 3
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        a /= a;
        a.roundZero();
        if(a.grade() != 0) return EXIT_FAILURE;
        if(fabs(a[c3ga::scalar] - 1.0) > epsilon) return EXIT_FAILURE;
    }

    // test 4
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        a[c3ga::E13] = 2.0;
        a[c3ga::E23] = 3.0;
        b[c3ga::scalar] = 1.0;
        b /= a;
//        std::cout << "inv: " << b << "    with norm : " << b.norm() << "   ";  // norm should be 5.1 (error with the inner product)
        if(b.grade() != 2) return EXIT_FAILURE;
        if(fabs(b[c3ga::E13] + 0.08) > 1e-2) return EXIT_FAILURE;
        if(fabs(b[c3ga::E23] + 0.12) > 1e-2) return EXIT_FAILURE;
        if(fabs(b[c3ga::E1]  - 0.08) > 1e-2) return EXIT_FAILURE;
        if(fabs(b[c3ga::E2]  - 0.12) > 1e-2) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int reverseTest(){

    double epsilon = 1.0e-7;

    //  scalar
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::scalar] = 42;
        c3ga::Mvec<double> mv2 = ~mv;
        if(mv2.grade() != 0) return EXIT_FAILURE;
        if(fabs(mv2[c3ga::scalar] - 42.0) > epsilon) return EXIT_FAILURE;
    }

    // grade 1
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 1;
        mv[c3ga::E2] = 2;
        c3ga::Mvec<double> mv2 = ~mv;
        if(mv2.grade() != 1) return EXIT_FAILURE;
        if(fabs(mv2[c3ga::E1] -1.0) > epsilon) return EXIT_FAILURE;
        if(fabs(mv2[c3ga::E2] -2.0) > epsilon) return EXIT_FAILURE;
    }

    // grade 2
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E12] = 2;
        mv[c3ga::E23] = 3;
        c3ga::Mvec<double> mv2 = ~mv;
        if(mv2.grade() != 2) return EXIT_FAILURE;
        if(fabs(mv2[c3ga::E12] + 2.0) > epsilon) return EXIT_FAILURE;
        if(fabs(mv2[c3ga::E23] + 3.0) > epsilon) return EXIT_FAILURE;
    }

    // grade 2 with "reverse()"
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E12] = 2;
        mv[c3ga::E23] = 3;
        c3ga::Mvec<double> mv2 = mv.reverse();
        if(mv2.grade() != 2) return EXIT_FAILURE;
        if(fabs(mv2[c3ga::E12] + 2.0) > epsilon) return EXIT_FAILURE;
        if(fabs(mv2[c3ga::E23] + 3.0) > epsilon) return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int dualTest(){

    double epsilon = 1.0e-7;

    // grade 1
//    {
//        c3ga::Mvec<double> mv;
//        mv[c3ga::E1] = 2;
//        mv[c3ga::E2] = 3;
//        c3ga::Mvec<double> mv2 = mv.dual();
//        if(mv2.grade() != 2) return EXIT_FAILURE;
//        if(fabs(mv2[c3ga::E23] + 2.0) > epsilon) return EXIT_FAILURE;
//        if(fabs(mv2[c3ga::E13] - 3.0) > epsilon) return EXIT_FAILURE;
//    }


    // c3ga point
    {
        c3ga::Mvec<float> mv;
        mv[c3ga::E0] = 1;
        mv[c3ga::E1] = 2;
        mv[c3ga::Ei] = 2;
        c3ga::Mvec<float> mv2 = mv.dual();
        if(mv2.grade() != 4) {std::cout << "dual: test " << 2 << " failed "<< ", mv3 is "<< mv2 <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv2[c3ga::E0123] + 1.0) > epsilon) {std::cout << "dual: test " << 2 << " failed "<< ", mv3 is "<< mv2 <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv2[c3ga::E123i] + 2.0) > epsilon) {std::cout << "dual: test " << 2 << " failed "<< ", mv3 is "<< mv2 <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv2[c3ga::E023i] + 2.0) > epsilon) {std::cout << "dual: test " << 2 << " failed "<< ", mv3 is "<< mv2 <<std::endl; return EXIT_FAILURE;}
    }

    // c3ga pair point
    {
        c3ga::Mvec<float> mv;
        mv[c3ga::E01] = -2;
        mv[c3ga::E12] =  4;
        mv[c3ga::E1i] =  4;
        mv[c3ga::E2i] = -4;
        c3ga::Mvec<float> mv2 = mv.dual();

        if(mv2.grade() != 3) {std::cout << "dual: test " << 2 << " failed "<< ", mv3 is "<< mv2 <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv2[c3ga::E23i] + 4.0) > epsilon) {std::cout << "dual: test " << 3 << " failed "<< ", mv3 is "<< mv2 <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv2[c3ga::E023] - 2.0) > epsilon) {std::cout << "dual: test " << 3 << " failed "<< ", mv3 is "<< mv2 <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv2[c3ga::E13i] + 4.0) > epsilon) {std::cout << "dual: test " << 3 << " failed "<< ", mv3 is "<< mv2 <<std::endl; return EXIT_FAILURE;}
        if(fabs(mv2[c3ga::E03i] + 4.0) > epsilon) {std::cout << "dual: test " << 3 << " failed "<< ", mv3 is "<< mv2 <<std::endl; return EXIT_FAILURE;}
    }
    return EXIT_SUCCESS;
}

int normTest(){

    double epsilon = 1.0e-7;

    // grade 1
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 2;
        mv[c3ga::E2] = 3;
        if(fabs(mv.norm() - 3.605551275) > epsilon) return EXIT_FAILURE;
    }

    // grade 2
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E13] = 2;
        mv[c3ga::E23] = 3;
        if(fabs(mv.norm() - 3.605551275) > epsilon) return EXIT_FAILURE;
    }

    // non homogeneous multivector
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 2;
        mv[c3ga::E2] = 3;
        mv[c3ga::E13] = 2;
        mv[c3ga::E23] = 3;
        if(fabs(mv.norm() - 5.099019513) > epsilon) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int quadraticNormTest(){

    double epsilon = 1.0e-7;

    // grade 1
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 2;
        mv[c3ga::E2] = 3;
        if(fabs(mv.quadraticNorm() - 13) > epsilon) return EXIT_FAILURE;
    }

    // grade 2
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E13] = 2;
        mv[c3ga::E23] = 3;
//        std::cout << "quadraticNorm : " << mv.quadraticNorm() << " ";
        if(fabs(mv.quadraticNorm() - 13) > epsilon) return EXIT_FAILURE;
    }

    // non homogeneous multivector
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 2;
        mv[c3ga::E2] = 3;
        mv[c3ga::E13] = 2;
        mv[c3ga::E23] = 3;
//        std::cout << mv.quadraticNorm() << " ";
        if(fabs(mv.quadraticNorm() - 26) > epsilon) return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int pseudoscalarTest(){

    c3ga::Mvec<double> mv = c3ga::I<double>();
    std::vector<unsigned int> g = mv.grades();
    if(g.size() != 1 )       return EXIT_FAILURE;
    if(g[0] != 5)            return EXIT_FAILURE;
    if(mv[c3ga::E0123i] != 1) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int operatorConstBracket(){

    double epsilon = 1.0e-7;

    const c3ga::Mvec<double> mvConst = 2.0*c3ga::e0<double>();
    c3ga::Mvec<double> mvOut;
    mvOut[c3ga::E0] = mvConst[c3ga::E0];

    if(fabs(mvOut[c3ga::E0] - mvConst[c3ga::E0]) > epsilon) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int pseudoscalarInverseTest(){
    c3ga::Mvec<double> mv = c3ga::Iinv<double>();
    std::vector<unsigned int> g = mv.grades();
    if(g.size() != 1 )        return EXIT_FAILURE;
    if(g[0] != 5)             return EXIT_FAILURE;
    if(mv[c3ga::E0123i] != -1) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

std::string success(const int value){
    if(value == EXIT_SUCCESS) return std::string("success");
    if(value == EXIT_FAILURE) return std::string("failure");
    return std::string("unknown");
}


int main(){

    std::cout << "metric              : \n" << c3ga::metric << std::endl;
    std::cout << "Mvec()              : " << success(constructorDefaultTest()) << std::endl;
    std::cout << "Mvec(mv)            : " << success(constructorCopyTest()) << std::endl;
    std::cout << "Mvec(scalar)        : " << success(constructorScalarTest()) << std::endl;
    std::cout << "Mvec<T>(mv<U>)      : " << success(constructorTemplateConverterTest()) << std::endl;
    std::cout << "double(mv)          : " << success(castTest()) << std::endl;
    std::cout << "operator =          : " << success(operatorEqualTest()) << std::endl;
    std::cout << "operator ==         : " << success(operatorEqualEqualTest()) << std::endl;
    std::cout << "operator !=         : " << success(operatorNotEqualTest()) << std::endl;
    std::cout << "operator +          : " << success(operatorPlusTest()) << std::endl;
    std::cout << "operator +=         : " << success(operatorPlusEqualTest()) << std::endl;
    std::cout << "operator -          : " << success(operatorMinusTest()) << std::endl;
    std::cout << "operator - (unary)  : " << success(operatorUnaryMinusTest()) << std::endl;
    std::cout << "operator -=         : " << success(operatorMinusEqualTest()) << std::endl;
    std::cout << "primal^dual (recurs): " << success(wedgePrimalDualTest()) << std::endl;
    std::cout << "operator ^          : " << success(wedgeTest()) << std::endl;
    std::cout << "operator ^=         : " << success(wedgeEqualTest()) << std::endl;
    std::cout << "operator |          : " << success(innerProductTest()) << std::endl;
    std::cout << "operator |=         : " << success(innerProductEqualTest()) << std::endl;
    std::cout << "operator <          : " << success(leftContractionTest()) << std::endl;
    std::cout << "operator < scalar   : " << success(leftContractionScalarTest()) << std::endl;
    std::cout << "operator >          : " << success(rightContractionTest()) << std::endl;
    std::cout << "operator > scalar   : " << success(rightContractionScalarTest()) << std::endl;
    std::cout << "operator *          : " << success(geometricProductTest()) << std::endl;
    std::cout << "operator *=         : " << success(geometricProductEqualTest()) << std::endl;
    std::cout << "operator /          : " << success(invertTest()) << std::endl;
    std::cout << "operator /=         : " << success(invertEqualTest()) << std::endl;
    std::cout << "operator ~          : " << success(reverseTest()) << std::endl;
    std::cout << "mv.dual()           : " << success(dualTest()) << std::endl;
    std::cout << "mv.norm()           : " << success(normTest()) << std::endl;
    std::cout << "quadraticNorm(mv)   : " << success(quadraticNormTest()) << std::endl;
    std::cout << "scalar product      : " << success(scalarProductTest()) << std::endl;
    std::cout << "Hestenes product    : " << success(hestenesProductTest()) << std::endl;
    std::cout << "mv.grade()          : " << success(gradeTest()) << std::endl;
    std::cout << "mv.isGrade()        : " << success(isGradeTest()) << std::endl;
    std::cout << "mv.grades()         : " << success(gradesTest()) << std::endl;
    std::cout << "mv1.sameGrade(mv2)  : " << success(sameGradesTest()) << std::endl;
    std::cout << "isHomogeneous()     : " << success(isHomogeneousTest()) << std::endl;
    std::cout << "mv.roundZero()      : " << success(roundZeroTest()) << std::endl;
    std::cout << "mv.isEmpty()        : " << success(isEmptyTest()) << std::endl;
    std::cout << "mv.clear()          : " << success(clearTest()) << std::endl;
    std::cout << "mv[basis]           : " << success(accessOperatorTest()) << std::endl;
    std::cout << "const mv[basis]     : " << success(operatorConstBracket()) << std::endl;
    std::cout << "I()                 : " << success(pseudoscalarTest()) << std::endl;
    std::cout << "Inv()               : " << success(pseudoscalarInverseTest()) << std::endl;
    std::cout << "e12()               : " << success(staticGetTest()) << std::endl;
    std::cout << "mv.e12()            : " << success(methodGetTest()) << std::endl;

    return EXIT_SUCCESS;
}
