// Garamon unit tests: c3ga, hand-computed expectations (formerly plugin/c3gaUnitTests.cpp)

#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <c3ga/Mvec.hpp>

using Catch::Approx;


TEST_CASE("c3ga legacy: constructorDefaultTest", "[c3ga][legacy]") {
    c3ga::Mvec<double> mv;
    CHECK(mv.isEmpty());
}

TEST_CASE("c3ga legacy: constructorCopyTest", "[c3ga][legacy]") {

    // empty vector
    {
        c3ga::Mvec<double> mv;
        c3ga::Mvec<double> mv2(mv);
        CHECK(mv2.isEmpty());
    }

    // any multivector
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 1;
        mv[c3ga::E12i] = 3;
        mv[c3ga::E0123] = 4;
        c3ga::Mvec<double> mv2(mv);
        CHECK(mv2[c3ga::E1] == 1);
        CHECK(mv2[c3ga::E12i] == 3);
        CHECK(mv2[c3ga::E0123] == 4);
    }
}

TEST_CASE("c3ga legacy: constructorTemplateConverterTest", "[c3ga][legacy]") {

    c3ga::Mvec<double> mv;
    mv[c3ga::E1] = 1;
    mv[c3ga::E12i] = 3;
    mv[c3ga::E0123] = 4;
    c3ga::Mvec<double> mv2(mv);
    CHECK(mv2[c3ga::E1] == 1);
    CHECK(mv2[c3ga::E12i] == 3);
    CHECK(mv2[c3ga::E0123] == 4);
}

TEST_CASE("c3ga legacy: constructorScalarTest", "[c3ga][legacy]") {
    c3ga::Mvec<double> mv(5);

    CHECK(mv[c3ga::scalar] == 5);
    CHECK(mv.grade() == 0);

    double a = 5;
    c3ga::Mvec<double> mv2(a);
    CHECK(mv2[c3ga::scalar] == 5);
    CHECK(mv2.grade() == 0);
}

TEST_CASE("c3ga legacy: castTest", "[c3ga][legacy]") {

    // empty vector
    {
        c3ga::Mvec<double> mv;
        CHECK(int(mv) == 0);
    }

    // scalar
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::scalar] = 42;
        CHECK(int(mv) == 42);
    }

    // vector
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 42;
        CHECK(int(mv) == 0);
    }

    // tri-vector
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E12i] = 42;
        CHECK(double(mv) == 0.0f);
    }

    // multivrctor with scalar
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::scalar] = 42;
        mv[c3ga::E12i] = 42;
        CHECK(int(mv) == 42);
    }

    // multivector without scalar
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E12]  = 42;
        mv[c3ga::E12i] = 42;
        CHECK(double(mv) == 0.0);
    }
}

TEST_CASE("c3ga legacy: operatorEqualTest", "[c3ga][legacy]") {

    // operator =
    c3ga::Mvec<double> a;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] = -6.0;
    c3ga::Mvec<double> b = a;
    CHECK(b[c3ga::E12] == 5.0);
    CHECK(b[c3ga::E13] == -6.0);

    // operator = (const Mvec)
    const c3ga::Mvec<double> c(a);
    c3ga::Mvec<double> d = c;
    CHECK(d[c3ga::E12] == 5.0);
    CHECK(d[c3ga::E13] == -6.0);

    // operator = (non ref Mvec)
    c3ga::Mvec<double> e = c3ga::Mvec<double>(c);
    CHECK(e[c3ga::E12] == 5.0);
    CHECK(e[c3ga::E13] == -6.0);
}

TEST_CASE("c3ga legacy: operatorEqualEqualTest", "[c3ga][legacy]") {

    c3ga::Mvec<double> a,b;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] = -6.0;
    b[c3ga::E12] =  a[c3ga::E12];
    b[c3ga::E13] =  a[c3ga::E13];
    CHECK(a==b);

    b[c3ga::E23] = 2.0;
    CHECK(a != b);
}

TEST_CASE("c3ga legacy: operatorNotEqualTest", "[c3ga][legacy]") {

    c3ga::Mvec<double> a,b;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] = -6.0;
    b[c3ga::E12] =  a[c3ga::E12];
    b[c3ga::E13] =  a[c3ga::E13];
    CHECK(a == b);

    b[c3ga::E23] = 2.0;
    CHECK(a!=b);
}

TEST_CASE("c3ga legacy: operatorPlusEqualTest", "[c3ga][legacy]") {
    c3ga::Mvec<double> a;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] = -6.0;
    c3ga::Mvec<double> b;
    b[c3ga::E13] = 1.0;
    b[c3ga::E1] =  2.0;
    b[c3ga::E3] =  1.0;
    a += b;
    CHECK(a[c3ga::E12] == 5.0);
    CHECK(a[c3ga::E13] == -5.0);
    CHECK(a[c3ga::E1] == 2.0);
    CHECK(a[c3ga::E3] == 1.0);
    std::vector<unsigned int> g = a.grades();
    CHECK(g[0] == 1);
    CHECK(g[1] == 2);
    CHECK(a.grade() == 2);


    a += 8.0; // from constructor with scalar

    CHECK(a[c3ga::scalar] == 8.0);
    g = a.grades();
    CHECK(g[0] == 0);
    CHECK(g[1] == 1);
    CHECK(g[2] == 2);
    CHECK(a.grade() == 2);
}


TEST_CASE("c3ga legacy: operatorPlusTest", "[c3ga][legacy]") {

    c3ga::Mvec<double> a;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] = -6.0;

    c3ga::Mvec<double> b;
    b[c3ga::E13] = 1.0;
    b[c3ga::E1] =  2.0;
    b[c3ga::E3] =  1.0;

    c3ga::Mvec<double> c;
    c = a + b;

    CHECK(c[c3ga::E12] == 5.0);
    CHECK(c[c3ga::E13] == -5.0);
    CHECK(c[c3ga::E1] == 2.0);
    CHECK(c[c3ga::E3] == 1.0);

    std::vector<unsigned int> g = c.grades();
    CHECK(g[0] == 1);
    CHECK(g[1] == 2);
    CHECK(c.grade() == 2);

    c3ga::Mvec<double> d;
    d = a + 42; // int
    CHECK(d[c3ga::E12] == 5.0);
    CHECK(d[c3ga::E13] == -6.0);
    CHECK(d[c3ga::scalar] == 42.0);
    g = d.grades();
    CHECK(g[0] == 0);
    CHECK(g[1] == 2);
    CHECK(d.grade() == 2);

    c3ga::Mvec<double> e;
    e = a + 42.0f; // double
    CHECK(e[c3ga::E12] == 5.0);
    CHECK(e[c3ga::E13] == -6.0);
    CHECK(e[c3ga::scalar] == 42.0);
    g = e.grades();
    CHECK(g[0] == 0);
    CHECK(g[1] == 2);
    CHECK(e.grade() == 2);

    c3ga::Mvec<double> f;
    f = 42 + a; // int
    CHECK(f[c3ga::E12] == 5.0);
    CHECK(f[c3ga::E13] == -6.0);
    CHECK(f[c3ga::scalar] == 42.0);
    g = f.grades();
    CHECK(g[0] == 0);
    CHECK(g[1] == 2);
    CHECK(f.grade() == 2);

    c3ga::Mvec<double> h;
    h = 42.0f + a;  // double
    CHECK(h[c3ga::E12] == 5.0);
    CHECK(h[c3ga::E13] == -6.0);
    CHECK(h[c3ga::scalar] == 42.0);
    g = h.grades();
    CHECK(g[0] == 0);
    CHECK(g[1] == 2);
    CHECK(h.grade() == 2);
}

TEST_CASE("c3ga legacy: operatorMinusTest", "[c3ga][legacy]") {

    c3ga::Mvec<double> a;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] = -6.0;

    c3ga::Mvec<double> b;
    b[c3ga::E13] = 1.0;
    b[c3ga::E1] =  2.0;
    b[c3ga::E3] =  1.0;

    c3ga::Mvec<double> c;
    c = a - b;

    CHECK(c[c3ga::E12] == 5.0);
    CHECK(c[c3ga::E13] == -7.0);
    CHECK(c[c3ga::E1] == -2.0);
    CHECK(c[c3ga::E3] == -1.0);

    std::vector<unsigned int> g = c.grades();
    CHECK(g[0] == 1);
    CHECK(g[1] == 2);
    CHECK(c.grade() == 2);

    c3ga::Mvec<double> d;
    d = a - 42; // int
    CHECK(d[c3ga::E12] == 5.0);
    CHECK(d[c3ga::E13] == -6.0);
    CHECK(d[c3ga::scalar] == -42.0);
    g = d.grades();
    CHECK(g[0] == 0);
    CHECK(g[1] == 2);
    CHECK(d.grade() == 2);

    c3ga::Mvec<double> e;
    e = a - 42.0f; // double
    CHECK(e[c3ga::E12] == 5.0);
    CHECK(e[c3ga::E13] == -6.0);
    CHECK(e[c3ga::scalar] == -42.0);
    g = e.grades();
    CHECK(g[0] == 0);
    CHECK(g[1] == 2);
    CHECK(e.grade() == 2);

    c3ga::Mvec<double> f;
    f = 42 - a; // int
    CHECK(f[c3ga::E12] == -5.0);
    CHECK(f[c3ga::E13] == 6.0);
    CHECK(f[c3ga::scalar] == 42.0);
    g = f.grades();
    CHECK(g[0] == 0);
    CHECK(g[1] == 2);
    CHECK(f.grade() == 2);

    c3ga::Mvec<double> h;
    h = 42.0f - a;  // double
    CHECK(h[c3ga::E12] == -5.0);
    CHECK(h[c3ga::E13] == 6.0);
    CHECK(h[c3ga::scalar] == 42.0);
    g = h.grades();
    CHECK(g[0] == 0);
    CHECK(g[1] == 2);
    CHECK(h.grade() == 2);
}

TEST_CASE("c3ga legacy: operatorMinusEqualTest", "[c3ga][legacy]") {
    c3ga::Mvec<double> a;
    a[c3ga::E12] =  5.0;
    a[c3ga::E13] = -6.0;
    c3ga::Mvec<double> b;
    b[c3ga::E1]  =  2.0;
    b[c3ga::E3]  =  1.0;
    b[c3ga::E13] =  1.0;

    a -= b;
    CHECK(a[c3ga::E12] == 5.0);
    CHECK(a[c3ga::E13] == -7.0);
    CHECK(a[c3ga::E1] == -2.0);
    CHECK(a[c3ga::E3] == -1.0);
    std::vector<unsigned int> g = a.grades();
    CHECK(g[0] == 1);
    CHECK(g[1] == 2);
    CHECK(a.grade() == 2);

    a-= 8.0;
    CHECK(a[c3ga::scalar] == -8.0);
    g = a.grades();
    CHECK(g[0] == 0);
    CHECK(g[1] == 1);
    CHECK(g[2] == 2);
}

TEST_CASE("c3ga legacy: operatorUnaryMinusTest", "[c3ga][legacy]") {

    c3ga::Mvec<double> a;
    a[c3ga::E1] =  2.0;
    a[c3ga::E3] =  1.0;
    a[c3ga::E12] =  5.0;
    a[c3ga::E03] = -6.0;
    a[c3ga::E23i] = -7.0;
    c3ga::Mvec<double> b;
    b = - a;

    CHECK(b[c3ga::E1] == -2.0);
    CHECK(b[c3ga::E3] == -1.0);
    CHECK(b[c3ga::E12] == -5.0);
    CHECK(b[c3ga::E03] == 6.0);
    CHECK(b[c3ga::E23i] == 7.0);
}

TEST_CASE("c3ga legacy: gradeTest", "[c3ga][legacy]") {

    c3ga::Mvec<double> a;
    CHECK(a.grade() == 0);

    a[c3ga::scalar] =  5.0;
    CHECK(a.grade() == 0);

    a[c3ga::E1] =  5.0;
    a[c3ga::E3] = -6.0;
    CHECK(a.grade() == 1);

    a[c3ga::E12] =  5.0;
    a[c3ga::E03] = -6.0;
    CHECK(a.grade() == 2);

    a[c3ga::E12i] =  5.0;
    a[c3ga::E013] = -6.0;
    CHECK(a.grade() == 3);

    a[c3ga::E012i] =  5.0;
    a[c3ga::E0123] = -6.0;
    CHECK(a.grade() == 4);

    a[c3ga::E0123i] = 5.0; // ne fonctionne pas
    CHECK(a.grade() == 5);
}

TEST_CASE("c3ga legacy: isGradeTest", "[c3ga][legacy]") {

    c3ga::Mvec<double> a;
    CHECK(a.isGrade(0) != false);
    CHECK(a.isGrade(1) != true);

    a[c3ga::scalar] =  5.0;
    CHECK(a.isGrade(0) != false);
    CHECK(a.isGrade(1) != true);

    a[c3ga::E1] =  5.0;
    a[c3ga::E3] = -6.0;
    CHECK(a.isGrade(1) != false);
    CHECK(a.isGrade(2) != true);

    a[c3ga::E12] =  5.0;
    a[c3ga::E03] = -6.0;
    CHECK(a.isGrade(2) != false);
    CHECK(a.isGrade(3) != true);

    a[c3ga::E12i] =  5.0;
    a[c3ga::E013] = -6.0;
    CHECK(a.isGrade(3) != false);
    CHECK(a.isGrade(4) != true);

    a[c3ga::E012i] =  5.0;
    a[c3ga::E0123] = -6.0;
    CHECK(a.isGrade(4) != false);
    CHECK(a.isGrade(5) != true);

    a[c3ga::E0123i] = 5.0;
    CHECK(a.isGrade(5) != false);
}

TEST_CASE("c3ga legacy: accessOperatorTest", "[c3ga][legacy]") {

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
}

TEST_CASE("c3ga legacy: staticGetTest", "[c3ga][legacy]") {

    // bivector
    {
        c3ga::Mvec<double> mv = c3ga::e12<double>();
        CHECK(mv.grade() == 2);
        CHECK(mv[c3ga::E12] == 1);
    }

    // trivector
    {
        c3ga::Mvec<double> mv = 2 * c3ga::e12i<double>();
        CHECK(mv.grade() == 3);
        CHECK(mv[c3ga::E12i] == 2);
    }

    // quadvector
    {
        c3ga::Mvec<double> mv =  c3ga::e1<double>() ^ c3ga::e2<double>() ^ c3ga::e3i<double>();
        CHECK(mv.grade() == 4);
        CHECK(mv[c3ga::E123i] == 1);
    }
}

TEST_CASE("c3ga legacy: methodGetTest", "[c3ga][legacy]") {

    // bivector
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E12] = 2;
        a[c3ga::E23] = 3;
        b = a.e12();
        CHECK(b.grade() == 2);
        CHECK(b[c3ga::E12] == 2);
        CHECK(b[c3ga::E23] == 0);
    }

    // empty
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E12i] = 2;
        a[c3ga::E23i] = 3;
        b = a.e12();
        CHECK(b.grade() == 0);
    }

    // trivector
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E012] = 2;
        a[c3ga::E023] = 3;
        b = 2 * a.e012() ^ c3ga::ei<double>();
        CHECK(b.grade() == 4);
        CHECK(b[c3ga::E012i] == 4);
        CHECK(b[c3ga::E023i] == 0);
    }
}


TEST_CASE("c3ga legacy: gradesTest", "[c3ga][legacy]") {

    c3ga::Mvec<double> a;

    // empty multivector
    std::vector<unsigned int> g = a.grades();
    CHECK(g.size() == 0);

    // grade 0
    a[c3ga::scalar] =  5.0;
    g = a.grades();
    CHECK(g.size() == 1);
    CHECK(g[0] == 0);

    // grade 1
    a[c3ga::E2] =  5.0;
    g = a.grades();
    CHECK(g.size() == 2);
    CHECK(g[0] == 0);
    CHECK(g[1] == 1);

    // grade 2
    a[c3ga::E13] =  5.0;
    g = a.grades();
    CHECK(g.size() == 3);
    CHECK(g[0] == 0);
    CHECK(g[1] == 1);
    CHECK(g[2] == 2);

    // grade 3
    a[c3ga::E13i] =  5.0;
    g = a.grades();
    CHECK(g.size() == 4);
    CHECK(g[0] == 0);
    CHECK(g[1] == 1);
    CHECK(g[2] == 2);
    CHECK(g[3] == 3);

    // grade 4
    a[c3ga::E123i] =  5.0;
    g = a.grades();
    CHECK(g.size() == 5);
    CHECK(g[0] == 0);
    CHECK(g[1] == 1);
    CHECK(g[2] == 2);
    CHECK(g[3] == 3);
    CHECK(g[4] == 4);

    // grade 4
    a[c3ga::E0123i] =  5.0;
    g = a.grades();
    CHECK(g.size() == 6);
    CHECK(g[0] == 0);
    CHECK(g[1] == 1);
    CHECK(g[2] == 2);
    CHECK(g[3] == 3);
    CHECK(g[4] == 4);
    CHECK(g[5] == 5);

    // non consecutives grades
    c3ga::Mvec<double> b;
    b[c3ga::E1]     =  5.0;
    b[c3ga::E123]   =  5.0;
    b[c3ga::E0123i] =  5.0;
    g = b.grades();
    CHECK(g.size() == 3);
    CHECK(g[0] == 1);
    CHECK(g[1] == 3);
    CHECK(g[2] == 5);
}

TEST_CASE("c3ga legacy: sameGradesTest", "[c3ga][legacy]") {

    // empty vectors
    {
        c3ga::Mvec<double> a,b;
        std::vector<unsigned int> g = a.grades();
        CHECK(g.size() == 0);
    }

    // scalars
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::scalar] = 42;
        b[c3ga::scalar] = 2;
        CHECK(a.sameGrade(b) != false);
    }

    // scalars
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::scalar] = 42;
        b[c3ga::E1] = 2;
        CHECK(a.sameGrade(b) != true);
    }

    // vectors
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E2] = 42;
        b[c3ga::E1] = 2;
        CHECK(a.sameGrade(b) != false);
    }

    // vectors
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E2] = 42;
        b[c3ga::E12] = 2;
        CHECK(a.sameGrade(b) != true);
    }

    // any k-vector
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E2] = 42;
        a[c3ga::E23] = 42;
        a[c3ga::E23i] = 42;
        b[c3ga::E123] = 2;
        CHECK(a.sameGrade(b) != false);
    }

    // any k-vector
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E2] = 42;
        a[c3ga::E23] = 42;
        a[c3ga::E23i] = 42;
        b[c3ga::E123i] = 2;
        CHECK(a.sameGrade(b) != true);
    }
}

TEST_CASE("c3ga legacy: isHomogeneousTest", "[c3ga][legacy]") {
    // grade 0
    {
        c3ga::Mvec<double> a;
        CHECK(a.isHomogeneous() != false);
    }

    // grade 0
    {
        c3ga::Mvec<double> a;
        a[c3ga::scalar] = 42;
        CHECK(a.isHomogeneous() != false);
    }

    // grade 1
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 42;
        CHECK(a.isHomogeneous() != false);
    }

    // grade 2
    {
        c3ga::Mvec<double> a;
        a[c3ga::E13] = 42;
        CHECK(a.isHomogeneous() != false);
    }

    // grade 3
    {
        c3ga::Mvec<double> a;
        a[c3ga::E13i] = 42;
        CHECK(a.isHomogeneous() != false);
    }

    // grade 4
    {
        c3ga::Mvec<double> a;
        a[c3ga::E013i] = 42;
        CHECK(a.isHomogeneous() != false);
    }

    // grade 5
    {
        c3ga::Mvec<double> a;
        a[c3ga::E0123i] = 42;
        CHECK(a.isHomogeneous() != false);
    }

    // non homogeneous 1
    {
        c3ga::Mvec<double> a;
        a[c3ga::E012] = 42;
        a[c3ga::E1]   = 42;
        CHECK(a.isHomogeneous() != true);
    }

    // non homogeneous 2
    {
        c3ga::Mvec<double> a;
        a[c3ga::E012] = 42;
        a[c3ga::scalar] = 42;
        a[c3ga::E12]  = 42;
        CHECK(a.isHomogeneous() != true);
    }
}

TEST_CASE("c3ga legacy: roundZeroTest", "[c3ga][legacy]") {

    c3ga::Mvec<double> mv;
    mv[c3ga::E12]  = 1.0e-10;
    mv[c3ga::E123] =  5.0;
    mv[c3ga::E12i] = -5.0;
    mv[c3ga::E13i] = -1.0e-10;

    mv.roundZero(1.0e-6);

    std::vector<unsigned int> g = mv.grades();
    CHECK(g.size() == 1);
    CHECK(g[0] == 3);
    CHECK(mv[c3ga::E12] == 0.0f); // when doing this, create a zero array for bivectors
    CHECK(mv[c3ga::E123] != 0.0f);
    CHECK(mv[c3ga::E12i] != 0.0f);
    CHECK(mv[c3ga::E13i] == 0.0f);
}

TEST_CASE("c3ga legacy: isEmptyTest", "[c3ga][legacy]") {
    c3ga::Mvec<double> mv;
    CHECK(mv.isEmpty());

    mv[c3ga::E3] = 42.0;
    CHECK_FALSE(mv.isEmpty());

    mv[c3ga::E123] = 42.0;
    CHECK_FALSE(mv.isEmpty());
}

TEST_CASE("c3ga legacy: clearTest", "[c3ga][legacy]") {

    c3ga::Mvec<double> mv;
    mv[c3ga::E3]   = 42.0;
    mv[c3ga::E123] = 42.0;


    // erase specific grade
    mv.clear(3);
    std::vector<unsigned  int> g = mv.grades();
    CHECK(g.size() == 1);
    CHECK(g[0] == 1);

    // erase non existing grade
    mv.clear(2);
    g = mv.grades();
    CHECK(g.size() == 1);
    CHECK(g[0] == 1);

    // erase all
    mv[c3ga::E123]  = 42.0;
    mv[c3ga::E0123] = 42.0;
    mv.clear();
    g = mv.grades();
    CHECK(g.size() == 0);
    CHECK(mv.isEmpty());
}

TEST_CASE("c3ga legacy: wedgeTest", "[c3ga][legacy]") {
    const double epsilon = 1.0e-7;

     // test1
     {
         c3ga::Mvec<double> a;
         a[c3ga::E1] = 1.0;
         a[c3ga::E2] = 1.0;
         c3ga::Mvec<double> b;
         b[c3ga::E1] = 1.0;
         b[c3ga::E3] = 1.0;
         c3ga::Mvec<double> c = a ^b;
         CHECK(c[c3ga::E13] == Approx(1.0).margin(epsilon));
         CHECK(c[c3ga::E12] == Approx(-1.0).margin(epsilon));
         CHECK(c[c3ga::E23] == Approx(1.0).margin(epsilon));
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
         CHECK(c[c3ga::E123] == Approx(1.0).margin(epsilon));
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
        CHECK(c[c3ga::E123] == Approx(10.0).margin(epsilon));
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
        CHECK(c[c3ga::E1] == Approx(2.0).margin(epsilon));
        CHECK(c[c3ga::E23] == Approx(2.0).margin(epsilon));
        CHECK(c[c3ga::E13] == Approx(-3.0).margin(epsilon));
    }
}

TEST_CASE("c3ga legacy: wedgeEqualTest", "[c3ga][legacy]") {
    const double epsilon = 1.0e-7;

    // test1
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 1.0;
        a[c3ga::E2] = 1.0;
        c3ga::Mvec<double> b;
        b[c3ga::E1] = 1.0;
        b[c3ga::E3] = 1.0;
        a ^= b;
        CHECK(a[c3ga::E13] == Approx(1.0).margin(epsilon));
        CHECK(a[c3ga::E12] == Approx(-1.0).margin(epsilon));
        CHECK(a[c3ga::E23] == Approx(1.0).margin(epsilon));
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
        CHECK(a[c3ga::E123] == Approx(1.0).margin(epsilon));
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
        CHECK(a[c3ga::E1] == Approx(2.0).margin(epsilon));
        CHECK(a[c3ga::E23] == Approx(2.0).margin(epsilon));
        CHECK(a[c3ga::E13] == Approx(-3.0).margin(epsilon));
    }

    // test4
    {
        c3ga::Mvec<double> a;
        a[c3ga::scalar] = 2.0;
        a[c3ga::E3] = 3.0;
        a ^= 42.0f; // double
        CHECK(a[c3ga::scalar] == Approx(84.0).margin(epsilon));
        CHECK(a[c3ga::E3] == Approx(126.0).margin(epsilon));
    }

    // test5
    {
        c3ga::Mvec<double> a;
        a[c3ga::scalar] = 2.0;
        a[c3ga::E3] = 3.0;
        a ^= 42; // int
        CHECK(a[c3ga::scalar] == Approx(84.0).margin(epsilon));
        CHECK(a[c3ga::E3] == Approx(126.0).margin(epsilon));
    }
}


/// test the outer product between the primal form of a multivector and the dual of another
TEST_CASE("c3ga legacy: wedgePrimalDualTest", "[c3ga][legacy]") {
    const double epsilon = 1.0e-7;

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
        CHECK(c[c3ga::E0123i] == Approx(8.0).margin(epsilon)); // first, all verifications are now wrong (dual instead of primal) I miss one verification here, I have to put it on the todo list
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
        // a ^ dual(b), with dual(b) = b _| I = -2 e03i - 5 e01i + 6 e13i
        // (the expectations of the former test were wrong: E0123 18, E013i -8, E023i -6)
        CHECK(c.grades() == std::vector<unsigned int>{4});
        CHECK(c[c3ga::E0123] == Approx(0.0).margin(epsilon));
        CHECK(c[c3ga::E012i] == Approx(-15.0).margin(epsilon));
        CHECK(c[c3ga::E013i] == Approx(8.0).margin(epsilon));
        CHECK(c[c3ga::E023i] == Approx(6.0).margin(epsilon));
        CHECK(c[c3ga::E123i] == Approx(-18.0).margin(epsilon));
        CHECK(c == (a ^ b.dual()));
    }
}


TEST_CASE("c3ga legacy: innerProductTest", "[c3ga][legacy]") {

    const double epsilon = 1.0e-7;
    // grade 1
    {
        c3ga::Mvec<double> mv1;
        mv1[c3ga::E1] = 1.2;
        mv1[c3ga::E2] = 2.2;
        c3ga::Mvec<double> mv2;
        mv2[c3ga::E1] = 1;
        mv2[c3ga::E3] = 2;
        c3ga::Mvec<double> mv3 = mv1 | mv2;
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx(1.2).margin(epsilon));
    }

    // grade 1 c3ga with eo and ei
    {
        c3ga::Mvec<double> mv1,mv2,mv3;
        mv1[c3ga::E0] = 2.0;
        mv1[c3ga::E2] = 3.0;
        mv1[c3ga::Ei] = 4.0;
        mv2[c3ga::E0] = 1;
        mv3 = mv1 | mv2;
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx(-4.0).margin(epsilon));
    }

    // grade 3 c3ga
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E012] = 2.0;
        mv1[c3ga::E023] = 3.0;
        mv1[c3ga::E23i] = 4.0;
        mv2[c3ga::E023] = 1;
        mv3 = mv1 | mv2;
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx(4.0).margin(epsilon));
    }

    // different grades
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E02] = 1;
        mv3 = mv1 | mv2;
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx(0).margin(epsilon));

    }

    // scalar
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::scalar] = 42;
        mv3 = mv1 | mv2;
        CHECK(mv3.isEmpty()); // the inner product with a scalar is 0
        mv3 = mv1.dotProduct(mv2);
        CHECK(mv3.grade() == 3);
        CHECK(mv3[c3ga::E123] == Approx(84.0).margin(epsilon));
        CHECK(mv3[c3ga::scalar] == Approx(0).margin(epsilon));
    }

    // scalar
    {
        c3ga::Mvec<double> mv1, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv3 = mv1 | 42;
        CHECK(mv3.isEmpty()); // the inner product with a scalar is 0

    }

    // scalar
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::scalar] = 42;
        mv3 = mv2 | mv1;
        CHECK(mv3.isEmpty()); // the inner product with a scalar is 0
        mv3 = mv2.dotProduct(mv1);
        CHECK(mv3.grade() == 3);
        CHECK(mv3[c3ga::E013] == Approx(126.0).margin(epsilon));
    }

    // scalar
    {
        c3ga::Mvec<double> mv1, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv3 = 42 | mv1;
        CHECK(mv3.isEmpty()); // the inner product with a scalar is 0
    }

    // scalar, vector and bivector
    {
        c3ga::Mvec<double> mv1, mv2;
        mv1[c3ga::E12]  = 2.0;
        mv1[c3ga::E2]   = 3.0;
        mv2[c3ga::E12]  = 4.0;
        mv2[c3ga::E23]  = 6.0;
        mv2[c3ga::E123] = 5.0;
        c3ga::Mvec<double> mv3 = mv1 | mv2; // -8.00 - 12.00*e1 + 8.00*e3 - 15.00*e1^e3
        CHECK(mv3.grade() == 2);
        CHECK(mv3[c3ga::scalar] == Approx(-8.0).margin(epsilon));
        CHECK(mv3[c3ga::E1] == Approx(-12.0).margin(epsilon));
        CHECK(mv3[c3ga::E3] == Approx(8.0).margin(epsilon));
        CHECK(mv3[c3ga::E13] == Approx(-15.0).margin(epsilon));
    }

    {
        c3ga::Mvec<double> mv1,mv2;
        mv1[c3ga::scalar]  = 5.0;
        mv1[c3ga::E1]  = 2.0;
        mv1[c3ga::Ei]  = 3.0;
        mv2[c3ga::E12] = 4.0;
        mv2[c3ga::E1i] = 5.0;
        // the scalar part of mv1 does not contribute to the inner product, but does contribute to the dot product
        c3ga::Mvec<double> mv3 = mv1 | mv2; // 8.00*e2 + 10.00*ni
        CHECK(mv3.grades() == std::vector<unsigned int>{1});
        CHECK(mv3[c3ga::E2] == Approx(8.0).margin(epsilon));
        CHECK(mv3[c3ga::Ei] == Approx(10.0).margin(epsilon));

        mv3 = mv1.dotProduct(mv2); // 8.00*e2 + 10.00*ni + 20.00*e1^e2 + 25.00*e1^ni
        CHECK(mv3.grade() == 2);
        CHECK(mv3[c3ga::E2] == Approx(8.0).margin(epsilon));
        CHECK(mv3[c3ga::Ei] == Approx(10.0).margin(epsilon));
        CHECK(mv3[c3ga::E12] == Approx(20.0).margin(epsilon));
        CHECK(mv3[c3ga::E1i] == Approx(25.0).margin(epsilon));
    }
}

TEST_CASE("c3ga legacy: innerProductEqualTest", "[c3ga][legacy]") {

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
        CHECK(mv1.grade() == 0);
        CHECK(mv1[c3ga::scalar] == Approx(1.2).margin(epsilon));
    }

    // grade 1 c3ga
    {
        c3ga::Mvec<double> mv1,mv2,mv3;
        mv1[c3ga::E0] = 2.0;
        mv1[c3ga::E2] = 3.0;
        mv1[c3ga::Ei] = 4.0;
        mv2[c3ga::E0] = 1;
        mv1 |= mv2;
        CHECK(mv1.grade() == 0);
        CHECK(mv1[c3ga::scalar] == Approx(-4.0).margin(epsilon));
    }


    // grade 3 c3ga
    {
        c3ga::Mvec<double> mv1, mv2;
        mv1[c3ga::E012] = 2.0;
        mv1[c3ga::E023] = 3.0;
        mv1[c3ga::E23i] = 4.0;
        mv2[c3ga::E023] = 1;
        mv1 |= mv2;
        CHECK(mv1.grade() == 0);
        CHECK(mv1[c3ga::scalar] == Approx(4.0).margin(epsilon));
    }


    // different grades
    {
        c3ga::Mvec<double> mv1, mv2;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E02] = 1;
        mv1 |= mv2;
        CHECK(mv1.grade() == 0);
        CHECK(mv1[c3ga::scalar] == Approx(0).margin(epsilon));
    }


    // scalar
    {
        c3ga::Mvec<double> mv1, mv2;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::scalar] = 42;
        mv1 |= mv2;
        CHECK(mv1.isEmpty()); // the inner product with a scalar is 0
    }

    // scalar
    {
        c3ga::Mvec<double> mv1;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv1 |= 42;
        CHECK(mv1.isEmpty()); // the inner product with a scalar is 0
    }

    // scalar
    {
        c3ga::Mvec<double> mv1, mv2;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::scalar] = 42;
        mv1 |= mv2;
        CHECK(mv1.isEmpty()); // the inner product with a scalar is 0

    }
}

TEST_CASE("c3ga legacy: leftContractionScalarTest", "[c3ga][legacy]") {

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

        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx((2.0*42.0)).margin(epsilon));
    }

    // mv < scalar = mv.scalar x scalar
    {
        c3ga::Mvec<double> mv1, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv1[c3ga::scalar] = 5.0;
        mv3 = mv1 < 42;

        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx((5.0*42.0)).margin(epsilon));
    }

    // scalar < mv = mv x scalar
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::scalar] = 2;
        mv2[c3ga::E123] = 2.0;
        mv2[c3ga::E013] = 3.0;
        mv2[c3ga::E023] = 4.0;
        mv3 = mv1 < mv2;

        CHECK(mv3.grade() == 3);
        CHECK(mv3 == (2 ^ mv2));
    }

    // scalar < mv
    {
        c3ga::Mvec<double> mv1, mv3;
        mv1[c3ga::scalar] = 2.0;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv3 = 2 < mv1;

        CHECK(mv3.grade() == 3);
        CHECK(mv3 == (2 ^ mv1));
    }
}

TEST_CASE("c3ga legacy: leftContractionTest", "[c3ga][legacy]") {

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

        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx(1.2).margin(epsilon));
    }

    // grade 1 c3ga
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E0] = 2.0;
        mv1[c3ga::E2] = 3.0;
        mv1[c3ga::Ei] = 4.0;
        mv2[c3ga::E0] = 1;
        mv3 = (mv2 < mv1);
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx(-4.0).margin(epsilon));
    }

    // grade 3 c3ga
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E012] = 2.0;
        mv1[c3ga::E023] = 3.0;
        mv1[c3ga::E23i] = 4.0;
        mv2[c3ga::E023] = 1;
        mv3 = mv2 < mv1;
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx(4.0).margin(epsilon));
    }

    // different grades
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E02] = 1;
        mv3 = mv1 < mv2;
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx(0).margin(epsilon));
    }

    // different grades
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E23] = 1;
        mv3 = mv2 < mv1;
        CHECK(mv3.grade() == 1);
        CHECK(mv3[c3ga::E0] == Approx(-4.0).margin(epsilon));
        CHECK(mv3[c3ga::E1] == Approx(-2.0).margin(epsilon));
    }

    // different grades with "leftContraction"
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E23] = 1;
        mv3 = c3ga::leftContraction(mv2,mv1);
        CHECK(mv3.grade() == 1);
        CHECK(mv3[c3ga::E0] == Approx(-4.0).margin(epsilon));
        CHECK(mv3[c3ga::E1] == Approx(-2.0).margin(epsilon));
    }
}


TEST_CASE("c3ga legacy: rightContractionScalarTest", "[c3ga][legacy]") {

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
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx((42*3.0)).margin(epsilon)); // grade(mv1) > grade(mv2)
    }

    // scalar > mv = mv.scalar x scalar
    {
        c3ga::Mvec<double> mv1, mv3;
        mv1[c3ga::scalar] = 5.0;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv3 = 42 > mv1;
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx((5.0*42)).margin(epsilon)); // grade(mv1) > grade(mv2)
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
        CHECK(mv3.grade() == 3);
        CHECK(mv3 == (2^mv1));
    }

    // mv > scalar = mv x scalar
    {
        c3ga::Mvec<double> mv1, mv3;
        mv1[c3ga::scalar] = 1.0;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv3 = mv1 > 2;
        CHECK(mv3.grade() == 3);
        CHECK(mv3 == (2^mv1));
    }
}


TEST_CASE("c3ga legacy: rightContractionTest", "[c3ga][legacy]") {

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
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx(1.2).margin(epsilon));
    }

    // grade 1 c3ga
    {
        c3ga::Mvec<double> mv1,mv2,mv3;
        mv1[c3ga::E0] = 2.0;
        mv1[c3ga::E2] = 3.0;
        mv1[c3ga::Ei] = 4.0;
        mv2[c3ga::E0] = 1;
        mv3 = (mv2 > mv1);
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx(-4.0).margin(epsilon));
    }

    // grade 3 c3ga
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E012] = 2.0;
        mv1[c3ga::E023] = 3.0;
        mv1[c3ga::E23i] = 4.0;
        mv2[c3ga::E023] = 1;
        mv3 = mv2 > mv1;
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx(4.0).margin(epsilon));
    }

    // different grades
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E02] = 1;
        mv3 = mv2 > mv1;
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx(0).margin(epsilon)); // grade(mv1) > grade(mv2)
    }

    // different grades
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E23] = 1;
        mv3 = mv1 > mv2;
        CHECK(mv3.grade() == 1);
        CHECK(mv3[c3ga::E0] == Approx(-4.0).margin(epsilon));
        CHECK(mv3[c3ga::E1] == Approx(-2.0).margin(epsilon));
    }

    // different grades with "rightContraction"
    {
        c3ga::Mvec<double> mv1, mv2, mv3;
        mv1[c3ga::E123] = 2.0;
        mv1[c3ga::E013] = 3.0;
        mv1[c3ga::E023] = 4.0;
        mv2[c3ga::E23] = 1;
        mv3 = c3ga::rightContraction(mv1, mv2);
        CHECK(mv3.grade() == 1);
        CHECK(mv3[c3ga::E0] == Approx(-4.0).margin(epsilon));
        CHECK(mv3[c3ga::E1] == Approx(-2.0).margin(epsilon));
    }
}

TEST_CASE("c3ga legacy: geometricProductTest", "[c3ga][legacy]") {

    const double epsilon = 1.0e-7;

    // test 1
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        b[c3ga::E2] = 4.0;
        b[c3ga::E3] = 5.0;
        c = a * b; // 12.00 + 15.00*e2^e3 + -10.00*e3^e1 + 8.00*e1^e2
        CHECK(c.grade() == 2);
        CHECK(c[c3ga::scalar] == Approx(12.0).margin(epsilon));
        CHECK(c[c3ga::E12] == Approx(8.0).margin(epsilon));
        CHECK(c[c3ga::E13] == Approx(10.0).margin(epsilon));
        CHECK(c[c3ga::E23] == Approx(15.0).margin(epsilon));
    }

    // test 2
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        b[c3ga::E12] = 4.0;
        b[c3ga::E13] = 5.0;
        c = a * b;  // -12.00*e1 + 8.00*e2 + 10.00*e3 + -15.00*e1^e2^e3
        CHECK(c.grade() == 3);
        CHECK(c[c3ga::E1] == Approx(-12.0).margin(epsilon));
        CHECK(c[c3ga::E2] == Approx(8.0).margin(epsilon));
        CHECK(c[c3ga::E3] == Approx(10.0).margin(epsilon));
        CHECK(c[c3ga::E123] == Approx(-15.0).margin(epsilon));
    }

    // test 3
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        c = a * 2.0f;  // double
        CHECK(c.grade() == 1);
        CHECK(c[c3ga::E1] == Approx(4.0).margin(epsilon));
        CHECK(c[c3ga::E2] == Approx(6.0).margin(epsilon));
    }

    // test 4
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        c = a * 2;  // int
        CHECK(c.grade() == 1);
        CHECK(c[c3ga::E1] == Approx(4.0).margin(epsilon));
        CHECK(c[c3ga::E2] == Approx(6.0).margin(epsilon));
    }

    // test 5
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        c = 2 * a;  // int
        CHECK(c.grade() == 1);
        CHECK(c[c3ga::E1] == Approx(4.0).margin(epsilon));
        CHECK(c[c3ga::E2] == Approx(6.0).margin(epsilon));
    }

    // test 6
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        c = 2 * a;  // int
        CHECK(c.grade() == 1);
        CHECK(c[c3ga::E1] == Approx(4.0).margin(epsilon));
        CHECK(c[c3ga::E2] == Approx(6.0).margin(epsilon));
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
        CHECK(c.grade() == 5);
        CHECK(c[c3ga::scalar] == Approx(-88.0).margin(epsilon));
    }

    // test 8  => very important test (do not restrict to outer and inner)
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1i] = 2.0;
        b[c3ga::E02] = 3.0;
        c = (a * b);  // - 6.00*e1^e2 - 6.00*eo^e1^e2^ei
        CHECK(c.grade() == 4);
        CHECK(c[c3ga::E12] == Approx(-6.0).margin(epsilon));
        CHECK(c[c3ga::E012i] == Approx(-6.0).margin(epsilon));
    }

    // test 9 => very important test (do not restrict to outer and inner)
    {
        c3ga::Mvec<double> a,b,c;
        a[c3ga::E1i] = 2.0;
        b[c3ga::E023] = 3.0;
        b[c3ga::E012] = 3.0;
        c = (a * b);  // - 6.00*e2 - 6.00*e1^e2^e3 + 6.00*e0^e2^ei + 6.00*eo^e1^e2^e3^ei
        CHECK(c.grade() == 5);
        CHECK(c[c3ga::E2] == Approx(-6.0).margin(epsilon));
        CHECK(c[c3ga::E123] == Approx(-6.0).margin(epsilon));
        CHECK(c[c3ga::E02i] == Approx(6.0).margin(epsilon));
        CHECK(c[c3ga::E0123i] == Approx(6.0).margin(epsilon));
    }
}

TEST_CASE("c3ga legacy: geometricProductEqualTest", "[c3ga][legacy]") {

    const double epsilon = 1.0e-7;

    // test 1
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        b[c3ga::E2] = 4.0;
        b[c3ga::E3] = 5.0;
        a *= b; // 12.00 + 15.00*e2^e3 + -10.00*e3^e1 + 8.00*e1^e2
        CHECK(a.grade() == 2);
        CHECK(a[c3ga::scalar] == Approx(12.0).margin(epsilon));
        CHECK(a[c3ga::E12] == Approx(8.0).margin(epsilon));
        CHECK(a[c3ga::E13] == Approx(10.0).margin(epsilon));
        CHECK(a[c3ga::E23] == Approx(15.0).margin(epsilon));
    }

    // test 2
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        b[c3ga::E12] = 4.0;
        b[c3ga::E13] = 5.0;
        a *= b;  // -12.00*e1 + 8.00*e2 + 10.00*e3 + -15.00*e1^e2^e3
        CHECK(a.grade() == 3);
        CHECK(a[c3ga::E1] == Approx(-12.0).margin(epsilon));
        CHECK(a[c3ga::E2] == Approx(8.0).margin(epsilon));
        CHECK(a[c3ga::E3] == Approx(10.0).margin(epsilon));
        CHECK(a[c3ga::E123] == Approx(-15.0).margin(epsilon));
    }

    // test 3
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        a *= 2.0f;  // double
        CHECK(a.grade() == 1);
        CHECK(a[c3ga::E1] == Approx(4.0).margin(epsilon));
        CHECK(a[c3ga::E2] == Approx(6.0).margin(epsilon));
    }

    // test 4
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        a *= 2;  // int
        CHECK(a.grade() == 1);
        CHECK(a[c3ga::E1] == Approx(4.0).margin(epsilon));
        CHECK(a[c3ga::E2] == Approx(6.0).margin(epsilon));
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
        CHECK(a.grade() == 5);
        CHECK(a[c3ga::scalar] == Approx(-88.0).margin(epsilon));
    }
}


TEST_CASE("c3ga legacy: scalarProductTest", "[c3ga][legacy]") {

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
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx(-2.0*4.0).margin(epsilon));
    }

    // test 2
    {
        c3ga::Mvec<double> mv1,mv2;
        mv1[c3ga::E1]  = 2.0;
        mv1[c3ga::Ei]  = 3.0;
        mv2[c3ga::E12] = 4.0;
        mv2[c3ga::E1i] = 5.0;
        c3ga::Mvec<double> mv3 = mv1.scalarProduct(mv2); //
        CHECK(mv3.grade() == 0);
        CHECK(mv3[c3ga::scalar] == Approx(0).margin(epsilon));
    }
}


// formerly hestenesProductTest: the Hestenes product (no contribution of the scalars) is now the inner product operator |,
// and dotProduct() is the inner product including the products with scalars (see commits a186ae2, 1496ed7)
TEST_CASE("c3ga legacy: dotProductTest", "[c3ga][legacy]") {

    const double epsilon = 1.0e-7;

    // test 1
    {
        c3ga::Mvec<double> mv1, mv2;
        mv1[c3ga::E12]  = 2.0;
        mv1[c3ga::E2]   = 3.0;
        mv2[c3ga::E12]  = 4.0;
        mv2[c3ga::E23]  = 6.0;
        mv2[c3ga::E123] = 5.0;
        c3ga::Mvec<double> mv3 = mv1.dotProduct(mv2); // -12.00*e1 + 8.00*e3 - 15.00*e1^e3
        CHECK(mv3.grade() == 2);
        CHECK(mv3[c3ga::E1] == Approx(-12.0).margin(epsilon));
        CHECK(mv3[c3ga::E3] == Approx(8.0).margin(epsilon));
        CHECK(mv3[c3ga::E13] == Approx(-15.0).margin(epsilon));
    }

    // test 2
    {
        c3ga::Mvec<double> mv1,mv2;
        mv1[c3ga::scalar]  = 5.0;
        mv1[c3ga::E1]  = 2.0;
        mv1[c3ga::Ei]  = 3.0;

        mv2[c3ga::E12] = 4.0;
        mv2[c3ga::E1i] = 5.0;
        c3ga::Mvec<double> mv3 = mv1 | mv2; // 8e2 + 10ei

        CHECK(mv3.grade() == 1);
        CHECK(mv3[c3ga::E2] == Approx(8.0).margin(epsilon));
        CHECK(mv3[c3ga::Ei] == Approx(10.0).margin(epsilon));
        CHECK(mv3[c3ga::E12] == Approx(0).margin(epsilon));
        CHECK(mv3[c3ga::E1i] == Approx(0).margin(epsilon));

        mv3 = mv1.dotProduct(mv2); // 8e2 + 10ei + 5 mv2
        CHECK(mv3.grade() == 2);
        CHECK(mv3[c3ga::E2] == Approx(8.0).margin(epsilon));
        CHECK(mv3[c3ga::Ei] == Approx(10.0).margin(epsilon));
        CHECK(mv3[c3ga::E12] == Approx(20.0).margin(epsilon));
        CHECK(mv3[c3ga::E1i] == Approx(25.0).margin(epsilon));
    }
}


TEST_CASE("c3ga legacy: invertTest", "[c3ga][legacy]") {

    const double epsilon = 1.0e-6;
    // test 1
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        b = 1.0 / a;
        CHECK(b.grade() == 1);
        CHECK(b[c3ga::E1] == Approx(0.153846).margin(epsilon));
        CHECK(b[c3ga::E2] == Approx(0.230769).margin(epsilon));
    }

    // test 2
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E13] = 2.0;
        a[c3ga::E23] = 3.0;
        b = 1.0 / a;

        CHECK(b.grade() == 2);
        CHECK(b[c3ga::E13] == Approx(-0.153846).margin(epsilon));
        CHECK(b[c3ga::E23] == Approx(-0.230769).margin(epsilon));
    }

    // test 2 with "inv(mv)"
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E13] = 2.0;
        a[c3ga::E23] = 3.0;
        b = a.inv();

        CHECK(b.grade() == 2);
        CHECK(b[c3ga::E13] == Approx(-0.153846).margin(epsilon));
        CHECK(b[c3ga::E23] == Approx(-0.230769).margin(epsilon));
    }

    // test 3
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E13] = 6.0;
        a[c3ga::E23] = 4.0;
        b = a / 2;
        CHECK(b.grade() == 2);
        CHECK(b[c3ga::E13] == Approx(3).margin(epsilon));
        CHECK(b[c3ga::E23] == Approx(2).margin(epsilon));
    }

    // test 4
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        a[c3ga::E13] = 2.0;
        a[c3ga::E23] = 3.0;
        b = 1.0 / a;
        CHECK(b.grade() == 2);
        CHECK(b[c3ga::E13] == Approx(-0.08).margin(1e-2)); // not sure about the result: gaviewer : inverse() or general_inverse() does not provide the same answer
        CHECK(b[c3ga::E23] == Approx(-0.12).margin(1e-2));
        CHECK(b[c3ga::E1] == Approx(0.08).margin(1e-2));
        CHECK(b[c3ga::E2] == Approx(0.12).margin(1e-2));
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
        CHECK(c.grade() == 1);
        CHECK(c[c3ga::E3] == Approx(-1.0).margin(epsilon));
    }
}

TEST_CASE("c3ga legacy: invertEqualTest", "[c3ga][legacy]") {

    const double epsilon = 1.0e-6;

    // test 1
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 6.0;
        a[c3ga::E2] = 4.0;
        a /= 2.0;
        CHECK(a.grade() == 1);
        CHECK(a[c3ga::E1] == Approx(3).margin(epsilon));
        CHECK(a[c3ga::E2] == Approx(2).margin(epsilon));
    }

    // test 2
    {
        c3ga::Mvec<double> a,b;
        a[c3ga::E1] = 6.0;
        a[c3ga::E2] = 4.0;
        b[c3ga::scalar] = 2;
        a /= b;
        CHECK(a.grade() == 1);
        CHECK(a[c3ga::E1] == Approx(3).margin(epsilon));
        CHECK(a[c3ga::E2] == Approx(2).margin(epsilon));
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
        CHECK(a.grade() == 1);
        CHECK(a[c3ga::E3] == Approx(-1.0).margin(epsilon));
    }

    // test 3
    {
        c3ga::Mvec<double> a;
        a[c3ga::E1] = 2.0;
        a[c3ga::E2] = 3.0;
        a /= a;
        a.roundZero();
        CHECK(a.grade() == 0);
        CHECK(a[c3ga::scalar] == Approx(1.0).margin(epsilon));
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
        CHECK(b.grade() == 2);
        CHECK(b[c3ga::E13] == Approx(-0.08).margin(1e-2));
        CHECK(b[c3ga::E23] == Approx(-0.12).margin(1e-2));
        CHECK(b[c3ga::E1] == Approx(0.08).margin(1e-2));
        CHECK(b[c3ga::E2] == Approx(0.12).margin(1e-2));
    }
}

TEST_CASE("c3ga legacy: reverseTest", "[c3ga][legacy]") {
    const double epsilon = 1.0e-7;

    //  scalar
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::scalar] = 42;
        c3ga::Mvec<double> mv2 = ~mv;
        CHECK(mv2.grade() == 0);
        CHECK(mv2[c3ga::scalar] == Approx(42.0).margin(epsilon));
    }

    // grade 1
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 1;
        mv[c3ga::E2] = 2;
        c3ga::Mvec<double> mv2 = ~mv;
        CHECK(mv2.grade() == 1);
        CHECK(mv2[c3ga::E1] == Approx(1.0).margin(epsilon));
        CHECK(mv2[c3ga::E2] == Approx(2.0).margin(epsilon));
    }

    // grade 2
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E12] = 2;
        mv[c3ga::E23] = 3;
        c3ga::Mvec<double> mv2 = ~mv;
        CHECK(mv2.grade() == 2);
        CHECK(mv2[c3ga::E12] == Approx(-2.0).margin(epsilon));
        CHECK(mv2[c3ga::E23] == Approx(-3.0).margin(epsilon));
    }

    // grade 2 with "reverse()"
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E12] = 2;
        mv[c3ga::E23] = 3;
        c3ga::Mvec<double> mv2 = mv.reverse();
        CHECK(mv2.grade() == 2);
        CHECK(mv2[c3ga::E12] == Approx(-2.0).margin(epsilon));
        CHECK(mv2[c3ga::E23] == Approx(-3.0).margin(epsilon));
    }
}

TEST_CASE("c3ga legacy: dualTest", "[c3ga][legacy]") {
    const double epsilon = 1.0e-7;

    // grade 1
    {
        // dual(A) = A _| I: e1 _| e0123i = -e023i, e2 _| e0123i = e013i
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 2;
        mv[c3ga::E2] = 3;
        c3ga::Mvec<double> mv2 = mv.dual();
        CHECK(mv2.grades() == std::vector<unsigned int>{4});
        CHECK(mv2[c3ga::E023i] == Approx(-2.0).margin(epsilon));
        CHECK(mv2[c3ga::E013i] == Approx(3.0).margin(epsilon));
    }


    // c3ga point
    {
        c3ga::Mvec<float> mv;
        mv[c3ga::E0] = 1;
        mv[c3ga::E1] = 2;
        mv[c3ga::Ei] = 2;
        c3ga::Mvec<float> mv2 = mv.dual();
        CHECK(mv2.grade() == 4);
        CHECK(mv2[c3ga::E0123] == Approx(-1.0).margin(epsilon));
        CHECK(mv2[c3ga::E123i] == Approx(-2.0).margin(epsilon));
        CHECK(mv2[c3ga::E023i] == Approx(-2.0).margin(epsilon));
    }

    // c3ga pair point
    {
        c3ga::Mvec<float> mv;
        mv[c3ga::E01] = -2;
        mv[c3ga::E12] =  4;
        mv[c3ga::E1i] =  4;
        mv[c3ga::E2i] = -4;
        c3ga::Mvec<float> mv2 = mv.dual();

        CHECK(mv2.grade() == 3);
        CHECK(mv2[c3ga::E23i] == Approx(-4.0).margin(epsilon));
        CHECK(mv2[c3ga::E023] == Approx(2.0).margin(epsilon));
        CHECK(mv2[c3ga::E13i] == Approx(-4.0).margin(epsilon));
        CHECK(mv2[c3ga::E03i] == Approx(-4.0).margin(epsilon));
    }
}

TEST_CASE("c3ga legacy: normTest", "[c3ga][legacy]") {
    const double epsilon = 1.0e-7;

    // grade 1
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 2;
        mv[c3ga::E2] = 3;
        CHECK(mv.norm() == Approx(3.605551275).margin(epsilon));
    }

    // grade 2
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E13] = 2;
        mv[c3ga::E23] = 3;
        CHECK(mv.norm() == Approx(3.605551275).margin(epsilon));
    }

    // non homogeneous multivector
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 2;
        mv[c3ga::E2] = 3;
        mv[c3ga::E13] = 2;
        mv[c3ga::E23] = 3;
        CHECK(mv.norm() == Approx(5.099019513).margin(epsilon));
    }
}

TEST_CASE("c3ga legacy: quadraticNormTest", "[c3ga][legacy]") {
    const double epsilon = 1.0e-7;

    // grade 1
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 2;
        mv[c3ga::E2] = 3;
        CHECK(mv.quadraticNorm() == Approx(13).margin(epsilon));
    }

    // grade 2
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E13] = 2;
        mv[c3ga::E23] = 3;
        CHECK(mv.quadraticNorm() == Approx(13).margin(epsilon));
    }

    // non homogeneous multivector
    {
        c3ga::Mvec<double> mv;
        mv[c3ga::E1] = 2;
        mv[c3ga::E2] = 3;
        mv[c3ga::E13] = 2;
        mv[c3ga::E23] = 3;
        CHECK(mv.quadraticNorm() == Approx(26).margin(epsilon));
    }
}

TEST_CASE("c3ga legacy: pseudoscalarTest", "[c3ga][legacy]") {

    c3ga::Mvec<double> mv = c3ga::I<double>();
    std::vector<unsigned int> g = mv.grades();
    CHECK(g.size() == 1);
    CHECK(g[0] == 5);
    CHECK(mv[c3ga::E0123i] == 1);
}

TEST_CASE("c3ga legacy: operatorConstBracket", "[c3ga][legacy]") {
    const double epsilon = 1.0e-7;

    const c3ga::Mvec<double> mvConst = 2.0*c3ga::e0<double>();
    c3ga::Mvec<double> mvOut;
    mvOut[c3ga::E0] = mvConst[c3ga::E0];

    CHECK(mvOut[c3ga::E0] == Approx(mvConst[c3ga::E0]).margin(epsilon));
}

TEST_CASE("c3ga legacy: pseudoscalarInverseTest", "[c3ga][legacy]") {
    c3ga::Mvec<double> mv = c3ga::Iinv<double>();
    std::vector<unsigned int> g = mv.grades();
    CHECK(g.size() == 1);
    CHECK(g[0] == 5);
    CHECK(mv[c3ga::E0123i] == -1);
}

