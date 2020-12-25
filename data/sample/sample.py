from project_namespace_py import *


mv1 = Mvec();
mv1[scalar] = 1.0;
mv1[Eproject_first_vector_basis] = 42.0;
print("mv1 : ", mv1 )

mv2 = Mvec()
mv2[Eproject_first_vector_basis] = 1.0;
mv2[Eproject_second_vector_basis] = 2.0;
mv2 += I() + 2*eproject_first_vector_basisproject_second_vector_basis();
print("mv2 : " , mv2 );

# some products
print("outer product     : ", (mv1 ^ mv2) )
print("inner product     : ", (mv1 | mv2) )
print("geometric product     : ", (mv1 * mv2) )


# some tools
print("grade : ", mv1.grade())
print("norm : ", mv1.norm())
print("grade of mv2 : ", mv2.grade())


