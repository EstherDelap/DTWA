using Test
using DTWA

@test DTWA.greet("hello") === "hello"

A = DTWA.Interactions(3,3,3)
#testing that the struct Interactions works
@test size(A.data)===(3,3,3,3,3,3)

@test DTWA.dims(A) === (3,3,3)

@test DTWA.norm(DTWA.CartesianIndices(A)[2]) === sqrt(6)

@test DTWA.Jz_Ising((3,3,3),3)[DTWA.CartesianIndices(A)[2],DTWA.CartesianIndices(A)[2]] === 0.0
#testing if the Jz_ising func works
@test DTWA.getindex(DTWA.Jz_Ising((3,3,3),2),[2,2,1],[2,2,1]) === 0.0
#testing if the J_XYZ func works, specifically the periodic boundary conditions
@test DTWA.getindex(DTWA.J_XYZ((3,3,3),1.0),[3,2,2],[1,2,2]) === 1.0 

@test DTWA.spin_array_3D([3,3,3], 1,1)[1,1,1][1] === 1.0

B = DTWA.euler_3D(5, (0.0,1.0), DTWA.spin_array_3D([3,3,3], 1,1),0, 0, 0,DTWA.Jx_Ising([3,3,3]),DTWA.Jx_Ising([3,3,3]), DTWA.Jx_Ising([3,3,3]))

@test size(B)[1] == 5

@test B[1] == B[2]

