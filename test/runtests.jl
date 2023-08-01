using Test
using DTWA

@test DTWA.greet("hello") === "hello"

A = DTWA.Interactions(3,3,3)

@test size(A.data)===(3,3,3,3,3,3)

@test DTWA.dims(A) === (3,3,3)

@test DTWA.norm(DTWA.CartesianIndices(A)[2]) === sqrt(6)