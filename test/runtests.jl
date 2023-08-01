using Test
using DTWA

@test DTWA.greet("hello") === "hello"

@test size(DTWA.Interactions(3,3,3).data)===(3,3,3,3,3,3)

A = DTWA.Interactions(3,3,3)

@test DTWA.dims(A) === (3,3,3)