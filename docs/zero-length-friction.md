# Coulomb friction in ZeroLength elements

Friction coupling is implemented in planar `ZeroLength` and spatial
`ZeroLength3` elements and exposed through both input-table profiles.
It couples shear response to another DOF's current normal force on the same
element; it is an element feature rather than a uniaxial material variant.

Compression is negative. Friction capacity is `mu * max(-normalForce, 0)`.
The shear response sticks with stiffness `k0`, then slides with post-slip
tangent `b * k0`. Committed slip persists between steps. The tangent includes
the derivative of shear resistance with respect to the normal response.

## JavaScript input

`zeroLengths.friction` entries are tuples:

```typescript
[elementRow, normalDof, shearDof, mu, k0, b]
```

`zeroLengths3.friction` entries are:

```typescript
[elementRow, normalDof, shearDof1, shearDof2, mu, k0, b]
```

Rows refer to the corresponding zero-length table. DOFs use zero-based global
indices (0–2 planar, 0–5 spatial). The normal DOF needs its own material in
the table's `materials` entries. `b` is explicit; a nonzero post-slip slope
can preserve a stiffness path while sliding.

Spatial shear directions slide independently against the shared normal
force, giving a square interaction surface. A circular biaxial friction
cone and velocity-dependent friction are not implemented. Friction follows an
element's `orient` frame like any other DOF, so its normal and shear
directions are the local axes.

Core tests in `core/src/model/elements/zero_length.rs` cover sticking,
sliding, normal-force dependence, tension, permanent slip, and the coupled
tangent. The planar and spatial `m10_*input_v1.rs` bridge tests cover
friction input decoding and analysis.

The [archived friction design](obsolete/zero-length-friction-design.md)
preserves the original rationale and derivations. Its planar and independent
spatial phases are implemented; the circular-cone extension remains a proposal.
