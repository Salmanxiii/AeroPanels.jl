# Theory

AeroPanels.jl implements potential flow aerodynamic models using the Vortex Lattice Method (VLM) for both steady and unsteady flows. This section outlines the theoretical foundation and the specific implementation details used in this package.

## Steady Aerodynamics

The steady-state solver is based on a standard 3D Vortex Lattice Method. Lifting surfaces are discretized into a grid of panels, each containing a closed vortex ring element.

### Boundary Condition (No Permeability)
The primary boundary condition is the non-penetration (or no permeability) condition, which dictates that there can be no flow through the lifting surface. This is enforced at the collocation points of each panel (typically at 75% chord). The normal velocity induced by all vortex rings must cancel out the normal component of the body's velocity $\vec{v}_b$ (which includes freestream and kinematic motion).

For a system of $N$ panels, this leads to a linear system of equations:

$$[AIC]\{\Gamma_p\} = \{b\}$$

where:
-  $[AIC]$ is the Aerodynamic Influence Coefficient matrix. Each element $AIC_{ij}$ represents the normal velocity induced at panel $i$ by a unit vortex strength at panel $j$.
-  $\{\Gamma_p\}$ are the unknown panel vortex strengths.
-  $\{b\}$ is the normal wash vector, $b_i = -\vec{v}_b \cdot \vec{n}_i$, where $\vec{n}_i$ is the normal vector of panel $i$.

### Kutta Condition
The Kutta condition ensures that the flow leaves the trailing edge smoothly. In AeroPanels.jl, this is satisfied implicitly by including a wake that trails from the trailing edge. In the steady case, the wake is modeled as a set of flat, semi-infinite vortex rings extending downstream. The influence of these wake panels is calculated and added directly to the $[AIC]$ matrix columns corresponding to the trailing edge (TE) panels:

$$AIC_{TE} = AIC_{ring, TE} + AIC_{wake}$$

Because the steady wake circulation equals the TE panel circulation, adding their influences directly enforces the Kutta condition and allows the system to be solved for $\{\Gamma_p\}$.

### Steady Force Computation
The aerodynamic forces on each panel segment are computed using the Kutta-Joukowski theorem. A key difference between some literature approaches and our implementation is that **we explicitly include the induced velocity effect**.

The force $\vec{F}_i$ on a segment $i$ is calculated as:

$$\vec{F}_i = \rho \Gamma_{s,i} (\vec{V}_i + \vec{v}_b) \times \vec{r}_i$$

where:
-  $\rho$ is the air density.
-  $\Gamma_{s,i}$ is the circulation of segment $i$.
-  $\vec{r}_i$ is the segment vector.
-  $\vec{v}_b$ is the body velocity.
-  $\vec{V}_i$ is the **induced velocity** at the segment, computed explicitly from the influence of all other rings and the wake:

$\vec{V}_i = [AIC3_{ring}]\{\Gamma_p\} + [AIC3_{wake}]\{\Gamma_w\}$.

---

## Unsteady Aerodynamics

The unsteady solver implements a **Continuous-Time Unsteady Vortex Lattice Method (UVLM)**. Unlike traditional discrete-time models that use a fixed time step, this approach formulates the aerodynamics as a continuous-time state-space system, which is well-suited for integration with standard ODE solvers and aeroelastic analysis.

### Governing Equations
The unsteady model is derived from the transport of vorticity in the wake. Based on the formulations in continuous UVLM literature (e.g., Binder 2017, Werter et al. 2018), the aerodynamic state can be represented by the wake circulation states $\{\Gamma_w\}$, which evolve according to the transport equation:

$$\dot{\{\Gamma_w\}} = [K_8]\{\Gamma_w\} + [K_9]\{b\}$$

where:
-  $\{\Gamma_w\}$ are the wake circulation states (the state vector of the system).
-  $\{b\}$ is the normal wash vector (the input vector).
-  $[K_8]$ represents the transport (convection) of vorticity through the wake.
-  $[K_9]$ represents the shedding of new vorticity from the trailing edge due to changes in the boundary condition.

The body circulations $\{\Gamma_b\}$ are then recovered from the wake states and the boundary condition using steady-like relations mapping matrix operators:

$$\{\Gamma_b\} = [L_3]\{\Gamma_w\} - [L_4]\{b\}$$

### Unsteady Force Computation
The total aerodynamic forces in the unsteady model are calculated using the unsteady Bernoulli equation, decomposed into two main components:

1. **Quasi-steady forces (Circulatory)**: These are computed identically to the steady case, using the Kutta-Joukowski theorem on the segments, but utilizing the instantaneous body and wake circulations ($\{\Gamma_b\}$ and $\{\Gamma_w\}$). Crucially, **the induced velocity effect is included** in this calculation as well:

   $$\vec{F}_{qs,i} = \rho \Gamma_{s,i} (\vec{V}_i + \vec{v}_b) \times \vec{r}_i$$

2. **Unsteady forces (Added Mass / Non-circulatory)**: These arise from the time rate of change of the potential field (or circulation) over the panels. For a panel $i$ with area $A_i$ and normal $\vec{n}_i$, the unsteady force is:

   $$\vec{F}_{u,i} = \rho \dot{\Gamma}_{b,i} A_i \vec{n}_i$$

The derivative of the body circulation $\dot{\Gamma}_b$ is computed analytically from the state-space matrices:

   $$\dot{\{\Gamma_b\}} = [L_5]\{\Gamma_w\} + [L_6]\{b\} - [L_4]\{\dot{b}\}$$

where $\{\dot{b}\}$ accounts for the fluid's acceleration or pitch-rate changes.

By combining the quasi-steady and unsteady force components, the complete dynamic loads on the lifting surfaces are obtained.

---

## Verification

The methods implemented in this package have been verified against classical analytical solutions.

### Wagner Problem (Sudden Acceleration)
The unsteady lift buildup on a flat plate following an impulsive start is compared against the Wagner function.

![Sudden Acceleration Verification](assets/SuddenAcceleration.png)

### Steady Sweep and Drag
The steady solver has been validated against data from Plotkin and Dimitriadis for swept wings and induced drag predictions.

![Drag Verification](assets/AeroPanels.png)
