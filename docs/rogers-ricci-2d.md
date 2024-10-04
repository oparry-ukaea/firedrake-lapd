# Rogers-Ricci 2D

Model based on the **2D** finite difference implementation described in "*Low-frequency turbulence in a linear magnetized plasma*", B.N. Rogers and P. Ricci, PRL **104**, 225002, 2010 ([link](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.225002)); see equations 5-7.

Script: [2Drogers-ricci.py](../scripts/2Drogers-ricci.py)

## Equations

In SI units:

$$
\begin{aligned}
\frac{d n}{dt} &= -\sigma\frac{n c_s}{R}\exp(\Lambda - e\phi/T_e) + S_n ~~~(1)\\
\frac{d T_e}{dt} &= -\sigma\frac{2}{3}\frac{T_e c_s}{R}\left[1.71\exp(\Lambda - e\phi/T_e)-0.71\right] + S_T ~~~(2)\\
\frac{d \nabla^2\phi}{dt} &= \sigma \frac{c_s m_i \Omega_{ci}^2}{eR}\left[1-\exp(\Lambda - e\phi/T_e)\right] ~~~(3)\\
\end{aligned}
$$

where

$$
\begin{aligned}
\sigma &= \frac{1.5 R}{L_z} \\
\frac{df}{dt} &= \frac{\partial f}{\partial t} - \frac{1}{B}\left[\phi,f\right] \\
\end{aligned}
$$

and the source terms have the form

$$
\begin{aligned}
S_n &= S_{0n}\frac{1-{\rm tanh[(r-r_s)/L_s]}}{2} \\
S_T &= S_{0T}\frac{1-{\rm tanh[(r-r_s)/L_s]}}{2} \\
\end{aligned}
$$

where $r = \sqrt{x^2 + y^2}$

## Parameter Choices

| Parameter     | Value               | Comment                                                                                                         |
| ------------- | ------------------- | --------------------------------------------------------------------------------------------------------------- |
| $T_{e0}$      | 6 eV                |                                                                                                                 |
| $L_z$         | 18 m                |                                                                                                                 |
| $n_0$         | 2e18 m<sup>-3</sup> |                                                                                                                 |
| $m_i$         | 6.67e-27 kg         | Inferred from the value of $c_{s0}$ quoted in the paper. Value is $\sim 4 m_p$, consistent with a Helium plasma |
| $\Omega_{ci}$ | $9.6e5$             |                                                                                                                 |
| $\Lambda$     | 3                   | Couloumb Logarithm                                                                                              |
| R             | 0.5 m               | Approx radius of the plasma column                                                                              |

Derived values
| Parameter   | Calculated as            | Value                               | Comment                                           |
| ----------- | ------------------------ | ----------------------------------- | ------------------------------------------------- |
| B           | $\Omega_{ci} m_i q_E$    | 40 mT                               |                                                   |
| $c_{s0}$    | $\sqrt{T_{e0}/m_i}$      | 1.2e4 ms<sup>-1</sup>               |                                                   |
| $\rho_{s0}$ | $c_{s0}/\Omega{ci}$      | 1.2e-2 m                            | Paper has 1.4e-2 m ... implies $m_i\sim 3 m_p$ !? |
| $S_{0n}$    | 0.03 $n_0 c_{s0}/R$      | 4.8e22 m<sup>-3</sup>s<sup>-1</sup> |                                                   |
| $S_{0T}$    | 0.03 $T_{e0} c_{s0} / R$ | 4318.4 Ks<sup>-1</sup>              |                                                   |
| $\sigma$    | $1.5 R/L_z$              | 1/24                                |                                                   |
| $L_s$       | $0.5\rho_{s0}$           | 6e-3 m                              |                                                   |
| $r_s$       | $20\rho_{s0}$            | 0.24 m                              | Approx radius of the LAPD plasma chamber          |

### Other implementation details

#### Boundary conditions
Homogeneous Dirichlet for the potential on all boundaries; all other BCs are homogeneous Neumann.


#### Domain and mesh

The mesh is a square with the origin at the centre and size $\sqrt{T_{e0}/m_i}/\Omega{ci} = 100\rho_{s0} = 1.2$ m.

By default, there are 64x64 quadrilateral (square) elements, giving sizes of 1.875 cm = 25/16 $\rho_{s0}$

Default res is substantially lower than that used in the finite difference model, which has 1024x1024 elements. (I think)

#### Normalisation

Normalisations follow those in Rogers & Ricci, that is:

|                       | Normalised to   |
| --------------------- | --------------- |
| Charge                | $e$             |
| Electric potential    | $e/T_{e0}$      |
| Energy                | $T_{e0}$        |
| Number densities      | $n_0$           |
| Perpendicular lengths | $100 \rho_{s0}$ |
| Parallel lengths      | $R$             |
| Time                  | $R/c_{S0}$      |


The normalised forms of the equations are:

$$
\begin{align}
\frac{\partial n}{\partial t} &= 40\left[\phi,n\right] -\frac{1}{24}\exp(3 - \phi/T_e)n + S_n  ~~~({\bf 4}) \\
\frac{\partial T_e}{\partial t} &= 40\left[\phi,T_e\right] -\frac{1}{36}\left[1.71\exp(3 - \phi/T_e)-0.71\right]T_e + S_T  ~~~({\bf 5}) \\
\frac{\partial  \nabla^2\phi}{\partial t} &= 40\left[\phi,\nabla^2\phi\right] + \frac{1}{24}\left[1-\exp(3 - \phi/T_e)\right] ~~~({\bf 6})\\
\nabla^2\phi &= \omega ~~({\bf 7}) \\
\end{align}
$$

with 

$$
\begin{equation}
S_n = S_T = 0.03\left\\{1-\tanh[(\rho_{s0}r-r_s)/L_s]\right\\}
\end{equation}
$$

where $\rho_{s0}$, $r_s$ and $Ls$ have the (SI) values listed in the tables above.
This system can be be obtained by applying the normalisation factors, then simplifying; see [here](./details/rogers-ricci-2d-normalised.md) for details. Note that the prime notation used in the derivations is dropped in the equations above for readability.

#### Simulation time

Based on Fig 4. of Rogers & Ricci, the simulation time for the 3D version might be $\sim 12$ in normalised units (= $500~{\rm{\mu}s}$), but it's not clear if the full duration is being shown. 
$500~{\rm{\mu}s}$ doesn't seem enough time for anything interesting to happen - we (arbitrarily) choose to run for ten times longer - $5~{\rm ms}$, or  $\sim 120$ in normalised units.

#### Decoupled density
Note that, since the density only features in equation 4, it is effectively decoupled from the rest of system. Implementing equations 5-7 only is therefore sufficient to capture most of the interesting behaviour.

## Example output

### CG version

<figure>
  <img src="media/rr2D_64x64_CG.gif" width="600">
  <figcaption>Temperature (normalised units) in the CG implementation, run on a 64x64 quad mesh for 5 ms (120 normalised time units).</figcaption>
</figure>

### DG version

<figure>
  <img src="media/rr2D_64x64_DG.gif" width="600">
  <figcaption>Temperature (normalised units) in the DG implementation, run on a 64x64 quad mesh for 5 ms (120 normalised time units).</figcaption>
</figure>

## Weak Forms

### CG version

Equations 4-7 are discretised over a domain $\Omega$ using a continuous Galerkin formulation.
The functions $n$, $T$ and $\omega$, and respective test functions $v_1$, $v_2$, $v_3$ live in separate CG function spaces ($V_1$ - $V_3$) which need not be of the same polynomial order.

The weak form of the equations, below are written using the shorthand $\left< f(n), v_1 \right> \equiv \int_\Omega f(n) v_1 d\mathbf{x}$ to indicate the standard process of multiplying terms by a test function and integrating over the domain. In practice we look for a solution in the combined function space $V=V_1\times V_2\times V_3$ where

ToDo: Add CG weak form.

### DG version

One difference in the default DG setup (see this [config file](../configs/2Drogers-ricciDG_config.yml)) is that the value of Ls is doubled in order to soften the sharp edges of the source profiles. It should be possible to drop this modification when using a sufficiently high resolution mesh. 

ToDo: Add DG weak form.