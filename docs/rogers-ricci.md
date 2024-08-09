# Rogers-Ricci

Model based on the finite difference implementation described in "*Low-frequency turbulence in a linear magnetized plasma*", B.N. Rogers and P. Ricci, PRL **104**, 225002, 2010 ([link](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.225002)).

## Equations

(Electrostatic Braginskii equations with $T_i \ll T_e$ and $\beta \ll 1$)
$$
\begin{aligned}
\dot{n} &= -\nabla_{\parallel}(n u_e) + S_n \\
\dot{u_i} &= -u_i\nabla_{\parallel}(u_i) - \frac{1}{n} \nabla_{\parallel}(p_e) \\
m_e\dot{u_e} &= -m_e u_e\nabla_{\parallel}(u_e) - \frac{T_e}{n}\nabla_{\parallel}(n) + e\nabla_{\parallel}(\phi) - 1.71\nabla_{\parallel}(T_e) + \frac{e j_\parallel}{\sigma_\parallel} \\
\dot{T_e} &= \frac{2}{3}\frac{T_e}{e n}0.71\nabla_\parallel j_\parallel - \frac{2}{3}T_e\nabla_{\parallel}(u_e) - u_e\nabla_{\parallel}(T_e) + S_T\\
\dot{\omega} &= -u_i\nabla_\parallel \omega + \frac{m_i \Omega_{ci}^2}{e^2 n}\nabla_\parallel j_\parallel\\
\nabla_{\perp}^2 \phi &= \omega
\end{aligned}
$$
where
$$
\begin{aligned}
u_x &\equiv V_{\parallel x} \\
p_e &= nT_e \\
j_\parallel &= e n (u_i-u_e) \\
\sigma_\parallel &= \frac{e^2 n_0 R}{m_{i}c_{s0}/\nu} \\
\Omega_{ci} &= e B / m_i c
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

Choices:

| Parameter     | Value       | Comment                                                                                                          |
| ------------- | ----------- | ---------------------------------------------------------------------------------------------------------------- |
| $T_{e0}$      | 6 eV        |                                                                                                                  |
| $L_z$         | 18 m        |                                                                                                                  |
| $n_0$         | 2e18 m^-3   |                                                                                                                  |
| $\nu$         | 0.03        |                                                                                                                  |
| $m_i$         | 6.67e-27 kg | N.B. This was inferred from the value of $c_{s0}$ quoted in the paper. Presumably they actually set $m_i=4 m_p$? |
|               |             | Later on they claim $m_i= 400 m_e$. We'll ignore this as it doesn't seem to fit anything else...                 |
| $\Omega_{ci}$ | $9.6e5$     | (Implies $B/m_i\sim1$)                                                                                           |
| $\Lambda$     | 3           | Couloumb Logarithm?                                                                                              |
| R             | 0.5 m       | Approx radius of the plasma column?                                                                              |

Derived values
| Parameter          | Calculated as                   | Value                 | Comment                                                 |
| ------------------ | ------------------------------- | --------------------- | ------------------------------------------------------- |
| $c_{s0}$           | $\sqrt{T_{e0}/m_i}$             | 1.2e4 m$s^{-1}$       | (If $m_i=400 m_e$, this would be $\sim5.1$e4 m$s^{-1}$) |
| $\rho_{s0}$        | $c_{s0}/\Omega{ci}$             | 1.2e-2 m              | Paper has 1.4e-2 m ... implies $m_i\sim 3 m_p$ !?       |
| $S_{0n}$           | 0.03 $n_0 c_{s0}/R$             | 4.8e22 $m^{-3}s^{-1}$ |                                                         |
| $S_{0T}$           | 0.03 $T_{e0} c_{s0} / R$        | 4318.4 K$s^{-1}$      |                                                         |
| $\sigma_\parallel$ | $e^2 n_0 R / (0.03 m_i c_{s0})$ | 10676.0               |                                                         |

## Weak Form

To be added

## Boundary Conditions
 Bohm BCs at the end walls ($z = \pm L_z/2$): $u_i= \pm c_s$, $u_e=\pm exp(\Lambda - e\phi/T_e)$


## CG version

Script: [rogers-ricci_CG.py](../scripts/rogers-ricci_CG.py)

<!-- ## DG version

Script: [rogers-ricci_DG.py](../scripts/rogers-ricci_DG.py) -->