## Normalisations

(Normalised quantities indicated by $'$)

$$
\begin{align}
n      &= n_0 n'                     \\
\omega &= n_0\omega'       \\
T_e    &= T_0 T_e'                   \\
\phi   &= \frac{T_0}{e}\phi'         \\
t      &= \frac{R}{c_{s}}t'          \\
\nabla &= \frac{1}{\rho_{s0}}\nabla' \\
[a,b] &= \left(\frac{1}{\rho_{s0}}\right)^2[a,b]' \\
L_\perp &= \rho_{s0}~L_\perp' \\
L_\parallel &= R~L_\parallel' \\
\end{align}
$$

## Normalised equations

The derivations below involve starting from the SI version of each equation, substituting in the normalised variables, then reducing the result to as simple a form as possible.
A normalised form of the Poisson brackets is derived in the final subsection, so those terms are ignored in the derivations for each time-evolved equation.

### Density

SI equation: $\frac{d n}{dt} = -\sigma\frac{n c_s}{R}\exp(\Lambda - e\phi/T_e) + S_n$

---
Let
$$
\begin{align}
S_n' &= S_n \frac{R}{n_0 c_s}\\
&= 0.03\left\{1-\tanh[(r-r_s)/L_s]\right\} \\
&= 0.03\left\{1-\tanh[(\rho_{s0}r'-r_s)/L_s]\right\}
\end{align}
$$

and using
$\sigma=1.5R/L_z$,

$$
\begin{align}
\frac{d n'}{dt'}n_0\frac{c_s}{R} &= -\sigma\frac{c_s}{R}\exp[\Lambda - e(\phi'T_0/e)/(T_0 T_e')]n'n_0 + S_n  \\
\frac{d n'}{dt'}n_0\frac{c_s}{R}  &= -1.5\frac{R}{L_z}\frac{c_s}{R}\exp(\Lambda - \phi'/T_e')n'n_0  + S_n\\
\frac{d n'}{dt'} &= -1.5\frac{R}{L_z}\exp(\Lambda - \phi'/T_e')n' + S_n' \\
&= -\frac{1}{24}\exp(3 - \phi'/T_e')n' + S_n'  ~~({\bf 1})
\end{align}
$$

### Temperature

SI equation: $\frac{d T_e}{dt} = -\sigma\frac{2}{3}\frac{T_e c_s}{R}\left[1.71\exp(\Lambda - e\phi/T_e)-0.71\right] + S_T$

---

Let
$$
\begin{align}
S_T' &= S_T \frac{R}{T_0 c_s}\\
&= 0.03\left\{1-\tanh[(r-r_s)/L_s]\right\} \\
&= 0.03\left\{1-\tanh[(\rho_{s0}r'-r_s)/L_s]\right\}
\end{align}
$$

and using $\sigma=1.5R/L_z$,

$$
\begin{align}
\frac{d T_e'}{dt'}\frac{T_0 c_s}{R} &= -\sigma\frac{2}{3}\frac{T_0 T_e'c_s}{R}\left[1.71\exp[\Lambda - e(\phi'T_0/e)/(T_0 T_e')]-0.71\right] + S_T \\
\frac{d T_e'}{dt'}\frac{T_0 c_s}{R} &= -\frac{R}{L_z}\frac{T_0 T_e'c_s}{R}\left[1.71\exp(\Lambda - \phi'/T_e')-0.71\right] + S_T \\
\frac{d T_e'}{dt'} &= -\frac{R}{L_z}T_e'\left[1.71\exp(\Lambda - \phi'/T_e')-0.71\right] + S_T'\\
&= -\frac{1}{36}\left[1.71\exp(3 - \phi'/T_e')-0.71\right]T_e' + S_T'  ~~({\bf 2})
\end{align}
$$

### Vorticity

SI equation: $\frac{d \nabla^2\phi}{dt} = \sigma \frac{c_s m_i \Omega_{ci}^2}{eR}\left[1-\exp(\Lambda - e\phi/T_e)\right]$

---

$$
\begin{align}
\frac{d \nabla'^2\phi'}{dt'}\left(\frac{1}{\rho_{s0}}\right)^2\frac{T_0}{e}\frac{c_s}{R} &= \sigma \frac{c_s m_i \Omega_{ci}^2}{eR}\left[1-\exp[\Lambda - e(\phi'T_0/e)\right] \\
\frac{d \nabla'^2\phi'}{dt'}\frac{T_0 c_s}{\rho_{s0}^2 e R} &= 1.5\frac{R}{L_z} \frac{c_s m_i \Omega_{ci}^2}{eR}\left[1-\exp(\Lambda - \phi'/T_e')\right] \\
\frac{d \nabla'^2\phi'}{dt'} &= 1.5\frac{R}{L_z} \frac{\rho_{s0}^2 m_i \Omega_{ci}^2}{T_0}\left[1-\exp(\Lambda - \phi'/T_e')\right] \\
\frac{d \nabla'^2\phi'}{dt'} &= 1.5\frac{R}{L_z} \left[1-\exp(\Lambda - \phi'/T_e')\right]\\
&= \frac{1}{24}\left[1-\exp(3 - \phi'/T_e')\right]~~({\bf 3})\\
\end{align}
$$

Where  $c_s^2=T_0/m_i=\rho_{s0}^2\Omega_{ci}^2$ was used to obtain the penultimate line

### Potential

SI equation: $\nabla^2\phi = \frac{eB^2}{m_i \bar{n}}\omega$

---

$$
\begin{align}
\nabla'^2\phi'\left(\frac{1}{\rho_{s0}}\right)^2\frac{T_0}{e} &= \frac{eB^2}{m_i \bar{n}}n_0\omega' \\
\nabla'^2\phi'\left(\frac{1}{\rho_{s0}}\right)^2\frac{T_0}{e} &= \frac{\Omega_{ci}^2 m_i}{e\bar{n}}n_0\omega \\
\nabla'^2\phi' &= \frac{\rho_{s0}^2\Omega_{ci}^2 m_i}{T_0\bar{n}}n_0\omega \\
\nabla'^2\phi' &= \frac{n_0}{\bar{n}}\omega \\
\nabla'^2\phi' &= \omega
~~({\bf 4}) \\
\end{align}
$$
Where  $c_s^2=T_0/m_i=\rho_{s0}^2\Omega_{ci}^2$ was used to obtain the penultimate line and the last line follows from the choice $\bar{n}=n_0$.
### Poisson Brackets

SI: $\frac{df}{dt} = \frac{\partial f}{\partial t} - \frac{1}{B}\left[\phi,f\right]$

---
$$
\begin{align}
\frac{df'}{dt'}\frac{c_s}{R} &= \frac{\partial f'}{\partial t'}\frac{c_s}{R} - \frac{1}{B}\frac{T_0}{e}\left(\frac{1}{\rho_{s0}}\right)^2\left[\phi',f'\right]' \\
\frac{df'}{dt'} &= \frac{\partial f'}{\partial t'} - \frac{T_0 R}{B \rho_{s0}^2 e c_s}\left[\phi',f'\right]' \\
\frac{df'}{dt'} &= \frac{\partial f'}{\partial t'} - \frac{R}{\rho_{s0}}\left[\phi',f'\right]' \\
\frac{df'}{dt'} &= \frac{\partial f'}{\partial t'} - 40\left[\phi',f'\right]' 
~~({\bf 5}) \\
\end{align}
$$
Where
$\rho_{s0}=\frac{c_s}{\Omega_{ci}}$, $c_s^2=\frac{T_0}{m_i}$ and $\Omega_{ci}=\frac{eB}{m_i}$ ($=> \frac{T_0}{B\rho_{s0}e c_s}=1$) were used to obtain the penultimate line  from the one above.
and all factors required to transform $f\mapsto f'$ on the first line cancel.
