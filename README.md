# The Fermi-Pasta-Ulam-Tsingou Recurrence Problem

## $\alpha$ Lattice

In 1953, Enrico Fermi, John Pasta, Stanislaw Ulam and Mary Tsingou implemented one of the first computer simulations of a physical problem on MANIAC I at Los Alamos. They considered a chain of oscillators with a non-linear nearest-neighbor interaction. Fermi hypothesized that after a substantial number of cycles, the influence of initially exicted normal modes would fade and eventually lead to theramalization. The team let the system run overnight, and to their surprise, the lattice had gone right back to its initial exicted state.

We consider a system of $N$ interacting particles with mass $m$, lattice spacing $\ell$ and linear spring constant $\kappa$. The displacement and momentum of the $\mathrm{n^{th}}$ particle are denoted $q_{n}$ and $p_{n}$ respectively. The Hamiltonian for the $\alpha$-FPUT lattice with nearest-neighbor interaction is:

$$\begin{align} 
  H_{\alpha} = \sum_{j=0}^{N} \frac{1}{2m}p_{j}^{2} + \frac{\kappa}{2}(q_{j+1} - q_{j})^2 + \frac{\alpha\kappa}{3}(q_{j+1} - q_{j})^3
\end{align}$$

Which leads to the canonical equations of motion for the $n^{\mathrm{th}}$ particle:

$$\begin{align}
    \dot{q}_{n} &= \frac{p_{n}}{m}\\
    \dot{p}_{n} &= \kappa\left(q_{n+1} + q_{n-1} - 2 q_{n}\right)
        \Big(1 + \alpha \left(q_{n+1} - q_{n-1}\right)\Big)
\end{align}$$

The frequency of the linearized system $(\alpha = 0)$ is $\displaystyle \omega^2 = \frac{\kappa}{m}$. We introduce dimensionless variables for ($q_n,p_n,t$) as follows:

$$\begin{align}
    t &= \tilde{t}\cdot \frac{1}{\omega}\\
    q_n &= \tilde{q}_n\cdot \ell \\
    p_n &= \tilde{p}_n\cdot m\omega\ell
\end{align}$$

The non-linear coefficient $\alpha$ has dimensions of reciprocal length. In the equations of motion, we see that the non-linearity extends over the region between the left and right nearest neighbors of the particle at $n$. Hence, we define the dimensionless parameter $\displaystyle \varepsilon = 2\alpha\ell$ to characterize the strength of the non-linear interaction between nearest neighbors.

The dimensionless equations of motion are then:

$$\begin{align}
    \frac{d\tilde{q}_n}{d\tilde{t}} &= \tilde{p}_n\\
    \frac{d\tilde{p}_n}{d\tilde{t}} &= \left( \tilde{q}_{n+1} + \tilde{q}_{n-1} - 2 \tilde{q}_{n} \right)
    \Big(1 + \frac{\varepsilon}{2}\left(\tilde{q}_{n+1} - \tilde{q}_{n-1}\right) \Big)
\end{align}$$

We scale the Hamiltonian as $\displaystyle H_{\alpha} = \tilde{H}_{\alpha}\cdot \frac{1}{2} m\omega^2\ell^2$, so the dimensionless Hamiltonian in terms of the dimensionless variables is:

$$\begin{align}
\tilde{H}_{\alpha} = \sum_{j=0}^{N} \tilde{p}_{j}^{2}+(\tilde{q}_{j+1}-\tilde{q}_{j})^2 +\frac{\varepsilon}{3}(\tilde{q}_{j+1}-\tilde{q}_{j})^3
\end{align}$$

For clarity, all dimensionless variables $\tilde{x}$ will be denoted $x$ unless indicated otherwise. Accordingly, the system is characterized by the dimensionless variables as follows:

$$\begin{align}
H_{\alpha} &= \sum_{j = 0}^{N} p_{j}^2 + (q_{j+1} - q_{j})^2 + \frac{\varepsilon}{3}(q_{j+1} - q_{j})^3\\
\dot{q}_n &= p_{n}\\
\dot{p}_n &= (q_{n+1} - 2q_{n} + q_{n-1})\Big(1 + \frac{\varepsilon}{2}\left(q_{n+1} - q_{n-1}\right)\Big)
\end{align}$$

## $\alpha$ String

We now consider the continuum approximation of the FPUT lattice with periodic boundary conditions. 

In terms of the dimension-full variables, the de-coupled equation of motion for the displacement $q_{n}$ is:

$$\begin{align}
    \ddot{q}_{n} = \omega^2 \left(q_{n+1} + q_{n-1} - 2 q_{n}\right)
        \Big(1 + \alpha \left(q_{n+1} - q_{n-1}\right)\Big)\quad \mathrm{where} \quad \omega^2 = \frac{\kappa}{m}
\end{align}$$

To approximate the lattice as a continuum, we take the limit where $\ell\to 0$ and $N \to \infty$ so as to keep the total length of the chain a constant:

$$\lim_{\ell\to 0}\lim_{N\to\infty} \Big(\ell N\Big) = L$$

The lattice vectors $\ell n$ will transform to the continuous position along the lattice, denoted $x$, and the displacements will transform as:

$$\begin{align}
    q_{n}(t) &\to u(x,t)\\
    q_{n\pm 1}(t) &\to u(x\pm \ell, t)
\end{align}$$

We perform $\mathcal{O}(\ell^4)$ Taylor expansions of the nearest neighbor displacements and substitute these into the equation of motion. This procedure yields:

$$\begin{align}
    u_{tt} &= (\omega\ell)^2 \left(u_{xx} + 2\alpha\ell \cdot u_{x}u_{xx} + \frac{\ell^2}{12}u_{xxxx} + \mathcal{O}(\ell^6)\right)
\end{align}$$

Defining the characteristic speed $\displaystyle c = \omega\ell$, and as before, the strength of the nonlinearity $\varepsilon = 2\alpha\ell$, we can write:

$$\begin{align}
    \frac{1}{c^2} u_{tt} = u_{xx} + \varepsilon \left(\frac{\ell}{24\alpha} u_{xxxx} + u_{x}u_{xx}\right)
\end{align}$$

In order to reduce this to an equation that is first order in time, we make the following change of variables:

$$\begin{align}
    z &= x - ct\\
    \tau &= \frac{\varepsilon}{2}ct
\end{align}$$

As $\ell\to 0$, the continuum limit is approached when $\alpha\ell \to 0$ and $\displaystyle \frac{\ell}{\alpha} \to \mathrm{const}$. Thus, we define the characteristic displacement $\displaystyle \delta = \lim_{\ell\to 0}\sqrt{\frac{\ell}{24\alpha}}$ and exclude any terms that are $\mathcal{O}(\varepsilon)$: 

$$\begin{align}
    u_{z\tau} + u_{z}u_{zz} + \delta^{2} u_{zzzz} = 0 
\end{align}$$

Denoting $f = u_{z}$, we have the following equation for $f(z,\tau)$ that is first order in $\tau$:

$$\begin{align}
    \boxed{f_{\tau} + f f_{z} + \delta^{2} f_{zzz} = 0}
\end{align}$$


## Finite Difference Schemes

### FTCS

We will follow the FTCS scheme that is first order accurate in time, and examine schemes that are second and sixth order accurate in z.

$$\begin{align}
    f(z,\tau) 
    &\longrightarrow 
        f^{n}_{m}\\
    f_{\tau}(z,\tau) 
    &\longrightarrow 
        \frac{f^{n+1}_{m} - f^{n}_{m}}{\Delta \tau}\\
    f_z(z,\tau) 
    &\longrightarrow 
        \frac{f^{n}_{m+1} - f^{n}_{n-1}}{2 \Delta z}\\
    f_{zzz}(z,\tau) 
    &\longrightarrow 
        \frac{(f^{n}_{m+2}- f^{n}_{m-2}) - 2 (f^{n}_{m+1} -  f^{n}_{m-1})}{2 \Delta z^{3}}\\
     f_{\tau} = -f f_{z} - \delta^{2} f_{zzz}
    &\longrightarrow 
    f^{n+1}_{m} = f^{n}_{m}
    - \frac{\Delta \tau}{2\Delta z}\left(f^{n}_{m}(f^{n}_{m+1} - f^{n}_{n-1}) 
    + \left(\frac{\delta}{\Delta z}\right)^{2}{(f^{n}_{m+2}- f^{n}_{m-2} - 2 f^{n}_{m+1} + 2f^{n}_{m-1})}\right)
\end{align}$$

### Lax Friedrichs

Considering the equation where $\delta^2 = 0$, we can express the non-linear contribution in the conservation law form that is treated by the Lax-Friedrichs method:

$$\begin{align}
    f_{t} + f f_{z} &= 0 \\ \Big\downarrow\\ f_{t} + \frac{1}{2}\left(f^2\right)_{z} &= 0 \longrightarrow 
    f^{n+1}_{m} 
    = \frac{1}{2}\left(f^{n}_{m+1} - f^{n}_{m-1}\right) 
    - \frac{\Delta \tau}{4\Delta z}\left( (f^{n}_{m+1})^2-(f^{n}_{m-1})^2 \right)
\end{align}$$

We will introduce the linear $\delta^2 f_{zzz}$ term to this equation with a centered finite difference approximation.

### Stability

The numerical stability is determined by the term that is third-order in $z$, so we consider the linearized equation:

$$f_{\tau} = -\delta^2 f_{zzz}$$

Conducting Von-Neumann stability analysis, we make the ansatz $\displaystyle f_{m}^{n}\propto e^{-i\omega t_{n}}e^{ik z_{m}}$. Substitution into the discretized linear equation (2nd order in z) yields:

$$\begin{align}
    e^{-i\omega\Delta t} &= 1 - \frac{\delta^2}{2}\frac{\Delta \tau}{\Delta z^3}
    ( e^{2ik \Delta z} - e^{-2ik \Delta z} - 2 (e^{ik \Delta z} - e^{-ik \Delta z}))\\
    \implies e^{-i\omega\Delta t} &= 1 - i\delta^2\frac{\Delta \tau}{\Delta z^3}
    ( \sin(2k \Delta z) - 2 \sin(k \Delta z))\\
    |e^{-i\omega\Delta t}| &= \sqrt{1 + \left(\delta^2\frac{\Delta \tau}{\Delta z^3}\right)^2
    \underbrace{\left(\sin(2k \Delta z) - 2 \sin(k \Delta z)\right)^2}_{\displaystyle 0\leq\dots\leq9\frac{3}{4}}}
\end{align}$$

Thus without the non-linear term, the KdV equation under a finite-difference scheme is conditionally stable. To minimize the instability from this term, we will choose

$$\begin{align}
\frac{\Delta \tau}{\Delta z^3} \leq \frac{3\sqrt{3}}{2\delta^2}
\end{align}$$
