# Summary

- What represents the `residues` of a transfer function on control systems?

In control systems, the residues of a transfer function provide crucial information about the weight or contribution of each pole to the system's dynamic response. They appear when performing a partial fraction expansion of a transfer function.

The physical interpretation, the `residue` determinates the amplitude of the mode $e^{p_{i}t}$ in the time domain.

If the residue $R_{i}$ islarge, the corresponding pole $p_{i}$​ has a dominant effect on the system dynamics.

If $R_{i}$​ is small, the pole has a negligible effect.

- What is the ODE solution for (homogeneous and non-homogeneous) ? what physical meaning has this?

The homogeneous solution describes the system's natural response (or free response) when no external input ($u(t)= 0$) is applied. It depends only on the initial condition $x(t_{0})$ and the system's dynamics (matrix $A$). This represents how the system evolves over time due to its internal dynamics, such as the decay of stored energy in electrical circuits or the oscillation of a mechanical system.

The non-homogeneous solution describes the system's total response, which includes:

1. The natural response (homogeneous part).
2. The forced response (integral term), which accounts for the effect of the external input $u(t)$. This represents how the system reacts to external forces or control signals, such as an applied voltage in an electrical circuit or a driving force in a mechanical system.