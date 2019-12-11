# note about code

## WeightResidual_female.py

- <https://machinelearningmastery.com/model-residual-errors-correct-time-series-forecasts-python/>

- dolfin
  - <https://pypi.org/project/DOLFIN/>
- FEniCS
  - <https://fenicsproject.org/fenics-project-roadmap-2019/>
- PyMor DEMO heat
  - $$ \begin{aligned}
        \partial_t z(x, y, t) &= \Delta z(x, y, t),      & 0 < x, y < 1,\ t > 0 \\
        -\nabla z(0, y, t) \cdot n &= z(0, y, t) - u(t), & 0 < y < 1, t > 0 \\
        -\nabla z(1, y, t) \cdot n &= z(1, y, t),        & 0 < y < 1, t > 0 \\
        -\nabla z(0, x, t) \cdot n &= z(0, x, t),        & 0 < x < 1, t > 0 \\
        -\nabla z(1, x, t) \cdot n &= z(1, x, t),        & 0 < x < 1, t > 0 \\
        z(x, y, 0) &= 0                                  & 0 < x, y < 1 \\
        y(t) &= \int_0^1 z(1, y, t) dy,                  & t > 0 \\
        \end{aligned}
    $$
