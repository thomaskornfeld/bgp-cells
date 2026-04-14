**Given:**
$$dq(t) = \gamma \Big( \mu(c(x(t),y(t))) - q(t) \Big) dt + \sqrt{2D}dW_t$$
$$\frac{dg(t)}{dt} = q(t) - \beta g(t)$$

**Solving for $q(t)$:**
$$dq(t) + \gamma q(t) dt = \gamma \mu(c(x(t),y(t))) dt + \sqrt{2D} dW_t$$
Multiplying by $e^{\gamma t}$ and applying Itô's Lemma:
$$d\left( q(t) e^{\gamma t} \right) = \gamma e^{\gamma t} \mu(c(x(t),y(t))) dt + e^{\gamma t} \sqrt{2D} dW_t$$
Integrating over $[0, t]$:
$$q(t) = q(0)e^{-\gamma t} + \gamma \int_0^t e^{-\gamma (t-s)} \mu(c(x(s),y(s))) ds + \sqrt{2D} \int_0^t e^{-\gamma (t-s)} dW_s$$

**Solving for $g(t)$:**
$$\frac{dg(t)}{dt} + \beta g(t) = q(t) \implies \frac{d}{dt}\left(g(t)e^{\beta t}\right) = q(t)e^{\beta t}$$
Integrating over $[0, t]$:
$$g(t) = g(0)e^{-\beta t} + \int_0^t q(s)e^{-\beta (t-s)} ds$$

**Substituting $q(s)$ into $g(t)$:**
Let $\mu(u) = \mu(c(x(u),y(u)))$ and assume $\beta \neq \gamma$.
$$g(t) = g(0)e^{-\beta t} + I_1 + I_2 + I_3$$

Where:
$$I_1 = \int_0^t q(0) e^{-\gamma s} e^{-\beta (t-s)} ds = q(0) \frac{e^{-\gamma t} - e^{-\beta t}}{\beta - \gamma}$$

By Fubini's Theorem:
$$I_2 = \int_0^t \left( \gamma \int_0^s e^{-\gamma (s-u)} \mu(u) du \right) e^{-\beta (t-s)} ds$$
$$I_2 = \gamma \int_0^t \mu(u) e^{\gamma u} e^{-\beta t} \left( \int_u^t e^{(\beta - \gamma) s} ds \right) du = \frac{\gamma}{\beta - \gamma} \int_0^t \mu(u) \left( e^{-\gamma(t-u)} - e^{-\beta(t-u)} \right) du$$

By Stochastic Fubini's Theorem:
$$I_3 = \int_0^t \left( \sqrt{2D} \int_0^s e^{-\gamma (s-u)} dW_u \right) e^{-\beta (t-s)} ds$$
$$I_3 = \frac{\sqrt{2D}}{\beta - \gamma} \int_0^t \left( e^{-\gamma(t-u)} - e^{-\beta(t-u)} \right) dW_u$$

**Evaluating $I_3$ via Itô Isometry:**
$$\mathbb{E}[I_3] = 0$$
$$\text{Var}(I_3) = \frac{2D}{(\beta - \gamma)^2} \int_0^t \left( e^{-\gamma(t-u)} - e^{-\beta(t-u)} \right)^2 du$$
Let $v = t - u$:
$$\text{Var}(I_3) = \frac{2D}{(\beta - \gamma)^2} \int_0^t \left( e^{-2\gamma v} - 2e^{-(\beta+\gamma)v} + e^{-2\beta v} \right) dv$$
$$\text{Var}(I_3) = \frac{2D}{(\beta - \gamma)^2} \left[ \frac{1 - e^{-2\gamma t}}{2\gamma} - \frac{2(1 - e^{-(\beta+\gamma)t})}{\beta + \gamma} + \frac{1 - e^{-2\beta t}}{2\beta} \right]$$

**Therefore:**
$$g(t) = g(0)e^{-\beta t} + q(0) \frac{e^{-\gamma t} - e^{-\beta t}}{\beta - \gamma} + \frac{\gamma}{\beta - \gamma} \int_0^t \mu(u) \left( e^{-\gamma(t-u)} - e^{-\beta(t-u)} \right) du + I_3$$

Thus:
$$I_3 \sim \mathcal{N}\left(0, \frac{2D}{(\beta - \gamma)^2} \left( \frac{1 - e^{-2\gamma t}}{2\gamma} - \frac{2(1 - e^{-(\beta+\gamma)t})}{\beta + \gamma} + \frac{1 - e^{-2\beta t}}{2\beta} \right)\right)$$

Therefore, the gene expression level at a specific time $t$, is
$$g(t) = g(0)e^{-\beta t} + q(0) \frac{e^{-\gamma t} - e^{-\beta t}}{\beta - \gamma} + \frac{\gamma}{\beta - \gamma} \int_0^t \mu(u) \left( e^{-\gamma(t-u)} - e^{-\beta(t-u)} \right) du + I_3\qquad I_3 \sim \mathcal{N}\left(0, \frac{2D}{(\beta - \gamma)^2} \left( \frac{1 - e^{-2\gamma t}}{2\gamma} - \frac{2(1 - e^{-(\beta+\gamma)t})}{\beta + \gamma} + \frac{1 - e^{-2\beta t}}{2\beta} \right)\right)$$

What if our production rates encouraged cells near by to produce similar production rates? I.e.

$$dq_i(t) + \gamma q_i(t) dt = [\gamma \mu(c(\vec{p_i})) + \eta \sum_{j\neq i}{\frac{q_j-q_i}{(\vec{p_j}-\vec{p_i})^2}}]dt + \sqrt{2D} dW_t$$

With our vector of product rates, we can redefine this equation to be
$$\[
L_{ij} =
\begin{cases}
\frac{1}{\lVert \vec{p}_j - \vec{p}_i \rVert^2} & i \neq j \\[6pt]
-\displaystyle \sum_{k \neq i} \frac{1}{\lVert \vec{p}_k - \vec{p}_i \rVert^2} & i = j
\end{cases}
\]$$
where
$$\boxed{d\mathbf{q} = \left[-\gamma \mathbf{q} + \gamma\boldsymbol{\mu} + \eta L\mathbf{q}\right]dt + \sqrt{2D}\, d\mathbf{W}_t}$$
