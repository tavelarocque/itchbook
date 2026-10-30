# Phase 3a: the tent oracle GT and why it should give a V

GT is a CHOICE, a power device for later phases, not a proposed mechanism of real markets. It gives the
price impact of a run a decay exponent that depends on the run's size, with an interior maximum. This note
writes out, line by line, the heuristic that such a kernel should produce a zeta_K(a) with an interior
minimum (a "V"). The heuristic is a CONJECTURE. The last section checks it against the simulated scan.

## Setup

Run t has sign eps_t and size q_t.
Sizes are iid with log q ~ N(mu, sigma^2) and independent of the signs (GT's default size model).
Signs are stationary with C(l) = E[eps_t eps_{t+l}].

The mid read just before run n is

    m_n = theta * sum_{k<n} G_{b(k)}(n - k) * eps_k * q_k^psi,

where b(k) is the log-size bin of run k and G_b(l) = (1 + l)^(-beta_b) for l >= 1, G_b(0) = 0.
So run k moves every later mid by an amount that decays at its own bin's rate beta_b.

The tent is

    beta(log q) = beta_max - c * |log q - l_star|,

clipped below at 0 and held constant within each of 12 bins.
Decay is fastest at log q = l_star and slower on both sides.

A window of T runs starting at n has

    Delta_T(n) = m_{n+T} - m_n,        I^a_T(n) = sum_{t=n}^{n+T-1} eps_t * q_t^a,

and the statistic is K_a(T) = E[Delta_T I^a_T].
Its exponent zeta_K(a) is the log-log slope of K_a(T) in T. The normalisation k_a(T) = K_a(T) / (T K_a(1)) in
moments.py divides by a constant in T, so it shifts the slope by exactly -1 and nothing else.

## Step 1: split K into a diagonal and a cross part

Insert the mid into Delta_T:

    Delta_T(n) = theta * sum_k [G_{b(k)}(n + T - k) - G_{b(k)}(n - k)] * eps_k * q_k^psi.

Each run k contributes the change of its own impact over the window.

Multiply by I^a_T and take expectations:

    K_a(T) = theta * sum_k sum_{t=n}^{n+T-1} E[ (G_{b(k)}(n+T-k) - G_{b(k)}(n-k)) eps_k q_k^psi eps_t q_t^a ].

Each term pairs an impacting run k with a flow run t in the window.

The diagonal part D_a(T) collects the terms with k = t.
The cross part X_a(T) collects the terms with k != t, including impacting runs k before the window.

## Step 2: the diagonal part

For k = t inside the window, n - t <= 0, so G(n - t) = 0, and eps_t^2 = 1:

    D_a(T) = theta * sum_{t=n}^{n+T-1} E[ G_{b(q)}(n + T - t) * q^(psi + a) ].

Each window run contributes its own impact, observed at the end of the window.

Re-index with l = n + T - t, which runs over 1..T:

    D_a(T) = theta * E[ q^(psi + a) * S_{b(q)}(T) ],    S_b(T) = sum_{l=1}^{T} (1 + l)^(-beta_b).

The sum S_b(T) grows like T^(1 - beta_b) / (1 - beta_b) for beta_b < 1 and large T.

## Step 3: the tilt

For a Gaussian log q, weighting by q^s is the same as shifting the mean:

    E[ q^s f(log q) ] = E[q^s] * E[ f(L_s) ],    L_s ~ N(mu + s sigma^2, sigma^2).

This is exponential tilting of a normal density: the weight e^(s x) multiplies the density, and completing
the square moves the mean by s sigma^2 and leaves the variance unchanged.

With s = psi + a:

    D_a(T) = theta * E[q^(psi + a)] * E[ S_{b(L_{psi+a})}(T) ].

So the diagonal part sees the kernel of runs whose log size sits near mu + (psi + a) sigma^2.

Approximate the tilted average by the value at the tilted mean:

    D_a(T) ~ T^(1 - beta_eff(a)),    beta_eff(a) ~ beta(mu + (psi + a) sigma^2).

Raising a by one moves the diagonal contribution up by sigma^2 in log q (about 2.2 at Q-T3S1, sigma_l =
1.48). This is APPROXIMATE (see the caveats below).

## Step 4: the cross part is not tilted

For k != t the sizes q_k and q_t are independent of each other and of the signs, so the expectation factors:

    E[ dG_{b(k)} eps_k q_k^psi eps_t q_t^a ] = C(t - k) * E[ q^psi dG_{b(q)} ] * E[q^a].

The flow size enters only through E[q^a], which does not depend on T.

Hence

    X_a(T) = E[q^a] * X(T),

with X(T) the same function of T for every a. The cross part has an a-independent exponent zeta_X, which
under sign memory with C(l) ~ l^(-gamma) is about 2 - gamma - beta_X. Here beta_X is a q^psi-weighted
decay exponent.

## Step 5: the weights shift toward the diagonal as a grows

The ratio of the two parts' prefactors is

    E[q^(psi + a)] / (E[q^psi] E[q^a]) = exp(a psi sigma^2).

It grows with a, so large a puts more weight on the diagonal.

zeta_K(a) is therefore a weighted mixture of two exponents:

    zeta_K(a) ~ w(a) * (1 - beta_eff(a)) + (1 - w(a)) * zeta_X,    w(a) increasing in a.

The weight moves from the cross exponent toward the diagonal exponent as a grows.

## Step 6: G0 and GT

G0 has beta_b = beta for every bin, so 1 - beta_eff(a) = 1 - beta for every a.
Since 1 - beta < 2 - gamma - beta whenever gamma < 1, the only a-dependence is the weight shift, and
zeta_K(a) falls monotonically (the Lemma 1 world: no V).

GT has beta_eff(a) rising while mu + (psi + a) sigma^2 < l_star and falling after:

    a* ~ (l_star - mu) / sigma^2 - psi.

The diagonal exponent 1 - beta_eff(a) has its minimum at a* and rises beyond it.

Beyond a*, two effects compete. The weight shift still pulls zeta_K(a) down, and the rising diagonal
exponent pushes it up. The CONJECTURE: if c is large enough that the diagonal exponent rises faster than the
weight shift lowers zeta_K, then zeta_K(a) has an interior minimum. Because the weight shift keeps pulling
down after a*, that minimum sits at or after a*.

## Where the argument is only approximate

1. Cross terms are not tilted. Step 4 holds exactly only for iid sizes independent of signs. With the size
   copula (gen3 copula=True), sizes of nearby runs are correlated, and cross terms between nearby runs are
   partly tilted too.
2. beta_eff(a) at the tilted mean ignores the tilted distribution's width sigma. The average of
   T^(1 - beta(L)) over L ~ N(mu + s sigma^2, sigma^2) is pulled toward the smallest beta within about a
   sigma of the mean, more so at large T. Near l_star this flattens the minimum and shifts it.
3. The exponents are local. S_b(T) is not yet a power law at T = 8, and the OLS slopes over [8, 64],
   [16, 256], [64, 1024] mix diagonal and cross parts whose weights change with T.
4. The tilt saturates in a finite day. A day of N runs has largest log q near mu + sigma sqrt(2 ln N) (about
   mu + 5.3 sigma at N = 10^6). Once (psi + a) sigma exceeds that, the empirical tilted mean is the
   sample maximum, not mu + (psi + a) sigma^2. At Q-T3S1 this starts near a = 3. The same few runs then
   carry V_a and K_a, which is why the standard errors explode at large a.
5. Binning and clipping. beta is constant over bins of 0.67 sigma, and clipped at 0 where the tent would go
   negative (c = 0.4 clips most bins, so the tent becomes a plateau).
6. The cross part also contains the reversal from runs before the window (G(n+T-k) - G(n-k) < 0 for k < n).
   Its share changes with T, not with a.

## Check against the simulation (Q-T3S1, 16 days of 10^6 runs per setting)

Grid: beta_max in {0.3, 0.5, 0.7}, c in {0.1, 0.2, 0.4}, l_star = mu + {1, 2, 3} sigma_l, psi in {0.25, 0.5}.
Depth = min_a zeta_K(a) - min(zeta_K(0), zeta_K(3)), a V when depth < -2 SE (delete-one-day jackknife).

- First version (bins spanning mu +- 3 sigma_l, top bin open): NO setting gave a V. Every zeta_K(a) fell
  with a. The top bin's beta was that of its centre near mu + 2.75 sigma, while the q^(psi+a) weight of a = 3
  sits near mu + 5 sigma. The tent therefore never fell on the side where the high-a weight lives, so this
  version did not test the mechanism. The bins now span [mu - 3 sigma, mu + 5 sigma]. No other change was
  made, and no new mechanism was added.
- Second version: V's on [8, 64] at 9 of 54 settings, depths -0.037 to -0.075 (SE 0.002 to 0.018). All have
  c = 0.1 or 0.2 and l_star = mu + 1 or 2 sigma; none has beta_max = 0.3 or c = 0.4. On [16, 256] and
  [64, 1024] no depth is significant at 16 days. The a = 3 end has SE 0.03 to 0.1 there.
- Location (MEASURED vs Step 6): for l_star = mu + sigma, a* predicts 0.18 (psi 0.5) and 0.43 (psi 0.25),
  while the measured argmin is 0.75 and 1.0. For l_star = mu + 2 sigma, a* predicts 0.86 and the measured
  argmin is 1.5. The minimum sits 0.5 to 0.65 later in a than a*, as Step 6 expects from the continuing
  weight shift, but a* is not a quantitative predictor.

The setting used for the all-cell ground truth is beta_max = 0.7, c = 0.1, l_star = mu + sigma_l, psi = 0.5
(depth -0.045 +- 0.002 on [8, 64] and -0.058 +- 0.026 on [16, 256] at Q-T3S1), expressed relative to each
cell's mu and sigma_l. This is a CHOICE.
