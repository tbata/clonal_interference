# Some notes on modelling selective interference

TB — original notes dated July 4, 2016; extended October 2026

---

## 1. Introduction

Widespread interest among evolutionary geneticists in making predictions about parallel (//) evolution. Theory available so far envisions both mutational and selection heterogeneity across the genome and its consequences for expected patterns of parallel evolution. A central quantity of interest is the probability $P_i$ that a specific locus $i$ will contribute the next beneficial mutation to be fixed in an evolving population. If one considers a set of $L$ loci, one can use the $P_i$'s to predict, among other things, the probability of "gene reuse" as an operational measure of the degree of parallel evolution between 2 populations, or the distribution of counts of mutations per gene across a larger set of populations.

Theory available so far to calculate $P_i$ as a function of the heterogeneity in mutation and selection relies on the so-called SSWM assumption. Here the goal is to formulate a simple model for selective interference between mutations and consider whether a heuristic scaling proposed by TB & FB can be recovered under that set of assumptions.

### 1.1 Notation

We consider a set of $L$ loci. Each locus is characterised by its mutation rate $\mu_i$ and its distribution of fitness effects $\varphi_i$. The population is assumed to be composed of $N$ haploid asexual individuals. It is useful to denote the total genome-wide mutation rate for fitness as $U = \sum_i \mu_i$.

For simplicity we collapse the DFE of each locus to a simple atom at $s_i$, but should try to relax that assumption later. We denote by $T_i$ the waiting time for a new beneficial mutation to arise and escape early stochastic loss at locus $i$, and by $F_i$ the time needed for a mutation to go to quasi-fixation (conditional on fixation). To further simplify the calculations, and without loss of generality, we first order the loci by increasing beneficial selection coefficient $s_i$.

### 1.2 Known results for the SSWM condition

When a new beneficial mutation arises it escapes stochastic loss with probability $p(s_i)$. If one assumes that the population evolves under the strong selection weak mutation limit (also called the "Gillespie domain"), we can ignore selective interference between mutations. Define $\lambda_i = \mu_i s_i$; then the probability that focal locus $i$ contributes the next beneficial mutation to go to fixation is

$$P_i = \frac{\lambda_i}{\sum_j \lambda_j}.$$

### 1.3 Perfect interference limit and a heuristic scaling

Nagel et al. (2011) show that the locus with $\max_i(s_i)$ will contribute the next beneficial mutation. FB & TB recently proposed a purely heuristic modelling of the probability of a new mutation reaching fixation, using an alternative weighting of selection through $\lambda_i = \mu_i s_i^{\gamma}$, and noted that in the limit $\gamma \to \infty$ we recover the max rule and in the limit $\gamma \to 0$ we recover the SSWM limit. *(See note in §3.1: SSWM actually corresponds to an exponent of 1.)*

Below we formulate a probabilistic model for deriving the conditions under which the heuristic arises as a limit of a stochastic process describing selective interference. We also explore, via stochastic simulations, the process and the accuracy of the proposed scaling. One important insight is to relate the scaling parameter $\gamma$, initially introduced as a fix for selective interference (see Lenormand et al. 2016), to the biology underlying selective interference. We (hope to) show that $\gamma$ has a straightforward interpretation as the expected number of contending mutations that can potentially displace a mutation on its way to fixation after it escapes stochastic loss.

---

## 2. Modelling selective interference in the low-interference limit

We define $\Pi_i$ as the probability that locus $i$ will contribute the next beneficial mutation to go to fixation despite potential contending mutations. To calculate $\Pi_i$, we first condition on the event that a beneficial mutation at locus $i$ has just reached an appreciable frequency $\epsilon$ and could, if no interference occurs, sweep to fixation.

One of 3 things can happen:

- **A:** no contending mutation arises in the time interval needed to reach fixation;
- **B:** one (or more) mutation(s) arise(s) during the fixation phase but cannot displace the original mutation;
- **C:** one mutation arises during the fixation phase, escapes drift and confers a superior selective advantage, so it displaces the original mutation.

We simplify the problem drastically by assuming that the number of interfering mutations that can arise and escape early stochastic loss is small, and restrict ourselves to a single interference. This greatly simplifies the probabilities of events A, B and C.

We denote by $F_i$ the time for a mutation at locus $i$ to reach quasi-fixation. A lot of theory is available to provide expressions that can be "plugged in" for $F_i$; we expect $F_i$ to scale as $1/s_i$. We denote by $\lambda_i^+$ the sum of $\lambda_j$ at loci exceeding $i$, and by $\gamma_i$ the expected number of contending mutations that can arise during the fixation phase:

$$\gamma_i = N F_i \lambda_i^+.$$

The number of potential contending mutations is Poisson distributed with mean $\gamma_i$. Here we effectively assume $\gamma_i \ll 1$ and consider either zero or at most one interfering mutation. Then

$$P[A + B] = e^{-\gamma_i}, \qquad P[C] = 1 - e^{-\gamma_i}.$$

In words, a beneficial mutation at locus $i$ will be the next beneficial mutation to be fixed if:

- $s_i$ is the first contender and can sweep without interference from a superior mutation;
- $s_{i-1}$ is the first contender, but $s_i$ arises and displaces $s_{i-1}$ without being displaced itself;
- …
- $s_1$ is the first contender, but $s_i$ arises and displaces $s_1$ without being displaced itself.

We can write $\Pi_i$ as a sum of mutually disjoint events:

$$\Pi_i = p_i e^{-\gamma_i} + p_{i-1}(1-e^{-\gamma_{i-1}})\,p_i e^{-\gamma_i} + p_{i-2}(1-e^{-\gamma_{i-2}})\,p_i e^{-\gamma_i} + \dots + p_1(1-e^{-\gamma_1})\,p_i e^{-\gamma_i}.$$

Collecting the $p_i e^{-\gamma_i}$ terms and making a Taylor approximation as $\gamma \to 0$, $e^{-\gamma} \approx 1 - \gamma$, the original notes give:

$$\Pi_i \approx p_i\gamma_i\,(1 + p_{i-1}\gamma_{i-1} + p_{i-2}\gamma_{i-2} + \dots + p_1\gamma_1). \tag{1}$$

*(Eq. 1 is corrected in §3.8, Step 7.)*

### Original to-do list

1. Add the probability of a displaced mutation "landing" on a superior background.
2. Incorporate an explicit expression for $F_i$ and see if things get nicer or uglier.
3. Check $\Pi$ against stochastic simulations of the type run by FB, to examine the quality of the $s_i^{\gamma}$ approximation.

---

## 3. Additions (October 2026): step-by-step derivation of $\Pi_i$ and of the $s^{\Gamma}$ scaling

> **Revision note.** An earlier version of this section stated $\Gamma \approx 1 + \bar\gamma$, obtained from the response of one locus to its *own* $s$ (Step 8 below). The quantity that corresponds to the heuristic is instead the *cross-locus* slope of $\ln(\Pi_i/\mu_i)$ against $\ln s_i$ (Step 9). That slope contains an additional, contender-density term $n(s)$, which can be as large as $\gamma$ itself. The conclusion "exponent minus one = expected number of contenders" survives, with a more precise accounting of which contenders count.

### 3.0 Assumptions

All steps below rely on the following, stated once so each can be relaxed later.

- **(A1) Haploid asexual population of constant size $N$**, time in generations. Loci are ordered so that $s_1 < s_2 < \dots < s_L$, each locus having a single beneficial effect $s_i$ (atom DFE) with $s_i \ll 1$ and $N s_i \gg 1$.
- **(A2) Establishment.** A new mutation with advantage $s$ over the population it arises in escapes drift with probability $p(s) \approx \kappa s$. The constant $\kappa$ depends on offspring-number variance ($\kappa = 2$ for Wright–Fisher with Poisson offspring; Haldane 1927). Establishment is decided in a time short compared with a sweep, and once established a mutant grows deterministically.
- **(A3) Gerrish–Lenski contender rule.** A mutation arising in the ancestral background while another mutation is sweeping establishes with probability $p(s_k)$, i.e. its advantage is measured against the ancestor, not against the current population mean. It displaces the sweeping mutation if and only if it is superior ($s_k > s_j$). Inferior contenders are lost. Mutations arising on the sweeping background (double mutants) are ignored; this is original to-do item 1.
- **(A4) Renewal at takeover.** When a superior contender displaces the current leader, (i) the new leader starts its own sweep afresh, and (ii) any further contenders that arose during the previous leader's sweep are forgotten.
- **(A5) Start condition.** The population starts monomorphic for the ancestral genotype; we ask which locus supplies the *next* fixation.

### 3.1 Bookkeeping corrections to §2

- **First-contender weight.** The probability that locus $j$ supplies the first *established* mutation is $\lambda_j/\Lambda$ (Step 1), not $p_j$. Here $\lambda_j = N\mu_j\,p(s_j)$ is an establishment rate, and $\Lambda = \sum_k \lambda_k$.
- **Taylor step.** Eq. (1) vanishes as $\gamma \to 0$ instead of recovering SSWM. The correct expansion is in Step 7.
- **Limits of the heuristic.** With $\lambda_i = \mu_i s_i^{\gamma}$, the limit $\gamma \to 0$ gives weights $\propto \mu_i$, i.e. mutation only. SSWM ($\propto \mu_i s_i$) corresponds to an exponent of 1. Below the heuristic exponent is written $\Gamma$ to distinguish it from the contender number $\gamma_i$.

### 3.2 Step 1 — Arrival of established mutations

**Claim.** Starting from the monomorphic ancestor, established mutations at locus $j$ arrive as a Poisson process of rate

$$\lambda_j = N\mu_j\,p(s_j) \approx N\mu_j\,\kappa s_j,$$

independently across loci, and the first established mutation comes from locus $j$ with probability $\lambda_j/\Lambda$.

**Justification.**

1. *Mutation input.* Each generation, $N$ individuals each mutate at locus $j$ with probability $\mu_j \ll 1$, so the number of new $j$-mutations per generation is Binomial$(N,\mu_j) \approx$ Poisson$(N\mu_j)$. Over many generations this is a Poisson process with rate $N\mu_j$.
2. *Thinning.* By (A2) each new mutation independently establishes with probability $p(s_j)$. Independently marking the points of a Poisson process and keeping those marked "established" gives a Poisson process with rate $N\mu_j p(s_j)$ (thinning theorem).
3. *Superposition.* Independent Poisson processes with rates $\lambda_1,\dots,\lambda_L$ superpose into one Poisson process of rate $\Lambda = \sum_j \lambda_j$.
4. *First arrival.* The first arrival from locus $j$ occurs at an exponential time $T_j \sim \text{Exp}(\lambda_j)$, independently across $j$. For independent exponentials,
$$P\Big(T_j = \min_k T_k\Big) = \int_0^\infty \lambda_j e^{-\lambda_j t}\prod_{k\neq j} e^{-\lambda_k t}\,dt = \lambda_j \int_0^\infty e^{-\Lambda t}\,dt = \frac{\lambda_j}{\Lambda}.$$

Note that $\kappa$ and $N$ cancel in $\lambda_j/\Lambda$, so in SSWM $P_j = \mu_j s_j/\sum_k \mu_k s_k$, recovering §1.2.

### 3.3 Step 2 — The sweep window $F_j$

**Claim.** The effective time window during which contenders can arise in the ancestral background while a $j$-mutation sweeps is

$$F_j \equiv \int_0^\infty \big(1-x_j(t)\big)\,dt = \frac{\ln(N\kappa s_j)}{s_j}.$$

**Justification.**

1. *Why $1-x$.* By (A3), contenders arise only in ancestral individuals, whose number at time $t$ is $N(1-x_j(t))$. So the contender input rate at time $t$ is proportional to $1-x_j(t)$; integrating over the sweep gives an "effective number of generations of a full-size ancestral population", which is $F_j$.
2. *Starting frequency.* A mutation is effectively established once its copy number reaches order $1/p(s_j) = 1/(\kappa s_j)$; above that size, loss by drift is unlikely. So the deterministic phase starts at
$$x_0 = \frac{1}{N\kappa s_j}.$$
3. *Deterministic sweep.* Above $x_0$, (A2) gives logistic growth $\dot x = s_j x(1-x)$. Separating variables,
$$x_j(t) = \frac{x_0 e^{s_j t}}{1 - x_0 + x_0 e^{s_j t}}, \qquad 1 - x_j(t) = \frac{A e^{-s_j t}}{1 + A e^{-s_j t}}, \quad A = \frac{1-x_0}{x_0}.$$
4. *Integral.* Since $\frac{d}{dt}\big[-\tfrac{1}{s_j}\ln(1+Ae^{-s_j t})\big] = \frac{Ae^{-s_j t}}{1+Ae^{-s_j t}}$,
$$F_j = \Big[-\tfrac{1}{s_j}\ln\big(1+Ae^{-s_j t}\big)\Big]_0^\infty = \frac{1}{s_j}\ln(1+A) = \frac{1}{s_j}\ln\frac{1}{x_0} = \frac{\ln(N\kappa s_j)}{s_j}.$$

This confirms the intuition in §2 that $F_j \propto 1/s_j$, with a logarithmic correction. Two properties are used later:

$$\frac{dF}{ds} = \frac{1 - \ln(N\kappa s)}{s^2} < 0 \quad (\text{as } N\kappa s > e), \qquad s\,\frac{dF}{ds} = -F\Big(1 - \frac{1}{\ln N\kappa s}\Big).$$

The first says higher-$s$ loci sweep faster; the second follows by dividing the first by $F/s = \ln(N\kappa s)/s^2$. For brevity write $\ell_s \equiv \ln(N\kappa s)$.

### 3.4 Step 3 — Number of superior contenders, $\gamma_j$

**Claim.** During a $j$-sweep, the number of established superior contenders is Poisson with mean

$$\gamma_j = F_j\,\Lambda_j^+, \qquad \Lambda_j^+ \equiv \sum_{k>j}\lambda_k,$$

so the leader is not displaced with probability $e^{-\gamma_j}$.

**Justification.**

1. At time $t$ into the sweep, established $k$-contenders arise at rate $N(1-x_j(t))\,\mu_k\,p(s_k) = \lambda_k(1-x_j(t))$, by Steps 1–2 and (A3).
2. This is an inhomogeneous Poisson process; its total count over $[0,\infty)$ is Poisson with mean equal to the integrated rate, $\lambda_k \int_0^\infty (1-x_j)\,dt = \lambda_k F_j$.
3. Only $k>j$ can displace $j$ (A3); superposing these independent processes gives a Poisson count with mean $F_j\sum_{k>j}\lambda_k = \gamma_j$.
4. The probability of zero superior contenders is the Poisson zero class, $e^{-\gamma_j}$. Inferior contenders (event B of §2) are lost and do not alter $j$'s fate, so

$$P[A+B] = e^{-\gamma_j}, \qquad P[C] = 1 - e^{-\gamma_j},$$

exactly as in §2, but now without needing $\gamma_j \ll 1$: the Poisson form already accounts for any number of contenders.

Consistency with §2: with $\lambda = N\mu\kappa s$, $\gamma_j = N F_j\,\kappa\sum_{k>j}\mu_k s_k$, which is the original $\gamma_i = N F_i \lambda_i^+$ up to the constant $\kappa$. The top locus has $\Lambda_L^+ = 0$, so $\gamma_L = 0$.

### 3.5 Step 4 — Identity of the displacer

**Claim.** Given that $j$ is displaced, the displacer comes from locus $k>j$ with probability $\lambda_k/\Lambda_j^+$.

**Justification.** All superior processes share the same time profile $1-x_j(t)$; their intensities are $\lambda_k(1-x_j(t))$, i.e. proportional to each other with constant ratios. An inhomogeneous Poisson process whose intensity splits into components in fixed proportions is equivalent to one process of intensity $\Lambda_j^+(1-x_j(t))$ whose points carry independent labels, label $k$ with probability $\lambda_k/\Lambda_j^+$ (the marking theorem). The first point therefore carries label $k$ with probability $\lambda_k/\Lambda_j^+$, independently of when it occurs. Under (A4) the first superior contender becomes the new leader.

### 3.6 Step 5 — The "current leader" Markov chain

Combining Steps 1–4 and (A4), the identity of the current leader is a Markov chain on $\{1,\dots,L\}$:

- **Initial distribution:** $P(X_0 = i) = \lambda_i/\Lambda$ (Step 1).
- **Absorption (fixation of $j$):** probability $e^{-\gamma_j}$ (Step 3).
- **Transition $j \to k$, $k>j$:** probability $(1-e^{-\gamma_j})\,\lambda_k/\Lambda_j^+$ (Steps 3 and 4, multiplied because "displaced" and "by $k$" are a conditional chain of events).

Rows sum to one: $e^{-\gamma_j} + (1-e^{-\gamma_j})\sum_{k>j}\lambda_k/\Lambda_j^+ = e^{-\gamma_j} + (1-e^{-\gamma_j}) = 1$. The chain moves only to higher indices and locus $L$ is always absorbed ($\gamma_L = 0$), so absorption occurs with probability 1 in at most $L$ steps.

**Recursion.** Let $a_i$ be the probability that locus $i$ is ever the leader. Then

$$a_i = \frac{\lambda_i}{\Lambda} + \sum_{j<i} a_j\,(1-e^{-\gamma_j})\,\frac{\lambda_i}{\Lambda_j^+}, \qquad \Pi_i = a_i\,e^{-\gamma_i}.$$

**Justification.**

1. Because indices strictly increase, a path visits $i$ at most once.
2. The event "$i$ is visited" is the union of "the path starts at $i$" and, for each $j<i$, "the path visits $j$ and moves next to $i$". These events are disjoint because the state immediately before $i$ (or the absence of one) is unique.
3. By the Markov property, $P(\text{visit } j \text{ and next is } i) = a_j \cdot P(j\to i)$.
4. Locus $i$ supplies the fixation iff it is visited and then absorbed, so $\Pi_i = a_i e^{-\gamma_i}$.

### 3.7 Step 6 — Closed-form solution

Define the per-locus quantity

$$r_j \equiv \frac{\lambda_j\,(1-e^{-\gamma_j})}{\Lambda_j^+}\quad (j<L).$$

**Claim.**

$$\boxed{\;\Pi_i = \frac{\lambda_i}{\Lambda}\;e^{-\gamma_i}\;B_i,\qquad B_i = \prod_{j<i}\big(1 + r_j\big),\quad B_1 = 1.\;}$$

**Proof.** Set $b_i = a_i/\lambda_i$. Dividing the recursion by $\lambda_i$ (which factors out of every term because $P(j\to i)$ is proportional to $\lambda_i$):

$$b_i = \frac{1}{\Lambda} + \sum_{j<i}\frac{a_j(1-e^{-\gamma_j})}{\Lambda_j^+}.$$

The summand does not depend on $i$; only the upper limit does. Hence

$$b_{i+1} - b_i = \frac{a_i(1-e^{-\gamma_i})}{\Lambda_i^+} = b_i\,\frac{\lambda_i(1-e^{-\gamma_i})}{\Lambda_i^+} = b_i\,r_i,$$

so $b_{i+1} = b_i(1+r_i)$ with $b_1 = 1/\Lambda$. Iterating, $b_i = \Lambda^{-1}\prod_{j<i}(1+r_j)$, and $\Pi_i = \lambda_i b_i e^{-\gamma_i}$. $\square$

**Reading the formula.** $\lambda_i/\Lambda$ is the SSWM weight. $e^{-\gamma_i}$ is the cost of being displaced by superior loci. $B_i \ge 1$ is the gain from displacing inferior leaders: it increases with $i$, because higher loci can enter as contenders during the sweeps of more inferior loci.

**Relation to classical Gerrish–Lenski.** GL's distribution of fixed effects corresponds to $\Pi_i \propto \lambda_i e^{-\gamma_i}$, i.e. $B_i \equiv 1$. That version lets every established mutation try to fix at the same rate, including inferior ones that arise while a superior mutation is already sweeping (and so would in fact be lost). The chain removes those and, in exchange, credits superior loci with the extra entries counted by $B_i$.

**Check, $L = 2$.** $\gamma_1 = F_1\lambda_2$, $\gamma_2 = 0$, $r_1 = \lambda_1(1-e^{-\gamma_1})/\lambda_2$. Then

$$\Pi_1 = \frac{\lambda_1 e^{-\gamma_1}}{\Lambda}, \qquad \Pi_2 = \frac{\lambda_2}{\Lambda}\Big(1 + \frac{\lambda_1(1-e^{-\gamma_1})}{\lambda_2}\Big) = \frac{\lambda_2 + \lambda_1 - \lambda_1 e^{-\gamma_1}}{\Lambda},$$

and $\Pi_1 + \Pi_2 = (\lambda_1+\lambda_2)/\Lambda = 1$. Directly: locus 2 fixes if it comes first, or if 1 comes first and is displaced, as expected.

### 3.8 Step 7 — Weak-interference expansion (replaces Eq. 1)

Let $\varepsilon \equiv \Lambda F_1$. Because $F$ decreases with $s$ (Step 2), $F_1$ is the longest window, so $\gamma_j = F_j\Lambda_j^+ \le \varepsilon$ for all $j$ and $\sum_j \lambda_j F_j \le \varepsilon$. Assume $\varepsilon \ll 1$.

1. $e^{-\gamma_i} = 1 - \gamma_i + O(\gamma_i^2)$ (Taylor series of $e^{-x}$).
2. $1 - e^{-\gamma_j} = \gamma_j + O(\gamma_j^2)$, so $r_j = \lambda_j\gamma_j/\Lambda_j^+ + O(\varepsilon^2) = \lambda_j F_j + O(\varepsilon^2)$.
3. $\prod_{j<i}(1+r_j) = 1 + \sum_{j<i} r_j + O\big((\sum_j r_j)^2\big) = 1 + \sum_{j<i}\lambda_j F_j + O(\varepsilon^2)$, expanding the product and bounding all cross-products by $(\sum r_j)^2 \le \varepsilon^2$.
4. Multiplying 1 and 3 and dropping products of two small terms:

$$\Pi_i = \frac{\lambda_i}{\Lambda}\Big[\,1 - F_i\Lambda_i^+ + \sum_{j<i}\lambda_j F_j\,\Big] + O(\varepsilon^2).$$

**Normalisation check.** Summing $\lambda_i\times$ the two correction terms over $i$:

$$\sum_i \lambda_i F_i \Lambda_i^+ = \sum_{i<k}\lambda_i\lambda_k F_i, \qquad \sum_i \lambda_i\sum_{j<i}\lambda_j F_j = \sum_{j<i}\lambda_j\lambda_i F_j,$$

which are the same double sum (rename $j \to i$, $i \to k$), so the corrections cancel and $\sum_i\Pi_i = 1 + O(\varepsilon^2)$. In words, every displacement is a loss for one locus and an equal gain for another.

Comparison with Eq. (1): the original single-interference sum has the right structure (a first contender, then at most one displacement), but uses $p_j$ in place of $\lambda_j/\Lambda$ and $\lambda_i/\Lambda_j^+$, and its expansion replaced $e^{-\gamma_i}$ by $\gamma_i$ instead of $1-\gamma_i$, which is why Eq. (1) vanished as $\gamma\to 0$. The expansion above is the corrected version.

### 3.9 Step 8 — Own-locus elasticity (what changes if one locus's $s$ changes)

Hold all $\mu$'s and all other $s$'s fixed, move $s_i$ within $(s_{i-1}, s_{i+1})$ so the ordering is unchanged, and compute $\partial \ln\Pi_i/\partial\ln s_i$ from the closed form

$$\ln\Pi_i = \ln\lambda_i - \ln\Lambda - \gamma_i + \sum_{j<i}\ln(1+r_j).$$

Term by term, using $\partial/\partial\ln s_i = s_i\,\partial/\partial s_i$ and $\partial\lambda_i/\partial\ln s_i = \lambda_i$ (since $\lambda_i \propto s_i$):

- (a) $\partial\ln\lambda_i/\partial\ln s_i = 1$.
- (b) $\partial(-\ln\Lambda)/\partial\ln s_i = -\lambda_i/\Lambda$.
- (c) $\Lambda_i^+$ contains only $k>i$, so it does not depend on $s_i$; by Step 2, $\partial(-\gamma_i)/\partial\ln s_i = -\Lambda_i^+ s_i F'(s_i) = \gamma_i(1 - 1/\ell_{s_i})$.
- (d) For $j<i$, $s_i$ enters $r_j$ through $\Lambda_j^+$ (which contains $\lambda_i$) and through $\gamma_j = F_j\Lambda_j^+$. With $\partial\Lambda_j^+/\partial\ln s_i = \lambda_i$ and $\partial\gamma_j/\partial\ln s_i = F_j\lambda_i$, the quotient rule gives
$$\frac{\partial r_j}{\partial\ln s_i} = \frac{\lambda_j\lambda_i}{(\Lambda_j^+)^2}\,g(\gamma_j),\qquad g(\gamma) = \gamma e^{-\gamma} - (1-e^{-\gamma}).$$
Expanding, $g(\gamma) = -\gamma^2/2 + O(\gamma^3)$, and $-1 \le g \le 0$.

Collecting,

$$\frac{\partial\ln\Pi_i}{\partial\ln s_i} = 1 + \gamma_i\Big(1 - \frac{1}{\ell_{s_i}}\Big) - \frac{\lambda_i}{\Lambda} + \sum_{j<i}\frac{\lambda_i\lambda_j\,g(\gamma_j)}{(\Lambda_j^+)^2(1+r_j)}.$$

The last two terms are small when locus $i$ is a small share of total supply, and the final sum is additionally $O(\gamma^2)$. So the own-locus response is $\approx 1 + \gamma_i$. This was the result reported in the earlier version. It answers "how much does locus $i$ gain if its own $s$ rises", but it is *not* what the heuristic describes: the heuristic is a statement about how $\Pi$ varies *across* loci with different $s$ (Step 9).

### 3.10 Step 9 — Cross-locus slope: the quantity that matches the heuristic

The heuristic says $\Pi_i = \mu_i s_i^{\Gamma}/Z$, with $Z$ a normalising constant common to all loci. So it predicts that, for any two loci $i<k$,

$$\ln\frac{\Pi_k/\mu_k}{\Pi_i/\mu_i} = \Gamma\,\ln\frac{s_k}{s_i}.$$

**Discrete (secant) exponent.** From the closed form, $\Pi_i/\mu_i = (N\kappa/\Lambda)\,s_i\,e^{-\gamma_i}B_i$. The factor $N\kappa/\Lambda$ is common and cancels in the ratio, and $B_k/B_i = \prod_{i\le j<k}(1+r_j)$. Hence

$$\ln\frac{\Pi_k/\mu_k}{\Pi_i/\mu_i} = \ln\frac{s_k}{s_i} + (\gamma_i - \gamma_k) + \sum_{i\le j<k}\ln(1+r_j),$$

and the exponent that makes the heuristic exact for this pair is

$$\Gamma_{ik} = 1 + \frac{(\gamma_i - \gamma_k) + \sum_{i\le j<k}\ln(1+r_j)}{\ln(s_k/s_i)}.$$

The difference in contender numbers splits naturally. Writing $\Lambda_i^+ = \sum_{i<j\le k}\lambda_j + \Lambda_k^+$,

$$\gamma_i - \gamma_k = \underbrace{(F_i - F_k)\,\Lambda_k^+}_{\text{common superiors; } i \text{ is exposed longer}} + \underbrace{F_i\sum_{i<j\le k}\lambda_j}_{\text{loci between } i \text{ and } k}.$$

All three contributions to $\Gamma_{ik} - 1$ are non-negative: $F_i > F_k$ (Step 2), $\lambda_j \ge 0$, and $r_j \ge 0$.

**Continuum limit.** Now suppose loci are dense in $s$. Let $u(s)\,ds$ be the total mutation rate of loci with effect in $(s, s+ds)$, and

$$\ell(s) = N\kappa s\,u(s), \qquad \Lambda^+(s) = \int_s^\infty \ell(s')\,ds', \qquad \gamma(s) = F(s)\Lambda^+(s).$$

Each $\lambda_j \to \ell(s_j)\Delta s$ is infinitesimal, so $r_j$ is small, $\ln(1+r_j) \approx r_j$, and the sum becomes an integral (Riemann sum):

$$\ln B(s) = \int_0^s \frac{\ell(s')\big(1-e^{-\gamma(s')}\big)}{\Lambda^+(s')}\,ds'.$$

The weight per unit mutation rate is $W(s) = s\,e^{-\gamma(s)}B(s)$, and the local exponent is its log-log slope, $\Gamma(s) = d\ln W/d\ln s = s\,d\ln W/ds$:

1. $s\,\dfrac{d\ln s}{ds} = 1$.
2. $-s\gamma'(s)$, with $\gamma' = F'\Lambda^+ + F\,\dfrac{d\Lambda^+}{ds}$ (product rule) and $\dfrac{d\Lambda^+}{ds} = -\ell(s)$ (Leibniz rule for a variable lower limit). Using Step 2 for $sF'$:
$$-s\gamma' = \gamma\Big(1-\frac{1}{\ell_s}\Big) + s\,\ell(s)F(s).$$
3. $s\,\dfrac{d\ln B}{ds} = s\,\ell(s)\,\dfrac{1-e^{-\gamma(s)}}{\Lambda^+(s)}$ (fundamental theorem of calculus).

Define

$$n(s) \equiv s\,\ell(s)\,F(s),$$

the expected number of established mutations arising during one sweep window whose effects lie within one unit of $\ln s$ (since $\ell\,ds = s\ell\,d\ln s$). Noting that $s\ell(1-e^{-\gamma})/\Lambda^+ = n\,(1-e^{-\gamma})/\gamma$, the three pieces add to

$$\boxed{\;\Gamma(s) = 1 + \gamma(s)\Big(1 - \frac{1}{\ell_s}\Big) + n(s)\Big[1 + \frac{1-e^{-\gamma(s)}}{\gamma(s)}\Big].\;}$$

Limits of the bracket: it tends to 2 as $\gamma \to 0$ and to 1 as $\gamma \to \infty$. So

$$\Gamma(s) \approx 1 + \gamma(s) + 2n(s)\ \ (\gamma \ll 1), \qquad \Gamma(s) \approx 1 + \gamma(s) + n(s)\ \ (\gamma \gg 1),$$

up to the $1/\ell_s$ correction.

### 3.11 Step 10 — Interpretation and a worked example

**Interpretation.** $\Gamma - 1$ is a count of contenders, made of two kinds:

- $\gamma(s)$: superior contenders expected during one sweep (the cost of being displaced);
- $n(s)$: contenders of *similar* effect (within a factor $e$ of $s$). Moving up the $s$ axis removes such neighbours from the set of superiors (first $n$) and adds them to the set of inferiors one can displace (second $n$, weighted by $(1-e^{-\gamma})/\gamma$).

This is the precise form of the conjecture in §1.3. In the classical GL weighting ($B \equiv 1$, Step 6) the second $n$ term is absent: $\Gamma_{GL} = 1 + \gamma(1-1/\ell_s) + n$.

**How big is $n$ relative to $\gamma$?** Their ratio is

$$\frac{n(s)}{\gamma(s)} = \frac{s\,\ell(s)}{\Lambda^+(s)},$$

which depends on the shape of the distribution of effects across loci, not on $N$ or $U$.

**Example: exponential spread of effects.** Take $u(s) = (U/\sigma)\,e^{-s/\sigma}$. Then $\ell(s) = N\kappa U (s/\sigma)e^{-s/\sigma}$ and, integrating by parts,

$$\Lambda^+(s) = N\kappa U\int_s^\infty \frac{s'}{\sigma}e^{-s'/\sigma}ds' = N\kappa U\,(s+\sigma)\,e^{-s/\sigma}.$$

Hence

$$\frac{n}{\gamma} = \frac{s^2}{\sigma(s+\sigma)}, \qquad \gamma(s) = \frac{\ell_s}{s}\,N\kappa U\,(s+\sigma)\,e^{-s/\sigma}.$$

For the effects that actually fix under interference ($s$ a few $\sigma$), $n/\gamma \approx s/\sigma > 1$, so in the weak-interference regime $\Gamma \approx 1 + \gamma(1 + 2s^2/[\sigma(s+\sigma)])$, e.g. $\Gamma \approx 1 + 2.5\gamma$ at $s = 3\sigma$. The density term is not a small correction for this DFE.

### 3.12 Step 11 — From local slopes to a single fitted $\Gamma$

$\Gamma(s)$ varies with $s$, while the heuristic uses one number. If $\Gamma$ is estimated by (weighted) least-squares regression of $y = \ln(\Pi_i/\mu_i)$ on $x = \ln s_i$, the fitted slope is a weighted average of the local slopes:

$$\hat\Gamma = \frac{\int \Gamma(x)\,\omega(x)\,dx}{\int\omega(x)\,dx},\qquad \omega(x) = \int_x^\infty (t - \bar x)\,dG(t) \ \ge 0,$$

where $G$ is the (weighted) distribution of $x$ across loci and $\bar x$ its mean. This identity follows from writing $\text{Cov}(x, y)$ with $y(x) = y(x_{\min}) + \int_{x_{\min}}^x y'(t)\,dt$ and exchanging the order of integration (Yitzhaki 1996). The weight $\omega$ is largest near the middle of the $\ln s$ range. So the single $\Gamma$ of the heuristic is

$$\hat\Gamma \approx 1 + \langle\gamma\rangle_\omega + c\,\langle n\rangle_\omega,\qquad c\in[1,2],$$

an average of contender numbers over the loci, and a single exponent describes the data well when $\Gamma(s)$ varies little over the $\ln s$ range that carries most weight. Curvature of $\ln(\Pi_i/\mu_i)$ against $\ln s_i$ is the diagnostic.

### 3.13 Limits and validity

- **SSWM:** $NU F \to 0$ makes $\gamma, n \to 0$ for all loci, so $\Gamma \to 1$ and $\Pi_i \to \mu_i s_i/\sum_k \mu_k s_k$. ✓
- **Perfect interference:** as $N$ grows with everything else fixed, $\gamma_i$ grows like $N\ln N$ for every $i<L$ (Steps 2–3), while $B_i \le \prod_j(1 + \lambda_j/\Lambda_j^+)$ stays bounded. So $\Pi_i \to 0$ for $i<L$ and $\Pi_L \to 1$: the max rule. Correspondingly $\Gamma \to \infty$. ✓
- **Range of validity:** (A3)–(A4) require that contenders are rare enough to be treated one sweep at a time, and that double mutants are negligible; this fails roughly when $NU\ln(Ns) \gtrsim 1$ with effects of similar size (Gomez et al. 2020; Good et al. 2012). Beyond that point, the multiple-mutations theory replaces Steps 3–5.
- **Refinement of (A3):** strictly, a contender arising at time $t$ has advantage $s_k - s_j x_j(t)$ over the population mean, so its establishment probability is $\kappa(s_k - s_j x_j(t))$ rather than $\kappa s_k$. This changes $\gamma_j$ to $N\kappa\sum_{k>j}\mu_k\int(1-x_j)(s_k - s_j x_j)\,dt$, which can be evaluated with the logistic $x_j(t)$ of Step 2.

### 3.14 Simulation test (original to-do item 3)

- Simulate Wright–Fisher populations with $L$ loci and heterogeneous $(\mu_i, s_i)$, started monomorphic.
- Estimate $\Pi_i$ as the frequency with which locus $i$ supplies the first fixed mutation across replicates.
- Compare the simulated $\Pi_i$ with the closed form of Step 6, and with the GL weighting ($B\equiv 1$), to measure the contribution of $B_i$.
- Regress $\ln(\Pi_i/\mu_i)$ on $\ln s_i$; compare the slope with $\hat\Gamma$ from Step 11 computed from the model parameters.
- Scan $NU\ln(Ns)$ to locate where the chain (and hence the $s^{\Gamma}$ scaling) breaks down.

### 3.15 Updated to-do list

1. Extend the chain to allow displacers to arise on the leader's background (multiple mutations), relaxing (A3).
2. Implement the refined establishment probability of §3.13 and check how much it changes $\gamma_j$.
3. Relax the atom DFE: replace $s_i$ by $\varphi_i$; $\Gamma$ should become an average elasticity over each locus's DFE.
4. Joint exponents: fit $\Pi_i \propto \mu_i^{\beta}s_i^{\Gamma}$ and check whether $\beta$ departs from 1 under interference (cf. Cano et al. 2022). In the chain, $\mu_i$ also enters other loci's weights (through $r_i$ in $B_k$ for $k>i$, and through $\Lambda_j^+$ for $j<i$), so the cross-locus $\beta$ need not equal 1 — worth deriving along the lines of Step 9.
5. Connect to the rigorous Poisson-contender limit of Hermann et al. (2025).

**Additional reference for this section:** Yitzhaki, S. (1996) On using linear regressions in welfare economics. *J. Bus. Econ. Stat.* 14: 478–486 (regression slope as a weighted average of derivatives). Haldane, J.B.S. (1927) A mathematical theory of natural and artificial selection, part V. *Proc. Camb. Phil. Soc.* 23: 838–844 (establishment probability).

---

## 4. Literature notes (search done Oct 2026)

No 2021–2026 paper in *Genetics* or *Theoretical Population Biology* was found that derives a locus-level weight like $\mu_i s_i^{\gamma}$ under interference. The search was not exhaustive; worth a check of recent *TPB* issues for "Gerrish–Lenski".

**In target journals, tangential**

- *Genetics* (2022) 220(1): iyab134. Polygenic adaptation under global vs local adaptation; repeatability of the loci involved depends on genetic redundancy. Sexual populations with migration, no clonal interference. <https://academic.oup.com/genetics/article/220/1/iyab134/6497714>
- *Genetics* (2024) 228(1): iyae111. Pleiotropy, migration and repeated genetic adaptation; large-effect loci can still contribute repeatedly despite pleiotropic trade-offs. No interference. <https://academic.oup.com/genetics/article/228/1/iyae111/7712980>
- *Genetics* (2019). Interference effects of deleterious and beneficial mutations in large asexual populations; fixation probability under combined interference. Just outside the 5-year window. <https://doi.org/10.1534/genetics.119.301960>

**Closest in substance, other journals**

- Hermann, González Casanova, Soares dos Santos, Tóbiás & Wakolbinger (2025) *Ann. Appl. Probab.* 35(4): 2823–2865. "From clonal interference to Poissonian interacting trajectories." Rigorous limit for the Gerrish–Lenski regime with a finite number of Poisson contenders; closest formalisation of the setup in §2. <https://doi.org/10.1214/25-AAP2188>
- Cano et al. (2022) *PNAS* 119: e2119720119. "Mutation bias shapes the spectrum of adaptive substitutions." Exponent $\beta$ for mutation rate, the mutational twin of $\Gamma$; clonal interference is expected to lower $\beta$. <https://www.pnas.org/doi/10.1073/pnas.2119720119>
- Gomez, Bertram & Masel (2020) *Proc. R. Soc. B*. "Mutation bias can shape adaptation in large asexual populations experiencing clonal interference." The ratio $E_s/E_U$ of elasticities of adaptation to $s$ vs $U$ is essentially $\Gamma$ in the traveling-wave regime. <https://royalsocietypublishing.org/doi/10.1098/rspb.2020.1503>
- Schenk et al. (2022) *Nat. Ecol. Evol.* 6: 439–447. "Population size mediates the contribution of high-rate and large-benefit mutations to parallel evolution." Empirical test: larger populations shift fixations from high-rate to large-benefit mutations. <https://www.nature.com/articles/s41559-022-01669-3>
- Gitschlag, Cano, Payne, McCandlish & Stoltzfus (2023) *Am. Nat.* 202(4): 534–557. Joint distributions of mutation rates and selection coefficients; SSWM only. <https://www.journals.uchicago.edu/doi/10.1086/726014>
- Good, Rouzine, Balick, Hallatschek & Desai (2012) *PNAS*. "Distribution of fixed beneficial mutations and the rate of adaptation in asexual populations." Needed beyond Gerrish–Lenski. <https://www.pnas.org/doi/10.1073/pnas.1119910109>

**Background**

- Gerrish & Lenski (1998) *Genetica* 102/103: 127–144.
- Bailey, Blanquart, Bataillon & Kassen (2017) *BioEssays*. "What drives parallel evolution?"
- Lenormand, Chevin & Bataillon (2016) *Phil. Trans. R. Soc. B*. Parallel evolution: what does it (not) tell us.
- Nagel et al. (2011).
