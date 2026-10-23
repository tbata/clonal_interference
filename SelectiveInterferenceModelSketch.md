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

Nagel et al. (2011) show that the locus with $\max_i(s_i)$ will contribute the next beneficial mutation. FB & TB recently proposed a purely heuristic modelling of the probability of a new mutation reaching fixation, using an alternative weighting of selection through $\lambda_i = \mu_i s_i^{\gamma}$, and noted that in the limit $\gamma \to \infty$ we recover the max rule and in the limit $\gamma \to 0$ we recover the SSWM limit. *(See note in §3.2: SSWM actually corresponds to $\gamma = 1$.)*

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

*(Eq. 1 is corrected in §3.1.)*

### Original to-do list

1. Add the probability of a displaced mutation "landing" on a superior background.
2. Incorporate an explicit expression for $F_i$ and see if things get nicer or uglier.
3. Check $\Pi$ against stochastic simulations of the type run by FB, to examine the quality of the $s_i^{\gamma}$ approximation.

---

## 3. Additions (October 2026): an exact Gerrish–Lenski formulation and a derivation of the $s^{\Gamma}$ scaling

### 3.1 Bookkeeping corrections to §2

- **First-contender weight.** The probability that locus $j$ supplies the first *established* mutation is $\lambda_j/\Lambda$, not $p_j$. Here $\lambda_j = N\mu_j\,p(s_j)$ is the establishment rate and $\Lambda = \sum_k \lambda_k$. Likewise, given that leader $j$ is displaced, the displacer comes from locus $i>j$ with probability $\lambda_i/\Lambda_j^+$, where $\Lambda_j^+ = \sum_{k>j}\lambda_k$.
- **Taylor step.** Eq. (1) vanishes as $\gamma \to 0$ instead of recovering SSWM. The correct first-order expansion is given in §3.3.
- **Limits of the heuristic.** With $\lambda_i = \mu_i s_i^{\gamma}$, the limit $\gamma \to 0$ gives weights $\propto \mu_i$, i.e. mutation only. SSWM ($\propto \mu_i s_i$) corresponds to $\gamma = 1$.

### 3.2 The "current leader" as an absorbing Markov chain

The state of the chain is the locus currently sweeping. From leader $j$:

- the chain is absorbed (locus $j$ fixes) with probability $e^{-\gamma_j}$, where $\gamma_j = F_j\Lambda_j^+$;
- otherwise it jumps to superior locus $k>j$ with probability $(1-e^{-\gamma_j})\,\lambda_k/\Lambda_j^+$.

Let $a_i$ be the probability that locus $i$ is ever the leader. Then

$$a_i = \frac{\lambda_i}{\Lambda} + \sum_{j<i} a_j\,(1-e^{-\gamma_j})\,\frac{\lambda_i}{\Lambda_j^+}, \qquad \Pi_i = a_i\,e^{-\gamma_i}.$$

This is exact within Gerrish–Lenski assumptions and does **not** require $\gamma \ll 1$ or a single interference event. To-do item 1 (a displaced mutation landing on a superior background) is exactly what this chain ignores: including it is the step into the multiple-mutations regime.

### 3.3 Small-$\gamma$ check

To first order in $\gamma$:

$$\Pi_i \approx \frac{\lambda_i}{\Lambda}\Big[\,1 - F_i\Lambda_i^+ + \sum_{j<i}\lambda_j F_j\,\Big].$$

The term $-F_i\Lambda_i^+$ is what locus $i$ loses by being displaced; $\sum_{j<i}\lambda_j F_j$ is what it gains by displacing inferior loci. The $\Pi_i$ sum to 1 exactly at this order. This replaces Eq. (1).

### 3.4 The exponent as an elasticity (to-do item 2)

Take $p(s)\approx s$, and take the window over which contenders can arise in the ancestral background as

$$F_i \approx \frac{\ln(N s_i)}{s_i},$$

i.e. the integral of $1-x(t)$ over the sweep. Hold $\mu_i$ fixed and move $s_i$ between its neighbours, so that $\Lambda_i^+$ is constant. Then

$$\Gamma_i \equiv \frac{\partial \ln \Pi_i}{\partial \ln s_i} = 1 + \gamma_i\Big(1 - \frac{1}{\ln N s_i}\Big) + O\!\left(\lambda_i/\Lambda_j^+\right),$$

where the correction comes from $\lambda_i$ also appearing in the normalisations $\Lambda$ and $\Lambda_j^+$. This holds for any $\gamma_i$, not only small ones.

### 3.5 Recovering the heuristic

The weight $\mu_i s_i^{\Gamma}$ is the log-linear approximation to the Gerrish–Lenski weight $s_i\,e^{-\gamma_i}$, with

$$\Gamma \approx 1 + \bar\gamma,$$

where $\bar\gamma$ is the expected number of established superior contenders during a sweep, averaged over the loci that actually fix (weights $\Pi_i$).

**Interpretation:** the exponent minus one equals the expected number of contenders. Limits: $\bar\gamma \to 0$ gives SSWM ($\mu s$); $\bar\gamma \to \infty$ gives the max rule.

### 3.6 Where a single exponent breaks down

$\gamma_i$ decreases with $i$ (fewer superior loci, shorter sweeps), so $\ln(\Pi_i/\mu_i)$ against $\ln s_i$ is curved rather than linear. A clean summary is the secant exponent between loci $i<k$:

$$\Gamma_{ik} = 1 + \frac{\gamma_i - \gamma_k}{\ln(s_k/s_i)}.$$

A single $\Gamma$ works when $\Gamma_{ik}$ varies little across the loci that carry most fixations.

### 3.7 Simulation test (to-do item 3)

- Simulate Wright–Fisher populations with $L$ loci and heterogeneous $(\mu_i, s_i)$.
- Estimate $\Pi_i$ from the identity of the first fixed mutation over many replicates.
- Regress $\ln(\Pi_i/\mu_i)$ on $\ln s_i$ (same spirit as the $\beta$ regression of Cano et al. 2022) and compare the slope with $1 + \bar\gamma$ (weighted by $\Pi_i$).
- Compare the Markov-chain $\Pi_i$ of §3.2 directly against simulated $\Pi_i$.
- Scan $N U \ln(Ns)$: the prediction should hold in the Gerrish–Lenski window and fail once multiple mutations matter (Gomez et al. 2020 place that boundary at roughly $UN\ln(Ns) \gtrsim 1$). Beyond it, the fixation probability of Good et al. (2012) implies a different, larger elasticity tied to how far the fitness "nose" leads the population mean.

### 3.8 Updated to-do list

1. Extend the chain of §3.2 to allow displacers to arise on the leader's background (multiple mutations).
2. Relax the atom DFE: replace $s_i$ by $\varphi_i$ and check whether $\Gamma$ becomes an average elasticity over each locus's DFE.
3. Joint exponents: fit $\Pi_i \propto \mu_i^{\beta} s_i^{\Gamma}$ and check whether $\beta$ departs from 1 under interference (cf. Cano et al. 2022).
4. Connect to the rigorous Poisson-contender limit of Hermann et al. (2025).

---

## 4. Literature notes (search done Oct 2026)

No 2021–2026 paper in *Genetics* or *Theoretical Population Biology* was found that derives a locus-level weight like $\mu_i s_i^{\gamma}$ under interference. The search was not exhaustive; worth a check of recent *TPB* issues for "Gerrish–Lenski".

**In target journals, tangential**

- *Genetics* (2022) 220(1): iyab134. Polygenic adaptation under global vs local adaptation; repeatability of the loci involved depends on genetic redundancy. Sexual populations with migration, no clonal interference. <https://academic.oup.com/genetics/article/220/1/iyab134/6497714>
- *Genetics* (2024) 228(1): iyae111. Pleiotropy, migration and repeated genetic adaptation; large-effect loci can still contribute repeatedly despite pleiotropic trade-offs. No interference. <https://academic.oup.com/genetics/article/228/1/iyae111/7712980>
- *Genetics* (2019). Interference effects of deleterious and beneficial mutations in large asexual populations; fixation probability under combined interference. Just outside the 5-year window. <https://doi.org/10.1534/genetics.119.301960>

**Closest in substance, other journals**

- Logarithmic scaling of selective sweep curves: from tents to houses, preprint (submitted to JOMB), (2026) with F. Boenkost,  Felix Hermann, A. Tóbiás, A. Wakolbinger  
https://arxiv.org/abs/2605.10617

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
