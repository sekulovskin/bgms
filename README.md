
# bgms: Bayesian Analysis of Networks of Binary and/or Ordinal Variables

The `R` package <strong>bgms</strong> provides tools for Bayesian
analysis of the ordinal Markov random field, a graphical model
describing a network of binary and/or ordinal variables (Marsman, van
den Bergh, and Haslbeck 2025). A pseudolikelihood is used to approximate
the likelihood of the graphical model, and Markov chain Monte Carlo
methods are used to simulate from the corresponding pseudoposterior
distribution of the graphical model parameters.

The <strong>bgm</strong> function can be used for a one-sample design
and the <strong>bgmCompare</strong> function can be used for an
independent-sample design (see Marsman et al. 2024). Both functions can
model the selection of effects. In one-sample designs, the
<strong>bgm</strong> function models the presence or absence of edges
between pairs of variables in the network. The estimated posterior
inclusion probability indicates how plausible it is that a network with
an edge between the two corresponding variables produced the observed
data, and can be converted into a Bayes factor test for conditional
independence. The <strong>bgm</strong> function can also model the
presence or absence of communities or clusters of variables in the
network. The estimated posterior probability distribution of the number
of clusters indicates how plausible it is that a network with the
corresponding number of clusters produced the observed data, and can be
converted into a Bayes factor test for clustering (see Sekulovski et al.
2025).

In an independent-sample design, the <strong>bgmCompare</strong>
function models the selection of group differences in edge weights and
possibly category thresholds. The estimated posterior inclusion
probability indicates how plausible it is that graphical models with a
difference in the corresponding edge weight or category threshold
generated the data at hand, and can be converted to a Bayes factor test
for parameter equivalence.

## Why use Markov Random Fields?

Multivariate analysis using graphical models has received much attention
in the recent psychological and psychometric literature (Robinaugh et
al. 2020; Marsman and Rhemtulla 2022; Contreras et al. 2019). Most of
these graphical models are Markov Random Field (MRF) models, whose graph
structure reflects the partial associations between variables
(Kindermann and Snell 1980). In these models, a missing edge between two
variables in the network implies that these variables are independent,
given the remaining variables (Lauritzen 2004). In other words, the
remaining variables of the network fully account for the potential
association between the unconnected variables.

## Why use a Bayesian approach to analyze the MRF?

Testing the structure of the MRF in a one-sample design requires us to
determine the plausibility of the opposing hypotheses of conditional
dependence and conditional independence. That is, how plausible is it
that the observed data come from a network with a structure that
includes the edge between two variables compared to a network structure
that excludes that edge? Similarly, testing for group differences in the
MRF in an independent-samples design requires us to determine the
plausibility of the opposing hypotheses of parameter difference and
parameter equivalence. That is, how plausible is it that the observed
data come from MRFs with differences in specific edge weights or
threshold parameters compared to MRFs that do not differ in these
parameter?

Frequentist approaches are limited in this respect because they can only
reject, not support, null hypotheses of conditional independence or
parameter equivalence. This leads to the problem that if an edge is
excluded, we do not know whether this is because the edge is absent in
the population or because we do not have enough data to reject the null
hypothesis of independence. Similarly, if a difference is excluded, we
do not know whether this is because there is no difference in the
parameter between the different groups or because we do not have enough
data to reject the null hypothesis of parameter equivalence.

To avoid this problem, we will advocate a Bayesian approach using Bayes
factors. In one-sample designs, the inclusion Bayes factor (Huth et al.
2023; Sekulovski et al. 2024) allows us to quantify how much the data
support both conditional dependence -<em>evidence of edge
presence</em> - or conditional independence -<em>evidence of edge
absence</em>. It also allows us to conclude that there is limited
support for either hypothesis - an <em>absence of evidence</em>. In
independent-sample designs, they can be used to quantify how much the
data support the hypotheses of parameter difference and equivalence. The
output of the <strong>bgm</strong> and <strong>bgmCompare</strong>
functions can be used to estimate these inclusion Bayes factors.

## Installation

The current developmental version can be installed with

``` r
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   
remotes::install_github("MaartenMarsman/bgms")
```

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-ContrerasEtAl_2019" class="csl-entry">

Contreras, A., I. Nieto, C. Valiente, R. Espinosa, and C. Vazquez. 2019.
“The Study of Psychopathology from the Network Analysis Perspective: A
Systematic Review.” *Psychotherapy and Psychosomatics* 88: 71–83.
<https://doi.org/10.1159/000497425>.

</div>

<div id="ref-HuthEtAl_2023_intro" class="csl-entry">

Huth, K., J. de Ron, A. E. Goudriaan, K. Luigjes, R. Mohammadi, R. J.
van Holst, E.-J. Wagenmakers, and M. Marsman. 2023. “Bayesian Analysis
of Cross-Sectional Networks: A Tutorial in R and JASP.” *Advances in
Methods and Practices in Psychological Science* 6: 1–18.
<https://doi.org/10.1177/25152459231193334>.

</div>

<div id="ref-KindermannSnell1980" class="csl-entry">

Kindermann, R., and J. L. Snell. 1980. *Markov Random Fields and Their
Applications*. Vol. 1. Contemporary Mathematics. Providence: American
Mathematical Society.

</div>

<div id="ref-Lauritzen2004" class="csl-entry">

Lauritzen, S. L.. 2004. *Graphical Models*. Oxford: Oxford University
Press.

</div>

<div id="ref-MarsmanRhemtulla_2022_SIintro" class="csl-entry">

Marsman, M., and M. Rhemtulla. 2022. “Guest Editors’ Introduction to the
Special Issue ‘Network Psychometrics in Action’: Methodological
Innovations Inspired by Empirical Problems.” *Psychometrika* 87: 1–11.
<https://doi.org/10.1007/s11336-022-09861-x>.

</div>

<div id="ref-MarsmanVandenBerghHaslbeck_2024" class="csl-entry">

Marsman, M., D. van den Bergh, and J. M. B. Haslbeck. 2025. “Bayesian
Analysis of the Ordinal Markov Random Field.” *Psychometrika* 90:
146--182.

</div>

<div id="ref-MarsmanWaldorpSekulovskiHaslbeck_2024" class="csl-entry">

Marsman, M., L. J. Waldorp, N. Sekulovski, and J. M. B. Haslbeck. 2024.
“A Bayesian Independent Samples $t$ Test for Parameter Differences in
Networks of Binary and Ordinal Variables.” *Retrieved from
Https://Osf.io/Preprints/Osf/F4pk9*.

</div>

<div id="ref-RobinaughEtAl_2020" class="csl-entry">

Robinaugh, D. J., R. H. A. Hoekstra, E. R. Toner, and D. Borsboom. 2020.
“The Network Approach to Psychopathology: A Review of the Literature
2008–2018 and an Agenda for Future Research.” *Psychological Medicine*
50: 353–66. <https://doi.org/10.1017/S0033291719003404>.

</div>

<div id="ref-SekulovskiEtAl_2025" class="csl-entry">

Sekulovski, N., G. Arena, J. M. B. Haslbeck, K. B. S. Huth, N. Friel,
and M. Marsman. 2025. “A Stochastic Block Prior for Clustering in
Graphical Models.” *Retrieved from
<a href="https://osf.io/preprints/psyarxiv/29p3m_v1"
class="uri">Https://Osf.io/Preprints/Psyarxiv/29p3m_v1</a>*.

</div>

<div id="ref-SekulovskiEtAl_2024" class="csl-entry">

Sekulovski, N., S. Keetelaar, K. B. S. Huth, Eric-Jan Wagenmakers, R.
van Bork, D. van den Bergh, and M. Marsman. 2024. “Testing Conditional
Independence in Psychometric Networks: An Analysis of Three Bayesian
Methods.” *Multivariate Behavioral Research* 59: 913–33.
<https://doi.org/10.1080/00273171.2024.2345915>.

</div>

</div>
