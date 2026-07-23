---
title: 'ESRIcascade.jl: Economic Systemic Risk Index cascade computation in Julia'
tags:
  - Julia
  - supply chains
  - systemic risk
  - input-output networks
  - economic networks
authors:
  - name: Mitja Devetak
    affiliation: "1, 2"
  - name: Peter Klimek
    affiliation: "1, 2, 3, 4"
affiliations:
  - name: Section for the Science of Complex Systems, CeMSIIS, Medical University of Vienna, Vienna, Austria
    index: 1
  - name: Complexity Science Hub Vienna, Vienna, Austria
    index: 2
  - name: Supply Chain Intelligence Institute Austria, Vienna, Austria
    index: 3
  - name: Division of Insurance Medicine, Karolinska Institutet, Department of Clinical Neuroscience, Karolinska, Sweden
    index: 4
date: 20 May 2026
bibliography: paper.bib
repository: https://github.com/Devetak/ESRIcascade.jl
---

# Summary

ESRIcascade.jl is a Julia [@bezanson2017julia] implementation of the Economic Systemic Risk Index (ESRI) algorithm introduced by Diem et al. [@diem2022esri]. Given the production network of an economy, the package computes for each firm the share of total output that is lost after a scenario-level shock has propagated through the supply chain. When the scenario involves the complete shutdown of a single firm, this is known as the ESRI score of that firm, a measure of systemic importance. The package is intended for researchers who work with firm-level supply networks, reconstructed production networks, or synthetic input-output networks. These data are increasingly used to study supply-chain risk [@pichler2023alliance]. ESRIcascade.jl provides a package interface for economy-wide ESRI computations via single-firm shock scenarios and more complex capacity shocks.

# Statement of need

The Economic Systemic Risk Index was introduced to measure how much economy-wide output depends on each firm in a production network [@diem2022esri]. The method has since been used in work on reconstructed supply networks, economic predictability, financial stability, decarbonization, food security, and network rewiring [@reisch2022mobile; @diem2024predictability; @tabachova2024financial; @stangl2024decarbonization; @fessina2026inferring; @reisch2026rewiring; @zelbi2026mitigation]. These applications require repeated cascade computations over many firms, many network variants, or many shock scenarios.

The original ESRI implementation is distributed as scripts around R/C++ cascade code [@fastcascade; @esri_tutorial]. ESRIcascade.jl packages the core ESRI cascade calculation as a Julia library with documented functions, tests, continuous integration, and examples. It is designed for use inside Julia analysis pipelines and for reproducible comparisons across empirical and reconstructed networks.

This package makes ESRI available through a tested Julia API and exposes scenario-level calls that can be reused in research workflows. The supported interface includes whole-economy ESRI scores, individual firm closures, explicit capacity vectors for modeling complex scenarios, optional upstream or downstream components, and custom final aggregation weights.

# State of the field

Firm-level production networks have become central objects for studying supply-chain fragility because aggregate input-output tables can hide the concentration of critical suppliers and customers inside broad sectors [@pichler2023alliance; @diem2024predictability]. ESRI is one of the established measures for this setting: it converts a weighted firm-to-firm network into firm-level loss estimates by propagating a shutdown through upstream and downstream production dependencies [@diem2022esri]. This makes it useful not only for ranking firms by systemic importance, but also for comparing reconstructed networks, alternative network models, and counterfactual shock scenarios.

The public ESRI software landscape has so far been centered on reference and tutorial code distributed around the original R/C++ workflow [@fastcascade; @esri_tutorial]. That code is valuable for reproducing the method, but it is not organized as a reusable Julia package with a stable typed interface, package tests, documentation, and scenario-level composition points. In practice, this makes extension and repeated use harder: researchers who want to vary networks, reuse a prepared economy, inspect upstream and downstream channels, or integrate ESRI into Julia-based simulation and reconstruction pipelines need to adapt scripts rather than call a package API.

Computationally, economy-wide ESRI is expensive because one fixed-point cascade must be solved for each candidate firm. Repeated scenario evaluation therefore becomes the main bottleneck in reconstruction studies, sensitivity analysis, and comparisons across alternative networks. ESRIcascade.jl addresses this bottleneck as software engineering work around an existing method, not as a change to the ESRI definition.

# Software design

ESRIcascade.jl is organized around the workflow used in ESRI studies. A user starts from a weighted firm-to-firm network and an industry label for each firm. `IndustryInfo` stores these labels together with either a global essentiality flag for each supplier industry or a customer-specific input-classification matrix; omitting a classification gives the all-essential default. This follows the ESRI production setup, where essential inputs enter through a Leontief part and inessential inputs enter through a linear part [@diem2022esri]. `ESRIEconomy` then stores the normalized upstream and downstream operators, output weights, and network totals. This setup step is separate from the scenario calculation, so the same economy object can be reused for many firm closures, counterfactual networks, or partial shock vectors.

The main user calls compute ESRI scores for all selected firms, solve one firm-closure scenario, or solve a scenario from an explicit capacity vector, where one means normal operation, zero means closure, and intermediate values mean partial capacity. The result can be a single ESRI value or an `ESRIResult` with the converged upstream and downstream states.

The implementation targets sparse supply-chain networks with broad, power-law-like firm degree distributions, where most firms have few links and a small number of firms have many. This is a common property of real-world production networks [@reisch2026rewiring; @bacilieri2023partial]. For full-economy sparse runs without threading, ESRIcascade.jl sorts firms by degree before solving and then maps the scores back to the original firm order. This puts high-degree firms early in the scenario sequence, which is useful for the heavy-tailed graphs the package is intended to handle. Each scenario is solved by iterating the upstream and downstream states together until both have converged. The upstream step traverses the stored customer links needed for each firm, computing the loss of sales. The downstream step computes the loss of production due to missing inputs. One workspace is allocated for a serial economy-wide run and reused across its scenarios.

Compared with the original script-oriented R/C++ workflow, ESRIcascade.jl exposes the cascade calculation as a Julia package with a reusable API. Both implementations can distribute the independent firm-closure scenarios across workers, since the economy-wide ESRI calculation is embarrassingly parallel. ESRIcascade.jl adds a package-level interface for the parts that researchers often need to vary: the prepared economy object, the shock vector, the selected firms, the final weights, and the returned upstream and downstream states.

The implementation also differs in the intermediate steps used inside one scenario solve. The R/C++ comparison path expresses downstream propagation through matrix batches and forms scenario-level supplier-industry intermediates during the iteration. ESRIcascade.jl implements the same update as direct accumulation over observed supplier-customer links into preallocated workspaces. It still uses an essential-loss workspace by firm and industry, but it avoids rebuilding the larger matrix intermediates for each scenario step. The upstream and downstream states are advanced in the same iteration loop until both have converged.

These implementation changes preserve the ESRI fixed-point updates and final loss definition while reorganizing the same operations into a package-oriented sparse Julia implementation.

# Benchmark

The benchmark compares ESRIcascade.jl with the Diem et al. R/C++ implementation [@fastcascade] on the same full-economy task. For each network, every firm is closed once and the complete ESRI score vector is computed. The tested networks are sparse synthetic power-law graphs, with the number of links growing in proportion to the number of firms. This gives the benchmark the degree imbalance expected in firm-level supply networks without using confidential firm data.

Both implementations used 12 threads or workers. The stopping tolerance was `0.01` for both implementations, and the same generated networks were passed to both solvers. The benchmark machine used an Intel Core Ultra 7 155H CPU with 16 cores, 22 logical CPUs, and 24 MiB L3 cache.

The two implementations matched one-to-one up to numerical precision for this benchmark. Across all reported runs, the largest absolute difference between paired ESRI scores was at most `3e-6`. No nonfinite scores were produced by the R/C++ implementation. On these full-firm sparse graphs, ESRIcascade.jl was about 5.6 to 6.5 times faster.

| Number of firms | ESRIcascade.jl time | Diem et al. implementation time | Speedup |
| ---: | ---: | ---: | ---: |
| 10,000 | 8.29 s | 47.17 s | 5.69x |
| 20,000 | 35.25 s | 197.73 s | 5.61x |
| 50,000 | 208.85 s | 1362.74 s | 6.53x |
| 100,000 | 1016.20 s | 6602.74 s | 6.50x |

## Customer-specific classification study

The package also supports the reference R/C++ convention in which input type depends on both the supplying and purchasing industries. The supplied 616 by 616 IHS classification matrix uses rows for supplier industries and columns for customer industries. Each entry is `2` for an essential input, `1` for a non-essential input, or `0` for an input with no short-term downstream production impact. Zero-classified links remain available to upstream shock propagation.

A compatibility study using this classification matrix evaluated 16 firm-closure scenarios on a network with 1,232 firms and 9,856 edges. The maximum absolute difference between the Julia and C++ ESRI results was `5.995e-15`. The Julia calculation took `0.287 s`, compared with `4.899 s` for C++, corresponding to a `17.1x` speedup.

A separate prepared-economy run on the same 16 scenarios reduced the median serial solve from about `39 ms` to `5 ms` (about `8x`) by reusing the scenario workspace and reducing only touched sparse essential-input cells. The scores were unchanged exactly. This isolates the improvement in the cascade solver from CSV parsing and economy construction.

# Limitations

ESRIcascade.jl does not change the ESRI algorithm itself. A full economy-wide ESRI vector still requires one cascade solve per candidate firm, and each cascade solve still iterates propagation over the network until convergence. The package improves the software interface, extensibility, memory behavior, and constant factors, but it does not remove the fundamental scaling bottleneck of repeated all-firm cascade evaluation. For sparse benchmark families where the number of links grows proportionally with the number of firms, this means the full all-firm calculation still has the same approximately quadratic scaling in the number of edge-level dependencies as the underlying ESRI workflow.

# Research impact statement

ESRIcascade.jl supports research workflows where ESRI must be computed repeatedly on firm-level production networks. This includes empirical stress tests, sensitivity analysis over reconstructed networks, comparison of network inference methods, and repeated scenario evaluation. The package is used in industry-aware firm-level network reconstruction work to compare empirical and reconstructed networks with ESRI-based outcomes [@devetak2026industry]. The package improves the reproducibility of ESRI-based analyses by providing a documented package interface, deterministic examples, public tests, and benchmark scripts. It also gives researchers access to the upstream and downstream components of a scenario, which helps inspect whether losses arise mainly through missing inputs or lost customers.

# Acknowledgements

No external funding was received for this work.

# AI usage disclosure

OpenAI Codex, using ChatGPT 5.4 and ChatGPT 5.5, was used during preparation of this JOSS submission for paper drafting and editing, bibliography scaffolding, and consistency checks against the repository, documentation, benchmark records, and cited papers. The first version of ESRIcascade.jl was written by the author without generative AI. After that first version, Codex was used to draft scaffolding for parts of the package code, tests, documentation, benchmark scripts, and paper files. Codex did not make the core design decisions. All AI-assisted outputs were extensively audited by human authors, who reviewed, edited, and validated the text, citations, code claims, benchmark claims, and licensing statements. The authors remain responsible for the accuracy, originality, licensing, and compliance of the submitted materials.

# References
