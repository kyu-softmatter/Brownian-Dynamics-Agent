<!-- 생성됨: docs/tools/wiki_index.py — 직접 고치지 말 것 -->
# 발표된 문헌 — 색인

우리 랩 발표 논문(`lab_authored: true`)과 외부 문헌. **발표 여부로 가른 폴더다** — 미발표 랩 자산은 `../lab/`(gitignore)에 있다.

| | |
|---|---|
| 항목 수 | **40** |
| 원문 확보 (`raw_file`) | 22 / 40 |
| 파라미터 추출 완료 | **11 / 40** |
| 우리가 재현함 (`reproduced: yes`) | **0 / 40** |

> `parameters_extracted: no` 인 항목은 **수집 기록**이지 근거가 아니다.
> `reproduced: no` 인 파라미터를 문헌 근거처럼 인용하지 않는다 — 계약: [`knowledge/wiki/CLAUDE.md`](../../wiki/CLAUDE.md)

## 관련도 `direct` — BD/MC 시뮬레이션 직결

| 연도 | 논문 | 엔진 | 원문 | 추출 | 재현 |
|---|---|---|---|---|---|
| 2025 | [Choi KH — Non-equilibrium Dynamics in Bio-inspired Soft Matter (박사학위논문)](2025-choi-phd-thesis-noneq-bioinspired-soft-matter.md) | HOOMD-blue (BD, GPU) | ✅ | ✅ | — |
| 2025 | [Quah T et al. — Learning continuum-level closures for control of interacting active p...](2025-quah-continuum-closures-active-control.md) | ABP BD (agent-based) + 신경 연산자 UDE | ✅ | ✅ | — |
| 2025 | [Takatori SC et al. — Feedback control of active matter](2025-takatori-feedback-control-active-matter.md) | 리뷰 | 🔒 | — | — |
| 2024 | [Barakat JM et al. — Surface topography induces and orients nematic swarms of active filam...](2024-barakat-surface-topography-nematic-swarms.md) | 실험 + 시뮬 | — | — | — |
| 2024 | [Cheon J et al. — Motility modulates the partitioning of bacteria in aqueous two-phase ...](2024-cheon-motility-partitioning-atps.md) | 실험 (광집게) + BD | ✅ | — | — |
| 2024 | [Modica KJ et al. — Soft confinement of self-propelled rods: simulation and theory](2024-modica-soft-confinement-self-propelled-rods.md) | BD | 🔒 | — | — |
| 2024 | [Quah T et al. — graybox_abp_mpc — Learning continuum-level closures for control of in...](2024-quah-graybox-abp-mpc-repo.md) | HOOMD-blue 3.8.1 (BD, ABP) | ✅ | ✅ | — |
| 2024 | [Quah T et al. — Model predictive control of non-interacting active Brownian particles](2024-quah-mpc-noninteracting-abp.md) | BD + MPC | 🔒 | — | — |
| 2024 | [Xu Y et al. — Dynamic surfactants drive anisotropic colloidal assembly](2024-xu-dynamic-surfactants-anisotropic-assembly.md) | HOOMD-blue (BD, GPU) | ✅ | — | — |
| 2023 | [Modica KJ et al. — Boundary design regulates the diffusion of active matter in heterogen...](2023-modica-boundary-design-active-diffusion.md) | BD + 분산이론 | 🔒 | — | — |
| 2023 | [Nagella SG et al. — Colloidal transport phenomena in dynamic, pulsating porous materials](2023-nagella-colloidal-transport-pulsating-porous.md) | 이론 + 시뮬 | 🔒 | — | — |
| 2023 | [Xu Y et al. — Dynamic interfaces for contact-time control of colloidal interactions](2023-xu-dynamic-interfaces-contact-time.md) | HOOMD-blue (BD, GPU) | ✅ | ✅ | — |
| 2023 | [Xu Y et al. — Nonequilibrium interactions between multi-scale colloids regulate the...](2023-xu-nonequilibrium-multiscale-colloids.md) | BD | 🔒 | — | — |
| 2022 | [Barakat JM et al. — Enhanced dispersion in an oscillating array of harmonic traps](2022-barakat-enhanced-dispersion-harmonic-traps.md) | HOOMD-blue (BD) + COMSOL (연속체) | ✅ | ✅ | — |
| 2022 | [Modica KJ et al. — Porous media microstructure determines the diffusion of active matter...](2022-modica-porous-media-active-diffusion.md) | BD (ABP) + Janus 입자 실험 | ✅ | — | — |
| 2020 | [Takatori SC et al. — Motility-induced buckling and glassy dynamics regulate three-dimensio...](2020-takatori-motility-induced-buckling.md) | BD (ABP, 2D 이분산) | ✅ | — | — |
| 1999 | [Zahn K et al. — Two-Stage Melting of Paramagnetic Colloidal Crystals in Two Dimensions](1999-zahn-two-stage-melting-2d.md) | 실험 (디지털 비디오 현미경) | ✅ | ✅ | — |

## 관련도 `method` — 검증 대조에 쓸 수 있는 측정

| 연도 | 논문 | 엔진 | 원문 | 추출 | 재현 |
|---|---|---|---|---|---|
| 2026 | [Choi KH et al. — Anisotropic diffusion in lyotropic chromonic liquid crystal using flu...](2026-choi-anisotropic-diffusion-lclc-frap.md) | 실험 (FRAP) | ✅ | — | — |
| 2025 | [Kim DY et al. — Colloidal hydrodynamic interactions in viscoelastic fluids](2025-kim-colloidal-hydrodynamics-viscoelastic.md) | 실험 (광집게) + 점탄성 이론 | ✅ | — | — |
| 2025 | [Kim DY et al. — Direct experimental measurement of many-body hydrodynamic interaction...](2025-kim-manybody-hydrodynamics-optical-tweezers.md) | 실험 (광집게) | ✅ | ✅ | — |
| 2023 | [Lee HM et al. — Direct measurements of the colloidal Debye force](2023-lee-colloidal-debye-force.md) | 실험 (광집게) | ✅ | — | — |
| 2020 | [Choi KH et al. — Interpretation of electrostatic self-potential measurements using int...](2020-choi-electrostatic-self-potential-heterogeneity.md) | 실험 (광집게) + MC 시뮬레이션 | ✅ | ✅ | — |
| 2019 | [Choi KH et al. — Direct measurement of electrostatic interactions between poly(methyl ...](2019-choi-electrostatic-pmma-optical-tweezers.md) | 실험 (광집게) | 🔒 | — | — |
| 2019 | [Kang DW et al. — Mapping anisotropic and heterogeneous colloidal interactions via opti...](2019-kang-mapping-anisotropic-colloidal-interactions.md) | 실험 (광집게) | 🔒 | — | — |
| 2008 | [Park BJ et al. — Direct measurements of the effects of salt and surfactant on interact...](2008-park-salt-surfactant-interface-forces.md) | 실험 (광집게) | ✅ | ✅ | — |

## 관련도 `adjacent` — 랩 산출물, 직접 관련은 약함

| 연도 | 논문 | 엔진 | 원문 | 추출 | 재현 |
|---|---|---|---|---|---|
| 2026 | [Choi KH et al. — Mechanics of heterogeneous fiber networks](2026-choi-heterogeneous-fiber-networks.md) | 실험 (능동 미세유변학 + 형광 이미징) | ✅ | — | — |
| 2025 | [Gubbala A et al. — Phase field model for viscous inclusions in anisotropic networks](2025-gubbala-phase-field-viscous-inclusions.md) | 연속체 (Cahn–Hilliard + Landau–de Gennes) | ✅ | ✅ | — |
| 2024 | [Arnold DP et al. — Lipid membrane domains control actin network viscoelasticity](2024-arnold-lipid-domains-actin-viscoelasticity.md) | 실험 (미세유변학 + 형광) | ✅ | — | — |
| 2024 | [Gubbala A et al. — Dynamic swarms regulate the morphology and distribution of soft membr...](2024-gubbala-dynamic-swarms-membrane-domains.md) | 연속체 (Toner–Tu + Cahn–Hilliard) | ✅ | — | — |
| 2024 | [Quah T et al. — Neural network augmented model predictive control: application to act...](2024-quah-nn-augmented-mpc.md) | BD + MPC | 🔒 | — | — |
| 2023 | [Arnold DP et al. — Active surface flows accelerate the coarsening of lipid membrane domains](2023-arnold-active-surface-flows-coarsening.md) | 실험 + 해석 (Smoluchowski 응집 + phase field) | ✅ | ✅ | — |
| 2023 | [Arnold DP et al. — Bio-enabled engineering of multifunctional "living" surfaces](2023-arnold-bioenabled-living-surfaces.md) | 리뷰 | 🔒 | — | — |
| 2023 | [Jung IH et al. — Quantification of polystyrene microsphere attachment probability at t...](2023-jung-microsphere-attachment-microfluidic.md) | 실험 | — | — | — |
| 2022 | [Arnold DP et al. — Antibody binding reports spatial heterogeneities in cell membrane org...](2022-arnold-antibody-binding-heterogeneities.md) | 실험 + 모델 | ✅ | — | — |
| 2022 | [Jeong HW et al. — Retardation of capillary force between Janus particles at the oil-wat...](2022-jeong-capillary-retardation-janus.md) | 실험 | 🔒 | — | — |
| 2020 | [Choi KH et al. — Interfacial configurations of lens-shaped particles](2020-choi-interfacial-configurations-lens-shaped.md) | 실험 | 🔒 | — | — |
| 2020 | [Choi KH et al. — Interpretation of interfacial interactions between lenticular particles](2020-choi-lenticular-particle-interactions.md) | 실험 | 🔒 | — | — |
| 2020 | [Lee HE et al. — Interactions between polystyrene particles with diameters of several ...](2020-lee-polystyrene-oil-water-interactions.md) | 실험 | 🔒 | — | — |
| 2018 | [Kang DW et al. — Geometric effects of colloidal particles on stochastic interface adso...](2018-kang-geometric-effects-stochastic-adsorption.md) | 실험 | 🔒 | — | — |
| 2018 | [Lim JH et al. — Heterogeneous capillary interactions of interface-trapped ellipsoid p...](2018-lim-heterogeneous-capillary-ellipsoid.md) | 실험 | 🔒 | — | — |

---

**범례** — 원문: ✅ 확보 · 🔒 paywall(기관 접근권으로 `knowledge/raw/`에 배치) · — 없음
  ·  추출/재현: ✅ 완료 · ◐ 부분 · — 아직
