# Numerical Analysis Repository Organization Design

> **Date:** 2026-03-26
> **Repository:** 26Spring_Numerical-Analysis
> **Status:** Approved

---

## Overview

Organize the `26Spring_Numerical-Analysis` repository with a two-level hierarchy: 6 Part directories (academic taxonomy from the textbook) containing 24 Chapter subdirectories, each with concept and lab markdown files.

## Source Materials

- **Lecture PDFs (N00-N28):** Located in iCloud `[강의자료]/`
- **Lab Notebooks (.ipynb):** Located in iCloud `[공통]/course305-main/`
- **Textbook:** "Applied Numerical Methods with Python for Engineers and Scientists" by Chapra & Clough (2022, McGraw Hill)

## Directory Structure

```
26Spring_Numerical-Analysis/
├── Part1_Modeling-Computers-and-Error-Analysis/
│   ├── Chapter01_Mathematical-Modeling/
│   │   ├── Concepts_Lecture.md
│   │   └── Concepts_Lab.md
│   ├── Chapter02_Programming-and-Software/
│   │   ├── Concepts_Lecture.md
│   │   └── Concepts_Lab.md
│   ├── Chapter03_Approximations-and-Round-Off-Errors/
│   │   ├── Concepts_Lecture.md
│   │   └── Concepts_Lab.md
│   └── Chapter04_Truncation-Errors/
│       └── Concepts_Lecture.md
├── Part2_Roots-and-Optimization/
│   ├── Chapter05_Roots-Bracketing-Methods/
│   ├── Chapter06_Roots-Open-Methods/
│   └── Chapter07_Optimization/
├── Part3_Linear-Algebraic-Equations/
│   ├── Chapter08_Linear-Algebraic-Equations/
│   │   ├── Concepts_Lecture.md
│   │   └── Concepts_Lab.md
│   ├── Chapter09_Gauss-Elimination/
│   │   ├── Concepts_Lecture.md
│   │   └── Concepts_Lab.md
│   ├── Chapter10_LU-Factorization/
│   │   ├── Concepts_Lecture.md
│   │   └── Concepts_Lab.md
│   ├── Chapter11_Matrix-Inverse-and-Condition/
│   ├── Chapter12_Eigenvalues-Power-Method/
│   └── Chapter13_Eigenvalues-Symmetric-Matrices/
├── Part4_Curve-Fitting/
│   ├── Chapter14_Linear-Regression/
│   ├── Chapter15_Nonlinear-Regression/
│   ├── Chapter16_Fourier-Analysis/
│   ├── Chapter17_Interpolation/
│   └── Chapter18_Splines/
├── Part5_Integration-and-Differentiation/
│   ├── Chapter19_Numerical-Integration-Formulas/
│   ├── Chapter20_Numerical-Integration-of-Equations/
│   └── Chapter21_Numerical-Differentiation/
├── Part6_Ordinary-Differential-Equations/
│   ├── Chapter22_Initial-Value-Problems/
│   ├── Chapter23_Stiff-ODEs/
│   └── Chapter24_Boundary-Value-Problems/
├── images/
├── LICENSE
├── README.ko.md
└── README.md
```

## File Conventions

### Concepts_Lecture.md
- Based on lecture PDFs
- Systematic theory summary in English
- Structure: Table of Contents, sectioned by topic, key formulas, examples

### Concepts_Lab.md
- Based on course305 .ipynb notebooks
- Concepts explained more thoroughly than the original notebook (step-by-step, detailed)
- All code included as markdown code blocks (not raw .ipynb)
- English first; Korean version created separately after English review

### Korean Versions
- `Concepts_Lecture.ko.md` and `Concepts_Lab.ko.md` created after English versions are reviewed and approved
- Content must be identical in substance to English versions

### Assignments
- Excluded for now

## Lab-to-Chapter Mapping

| course305 Notebook | Target Chapter |
|---|---|
| `01-PythonBasic/Lab01-python-fundamentals.ipynb` | Chapter01 |
| `01-PythonBasic/Lab01-python-programming.ipynb` | Chapter02 |
| `02-SourceOfErrors/Lab02-SourceOfErrors.ipynb` | Chapter03 |
| `03-LinearSystemOfEquations/Lab03-LinearSystem.ipynb` | Chapter08 |
| `04-LUfactorization/Lab04-LU.ipynb` | Chapter10 |

## Week-to-Chapter Mapping (for README tables)

| Week | Lecture Notes | Part / Chapter |
|:---:|:---|:---|
| 1 | N00-N03 | Part1: Ch01, Ch02, Ch03 |
| 2 | N04 | Part1: Ch04 |
| 3 | N05-N06 | Part3: Ch08, Ch09 |
| 4 | N07-N08 | Part3: Ch10, Ch11 |
| 5 | N09-N10 | Part3: Ch12, Ch13 |
| 6 | N11-N12 | Part2: Ch05, Ch06 |
| 7 | N13-N14 | Part2: Ch07 |
| 8 | Midterm | — |
| 9 | N15-N16 | Part4: Ch14, Ch15 |
| 10 | N17-N18 | Part4: Ch16 |
| 11 | N19-N20 | Part4: Ch17, Ch18 |
| 12 | N21-N22 | Part5: Ch19, Ch20 |
| 13 | N23-N24 | Part5: Ch21 / Part6: Ch22 |
| 14 | N25-N26 | Part6: Ch22 |
| 15 | N27-N28 | Part6: Ch23, Ch24 |
| 16 | Final | — |

## Scope of Initial Work

**Directories:** All 6 Parts, all 24 Chapters, images/ — created immediately.

**Content (weeks 1-4 covered so far, as of 2026-03-26):**
- Part1 Ch01-Ch03: Concepts_Lecture.md + Concepts_Lab.md
- Part1 Ch04: Concepts_Lecture.md only (no dedicated lab notebook for this chapter)
- Part3 Ch08-Ch11: Concepts_Lecture.md + Concepts_Lab.md

**Note:** Chapters with multiple lecture sessions (e.g., Ch07 spans N13-N14, Ch22 spans N24-N26) produce a single unified Concepts_Lecture.md per chapter.

**Deferred:** All other chapters — directory only, content added as course progresses.

**README:** Updated with new structure and week-to-chapter mapping table.
