# [Spring 2026] Numerical Analysis

![Last Commit](https://img.shields.io/github/last-commit/Choroning/26Spring_Numerical-Analysis)
![Languages](https://img.shields.io/github/languages/top/Choroning/26Spring_Numerical-Analysis)

This repository organizes and stores numerical method implementations written for university lectures and assignments.

*Author: Cheolwon Park (Korea University Sejong, CSE) – Year 3 (Junior) as of 2026*
<br><br>

## 📑 Table of Contents

- [About This Repository](#about-this-repository)
- [Course Information](#course-information)
- [Prerequisites](#prerequisites)
- [Weekly Schedule](#weekly-schedule)
- [Repository Structure](#repository-structure)
- [License](#license)

---


<br><a name="about-this-repository"></a>
## 📝 About This Repository

This repository contains bilingual study materials and code developed for a university-level Numerical Analysis course, including:

- Bilingual Concepts notes (Korean `.ko.md` + English `.md`) for every lecture and lab session
- Python implementations of numerical algorithms from the textbook
- Part/Chapter directory structure covering the full Chapra & Canale curriculum

> **🤖 AI-Assisted Development**
> This course encourages the use of AI agents.
> [Claude Code](https://claude.ai/download) and [Gemini CLI](https://github.com/google-gemini/gemini-cli) were used as coding assistants throughout the course.

<br><a name="course-information"></a>
## 📚 Course Information

- **Semester:** Spring 2026 (March - June)
- **Affiliation:** Korea University Sejong

| Course&nbsp;Code| Course            | Type          | Instructor      | Department                              |
|:----------:|:------------------|:-------------:|:---------------:|:----------------------------------------|
|`DCSS305-00`|NUMERICAL ANALYSIS (English)|Major Elective|Prof. Shinhoo&nbsp;Kang|Department of Computer Software|

- **📖 References**

| Type | Contents |
|:----:|:---------|
|Textbook|"Numerical Methods for Engineers" by Steven C. Chapra and Raymond P. Canale (8th Edition, McGraw Hill)|
|Lecture Notes|Slides and online video lectures provided by instructor|

<br><a name="prerequisites"></a>
## ✅ Prerequisites

- Understanding of Linear Algebra and Discrete Mathematics
- Python interpreter installed
- Familiarity with scientific computing libraries

- **💻 Development Environment**

| Tool | Company |  OS  | Notes |
|:-----|:-------:|:----:|:------|
|Python 3|Python Software Foundation|macOS|    |
|NumPy|NumFOCUS|macOS|Numerical computing|
|SciPy|NumFOCUS|macOS|Scientific computing|
|Matplotlib|NumFOCUS|macOS|Visualization|
|JupyterLab|Project Jupyter|macOS|Interactive notebooks|

<br><a name="weekly-schedule"></a>
## 📅 Weekly Schedule

| Week | Lecture Notes | Part / Chapter | Topic |
|:---:|:---|:---|:---|
| 1 | N00-N03 | Part1: Ch01, Ch02, Ch03 | Modeling, Python, Round-Off Errors |
| 2 | N04 | Part1: Ch04 | Truncation Errors |
| 3 | N05-N06 | Part3: Ch08, Ch09 | Linear Systems, Gauss Elimination |
| 4 | N07-N08 | Part3: Ch10, Ch11 | LU Factorization, Matrix Inverse |
| 5 | N09-N10 | Part3: Ch12, Ch13 | Eigenvalues |
| 6 | N11-N12 | Part2: Ch05, Ch06 | Root Finding |
| 7 | N13-N14 | Part2: Ch07 | Optimization |
| 8 | Midterm | — | — |
| 9 | N15-N16 | Part4: Ch14, Ch15 | Regression |
| 10 | N17-N18 | Part4: Ch16 | Fourier Analysis |
| 11 | N19-N20 | Part4: Ch17, Ch18 | Interpolation, Splines |
| 12 | N21-N22 | Part5: Ch19, Ch20 | Numerical Integration |
| 13 | N23-N24 | Part5: Ch21 / Part6: Ch22 | Differentiation, IVP |
| 14 | N25-N26 | Part6: Ch22 | IVP (continued) |
| 15 | N27-N28 | Part6: Ch23, Ch24 | Stiff ODEs, BVP |
| 16 | Final | — | — |

<br><a name="repository-structure"></a>
## 🗂 Repository Structure

```plaintext
26Spring_Numerical-Analysis
├── Part1_Modeling-Computers-and-Error-Analysis
│   ├── Chapter01_Mathematical-Modeling
│   │   ├── Concepts_Lab.ko.md
│   │   ├── Concepts_Lab.md
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter02_Programming-and-Software
│   │   ├── Concepts_Lab.ko.md
│   │   ├── Concepts_Lab.md
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter03_Approximations-and-Round-Off-Errors
│   │   ├── Concepts_Lab.ko.md
│   │   ├── Concepts_Lab.md
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   └── Chapter04_Truncation-Errors
│       ├── Concepts_Lecture.ko.md
│       └── Concepts_Lecture.md
├── Part2_Roots-and-Optimization
│   ├── Chapter05_Roots-Bracketing-Methods
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter06_Roots-Open-Methods
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   └── Chapter07_Optimization
│       ├── Concepts_Lecture.ko.md
│       └── Concepts_Lecture.md
├── Part3_Linear-Algebraic-Equations
│   ├── Chapter08_Linear-Algebraic-Equations
│   │   ├── Concepts_Lab.ko.md
│   │   ├── Concepts_Lab.md
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter09_Gauss-Elimination
│   │   ├── Concepts_Lab.ko.md
│   │   ├── Concepts_Lab.md
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter10_LU-Factorization
│   │   ├── Concepts_Lab.ko.md
│   │   ├── Concepts_Lab.md
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter11_Matrix-Inverse-and-Condition
│   │   ├── Concepts_Lab.ko.md
│   │   ├── Concepts_Lab.md
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter12_Eigenvalues-Power-Method
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   └── Chapter13_Eigenvalues-Symmetric-Matrices
│       ├── Concepts_Lecture.ko.md
│       └── Concepts_Lecture.md
├── Part4_Curve-Fitting
│   ├── Chapter14_Linear-Regression
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter15_Nonlinear-Regression
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter16_Fourier-Analysis
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter17_Interpolation
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   └── Chapter18_Splines
│       ├── Concepts_Lecture.ko.md
│       └── Concepts_Lecture.md
├── Part5_Integration-and-Differentiation
│   ├── Chapter19_Numerical-Integration-Formulas
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter20_Numerical-Integration-of-Equations
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   └── Chapter21_Numerical-Differentiation
│       ├── Concepts_Lecture.ko.md
│       └── Concepts_Lecture.md
├── Part6_Ordinary-Differential-Equations
│   ├── Chapter22_Initial-Value-Problems
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter23_Stiff-ODEs
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   └── Chapter24_Boundary-Value-Problems
│       ├── Concepts_Lecture.ko.md
│       └── Concepts_Lecture.md
├── images/
├── .gitignore
├── LICENSE
├── README.ko.md
└── README.md
```

<br><a name="license"></a>
## 🤝 License

This repository is released under the [MIT License](LICENSE).

---
