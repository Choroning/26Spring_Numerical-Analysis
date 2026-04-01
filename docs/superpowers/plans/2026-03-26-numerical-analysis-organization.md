# Numerical Analysis Repository Organization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Organize 26Spring_Numerical-Analysis with Part/Chapter directory hierarchy and create Concepts_Lecture.md + Concepts_Lab.md files for all covered chapters (weeks 1-4).

**Architecture:** Two-level directory structure: 6 Parts (textbook taxonomy) containing 24 Chapter subdirectories. Each chapter with content gets Concepts_Lecture.md (from lecture PDFs) and Concepts_Lab.md (from .ipynb notebooks) in English.

**Tech Stack:** Markdown, Python code blocks

**Spec:** `docs/superpowers/specs/2026-03-26-numerical-analysis-organization-design.md`

---

## Corrections to Spec

The following corrections apply (lecture PDF content revealed misnamed chapters):

### Part 3 Chapter Titles

| Chapter | Spec (incorrect) | Corrected | Source |
|:---:|:---|:---|:---|
| 08 | Gauss-Elimination | **Linear-Algebraic-Equations** | N05 |
| 09 | LU-Factorization | **Gauss-Elimination** | N06 |
| 10 | Matrix-Inverse-and-Condition | **LU-Factorization** | N07 |
| 11 | Special-Matrices-and-Iterative-Methods | **Matrix-Inverse-and-Condition** | N08 |

> **Note:** The spec's original Ch11 title "Special-Matrices-and-Iterative-Methods" does not correspond to any chapter in the textbook used (Chapra). The professor's Chapter 11 covers Matrix Inverse and Condition Number. The "iterative methods" topic (Gauss-Seidel, Jacobi) is not a standalone chapter in this course.

### Lab Mapping (corrected)

| Notebook | Chapters Covered |
|:---|:---|
| `Lab03-LinearSystem.ipynb` | Ch08 (naive Gauss, setup) + Ch09 (pivoting, Thomas) |
| `Lab04-LU.ipynb` | Ch10 (LU, Cholesky) + Ch11 (inverse, condition) |

### Scope Update

All 4 Part 3 chapters (Ch08-Ch11) now get **both** Concepts_Lecture.md and Concepts_Lab.md.

---

## File Structure

### Directories to Create (all 6 Parts, 24 Chapters)

```
26Spring_Numerical-Analysis/
├── Part1_Modeling-Computers-and-Error-Analysis/
│   ├── Chapter01_Mathematical-Modeling/
│   ├── Chapter02_Programming-and-Software/
│   ├── Chapter03_Approximations-and-Round-Off-Errors/
│   └── Chapter04_Truncation-Errors/
├── Part2_Roots-and-Optimization/
│   ├── Chapter05_Roots-Bracketing-Methods/
│   ├── Chapter06_Roots-Open-Methods/
│   └── Chapter07_Optimization/
├── Part3_Linear-Algebraic-Equations/
│   ├── Chapter08_Linear-Algebraic-Equations/
│   ├── Chapter09_Gauss-Elimination/
│   ├── Chapter10_LU-Factorization/
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
└── images/
```

### Files to Create (content scope)

**Part 1 — 7 files:**
- `Part1/.../Chapter01/Concepts_Lecture.md` (from N01)
- `Part1/.../Chapter01/Concepts_Lab.md` (from Lab01-python-fundamentals)
- `Part1/.../Chapter02/Concepts_Lecture.md` (from N02)
- `Part1/.../Chapter02/Concepts_Lab.md` (from Lab01-python-programming)
- `Part1/.../Chapter03/Concepts_Lecture.md` (from N03)
- `Part1/.../Chapter03/Concepts_Lab.md` (from Lab02-SourceOfErrors)
- `Part1/.../Chapter04/Concepts_Lecture.md` (from N04 + Lab02 Section 4: finite differences code)

**Part 3 — 8 files:**
- `Part3/.../Chapter08/Concepts_Lecture.md` (from N05)
- `Part3/.../Chapter08/Concepts_Lab.md` (from Lab03 — Sections 1-3: setup, naive Gauss, cost)
- `Part3/.../Chapter09/Concepts_Lecture.md` (from N06)
- `Part3/.../Chapter09/Concepts_Lab.md` (from Lab03 — Sections 4-6: pivoting, condition, Thomas)
- `Part3/.../Chapter10/Concepts_Lecture.md` (from N07)
- `Part3/.../Chapter10/Concepts_Lab.md` (from Lab04 — Sections 1-4: LU, solving, thermal, multi-RHS, Cholesky)
- `Part3/.../Chapter11/Concepts_Lecture.md` (from N08)
- `Part3/.../Chapter11/Concepts_Lab.md` (from Lab04 — Sections 6-7: inverse, condition number)

**Root — 3 files:**
- `README.md`
- `README.ko.md`
- `LICENSE`

---

## Markdown Conventions

### Header Template (Concepts_Lecture.md)

```markdown
# Chapter XX Lecture — Topic Name

> **Last Updated:** YYYY-MM-DD

---

<br>

## Table of Contents

- [1. First Section](#1-first-section)
  - [1.1 Subsection](#11-subsection)
- [Summary](#summary)

---

<br>

## 1. First Section

### 1.1 Subsection

Content...

---

<br>

## Summary

| Topic | Key Point |
|:------|:----------|
| ... | ... |

---
```

### Header Template (Concepts_Lab.md)

```markdown
# Chapter XX Lab — Topic Name

> **Last Updated:** YYYY-MM-DD

---

<br>

## Table of Contents

...same structure as Lecture...
```

### Additional Explanations (prerequisite / cross-course references)

```markdown
> **[Linear Algebra]** This is a special case of matrix decomposition...

> **[Calculus]** The Taylor series expansion is derived from...

> **[Python]** NumPy's broadcasting rules allow element-wise operations...

> **Note:** Important clarification about a tricky concept.
```

### Code Blocks

Use language-tagged fenced code blocks:

~~~markdown
```python
import numpy as np
x = np.linalg.solve(A, b)
```
~~~

### Math Formulas

Use inline `$...$` and display `$$...$$` LaTeX notation.

---

## Source Material Paths

```
LECTURE_BASE="/Users/choroning/Library/Mobile Documents/com~apple~CloudDocs/1. 교육 마이데이터 (동기화)/0202401 고려대학교 세종캠퍼스 컴퓨터융합소프트웨어학과/2026학년도 1학기 (3학년 1학기)/전선 DCSS305-00 수치해석(영강)(강신후)/[강의자료]"

LAB_BASE="/Users/choroning/Library/Mobile Documents/com~apple~CloudDocs/1. 교육 마이데이터 (동기화)/0202401 고려대학교 세종캠퍼스 컴퓨터융합소프트웨어학과/2026학년도 1학기 (3학년 1학기)/전선 DCSS305-00 수치해석(영강)(강신후)/[공통]/course305-main"

REPO="/Users/choroning/Desktop/RepoSync/26Spring_Numerical-Analysis"
```

---

## Tasks

### Task 0: Infrastructure — Directories, README, LICENSE

**Files:**
- Create: All 31 directories (6 Parts + 24 Chapters + images/)
- Create: `README.md`, `README.ko.md`
- Verify: `LICENSE` (already exists with MIT, copyright Cheolwon Park)
- Modify: `docs/superpowers/specs/2026-03-26-numerical-analysis-organization-design.md` (apply corrections)

- [ ] **Step 1: Create all Part and Chapter directories**

```bash
BASE="/Users/choroning/Desktop/RepoSync/26Spring_Numerical-Analysis"

# Part 1
mkdir -p "$BASE/Part1_Modeling-Computers-and-Error-Analysis/Chapter01_Mathematical-Modeling"
mkdir -p "$BASE/Part1_Modeling-Computers-and-Error-Analysis/Chapter02_Programming-and-Software"
mkdir -p "$BASE/Part1_Modeling-Computers-and-Error-Analysis/Chapter03_Approximations-and-Round-Off-Errors"
mkdir -p "$BASE/Part1_Modeling-Computers-and-Error-Analysis/Chapter04_Truncation-Errors"

# Part 2
mkdir -p "$BASE/Part2_Roots-and-Optimization/Chapter05_Roots-Bracketing-Methods"
mkdir -p "$BASE/Part2_Roots-and-Optimization/Chapter06_Roots-Open-Methods"
mkdir -p "$BASE/Part2_Roots-and-Optimization/Chapter07_Optimization"

# Part 3
mkdir -p "$BASE/Part3_Linear-Algebraic-Equations/Chapter08_Linear-Algebraic-Equations"
mkdir -p "$BASE/Part3_Linear-Algebraic-Equations/Chapter09_Gauss-Elimination"
mkdir -p "$BASE/Part3_Linear-Algebraic-Equations/Chapter10_LU-Factorization"
mkdir -p "$BASE/Part3_Linear-Algebraic-Equations/Chapter11_Matrix-Inverse-and-Condition"
mkdir -p "$BASE/Part3_Linear-Algebraic-Equations/Chapter12_Eigenvalues-Power-Method"
mkdir -p "$BASE/Part3_Linear-Algebraic-Equations/Chapter13_Eigenvalues-Symmetric-Matrices"

# Part 4
mkdir -p "$BASE/Part4_Curve-Fitting/Chapter14_Linear-Regression"
mkdir -p "$BASE/Part4_Curve-Fitting/Chapter15_Nonlinear-Regression"
mkdir -p "$BASE/Part4_Curve-Fitting/Chapter16_Fourier-Analysis"
mkdir -p "$BASE/Part4_Curve-Fitting/Chapter17_Interpolation"
mkdir -p "$BASE/Part4_Curve-Fitting/Chapter18_Splines"

# Part 5
mkdir -p "$BASE/Part5_Integration-and-Differentiation/Chapter19_Numerical-Integration-Formulas"
mkdir -p "$BASE/Part5_Integration-and-Differentiation/Chapter20_Numerical-Integration-of-Equations"
mkdir -p "$BASE/Part5_Integration-and-Differentiation/Chapter21_Numerical-Differentiation"

# Part 6
mkdir -p "$BASE/Part6_Ordinary-Differential-Equations/Chapter22_Initial-Value-Problems"
mkdir -p "$BASE/Part6_Ordinary-Differential-Equations/Chapter23_Stiff-ODEs"
mkdir -p "$BASE/Part6_Ordinary-Differential-Equations/Chapter24_Boundary-Value-Problems"

# images
mkdir -p "$BASE/images"
```

- [ ] **Step 2: Verify LICENSE**

LICENSE already exists (MIT, copyright Cheolwon Park). Verify it is present and correct; no changes needed.

- [ ] **Step 3: Create README.md**

Structure:
- Repository title and badges
- Course info table (DCSS305-00, Numerical Analysis, Prof. Shinhoo Kang, etc.)
- Textbook reference
- Tools/environment table
- **Week-to-Chapter mapping table** (16 weeks including Midterm/Final)
- Repository structure overview (Part/Chapter tree)

Week-to-Chapter table:

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

- [ ] **Step 4: Create README.ko.md**

Korean translation of README.md. Identical content, translated to Korean.

- [ ] **Step 5: Update spec document**

Apply the corrections listed in the "Corrections to Spec" section above to the spec file.

---

### Task 1: Chapter 01 — Mathematical Modeling (Lecture + Lab)

**Files:**
- Create: `Part1_Modeling-Computers-and-Error-Analysis/Chapter01_Mathematical-Modeling/Concepts_Lecture.md`
- Create: `Part1_Modeling-Computers-and-Error-Analysis/Chapter01_Mathematical-Modeling/Concepts_Lab.md`

**Sources:**
- Lecture: N01 PDF (Chapter 1 — Mathematical Modeling and Engineering Problem Solving)
- Lab: `01-PythonBasic/Lab01-python-fundamentals.ipynb`

- [ ] **Step 1: Write Concepts_Lecture.md for Chapter 01**

Sections to cover (from N01):
1. What is a Mathematical Model?
   - Definition, general form: dependent = f(independent, parameters, forcing)
2. Newton's 2nd Law — Bungee Jumper Example
   - Force balance: $F_{net} = mg - cv$
   - ODE: $dv/dt = g - (c/m)v$
   - Analytical solution: $v(t) = \frac{gm}{c}(1 - e^{-(c/m)t})$
   - Terminal velocity derivation
   > **[Physics]** Newton's 2nd law and force balance
   > **[Calculus]** First-order ODE and its analytical solution via separation of variables
3. Euler's Method — Numerical Solution
   - Formula: $v(t_{i+1}) = v(t_i) + \frac{dv}{dt}\bigg|_{t_i} \cdot h$
   - Step-by-step example: m=68.1kg, c=12.5 kg/s
   - Convergence as step size decreases (h=2, h=1, h=0.5)
   > **[Calculus]** Euler's method approximates the derivative as a finite difference
4. Types of Engineering Problems Solved Numerically
   - Roots, linear systems, curve fitting, integration, ODEs, PDEs, optimization
5. Summary table

- [ ] **Step 2: Write Concepts_Lab.md for Chapter 01**

Content from Lab01-python-fundamentals.ipynb:
1. Variables and Assignment
   - Dynamic typing, `type()`, complex numbers, `math.pi`, `math.e`
   - Format strings, Planck's constant from `scipy.constants`
2. NumPy Arrays, Vectors, and Matrices
   - Array creation: `np.array`, `np.zeros`, `np.ones`, `np.eye`, `np.arange`, `np.linspace`, `np.logspace`
   - Indexing and slicing: `A[0,:]`, `A[:,1]`, negative indices
   - `np.matrix` (deprecated note)
3. Mathematical Operations
   - Operators table: +, -, *, /, //, %, **, @
   - Precedence (PEMDAS), element-wise vs matrix multiply
   - Dot product, cross product, transpose
4. Built-in Functions
   - `abs`, `round`, `min`, `max`, `len`, `type`
   - **Bungee Jumper Velocity** application with vectorized computation
5. Visualization with Matplotlib
   - `plt.plot`, `plt.scatter`, labels, legend, grid, `savefig`
   - Multi-curve example with different drag coefficients
6. Views vs. Copies in NumPy
   - Slicing creates views, fancy indexing creates copies, `y.base`

All code cells included as Python code blocks. Concepts explained step-by-step, more thoroughly than the notebook.

---

### Task 2: Chapter 02 — Programming and Software (Lecture + Lab)

**Files:**
- Create: `Part1_.../Chapter02_Programming-and-Software/Concepts_Lecture.md`
- Create: `Part1_.../Chapter02_Programming-and-Software/Concepts_Lab.md`

**Sources:**
- Lecture: N02 PDF (Chapter 2 — MATLAB/Python Fundamentals)
- Lab: `01-PythonBasic/Lab01-python-programming.ipynb`

- [ ] **Step 1: Write Concepts_Lecture.md for Chapter 02**

Sections (from N02):
1. Programming Environment — Python with NumPy/SciPy/Matplotlib
2. Python Basics — Variables, data types, NumPy arrays
3. Matrix Operations — `np.dot`, `@`, `.T`, `np.linalg.solve`, `np.linalg.inv`, `np.linalg.det`
4. Plotting with Matplotlib — `plt.plot`, labels, titles
5. Control Structures — for, while, if/elif/else
6. Functions — `def`, return values
7. Example: Bungee Jumper in Python — Euler's method implementation
   > **[Python]** Using scipy.integrate.solve_ivp for ODE solving
8. Summary table

- [ ] **Step 2: Write Concepts_Lab.md for Chapter 02**

Content from Lab01-python-programming.ipynb:
1. Script Files and Functions
   - Bungee jumper inline script → function refactoring
   - Docstrings, multiple return values (tuples)
2. Variable Scope (LEGB Rule)
   - Local, Enclosing, Global, Built-in
   - `global` keyword (with caution), nested functions
   > **[Programming Languages]** Scope rules are universal across languages; Python's LEGB is analogous to other lexical scoping rules
3. Arguments
   - Positional vs keyword, `*args`, `**kwargs`
4. File I/O
   - `np.loadtxt`, `np.save`/`np.load`, `np.savez`
   - KMA weather data example
5. Structured Programming
   - if/elif/else, for loops, while loops, break/continue
   - **Vectorization**: for-loop vs NumPy vectorized `np.cos(10*t)` with timing comparison
6. Functions as First-Class Objects
   - Lambda, closures, `functools.partial`, higher-order functions
   - `zip`, `filter`, comprehensions
7. Exception Handling — `try/except`
8. Decorators — function and method decorators
9. OOP Basics — `__init__`, `__repr__`, `@dataclass`, `@property`
10. Type Hints — `from typing import List, Optional, Tuple`
11. Generators — `yield`, lazy evaluation, memory efficiency

---

### Task 3: Chapter 03 — Approximations and Round-Off Errors (Lecture + Lab)

**Files:**
- Create: `Part1_.../Chapter03_Approximations-and-Round-Off-Errors/Concepts_Lecture.md`
- Create: `Part1_.../Chapter03_Approximations-and-Round-Off-Errors/Concepts_Lab.md`

**Sources:**
- Lecture: N03 PDF (Chapter 3)
- Lab: `02-SourceOfErrors/Lab02-SourceOfErrors.ipynb`

- [ ] **Step 1: Write Concepts_Lecture.md for Chapter 03**

Sections (from N03):
1. Significant Figures — definition, examples
2. Accuracy vs. Precision — bias, systematic error
3. Error Definitions
   - True error: $E_t = \text{true} - \text{approx}$
   - True percent relative error: $\varepsilon_t = \frac{E_t}{\text{true}} \times 100\%$
   - Approximate percent relative error: $\varepsilon_a = \frac{\text{current} - \text{previous}}{\text{current}} \times 100\%$
   > **[Calculus]** The iterative error concept connects to convergence of sequences
4. Scarborough Criterion
   - $|\varepsilon_a| < 0.5 \times 10^{2-n}\%$ guarantees n significant figures
5. Round-Off Errors
   - Floating-point representation: $\pm m \times b^e$
   - IEEE 754 double precision: 1 + 11 + 52 bits
   - Machine epsilon: $\varepsilon_{mach} = 2^{-52} \approx 2.22 \times 10^{-16}$
   > **[Computer Architecture]** IEEE 754 floating-point standard and how CPUs store real numbers
6. Subtractive Cancellation — when subtracting nearly equal numbers
7. Taylor Series Preview — connects to Chapter 4
8. Summary table

- [ ] **Step 2: Write Concepts_Lab.md for Chapter 03**

Content from Lab02-SourceOfErrors.ipynb (Sections 1-3):
1. Error Definitions and Metrics
   - Scarborough criterion lambda function
   - **sqrt(2) Approximation** (Babylonian/Heron's method) showing quadratic convergence
   > **[Linear Algebra]** The Babylonian method is a special case of Newton's method for root finding
2. Truncation Error and Taylor Series
   - Maclaurin series for $e^x$: term-by-term approximation
   - Iterative `exp_series()` function with stopping criterion
   - Visualization: `plt.semilogy` of error vs terms
3. Round-Off Error and Floating-Point
   - Machine epsilon demo: `np.finfo(float).eps`
   - Binary representation function
   - Accumulation of round-off: summing 0.0001 ten thousand times
   - **Kahan Summation Algorithm** — compensation technique
   - **Catastrophic Cancellation**: quadratic formula example with stable alternative

---

### Task 4: Chapter 04 — Truncation Errors (Lecture only)

**Files:**
- Create: `Part1_.../Chapter04_Truncation-Errors/Concepts_Lecture.md`

**Source:** N04 PDF (55 pages — the largest lecture, requires thorough coverage)

- [ ] **Step 1: Write Concepts_Lecture.md for Chapter 04**

Sections (from N04):
1. Taylor Series Expansion — Complete Treatment
   - Definition, step size form, remainder term
   > **[Calculus]** Taylor's theorem with Lagrange remainder
2. Order of Approximation
   - Zero through nth-order, O(h^{n+1}) error
   - Example: $e^x$ at $x=1$ about $x=0$ (0th through 6th order with error percentages)
3. Numerical Differentiation Using Taylor Series
   - **Forward Difference**: $f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h} + O(h)$
   - **Backward Difference**: $f'(x_i) = \frac{f(x_i) - f(x_{i-1})}{h} + O(h)$
   - **Central Difference**: $f'(x_i) = \frac{f(x_{i+1}) - f(x_{i-1})}{2h} + O(h^2)$
   - **Second Derivative (Central)**: $f''(x_i) = \frac{f(x_{i+1}) - 2f(x_i) + f(x_{i-1})}{h^2}$
   - Higher-order forward difference formula (3-point O(h^2))
   > **[Calculus]** Connection between finite differences and derivative definitions via limits
4. Detailed Example
   - $f(x) = -0.1x^4 - 0.15x^3 - 0.5x^2 - 0.25x + 1.2$
   - Forward, backward, central differences at x=0.5 with h=0.5 and h=0.25
   - Error comparison showing central difference superiority
5. Total Numerical Error
   - Truncation error + round-off error
   - Optimal step size $h^*$ minimizing total error
   - V-shaped error curve on log-log plot
6. Stability and Condition Number (of a function)
   - $\text{Cond} = |x \cdot f'(x) / f(x)|$
   - Ill-conditioned vs well-conditioned problems
7. Summary table

Also include content from Lab02 Section 4 (Finite Differences):
- Step size analysis code showing optimal h
- Forward vs central difference comparison exercise

---

### Task 5: Chapter 08 — Linear Algebraic Equations (Lecture + Lab)

**Files:**
- Create: `Part3_Linear-Algebraic-Equations/Chapter08_Linear-Algebraic-Equations/Concepts_Lecture.md`
- Create: `Part3_Linear-Algebraic-Equations/Chapter08_Linear-Algebraic-Equations/Concepts_Lab.md`

**Sources:**
- Lecture: N05 PDF (Chapter 8)
- Lab: Lab03-LinearSystem.ipynb (Sections 1-3)

- [ ] **Step 1: Write Concepts_Lecture.md for Chapter 08**

Sections (from N05):
1. Motivation — Mass-spring, truss, circuit examples
   > **[Physics]** Hooke's law and force equilibrium in structural systems
2. Small Systems (2x2) — Graphical interpretation
   - Three cases: unique, no solution, infinite solutions
3. Determinant and Solvability
   - 2x2 determinant, singular/ill-conditioned systems
   > **[Linear Algebra]** Determinant properties and their geometric interpretation
4. Cramer's Rule
   - Derivation by elimination, determinant form
   - Limitation: O(n!) computational cost
5. Naive Gauss Elimination
   - **Forward Elimination**: pivot element, factor computation, row operations
   - Pseudocode with nested loops
   - **Back Substitution**: formula $x_i = \frac{b_i^{(i-1)} - \sum_{j=i+1}^{n} a_{ij}^{(i-1)} x_j}{a_{ii}^{(i-1)}}$
6. Computational Cost
   - Forward elimination: $(2/3)n^3 + O(n^2)$
   - Back substitution: $(1/2)n^2 + O(n)$
   - Total dominated by $n^3$
7. Example (3x3 system)
   - Complete step-by-step: 3x₁ - 0.1x₂ - 0.2x₃ = 7.85, etc.
8. Pitfalls of Naive Gauss Elimination
   - Division by zero, near-zero pivots, large factor values
9. Summary table

- [ ] **Step 2: Write Concepts_Lab.md for Chapter 08**

Content from Lab03 (Sections 1-3):
1. Setting Up Ax = b
   - Matrix form, coefficient matrix, unknown vector, RHS
   - Geometric interpretation code (3-subplot figure)
2. Spring-Mass System Application
   - Engineering example with 3 masses, 4 springs
   - `np.linalg.solve(A, b)` usage, residual verification
3. Naive Gaussian Elimination Implementation
   - Augmented matrix, forward elimination, back substitution
   - Complete `gauss_naive()` implementation with code
4. Computational Cost Analysis
   - Timing benchmark code: `np.linalg.solve` for n=50 to 800
   - Log-log plot confirming O(n^3)

---

### Task 6: Chapter 09 — Gauss Elimination (Lecture + Lab)

**Files:**
- Create: `Part3_.../Chapter09_Gauss-Elimination/Concepts_Lecture.md`
- Create: `Part3_.../Chapter09_Gauss-Elimination/Concepts_Lab.md`

**Sources:**
- Lecture: N06 PDF (Chapter 9)
- Lab: Lab03-LinearSystem.ipynb (Sections 4-6)

- [ ] **Step 1: Write Concepts_Lecture.md for Chapter 09**

Sections (from N06):
1. Partial Pivoting
   - Problem: zero/near-zero pivots
   - Solution: swap rows to use largest absolute value as pivot
   - Algorithm modification with pseudocode
2. Scaled Partial Pivoting
   - Normalize rows before comparing pivot candidates
   - Scale factor: $s_i = \max_j |a_{ij}|$
3. Example Demonstrating Pivoting Necessity
   - 2x2 system with 0.0003 pivot: wrong without pivoting, correct with pivoting
   - Shows how limited precision causes errors
4. Gauss-Jordan Method
   - Eliminate above AND below diagonal → diagonal matrix
   - No back substitution needed, useful for matrix inverse
   - More expensive than Gauss elimination
5. Tridiagonal Systems — Thomas Algorithm
   - Structure: only main diagonal, one above, one below
   - Forward sweep + back substitution
   - **O(n) cost** — extremely efficient!
   > **[Data Structures]** The Thomas algorithm exploits the banded structure, similar to how sparse data structures save memory
6. Python Implementation
   - `np.linalg.solve()` for general systems
7. Summary table

- [ ] **Step 2: Write Concepts_Lab.md for Chapter 09**

Content from Lab03 (Sections 4-6):
1. Partial Pivoting Implementation
   - Zero-pivot problem demo
   - Complete `gauss_pivot()` implementation
   - `np.argmax` for finding max element, row swap
2. Condition Number
   - $\kappa(A) = \|A\|\|A^{-1}\|$ interpretation
   - Well-conditioned vs Hilbert matrix examples
   - Perturbation sensitivity demo
3. Tridiagonal Systems and Thomas Algorithm
   - Complete `thomas()` implementation
   - **1D Heat Conduction Application**: `heat_conduction_1d()` with visualization
   - Performance comparison: Thomas O(n) vs numpy O(n^3) with speedup table

---

### Task 7: Chapter 10 — LU Factorization (Lecture + Lab)

**Files:**
- Create: `Part3_.../Chapter10_LU-Factorization/Concepts_Lecture.md`
- Create: `Part3_.../Chapter10_LU-Factorization/Concepts_Lab.md`

**Sources:**
- Lecture: N07 PDF (Chapter 10)
- Lab: Lab04-LU.ipynb (Sections 1-5)

- [ ] **Step 1: Write Concepts_Lecture.md for Chapter 10**

Sections (from N07):
1. Overview of LU Factorization
   - Goal: $A = LU$, then solve $Ld = b$ and $Ux = d$
   - Key advantage: decompose once, solve many RHS
   - Cost: LU O(n^3), each solve O(n^2)
   > **[Linear Algebra]** LU is related to Gaussian elimination; the factors used during elimination form L
2. Gauss Elimination as LU Factorization
   - Factors $f_{ij}$ become L entries
   - U is the upper triangular result
   - Verification: LU = A
3. Storage Efficiency — Store L and U in single matrix
4. Forward and Back Substitution Formulas
5. Example 10.1 — 3x3 LU Decomposition step-by-step
6. Example 10.2 — Solving with LU (forward sub → back sub)
7. LU Factorization with Pivoting
   - $PA = LU$ where P is permutation matrix
   - Procedure: apply P, decompose, solve in two steps
   - Example 10.3 with 2x2 system
8. Cholesky Factorization
   - For symmetric positive definite matrices: $A = U^T U$
   - Saves memory and computation
   - General formulas for diagonal and off-diagonal elements
   - Example 10.5 (3x3 Cholesky)
   > **[Linear Algebra]** SPD matrices arise naturally in physics (stiffness matrices, covariance matrices)
9. Python Code — `scipy.linalg.lu`, `scipy.linalg.cholesky`
10. Summary table

- [ ] **Step 2: Write Concepts_Lab.md for Chapter 10**

Content from Lab04 (Sections 1-5):
1. LU Factorization — Doolittle's Algorithm
   - Formulas and step-by-step 3x3 example
   - Complete `lu_doolittle()` implementation
   - Verification against `sc.linalg.lu`
   - `scipy.linalg.lu` with permutation: $A = PLU$
2. Solving via L and U
   - Complete `forward_sub()` and `backward_sub()` implementations
   - Complete `lu_solve()` combining P, L, U
   - Test with 3x3 system, residual verification
3. Thermal Resistor Network Application
   - 4-node network problem setup and solution
4. Multiple Right-Hand Sides
   - Why LU shines: factor once + m solves
   - Timing comparison code
   - Spring-mass system with 3 load cases
   - **Superposition verification**: ||x_comb - (x_grav + x_wind)|| ≈ 0
5. Cholesky Factorization
   - SPD definition, `is_spd()` function
   - Complete `cholesky_solve()` implementation
   - FEM spring network application

---

### Task 8: Chapter 11 — Matrix Inverse and Condition (Lecture + Lab)

**Files:**
- Create: `Part3_.../Chapter11_Matrix-Inverse-and-Condition/Concepts_Lecture.md`
- Create: `Part3_.../Chapter11_Matrix-Inverse-and-Condition/Concepts_Lab.md`

**Sources:**
- Lecture: N08 PDF (Chapter 11)
- Lab: Lab04-LU.ipynb (Sections 6-7)

- [ ] **Step 1: Write Concepts_Lecture.md for Chapter 11**

Sections (from N08):
1. Matrix Inverse Definition — $AA^{-1} = A^{-1}A = I$
2. Computing the Inverse via LU Factorization
   - Solve $Ax_k = e_k$ for each column of identity
   - Example 11.1: complete 3x3 inverse calculation
3. Physical Interpretation — Mass-spring system
   - $a_{ij}^{-1}$ = contribution of force $b_j$ to displacement $x_i$
4. Error Analysis Methods
   - Three checks for ill-conditioning
5. Vector Norms
   - General Lp norm, L1 (Manhattan), L2 (Euclidean), L∞ (max)
   > **[Linear Algebra]** Norms generalize the concept of "length" to arbitrary vector spaces
6. Matrix Norms
   - Frobenius, column-sum (L1), row-sum (L∞), spectral (L2)
   - Formal norm properties (non-negativity, definiteness, homogeneity, triangle inequality)
   - Induced matrix norm definition
   - Submultiplicativity: $\|Ax\| \leq \|A\| \cdot \|x\|$
7. Matrix Condition Number
   - $\text{Cond}[A] = \|A\| \cdot \|A^{-1}\|$
   - Well-conditioned (≈1) vs ill-conditioned (>>1)
8. Error Bound Theorem
   - $\frac{\|\delta x\|}{\|x\|} \leq \text{Cond}[A] \cdot \frac{\|\delta A\|}{\|A\|}$
9. Precision Rule
   - $t$-digit coefficients → $(t - \log_{10}(\text{Cond}[A]))$-digit solution
10. Example 11.3 — Hilbert Matrix
    - 3x3 Hilbert matrix: Cond = 451.2
    > **[Linear Algebra]** Hilbert matrices are a classic example of ill-conditioning in numerical linear algebra
11. Python — `np.linalg.inv`, `np.linalg.cond`, `np.linalg.norm`
12. Summary table

- [ ] **Step 2: Write Concepts_Lab.md for Chapter 11**

Content from Lab04 (Sections 6-7):
1. Matrix Inverse
   - Definition, useful identities (inverse of product, transpose, 2x2 formula)
   - Complete `lu_inverse()` implementation
   - **Golden Rule: Solve, Don't Invert** — timing demo showing solve is ~3x faster
2. Condition Number and Ill-Conditioning
   - $\kappa(A) = \|A\|\|A^{-1}\|$, error bound formula
   - Computing: `np.linalg.cond(A)` with different norms
   - 2-norm = ratio of largest/smallest singular value
   - Geometric interpretation (well-conditioned vs ill-conditioned plots)
   - **Hilbert Matrix**: condition number growth table (n=2..15)
   - Error amplification demo with perturbation
   - Diagnostic exercise: analyze suspicious stiffness matrix

---

### Task 9: Final — Verify and Commit

- [ ] **Step 1: Verify all files exist and are non-empty**

```bash
find "$BASE" -name "*.md" -not -empty | sort
```

- [ ] **Step 2: Verify directory structure matches plan**

```bash
find "$BASE" -type d | sort
```

- [ ] **Step 3: Commit**

```bash
git add 26Spring_Numerical-Analysis/
git commit -m "$(cat <<'EOF'
Organize Numerical Analysis repo: Part/Chapter structure with concept notes (weeks 1-4)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Execution Notes

### Parallelization

Tasks 1-4 (Part 1, Ch01-Ch04) are independent and can run in parallel.
Tasks 5-8 (Part 3, Ch08-Ch11) are independent and can run in parallel.
Task 0 must run first. Task 9 must run last.

### Content Guidelines

1. **No content from lectures should be omitted.** Every formula, example, derivation, and remark must appear in the corresponding Concepts_Lecture.md.
2. **Lab explanations should be more detailed than the original notebooks.** Add step-by-step explanations, clarify why each step works, explain the intuition.
3. **Additional explanations** use blockquote format: `> **[Subject]** explanation`
   - Relevant subjects: Linear Algebra, Calculus, Physics, Python, Programming Languages, Computer Architecture, Data Structures
4. **Code blocks** include all code from notebooks with language tags.
5. **LaTeX math** for all formulas: inline `$...$`, display `$$...$$`.
