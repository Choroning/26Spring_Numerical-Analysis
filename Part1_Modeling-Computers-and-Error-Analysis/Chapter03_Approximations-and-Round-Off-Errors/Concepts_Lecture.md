# Chapter 3 Lecture — Approximations and Round-Off Errors

> **Last Updated:** 2026-03-26

---

<br>

## Table of Contents

- [1. Significant Figures](#1-significant-figures)
  - [1.1 Definition](#11-definition)
  - [1.2 Counting Rules](#12-counting-rules)
- [2. Accuracy vs. Precision](#2-accuracy-vs-precision)
  - [2.1 Definitions and Distinction](#21-definitions-and-distinction)
  - [2.2 The Dartboard Analogy](#22-the-dartboard-analogy)
- [3. Error Definitions](#3-error-definitions)
  - [3.1 True Error and True Percent Relative Error](#31-true-error-and-true-percent-relative-error)
  - [3.2 Approximate Percent Relative Error](#32-approximate-percent-relative-error)
  - [3.3 Stopping Criterion](#33-stopping-criterion)
- [4. Scarborough Criterion](#4-scarborough-criterion)
- [5. Round-Off Errors](#5-round-off-errors)
  - [5.1 Source of Round-Off Errors](#51-source-of-round-off-errors)
  - [5.2 Floating-Point Representation](#52-floating-point-representation)
  - [5.3 IEEE 754 Double Precision](#53-ieee-754-double-precision)
  - [5.4 Machine Epsilon](#54-machine-epsilon)
  - [5.5 Truncation Error vs. Round-Off Error](#55-truncation-error-vs-round-off-error)
- [6. Problems with Floating-Point Arithmetic](#6-problems-with-floating-point-arithmetic)
  - [6.1 Subtractive Cancellation](#61-subtractive-cancellation)
  - [6.2 Adding Large and Small Numbers](#62-adding-large-and-small-numbers)
  - [6.3 Smearing](#63-smearing)
- [7. Taylor Series Preview](#7-taylor-series-preview)
- [Summary](#summary)

---

<br>

## 1. Significant Figures

### 1.1 Definition

**Significant figures** (also called significant digits) are the number of **reliable digits** in a number. They indicate the precision with which a quantity is known. In numerical methods, significant figures provide a practical way to express how much confidence we have in a computed or measured result.

For example, the number $0.00456$ has **3 significant figures**: the digits 4, 5, and 6. The leading zeros (0.00...) serve only as placeholders to indicate the magnitude of the number — they do not convey any information about precision.

### 1.2 Counting Rules

The rules for counting significant figures are:

1. **All nonzero digits are significant**: $3.285$ has 4 significant figures
2. **Leading zeros are NOT significant**: $0.00456$ has 3 significant figures (the zeros before 4 are just placeholders)
3. **Trailing zeros after a decimal point ARE significant**: $2.3400$ has 5 significant figures (the trailing zeros indicate measured precision)
4. **Captive zeros (between nonzero digits) ARE significant**: $3.0025$ has 5 significant figures
5. **Trailing zeros in an integer are ambiguous**: $1200$ could have 2, 3, or 4 significant figures. Scientific notation resolves this: $1.2 \times 10^3$ (2 sig. figs.) vs. $1.200 \times 10^3$ (4 sig. figs.)

> **Note:** Significant figures become critically important when we report the results of numerical computations. Reporting more digits than are justified by the accuracy of the method gives a false sense of precision. Conversely, reporting too few digits throws away information that the computation actually provides.

---

<br>

## 2. Accuracy vs. Precision

### 2.1 Definitions and Distinction

These two terms are often confused in everyday language, but they have distinct and important meanings in numerical methods and science:

- **Accuracy**: How close a computed or measured value is to the **true value**. An accurate result is one that is near the correct answer.
- **Precision**: How close **repeated** computations or measurements are to **each other**. A precise result is one that is reproducible and consistent, regardless of whether it is correct.

The key insight is that **you can have precision without accuracy**. If a measurement instrument has a systematic error (bias) — for example, a scale that always reads 0.5 kg too high — then repeated measurements will be very close to each other (high precision) but consistently wrong (low accuracy).

> **[Physics]** In experimental sciences, accuracy relates to systematic errors (calibration), while precision relates to random errors (noise). Both must be considered when evaluating measurements.

### 2.2 The Dartboard Analogy

The relationship between accuracy and precision is best visualized with a dartboard:

| Scenario | Accuracy | Precision | Dart Pattern |
|:---------|:---------|:----------|:-------------|
| Accurate and precise | High | High | Tightly clustered around the bullseye |
| Accurate but not precise | High | Low | Scattered around the bullseye (average is correct) |
| Precise but not accurate | Low | High | Tightly clustered but away from the bullseye (systematic bias) |
| Neither accurate nor precise | Low | Low | Scattered all over the board |

The "precise but not accurate" case is particularly dangerous because the tight clustering gives a false sense of confidence. In numerical methods, this corresponds to a computation that converges consistently to the **wrong answer** — for example, due to a flawed algorithm or an incorrect model.

---

<br>

## 3. Error Definitions

### 3.1 True Error and True Percent Relative Error

When we know the exact (true) value, we can compute the **true error**:

$$E_t = \text{true value} - \text{approximate value}$$

The sign of $E_t$ tells us the direction of the error: positive means the approximation is an underestimate, negative means it is an overestimate. However, the magnitude of $E_t$ alone does not tell us how significant the error is — an error of 1 cm is negligible when measuring a bridge but catastrophic when measuring a microchip.

To normalize the error, we compute the **true percent relative error**:

$$\varepsilon_t = \frac{\text{true value} - \text{approximate value}}{\text{true value}} \times 100\%$$

This expresses the error as a percentage of the true value, making it dimensionless and comparable across different scales. For example:

- Measuring a bridge: true = 10,000 m, approx = 10,001 m → $\varepsilon_t = 0.01\%$ (negligible)
- Measuring a chip: true = 0.001 m, approx = 0.002 m → $\varepsilon_t = 100\%$ (catastrophic)

### 3.2 Approximate Percent Relative Error

In most practical situations, **we do not know the true value** — that is precisely why we are computing a numerical approximation in the first place. In these cases, we use the **approximate percent relative error**, which compares successive iterations of an approximation:

$$\varepsilon_a = \frac{\text{current approximation} - \text{previous approximation}}{\text{current approximation}} \times 100\%$$

The idea is simple: if the answer is not changing much between iterations, then the method has likely converged to a stable result.

> **[Calculus]** The approximate error $\varepsilon_a$ is fundamental to iterative methods. When we can't know the true value, we measure how much our answer is still changing. If it stops changing (small $\varepsilon_a$), we assume convergence.

### 3.3 Stopping Criterion

Iterative methods require a rule for when to stop computing. The standard approach is to pre-specify a **tolerance** $\varepsilon_s$ and iterate until:

$$|\varepsilon_a| < \varepsilon_s$$

This means: keep iterating until the approximate percent relative error drops below the specified tolerance. The choice of $\varepsilon_s$ depends on the required precision of the application — engineering design calculations may tolerate $\varepsilon_s = 0.1\%$, while financial computations may require $\varepsilon_s = 0.0001\%$.

---

<br>

## 4. Scarborough Criterion

The **Scarborough criterion** provides a direct link between the stopping tolerance and the number of correct significant figures in the result. Specifically:

$$\text{If } |\varepsilon_a| < (0.5 \times 10^{2-n})\%, \text{ then the result is correct to at least } n \text{ significant figures.}$$

This gives us a practical way to choose the tolerance $\varepsilon_s$:

| Desired Significant Figures ($n$) | Required $\|\varepsilon_a\|$ |
|:----------------------------------|:----------------------------|
| 1 | $< 5\%$ |
| 2 | $< 0.5\%$ |
| 3 | $< 0.05\%$ |
| 4 | $< 0.005\%$ |
| 5 | $< 0.0005\%$ |
| 6 | $< 0.00005\%$ |
| 7 | $< 0.000005\%$ |

For example, if we need our result to be correct to **3 significant figures**, we set $\varepsilon_s = 0.05\%$ and iterate until $|\varepsilon_a| < 0.05\%$.

> **Note:** The Scarborough criterion provides a sufficient (but not necessary) condition. Meeting the criterion guarantees at least $n$ significant figures, but the result may actually be more accurate than the criterion suggests.

---

<br>

## 5. Round-Off Errors

### 5.1 Source of Round-Off Errors

Round-off errors arise because computers can only store a **finite number of digits**. Since most real numbers have infinitely many digits (e.g., $\pi = 3.14159265...$, $1/3 = 0.33333...$), the computer must truncate or round the representation, introducing a small error every time a number is stored or a calculation is performed.

### 5.2 Floating-Point Representation

Computers represent real numbers in **floating-point** format:

$$\pm m \times b^e$$

where:
- $m$ = **mantissa** (also called significand): contains the significant digits of the number
- $b$ = **base**: the number system used (base 2 for binary computers)
- $e$ = **exponent**: determines the magnitude (position of the decimal/binary point)

This is analogous to scientific notation in base 10: $6.022 \times 10^{23}$ has mantissa $6.022$, base $10$, and exponent $23$.

### 5.3 IEEE 754 Double Precision

The **IEEE 754** standard defines the universal format for floating-point arithmetic in modern computers. The **double precision** (64-bit) format allocates the 64 bits as follows:

| Component | Bits | Purpose |
|:----------|:-----|:--------|
| Sign | 1 bit | $0$ = positive, $1$ = negative |
| Exponent | 11 bits | Range: $-1022$ to $+1023$ (biased by 1023) |
| Mantissa | 52 bits | Fractional part of the significand |

Because binary floating-point numbers are normalized to have a leading 1 (i.e., $1.xxxxx \times 2^e$), the leading 1 is **implicit** and does not need to be stored. This gives an effective 53 bits of precision in the mantissa.

> **[Computer Architecture]** IEEE 754 is the universal standard for floating-point arithmetic in modern CPUs. The 64-bit double precision format stores approximately 15-17 significant decimal digits. Understanding this limitation is essential for all numerical computing.

### 5.4 Machine Epsilon

**Machine epsilon** ($\varepsilon_{\text{mach}}$) is defined as the smallest number such that:

$$1 + \varepsilon_{\text{mach}} \neq 1$$

in floating-point arithmetic. For IEEE 754 double precision:

$$\varepsilon_{\text{mach}} = 2^{-52} \approx 2.22 \times 10^{-16}$$

Machine epsilon represents the **relative precision** of the floating-point system. It tells us the smallest relative perturbation that a floating-point number can "notice." Any number smaller than $\varepsilon_{\text{mach}}$ relative to 1 is effectively zero when added to 1.

This means that double precision provides approximately **15 to 17 significant decimal digits** of precision. Beyond that, the result is unreliable due to round-off.

### 5.5 Truncation Error vs. Round-Off Error

These two types of error have fundamentally different origins:

| Property | Truncation Error | Round-Off Error |
|:---------|:----------------|:----------------|
| **Source** | Approximating a mathematical procedure (e.g., truncating an infinite Taylor series after finitely many terms) | Representing numbers with finite digits in computer memory |
| **Depends on** | Algorithm and step size | Computer hardware and number format |
| **Behavior with more terms/steps** | Generally decreases | May accumulate |
| **Example** | Using $e^x \approx 1 + x + x^2/2$ instead of the full infinite series | Storing $1/3 = 0.333...$ as $0.333333333333333$ |

Understanding both types of error is essential because they often interact: reducing step size to decrease truncation error requires more arithmetic operations, which can increase accumulated round-off error.

---

<br>

## 6. Problems with Floating-Point Arithmetic

### 6.1 Subtractive Cancellation

**Subtractive cancellation** (also called catastrophic cancellation) occurs when two nearly equal numbers are subtracted. The result loses most of its significant digits:

$$0.12345678 - 0.12345671 = 0.00000007$$

The two operands each have 8 significant digits, but the result has **only 1 significant digit**. The 7 leading digits cancel, and all that remains is the least reliable digit of each operand. This is one of the most dangerous sources of numerical error.

Subtractive cancellation is particularly insidious because it can occur in the middle of a longer computation, silently corrupting the result. Any algorithm that involves subtracting quantities that are expected to be close in value should be carefully examined for this issue.

### 6.2 Adding Large and Small Numbers

When a very large number and a very small number are added, the small number may be **completely lost** due to the finite precision of floating-point arithmetic:

$$10^{15} + 1 = 10^{15} \quad \text{(in double precision)}$$

This happens because the number $1$ is below machine epsilon relative to $10^{15}$. When the computer aligns the two numbers to the same exponent for addition, all the digits of the smaller number fall outside the 52-bit mantissa and are discarded.

More precisely, machine epsilon is about $2.22 \times 10^{-16}$, so the smallest number that can be meaningfully added to $10^{15}$ is approximately $10^{15} \times 2.22 \times 10^{-16} \approx 0.222$. The value $1$ is above this threshold and would actually be preserved, but for illustration, $10^{16} + 1 = 10^{16}$ demonstrates the issue clearly — the 1 is lost entirely because it falls below the relative precision.

### 6.3 Smearing

**Smearing** occurs when an alternating series is summed and the intermediate partial sums are much larger in magnitude than the final answer. Each intermediate addition/subtraction introduces a round-off error that is proportional to the magnitude of the intermediate result. Since the intermediate results are much larger than the final answer, these round-off errors can be comparable to or larger than the final result itself.

For example, computing $e^{-20}$ using the Maclaurin series directly:

$$e^{-20} = 1 - 20 + \frac{400}{2!} - \frac{8000}{3!} + \cdots$$

The terms grow very large before eventually decreasing, and the alternating signs mean that massive cancellation occurs. The intermediate sums can reach values on the order of $10^8$, while the final answer is approximately $2.06 \times 10^{-9}$. The round-off errors accumulated during the intermediate sums overwhelm the tiny final result.

A better approach is to compute $e^{20}$ (no smearing because all terms are positive) and then take the reciprocal: $e^{-20} = 1/e^{20}$.

---

<br>

## 7. Taylor Series Preview

The **Taylor series** expansion of a function $f(x)$ about the point $x_i$ is:

$$f(x_{i+1}) = f(x_i) + f'(x_i)h + \frac{f''(x_i)}{2!}h^2 + \frac{f'''(x_i)}{3!}h^3 + \cdots$$

where $h = x_{i+1} - x_i$ is the step size.

This infinite series is exact, but in practice we must truncate it after a finite number of terms, introducing **truncation error**. The Taylor series is the foundation for understanding truncation errors and will be treated in detail in Chapter 4.

Key points to anticipate:

- The **remainder term** (truncation error) depends on the step size $h$ raised to some power: keeping more terms in the series means the error decreases faster as $h$ shrinks
- Taylor series provides the theoretical basis for most numerical methods: finite differences, interpolation, integration, and ODE solvers all derive from Taylor expansions
- The trade-off between truncation error (fewer terms or larger $h$) and round-off error (more computation with smaller $h$) is a central theme in numerical analysis

---

<br>

## Summary

| Topic | Key Point |
|:------|:----------|
| Significant Figures | Number of reliable digits; leading zeros don't count |
| Accuracy vs Precision | Accuracy = closeness to true; Precision = repeatability |
| True Error | $E_t = \text{true} - \text{approx}$ |
| Relative Error | $\varepsilon_t = (E_t / \text{true}) \times 100\%$ |
| Approximate Error | $\varepsilon_a$: used when true value is unknown |
| Scarborough | $|\varepsilon_a| < 0.5 \times 10^{2-n}\%$ → n sig. figs |
| IEEE 754 | 64-bit: 1+11+52 bits, $\varepsilon_{\text{mach}} \approx 2.22 \times 10^{-16}$ |
| Subtractive Cancellation | Nearly equal subtraction → catastrophic digit loss |
| Smearing | Alternating series with large intermediates → accumulated round-off |
| Taylor Series | Foundation for truncation error analysis (Chapter 4) |

---
