# Chapter 3 — Roundoff Errors

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-10-27
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 4

> **Learning Objectives:**
> 1. Count significant figures and distinguish accuracy from precision
> 2. Analyze roundoff errors in floating-point arithmetic
> 3. Apply the IEEE 754 floating-point representation

---

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
- [5. Roundoff Errors](#5-roundoff-errors)
  - [5.1 Source of Roundoff Errors](#51-source-of-roundoff-errors)
  - [5.2 Floating-Point Representation](#52-floating-point-representation)
  - [5.3 IEEE 754 Double Precision](#53-ieee-754-double-precision)
  - [5.4 Machine Epsilon](#54-machine-epsilon)
  - [5.5 Truncation Error vs. Roundoff Error](#55-truncation-error-vs-roundoff-error)
- [6. Problems with Floating-Point Arithmetic](#6-problems-with-floating-point-arithmetic)
  - [6.1 Subtractive Cancellation](#61-subtractive-cancellation)
  - [6.2 Adding Large and Small Numbers](#62-adding-large-and-small-numbers)
  - [6.3 Smearing](#63-smearing)
- [7. Taylor Series Preview](#7-taylor-series-preview)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Significant Figures

### 1.1 Definition

**Significant figures** (also called significant digits) are the number of **reliable digits** in a number. They indicate the precision with which a quantity is known. In numerical methods, significant figures provide a practical way to express how much confidence we have in a computed or measured result.

For example, the number $`0.00456`$ has **3 significant figures**: the digits 4, 5, and 6. The leading zeros (0.00...) serve only as placeholders to indicate the magnitude of the number and convey no information about precision.

### 1.2 Counting Rules

The rules for counting significant figures are as follows.

1. **All nonzero digits are significant.** $`3.285`$ has 4 significant figures.
2. **Leading zeros are not significant.** $`0.00456`$ has 3 significant figures; the zeros before 4 are only placeholders.
3. **Trailing zeros after a decimal point are significant.** $`2.3400`$ has 5 significant figures; the trailing zeros indicate measured precision.
4. **Captive zeros (between nonzero digits) are significant.** $`3.0025`$ has 5 significant figures.
5. **Trailing zeros in an integer are ambiguous.** $`1200`$ could have 2, 3, or 4 significant figures. Scientific notation resolves this: $`1.2 \times 10^3`$ (2 significant figures) versus $`1.200 \times 10^3`$ (4 significant figures).

> **Note:** Significant figures become critically important when we report the results of numerical computations. Reporting more digits than are justified by the accuracy of the method gives a false sense of precision. Conversely, reporting too few digits throws away information that the computation actually provides.

---

<br>

## 2. Accuracy vs. Precision

### 2.1 Definitions and Distinction

These two terms are often confused in everyday language, but they have distinct and important meanings in numerical methods and science.

- **Accuracy** is how close a computed or measured value is to the **true value**. An accurate result is one that is near the correct answer.
- **Precision** is how close **repeated** computations or measurements are to **each other**. A precise result is one that is reproducible and consistent, regardless of whether it is correct.

The key insight is that **precision is possible without accuracy**. If a measurement instrument has a systematic error (bias), for example a scale that always reads 0.5 kg too high, then repeated measurements will be very close to each other (high precision) but consistently wrong (low accuracy).

> **Note:** **The experimental view.** In experimental sciences, accuracy relates to systematic errors (calibration), while precision relates to random errors (noise). Both must be considered when evaluating measurements.

### 2.2 The Dartboard Analogy

The relationship between accuracy and precision is easiest to visualize with a dartboard.

| Scenario | Accuracy | Precision | Dart Pattern |
|:---------|:---------|:----------|:-------------|
| Accurate and precise | High | High | Tightly clustered around the bullseye |
| Accurate but not precise | High | Low | Scattered around the bullseye (average is correct) |
| Precise but not accurate | Low | High | Tightly clustered but away from the bullseye (systematic bias) |
| Neither accurate nor precise | Low | Low | Scattered all over the board |

*Table 1. Four combinations of accuracy and precision in the dartboard analogy*

The precise but not accurate case is particularly dangerous because the tight clustering gives a false sense of confidence. In numerical methods, this corresponds to a computation that converges consistently to the **wrong answer**, for example because of a flawed algorithm or an incorrect model.

---

<br>

## 3. Error Definitions

### 3.1 True Error and True Percent Relative Error

When we know the exact (true) value, we can compute the **true error**.

```math
E_t = \text{true value} - \text{approximate value}
```

The sign of $`E_t`$ tells us the direction of the error: positive means the approximation is an underestimate, negative means it is an overestimate. However, the magnitude of $`E_t`$ alone does not tell us how significant the error is. An error of 1 cm is negligible when measuring a bridge but catastrophic when measuring a microchip.

To normalize the error, we compute the **true percent relative error**.

```math
\varepsilon_t = \frac{\text{true value} - \text{approximate value}}{\text{true value}} \times 100\%
```

This expresses the error as a percentage of the true value, making it dimensionless and comparable across different scales. For example:

- Measuring a bridge with true value 10,000 m and approximation 10,001 m gives $`|\varepsilon_t| = 0.01\%`$, which is negligible.
- Measuring a chip with true value 0.001 m and approximation 0.002 m gives $`|\varepsilon_t| = 100\%`$, which is catastrophic.

### 3.2 Approximate Percent Relative Error

In most practical situations, **we do not know the true value**; that is precisely why we are computing a numerical approximation in the first place. In these cases, we use the **approximate percent relative error**, which compares successive iterations of an approximation.

```math
\varepsilon_a = \frac{\text{current approximation} - \text{previous approximation}}{\text{current approximation}} \times 100\%
```

The idea is simple: if the answer is not changing much between iterations, then the method has likely converged to a stable result.

> **Note:** **Judging convergence of iterative methods.** The approximate error $`\varepsilon_a`$ is fundamental to iterative methods. When we cannot know the true value, we measure how much our answer is still changing. If it stops changing (small $`\varepsilon_a`$), we assume convergence.

### 3.3 Stopping Criterion

Iterative methods require a rule for when to stop computing. The standard approach is to pre-specify a **tolerance** $`\varepsilon_s`$ and iterate until the following holds.

```math
|\varepsilon_a| < \varepsilon_s
```

This means that we keep iterating until the approximate percent relative error drops below the specified tolerance. The choice of $`\varepsilon_s`$ depends on the required precision of the application. Engineering design calculations may tolerate $`\varepsilon_s = 0.1\%`$, while financial computations may require $`\varepsilon_s = 0.0001\%`$.

---

<br>

## 4. Scarborough Criterion

The **Scarborough criterion** provides a direct link between the stopping tolerance and the number of correct significant figures in the result. Specifically:

```math
\text{If } |\varepsilon_a| < (0.5 \times 10^{2-n})\%, \text{ then the result is correct to at least } n \text{ significant figures.}
```

This gives us a practical way to choose the tolerance $`\varepsilon_s`$.

| Desired Significant Figures ($`n`$) | Required $`\vert\varepsilon_a\vert`$ |
|:----------------------------------|:----------------------------|
| 1 | $`< 5\%`$ |
| 2 | $`< 0.5\%`$ |
| 3 | $`< 0.05\%`$ |
| 4 | $`< 0.005\%`$ |
| 5 | $`< 0.0005\%`$ |
| 6 | $`< 0.00005\%`$ |
| 7 | $`< 0.000005\%`$ |

*Table 2. Scarborough tolerance for each desired number of significant figures*

For example, if we need our result to be correct to **3 significant figures**, we set $`\varepsilon_s = 0.05\%`$ and iterate until $`|\varepsilon_a| < 0.05\%`$.

> **Note:** The Scarborough criterion provides a sufficient (but not necessary) condition. Meeting the criterion guarantees at least $`n`$ significant figures, but the result may actually be more accurate than the criterion suggests.

---

<br>

## 5. Roundoff Errors

### 5.1 Source of Roundoff Errors

Roundoff errors arise because computers can only store a **finite number of digits**. Since most real numbers have infinitely many digits (e.g., $`\pi = 3.14159265...`$, $`1/3 = 0.33333...`$), the computer must truncate or round the representation, introducing a small error every time a number is stored or a calculation is performed.

### 5.2 Floating-Point Representation

Computers represent real numbers in **floating-point** format.

```math
\pm m \times b^e
```

where the symbols mean the following:
- $`m`$ is the **mantissa**, which contains the significant digits of the number. It is also called the significand.
- $`b`$ is the **base**, the number system used. It is 2 for binary computers.
- $`e`$ is the **exponent**, which determines the magnitude, that is, the position of the decimal or binary point.

This is analogous to scientific notation in base 10: $`6.022 \times 10^{23}`$ has mantissa $`6.022`$, base $`10`$, and exponent $`23`$.

### 5.3 IEEE 754 Double Precision

The **IEEE 754** standard defines the universal format for floating-point arithmetic in modern computers. The **double precision** (64-bit) format allocates the 64 bits as follows.

| Component | Bits | Purpose |
|:----------|:-----|:--------|
| Sign | 1 bit | $`0`$ = positive, $`1`$ = negative |
| Exponent | 11 bits | Range for normalized numbers: $`-1022`$ to $`+1023`$ (stored with a bias of 1023) |
| Mantissa | 52 bits | Fractional part of the significand |

*Table 3. Bit allocation of the IEEE 754 double-precision format*

```text
  63    62            52 51                              0
┌──────┬────────────────┬─────────────────────────────────┐
│ sign │ exponent (11)  │          fraction (52)          │
└──────┴────────────────┴─────────────────────────────────┘
```

*Figure 1. The 64-bit layout of a double-precision number*

Because binary floating-point numbers are normalized to have a leading 1 (i.e., $`1.xxxxx \times 2^e`$), the leading 1 is **implicit** and does not need to be stored. This gives an effective 53 bits of precision in the mantissa.

> **[Computer Architecture]** The floating-point unit of a CPU adds and multiplies IEEE 754 bit patterns directly in hardware and rounds each result to the 53-bit mantissa before storing it. This hardware rounding is where the roundoff error studied in this chapter originates, and the 53-bit mantissa corresponds to about 15-17 significant decimal digits.

### 5.4 Machine Epsilon

**Machine epsilon** ($`\varepsilon_{\text{mach}}`$) is the gap between 1 and the next larger floating-point number. Equivalently, it is the smallest power of 2 that satisfies

```math
1 + \varepsilon_{\text{mach}} \neq 1
```

in floating-point arithmetic. For IEEE 754 double precision:

```math
\varepsilon_{\text{mach}} = 2^{-52} \approx 2.22 \times 10^{-16}
```

Machine epsilon represents the **relative precision** of the floating-point system. It tells us the smallest relative change that a floating-point number can resolve. Under the default round-to-nearest mode, adding to 1 any number whose size relative to 1 is at most $`\varepsilon_{\text{mach}}/2`$ gives a result that rounds back to 1.

This means that double precision provides approximately **15-17 significant decimal digits** of precision. Beyond that, the digits are unreliable due to roundoff.

### 5.5 Truncation Error vs. Roundoff Error

These two types of error have fundamentally different origins.

| Property | Truncation Error | Roundoff Error |
|:---------|:----------------|:----------------|
| **Source** | Approximating a mathematical procedure (e.g., truncating an infinite Taylor series after finitely many terms) | Representing numbers with finite digits in computer memory |
| **Depends on** | Algorithm and step size | Computer hardware and number format |
| **Behavior with more terms/steps** | Generally decreases | May accumulate |
| **Example** | Using $`e^x \approx 1 + x + x^2/2`$ instead of the full infinite series | Storing $`1/3 = 0.333...`$ as $`0.333333333333333`$ |

*Table 4. Comparison of truncation error and roundoff error*

Understanding both types of error is essential because they often interact. Reducing the step size to decrease truncation error requires more arithmetic operations, which can increase accumulated roundoff error.

---

<br>

## 6. Problems with Floating-Point Arithmetic

### 6.1 Subtractive Cancellation

**Subtractive cancellation** (also called catastrophic cancellation) occurs when two nearly equal numbers are subtracted. The result loses most of its significant digits.

```math
0.12345678 - 0.12345671 = 0.00000007
```

The two operands each have 8 significant digits, but the result has **only 1 significant digit**. The 7 leading digits cancel, and all that remains is the least reliable last digit of each operand. This is one of the most dangerous sources of numerical error.

Subtractive cancellation is particularly insidious because it can occur in the middle of a longer computation, silently corrupting the result. Any algorithm that involves subtracting quantities that are expected to be close in value should be carefully examined for this issue.

### 6.2 Adding Large and Small Numbers

When a very large number and a very small number are added, the small number may be **completely lost** due to the finite precision of floating-point arithmetic.

```math
10^{16} + 1 = 10^{16} \quad \text{(in double precision)}
```

This happens because the number $`1`$ is $`10^{-16}`$ relative to $`10^{16}`$, which is less than half of machine epsilon (about $`1.11 \times 10^{-16}`$). When the computer aligns the two numbers to the same exponent for addition, the digits of the smaller number fall outside the range the mantissa can hold and are discarded by rounding. Near $`10^{16}`$, adjacent double-precision numbers are 2 apart.

Applying the same estimate to $`10^{15}`$, the size of a number that can change the result when added to $`10^{15}`$ is roughly $`10^{15} \times 1.11 \times 10^{-16} \approx 0.111`$. Therefore, in $`10^{15} + 1`$ the 1 is above this threshold and is preserved ($`10^{15} + 1 \neq 10^{15}`$), whereas in $`10^{16} + 1`$ the 1 falls below the relative precision and is lost entirely.

### 6.3 Smearing

**Smearing** occurs when an alternating series is summed and the intermediate partial sums are much larger in magnitude than the final answer. Each intermediate addition or subtraction introduces a roundoff error that is proportional to the magnitude of the intermediate result. Since the intermediate results are much larger than the final answer, these roundoff errors can be comparable to or larger than the final result itself.

For example, computing $`e^{-20}`$ directly with the Maclaurin series gives the following.

```math
e^{-20} = 1 - 20 + \frac{400}{2!} - \frac{8000}{3!} + \cdots
```

The terms grow very large before eventually decreasing (the largest term is about $`4.3 \times 10^7`$), and the alternating signs mean that massive cancellation occurs. The intermediate partial sums reach magnitudes of about $`2.2 \times 10^7`$, while the final answer is approximately $`2.06 \times 10^{-9}`$. The roundoff errors accumulated during the intermediate sums overwhelm the tiny final result, so adding the terms in order in double precision gives a value far from the true one.

A better approach is to compute $`e^{20}`$ (no smearing because all terms are positive) and then take the reciprocal: $`e^{-20} = 1/e^{20}`$.

---

<br>

## 7. Taylor Series Preview

The **Taylor series** expansion of a function $`f(x)`$ about the point $`x_i`$ is:

```math
f(x_{i+1}) = f(x_i) + f'(x_i)h + \frac{f''(x_i)}{2!}h^2 + \frac{f'''(x_i)}{3!}h^3 + \cdots
```

where $`h = x_{i+1} - x_i`$ is the step size.

This infinite series is exact, but in practice we must truncate it after a finite number of terms, introducing **truncation error**. The Taylor series is the foundation for understanding truncation errors and is treated in detail in Chapter 4.

The key points that Chapter 4 develops are the following.

- The **remainder term** (truncation error) depends on the step size $`h`$ raised to some power. Keeping more terms in the series means the error decreases faster as $`h`$ shrinks.
- Taylor series provides the theoretical basis for most numerical methods. Taylor expansions are used to derive or analyze the error of finite differences, interpolation, integration, and ODE solvers.
- The trade-off between truncation error (fewer terms or larger $`h`$) and roundoff error (more computation with smaller $`h`$) is a central theme in numerical analysis.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Significant Figures | They are the number of reliable digits in a number, and leading zeros do not count. |
| Accuracy and Precision | Accuracy is closeness to the true value, while precision is closeness of repeated results to each other (repeatability). |
| True Error | It is defined as $`E_t = \text{true} - \text{approx}`$, and its sign gives the direction of the error. |
| Relative Error | $`\varepsilon_t = (E_t / \text{true}) \times 100\%`$ expresses the error as a percentage of the true value. |
| Approximate Error | $`\varepsilon_a`$ is computed from the change between successive approximations when the true value is unknown. |
| Scarborough Criterion | If $`\vert\varepsilon_a\vert < (0.5 \times 10^{2-n})\%`$, the result is correct to at least $`n`$ significant figures. |
| IEEE 754 Double Precision | It splits 64 bits into 1 sign bit, 11 exponent bits, and 52 mantissa bits, and $`\varepsilon_{\text{mach}} = 2^{-52} \approx 2.22 \times 10^{-16}`$. |
| Subtractive Cancellation | Subtracting two nearly equal numbers cancels the leading digits and loses most significant digits. |
| Smearing | In an alternating series whose partial sums are much larger than the final answer, accumulated roundoff overwhelms the result. |
| Taylor Series | It is the foundation for truncation error analysis and is treated in detail in Chapter 4. |

---

<br>

## Self-Check Questions

1. **Counting significant figures:** Find the number of significant figures in $`0.004050`$, $`30.0400`$, and $`7.10 \times 10^4`$.

   <details>
   <summary>Answer</summary>

   $`0.004050`$ has 4 (the digits 4, 0, 5, 0), because leading zeros do not count and the trailing zero after the decimal point does. $`30.0400`$ has 6 (3, 0, 0, 4, 0, 0), because captive zeros and trailing zeros after the decimal point are all significant. $`7.10 \times 10^4`$ has 3, from the mantissa $`7.10`$.

   </details>

2. **Scarborough criterion:** What $`\varepsilon_s`$ should we choose to obtain a result correct to 4 significant figures? If we stop at $`|\varepsilon_a| = 0.003\%`$, how many figures are guaranteed?

   <details>
   <summary>Answer</summary>

   For $`n = 4`$, $`\varepsilon_s = 0.5 \times 10^{2-4}\% = 0.005\%`$. Since $`0.003\%`$ is below the 4-figure threshold $`0.005\%`$ but above the 5-figure threshold $`0.0005\%`$, the Scarborough criterion guarantees 4 figures.

   </details>

3. **True and relative error:** Approximating $`\pi \approx 3.14159265`$ by $`22/7 \approx 3.14285714`$, find $`E_t`$ and $`\varepsilon_t`$.

   <details>
   <summary>Answer</summary>

   $`E_t = 3.14159265 - 3.14285714 = -0.00126449`$. $`\varepsilon_t = -0.00126449 / 3.14159265 \times 100\% \approx -0.0402\%`$, and the negative sign means that $`22/7`$ overestimates $`\pi`$.

   </details>

4. **Machine epsilon and addition:** In double precision, what are the results of $`1 + 2^{-53}`$ and $`1 + 2^{-52}`$? Explain why $`10^{16} + 1 = 10^{16}`$.

   <details>
   <summary>Answer</summary>

   $`1 + 2^{-52}`$ is the next floating-point number after 1, so it is represented exactly and differs from 1. $`1 + 2^{-53}`$ lies exactly halfway between 1 and $`1 + 2^{-52}`$, so round-half-to-even gives 1. Near $`10^{16}`$, adjacent double-precision numbers are 2 apart, so $`10^{16} + 1`$ also lies exactly halfway between two numbers and is rounded to $`10^{16}`$ by the same rule.

   </details>

5. **Avoiding smearing:** Explain why summing the Maclaurin series of $`e^{-20}`$ directly is inaccurate, and give a better way to compute it.

   <details>
   <summary>Answer</summary>

   The intermediate partial sums reach about $`2.2 \times 10^7`$, so the roundoff error of each addition is about $`2.2 \times 10^7 \times 1.11 \times 10^{-16} \approx 2.4 \times 10^{-9}`$, comparable to the final answer $`2.06 \times 10^{-9}`$. Computing $`e^{20}`$, whose terms are all positive, by the series and then taking the reciprocal $`1/e^{20}`$ avoids smearing.

   </details>
