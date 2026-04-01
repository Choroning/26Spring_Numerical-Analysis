# Chapter 16 Lecture -- Fourier Analysis

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 16

> **Prerequisites**: [Calculus] Trigonometric functions, integrals (Ch 1-15).
>
> **Learning Objectives**:
> 1. Apply discrete Fourier transform (DFT) for signal analysis
> 2. Implement Fast Fourier Transform (FFT) algorithm
> 3. Interpret frequency spectra and power spectral density

---

<br>

## Table of Contents

- [1. Overview and Goals](#1-overview-and-goals)
- [2. Even and Odd Functions](#2-even-and-odd-functions)
  - [2.1 Even Functions](#21-even-functions)
  - [2.2 Odd Functions](#22-odd-functions)
  - [2.3 Decomposition of Any Function](#23-decomposition-of-any-function)
- [3. Periodic Functions](#3-periodic-functions)
  - [3.1 Definition](#31-definition)
  - [3.2 Example: Square Wave](#32-example-square-wave)
- [4. Fourier Series (Period $2\pi$)](#4-fourier-series-period-2pi)
  - [4.1 Fourier's Claim](#41-fouriers-claim)
  - [4.2 Fourier Coefficients -- Derivation](#42-fourier-coefficients----derivation)
  - [4.3 Orthogonality Relations](#43-orthogonality-relations)
  - [4.4 Coefficient Formulas](#44-coefficient-formulas)
- [5. Fourier Series of the Square Wave](#5-fourier-series-of-the-square-wave)
- [6. Period Other Than $2\pi$](#6-period-other-than-2pi)
- [7. Fourier Integral (Non-Periodic Functions)](#7-fourier-integral-non-periodic-functions)
- [8. Fourier Transform](#8-fourier-transform)
  - [8.1 Definition](#81-definition)
  - [8.2 Properties of the Fourier Transform](#82-properties-of-the-fourier-transform)
- [9. Curve Fitting with Sinusoidal Functions](#9-curve-fitting-with-sinusoidal-functions)
  - [9.1 Sinusoidal Function Parameters](#91-sinusoidal-function-parameters)
  - [9.2 Phase Angle](#92-phase-angle)
  - [9.3 Alternative Form](#93-alternative-form)
- [10. Least-Squares Fit of a Sinusoid](#10-least-squares-fit-of-a-sinusoid)
  - [10.1 Problem Setup](#101-problem-setup)
  - [10.2 Normal Equations](#102-normal-equations)
  - [10.3 Simplification with Equispaced Points](#103-simplification-with-equispaced-points)
  - [10.4 Closed-Form Solution](#104-closed-form-solution)
- [11. Continuous Fourier Series (General Model)](#11-continuous-fourier-series-general-model)
  - [11.1 General Fourier Series](#111-general-fourier-series)
  - [11.2 Euler's Formula and Complex Form](#112-eulers-formula-and-complex-form)
- [12. Frequency and Time Domain](#12-frequency-and-time-domain)
- [13. Fourier Integral and Fourier Transformation (Revisited)](#13-fourier-integral-and-fourier-transformation-revisited)
- [14. Discrete Fourier Transform (DFT)](#14-discrete-fourier-transform-dft)
- [15. Nyquist Frequency](#15-nyquist-frequency)
- [16. Fast Fourier Transform (FFT)](#16-fast-fourier-transform-fft)
- [17. Power Spectrum](#17-power-spectrum)
- [18. Summary Table](#18-summary-table)

---

<br>

## 1. Overview and Goals

**From N17 (Fourier Series):**

- Understand Fourier Series and Fourier Integral

**From N18 (Fourier Analysis):**

- Understand sinusoids (sine and cosine waves)
- Use least-square regression to fit a sinusoid to data
- Understand Euler's formula
- Analyze signals in frequency domain
- Discrete Fourier transform
- Nyquist frequency
- Fast Fourier transform
- Power spectrum

> **[Context]** Fourier analysis is one of the most powerful tools in numerical methods. It allows us to decompose any periodic (or even non-periodic) function into a sum of simple sinusoidal components. This has applications in signal processing, image compression, solving differential equations, and more.

---

<br>

## 2. Even and Odd Functions

### 2.1 Even Functions

A function $f(x)$ is **even** if:

$$f(-x) = f(x)$$

The graph of $f(x)$ is **symmetric about** $x = 0$ (the y-axis).

**Examples:** $\cos(x)$, $x^2$, $|x|$

### 2.2 Odd Functions

A function $f(x)$ is **odd** if:

$$f(-x) = -f(x)$$

The graph of $f(x)$ is **antisymmetric about** $x = 0$ (origin symmetry).

**Examples:** $\sin(x)$, $x^3$, $x$

### 2.3 Decomposition of Any Function

Any function can be decomposed into the sum of an even function and an odd function:

$$f(x) = \frac{1}{2}f(x) + \frac{1}{2}f(x)$$

$$= \underbrace{\frac{1}{2}\bigl[f(x) + f(-x)\bigr]}_{\text{Even part}} + \underbrace{\frac{1}{2}\bigl[f(x) - f(-x)\bigr]}_{\text{Odd part}}$$

> **[Key Idea]** This decomposition is fundamental to understanding Fourier series, where cosine terms capture the even part and sine terms capture the odd part of a function.

---

<br>

## 3. Periodic Functions

### 3.1 Definition

A function $f(x)$ is **periodic** with period $T$ if:

$$f(x) = f(x + T) \quad \forall x$$

We always choose the **smallest possible** positive $T$.

**Example:** $\sin(x)$ is periodic with period $T = 2\pi$, since:

$$\sin(x + 2\pi) = \sin(x)\cos(2\pi) + \cos(x)\sin(2\pi) = \sin(x)$$

Note that $\sin(x) = \sin(x + 2\pi) = \sin(x + 4\pi) = \cdots$, but we choose $T = 2\pi$ as the fundamental period.

### 3.2 Example: Square Wave

Consider the square wave with $2\pi$ period:

$$f(x) = \begin{cases} 0, & -\pi < x \le 0 \\ 1, & 0 < x \le \pi \end{cases}$$

This function repeats with period $2\pi$ and is a classic example used to demonstrate Fourier series.

---

<br>

## 4. Fourier Series (Period $2\pi$)

### 4.1 Fourier's Claim

**Fourier (1768--1830)** claimed that an arbitrary function defined over $[-\pi, \pi]$ could be represented in the form:

$$f(x) = \sum_{n=0}^{\infty} (a_n \cos nx + b_n \sin nx)$$

Rewriting with the conventional $\frac{a_0}{2}$ notation:

$$f(x) = \frac{a_0}{2} + \sum_{n=1}^{\infty} (a_n \cos nx + b_n \sin nx)$$

> **[Why $a_0/2$?]** When $n = 0$, we have $a_0 \cos(0) + b_0 \sin(0) = a_0 \cdot 1 + 0 = a_0$. We introduce the factor $\frac{1}{2}$ for convenience so that the formula for $a_0$ has the same form as $a_n$ for $n \ge 1$.

Notice the structure:

$$f(x) = \underbrace{\left(\frac{a_0}{2} + \sum_{n=1}^{\infty} a_n \cos nx \right)}_{\text{Even part}} + \underbrace{\left(\sum_{n=1}^{\infty} b_n \sin nx \right)}_{\text{Odd part}}$$

### 4.2 Fourier Coefficients -- Derivation

**Finding $a_0$:** Integrate both sides over $[-\pi, \pi]$:

$$\int_{-\pi}^{\pi} f(x)\,dx = \int_{-\pi}^{\pi} \frac{a_0}{2}\,dx + \sum_{n=1}^{\infty} a_n \underbrace{\int_{-\pi}^{\pi} \cos nx\,dx}_{= 0} + \sum_{n=1}^{\infty} b_n \underbrace{\int_{-\pi}^{\pi} \sin nx\,dx}_{= 0}$$

$$= \frac{a_0}{2} \cdot 2\pi = a_0 \pi$$

$$\therefore \quad a_0 = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x)\,dx$$

**Finding $a_n$:** Multiply both sides by $\cos mx$ and integrate:

$$\int_{-\pi}^{\pi} f(x)\cos mx\,dx = \underbrace{\int_{-\pi}^{\pi} \frac{a_0}{2}\cos mx\,dx}_{= 0} + \sum_{n=1}^{\infty} a_n \int_{-\pi}^{\pi} \cos nx \cos mx\,dx + \sum_{n=1}^{\infty} b_n \underbrace{\int_{-\pi}^{\pi} \sin nx \cos mx\,dx}_{= 0}$$

**Finding $b_n$:** Multiply both sides by $\sin mx$ and integrate:

$$\int_{-\pi}^{\pi} f(x)\sin mx\,dx = \underbrace{\int_{-\pi}^{\pi} \frac{a_0}{2}\sin mx\,dx}_{= 0} + \sum_{n=1}^{\infty} a_n \underbrace{\int_{-\pi}^{\pi} \cos nx \sin mx\,dx}_{= 0} + \sum_{n=1}^{\infty} b_n \int_{-\pi}^{\pi} \sin nx \sin mx\,dx$$

### 4.3 Orthogonality Relations

The key to deriving the coefficients relies on **orthogonality** of trigonometric functions.

**Product-to-sum identities:**

$$\cos a \cos b = \frac{1}{2}[\cos(a+b) + \cos(a-b)]$$

$$\sin a \cos b = \frac{1}{2}[\sin(a+b) + \sin(a-b)]$$

$$\sin a \sin b = -\frac{1}{2}[\cos(a+b) - \cos(a-b)]$$

**Orthogonality results on $[-\pi, \pi]$:**

$$\int_{-\pi}^{\pi} \cos nx \cos mx\,dx = \begin{cases} 0 & \text{if } n \neq m \\ 2\pi & \text{if } n = m = 0 \\ \pi & \text{if } n = m \ge 1 \end{cases}$$

$$\int_{-\pi}^{\pi} \sin nx \cos mx\,dx = 0 \quad \text{for all } n, m$$

$$\int_{-\pi}^{\pi} \sin nx \sin mx\,dx = \begin{cases} 0 & \text{if } n \neq m \\ \pi & \text{if } n = m \ge 1 \end{cases}$$

### 4.4 Coefficient Formulas

Using the orthogonality relations:

$$\boxed{a_n = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x) \cos nx\,dx, \quad n = 0, 1, 2, \ldots}$$

$$\boxed{b_n = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x) \sin nx\,dx, \quad n = 1, 2, 3, \ldots}$$

---

<br>

## 5. Fourier Series of the Square Wave

For the square wave $f(x) = \begin{cases} 0, & -\pi < x \le 0 \\ 1, & 0 < x \le \pi \end{cases}$:

**Compute $a_0$:**

$$a_0 = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x)\,dx = \frac{1}{\pi} \cdot \pi = 1$$

**Compute $a_n$:**

$$a_n = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x)\cos nx\,dx = \frac{1}{\pi} \int_{0}^{\pi} \cos nx\,dx = \frac{1}{\pi} \left[\frac{1}{n}\sin nx\right]_0^{\pi} = 0$$

**Compute $b_n$:**

$$b_n = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x)\sin nx\,dx = \frac{1}{\pi} \int_0^{\pi} \sin nx\,dx = \frac{1}{\pi}\left[-\frac{1}{n}\cos nx\right]_0^{\pi}$$

$$= -\frac{1}{n\pi}(\cos n\pi - 1) = \begin{cases} \frac{2}{n\pi} & \text{for odd } n \\ 0 & \text{for even } n \end{cases}$$

**Result:**

$$f(x) = \frac{1}{2} + \sum_{\substack{n=1,3,5,\ldots}} \frac{2}{n\pi} \sin nx = \frac{1}{2} + \frac{2}{\pi} \sum_{n=1}^{\infty} \frac{\sin(2n-1)x}{2n-1}$$

As $n \to \infty$, $\frac{1}{n} \to 0$, so the Fourier series **converges**.

---

<br>

## 6. Period Other Than $2\pi$

For a function with period $2\ell$ (instead of $2\pi$), use the change of variable:

$$x' = \pi \cdot \frac{x}{\ell}$$

This maps the interval $[-\ell, \ell]$ to $[-\pi, \pi]$. The Fourier series becomes:

$$\boxed{f(x) = \frac{a_0}{2} + \sum_{n=1}^{\infty} a_n \cos\frac{n\pi}{\ell}x + \sum_{n=1}^{\infty} b_n \sin\frac{n\pi}{\ell}x}$$

with coefficients:

$$a_n = \frac{1}{\ell} \int_{-\ell}^{\ell} f(x) \cos\frac{n\pi}{\ell}x\,dx, \quad n = 0, 1, 2, \ldots$$

$$b_n = \frac{1}{\ell} \int_{-\ell}^{\ell} f(x) \sin\frac{n\pi}{\ell}x\,dx, \quad n = 1, 2, 3, \ldots$$

---

<br>

## 7. Fourier Integral (Non-Periodic Functions)

**Question:** Can $f(x) = e^{-x^2}$ be expanded in a Fourier series?

**Answer:** No -- $f(x)$ is NOT periodic.

**Solution:** Take $\ell \to \infty$. As $\ell$ increases, the frequency spectrum $\frac{n\pi}{\ell}$ becomes finer and finer, tending to a **continuous** spectrum.

| $\ell$ | Frequencies $\frac{n\pi}{\ell}$ |
|:---:|:---|
| $\pi$ | $1, 2, 3, 4, \ldots$ |
| $2\pi$ | $\frac{1}{2}, 1, \frac{3}{2}, 2, \ldots$ |
| $10\pi$ | $0.1, 0.2, 0.3, 0.4, \ldots$ |

In the limit, the discrete sum becomes a continuous integral:

$$f(x) = \frac{a_0}{2} + \sum_{n=1}^{\infty} (a_n \cos nx + b_n \sin nx)$$

$$\Downarrow$$

$$\boxed{f(x) = \int_0^{\infty} \bigl[a(\omega)\cos\omega x + b(\omega)\sin\omega x\bigr]\,d\omega}$$

This is the **Fourier Integral** representation.

---

<br>

## 8. Fourier Transform

### 8.1 Definition

Restating the Fourier Integral using complex exponentials:

$$f(x) = \frac{1}{2\pi} \int_{-\infty}^{\infty} e^{i\omega x} \underbrace{\int_{-\infty}^{\infty} f(\xi) e^{-i\omega\xi}\,d\xi}_{\hat{f}(\omega)}\,d\omega = \frac{1}{2\pi} \int_{-\infty}^{\infty} \hat{f}(\omega)\,e^{i\omega x}\,d\omega$$

The **Fourier Transform** and its **inverse**:

$$\boxed{\hat{f}(\omega) = \mathcal{F}\{f(x)\} = \int_{-\infty}^{\infty} f(\xi)\,e^{-i\omega\xi}\,d\xi}$$

$$\boxed{f(x) = \mathcal{F}^{-1}\{\hat{f}(\omega)\} = \frac{1}{2\pi} \int_{-\infty}^{\infty} \hat{f}(\omega)\,e^{i\omega x}\,d\omega}$$

### 8.2 Properties of the Fourier Transform

**1) Linearity:**

$$\mathcal{F}\{\alpha f(x) + \beta g(x)\} = \alpha\,\mathcal{F}\{f\} + \beta\,\mathcal{F}\{g\}$$

$$\mathcal{F}^{-1}\{\alpha\hat{f} + \beta\hat{g}\} = \alpha\,\mathcal{F}^{-1}\{\hat{f}\} + \beta\,\mathcal{F}^{-1}\{\hat{g}\}$$

**2) Derivative Property:**

$$\mathcal{F}\left\{\frac{d}{dx}f(x)\right\} = i\omega\,\mathcal{F}\{f\}$$

provided $f(x) \to 0$ as $x \to \pm\infty$.

*Proof (by integration by parts):*

$$\int_{-\infty}^{\infty} \frac{d}{d\xi}f(\xi)\,e^{-i\omega\xi}\,d\xi = \underbrace{\left[f(\xi)e^{-i\omega\xi}\right]_{-\infty}^{\infty}}_{= 0} - (-i\omega)\int_{-\infty}^{\infty} f(\xi)\,e^{-i\omega\xi}\,d\xi = i\omega\,\mathcal{F}\{f\}$$

**General $n$-th derivative:**

$$\mathcal{F}\left\{\frac{d^n}{dx^n}f(x)\right\} = (i\omega)^n\,\mathcal{F}\{f\}$$

> **[Significance]** The derivative property converts differentiation in the spatial/time domain into multiplication in the frequency domain. This is why Fourier transforms are extremely useful for solving differential equations.

**3) Convolution:**

$$f * g := \int_{-\infty}^{\infty} f(x - \xi)\,g(\xi)\,d\xi = \mathcal{F}^{-1}\{\hat{g}(\omega)\,\hat{f}(\omega)\}$$

That is, **convolution in the spatial domain** corresponds to **multiplication in the frequency domain**.

**4) $x$-shift (Spatial Shift):**

$$\mathcal{F}^{-1}\{e^{-ia\omega}\,\hat{f}(\omega)\} = f(x - a)$$

Multiplying by $e^{-ia\omega}$ in frequency domain shifts the function by $a$ in the spatial domain.

**5) $\omega$-shift (Frequency Shift):**

$$\mathcal{F}^{-1}\{\hat{f}(\omega - a)\} = e^{iax}\,f(x)$$

Shifting in frequency domain corresponds to modulation by $e^{iax}$ in the spatial domain.

---

<br>

## 9. Curve Fitting with Sinusoidal Functions

### 9.1 Sinusoidal Function Parameters

Consider a periodic function $f(t) = f(t + T)$, where $T$ is a constant (the period).

A **sinusoidal function** is any waveform that can be described as a sine or cosine:

$$f(t) = A_0 + C_1 \cos(\omega_0 t + \theta)$$

where:

| Parameter | Meaning |
|:---:|:---|
| $A_0$ | **Mean value** (vertical offset) |
| $C_1$ | **Amplitude** (height of oscillation) |
| $\omega_0$ | **Angular frequency** (radians/time) |
| $\theta$ | **Phase angle** (radians) |

The relationships between frequency quantities:

$$\omega_0 = 2\pi f, \qquad f = \frac{1}{T}$$

where $f$ is the **ordinary frequency** in cycles/time $[\text{Hz}]$.

### 9.2 Phase Angle

The phase angle $\theta$ represents the distance (in radians) from $t = 0$ to the point at which the cosine function begins a new cycle.

For $\theta > 0$:

- $\cos(\omega_0 t - \theta)$: **lagging** phase angle (shifted right)
- $\cos(\omega_0 t + \theta)$: **leading** phase angle (shifted left)

Also note the relationship between sine and cosine:

$$\sin(t + \pi/2) = \cos(t), \qquad \cos(t - \pi/2) = \sin(t)$$

### 9.3 Alternative Form

We can pull $\theta$ out of the cosine function using the angle addition formula:

$$\cos(\omega_0 t + \theta) = \cos(\omega_0 t)\cos(\theta) - \sin(\omega_0 t)\sin(\theta)$$

Therefore:

$$f(t) = A_0 + C_1\cos(\omega_0 t + \theta)$$

$$= A_0 + C_1\cos\theta\,\cos(\omega_0 t) - C_1\sin\theta\,\sin(\omega_0 t)$$

$$\boxed{= A_0 + A_1\cos(\omega_0 t) + B_1\sin(\omega_0 t)}$$

where:

$$A_1 = C_1\cos\theta, \qquad B_1 = -C_1\sin\theta$$

From these we can recover:

$$\theta = \arctan\left(-\frac{B_1}{A_1}\right), \qquad C_1 = \sqrt{A_1^2 + B_1^2}$$

---

<br>

## 10. Least-Squares Fit of a Sinusoid

### 10.1 Problem Setup

Given data points $(t_i, y_i)$, fit the model:

$$y = A_0 + A_1\cos(\omega_0 t) + B_1\sin(\omega_0 t) = \hat{y}(t)$$

Interpret as a linear regression problem with basis functions:

$$f_0 = 1, \quad f_1 = \cos(\omega_0 t), \quad f_2 = \sin(\omega_0 t)$$

$$\hat{y}(t) = \beta_0 f_0 + \beta_1 f_1 + \beta_2 f_2$$

### 10.2 Normal Equations

Minimize the sum of squared errors:

$$SSE = \sum (y_i - A_0 - A_1\cos(\omega_0 t_i) - B_1\sin(\omega_0 t_i))^2$$

Taking partial derivatives and setting them to zero:

$$\frac{\partial SSE}{\partial A_0} = 0, \quad \frac{\partial SSE}{\partial A_1} = 0, \quad \frac{\partial SSE}{\partial B_1} = 0$$

This yields the normal equation system:

$$\begin{pmatrix} \sum y_i \\ \sum y_i \cos(\omega_0 t) \\ \sum y_i \sin(\omega_0 t) \end{pmatrix} = \begin{pmatrix} n & \sum\cos(\omega_0 t) & \sum\sin(\omega_0 t) \\ \sum\cos(\omega_0 t) & \sum\cos^2(\omega_0 t) & \sum\cos(\omega_0 t)\sin(\omega_0 t) \\ \sum\sin(\omega_0 t) & \sum\sin(\omega_0 t)\cos(\omega_0 t) & \sum\sin^2(\omega_0 t) \end{pmatrix} \begin{pmatrix} A_0 \\ A_1 \\ B_1 \end{pmatrix}$$

### 10.3 Simplification with Equispaced Points

Collect $n$ observations within a period $T$ at equally spaced points: $T = n\,\Delta t$.

For **equispaced** points, the following orthogonality properties hold:

$$\sum_{i=1}^{n} \sin(\omega_0 t_i) = 0, \qquad \sum_{i=1}^{n} \cos(\omega_0 t_i) = 0$$

$$\sum_{i=1}^{n} \sin^2(\omega_0 t_i) = \frac{n}{2}, \qquad \sum_{i=1}^{n} \cos^2(\omega_0 t_i) = \frac{n}{2}$$

$$\sum_{i=1}^{n} \cos(\omega_0 t_i)\sin(\omega_0 t_i) = \frac{1}{2}\sum_{i=1}^{n}[\sin(2\omega_0 t_i) - \sin 0] = 0$$

The normal equations simplify to a **diagonal** system:

$$\begin{pmatrix} \sum y_i \\ \sum y_i\cos(\omega_0 t) \\ \sum y_i\sin(\omega_0 t) \end{pmatrix} = \begin{pmatrix} n & 0 & 0 \\ 0 & \frac{n}{2} & 0 \\ 0 & 0 & \frac{n}{2} \end{pmatrix} \begin{pmatrix} A_0 \\ A_1 \\ B_1 \end{pmatrix}$$

### 10.4 Closed-Form Solution

$$\boxed{A_0 = \frac{1}{n}\sum y_i = \bar{y}}$$

$$\boxed{A_1 = \frac{2}{n}\sum y_i \cos(\omega_0 t_i)}$$

$$\boxed{B_1 = \frac{2}{n}\sum y_i \sin(\omega_0 t_i)}$$

> **[Key Insight]** With equispaced data points, the normal equations decouple completely, giving simple closed-form expressions for the Fourier coefficients. This is a major computational advantage.

---

<br>

## 11. Continuous Fourier Series (General Model)

### 11.1 General Fourier Series

Extending to the general model with $m$ harmonics ($2m + 1$ coefficients):

$$f(t) = a_0 + \sum_{k=1}^{\infty} \bigl(a_k \cos(k\omega_0 t) + b_k \sin(k\omega_0 t)\bigr)$$

where:

- $\omega_0 = \frac{2\pi}{T}$ is the **fundamental frequency**
- $2\omega_0, 3\omega_0, \ldots, m\omega_0, \ldots$ are the **harmonics** (constant multiples of $\omega_0$)

The coefficients for the general model:

$$A_0 = \frac{1}{n}\sum y_i = \bar{y}$$

$$A_j = \frac{2}{n}\sum y_i \cos(j\omega_0 t_i), \quad j = 1, 2, \ldots, m$$

$$B_j = \frac{2}{n}\sum y_i \sin(j\omega_0 t_i), \quad j = 1, 2, \ldots, m$$

> **[Regression vs. Collocation]** If $n > 2m + 1$ (more data points than coefficients), we compute the coefficients to fit data in the **regression** (least-squares) sense. If $n = 2m + 1$ (the number of data points equals the number of Fourier coefficients), this is called **collocation** -- the Fourier series interpolates the data exactly.

### 11.2 Euler's Formula and Complex Form

**Euler's formula:**

$$e^{i\theta} = \cos\theta + i\sin\theta$$

This allows us to express the Fourier series in **complex exponential form**:

$$f(t) = a_0 + \sum_{k=1}^{\infty}\bigl(a_k\cos(k\omega_0 t) + b_k\sin(k\omega_0 t)\bigr) = \sum_{k=-\infty}^{\infty} \tilde{C}_k\,e^{ik\omega_0 t}$$

where the complex Fourier coefficients are:

$$\tilde{C}_k = \frac{1}{T}\int_{-T/2}^{T/2} f(t)\,e^{-ik\omega_0 t}\,dt$$

---

<br>

## 12. Frequency and Time Domain

The **frequency domain** provides an alternative perspective for characterizing the behavior of oscillating functions.

- **Time domain**: the usual representation $f(t)$ -- amplitude vs. time
- **Frequency domain**: the Fourier coefficients $\hat{f}(\omega)$ -- amplitude vs. frequency

Fourier analysis transforms a signal from the time domain to the frequency domain, revealing which frequencies are present and their relative strengths.

---

<br>

## 13. Fourier Integral and Fourier Transformation (Revisited)

Transition from a periodic function to a non-periodic function by taking $T \to \infty$:

**Fourier Integral (inverse transform):**

$$f(t) = \frac{1}{2\pi} \int_{-\infty}^{\infty} F(\omega)\,e^{i\omega t}\,d\omega$$

**Fourier Transform:**

$$F(\omega) = \int_{-\infty}^{\infty} f(t)\,e^{-i\omega t}\,dt$$

---

<br>

## 14. Discrete Fourier Transform (DFT)

Collect a **finite set** of function values. An interval from $0$ to $T$ can be divided into $n$ equispaced subintervals with $\Delta t = \frac{T}{n}$.

Discrete samples: $f_j = f(t = t_j)$, where $t_j = t_0 + j\,\Delta t$.

**DFT (Forward):**

$$\hat{f}_k = \mathcal{F}\{f\}_k \approx \sum_{j=0}^{n-1} f_j\,e^{-ik\omega_0 j}, \quad k = 0, 1, 2, \ldots, n-1$$

**Inverse DFT:**

$$f_j = \mathcal{F}^{-1}\{\hat{f}\}_j \approx \frac{1}{n}\sum_{k=0}^{n-1} \hat{f}_k\,e^{ik\omega_0 j}, \quad j = 0, 1, 2, \ldots, n-1$$

where $\omega_0 = \frac{2\pi}{n}$.

```python
import numpy as np

def dft(f):
    """Compute the Discrete Fourier Transform of a 1D array f."""
    n = len(f)
    F = np.zeros(n, dtype=complex)
    for k in range(n):
        for j in range(n):
            F[k] += f[j] * np.exp(-2j * np.pi * k * j / n)
    return F

def idft(F):
    """Compute the Inverse Discrete Fourier Transform of a 1D array F."""
    n = len(F)
    f = np.zeros(n, dtype=complex)
    for j in range(n):
        for k in range(n):
            f[j] += F[k] * np.exp(2j * np.pi * k * j / n)
    return f / n
```

> **[Complexity]** The naive DFT has time complexity $O(n^2)$ because each of the $n$ output values requires a sum over $n$ input values.

---

<br>

## 15. Nyquist Frequency

The **Nyquist frequency** is the highest frequency that can be measured in a signal. It equals **half the sampling frequency**:

$$f_{\max} = \frac{1}{2} f_s = \frac{1}{2\,\Delta t}$$

**Example:** Collect 100 samples ($n = 100$) at a sampling frequency of $f_s = 1000$ Hz:

- Sampling interval: $\Delta t = \frac{1}{f_s} = \frac{1}{1000} = 10^{-3}$ sec/sample
- Total time: $T = n \cdot \Delta t = 100 \times 10^{-3} = 0.1$ s
- **Nyquist frequency:** $f_{\max} = \frac{1}{2}f_s = 500$ Hz
- **Lowest detectable frequency:** $f_{\min} = \frac{1}{T} = \frac{1}{n\,\Delta t} = 10$ Hz

> **[Aliasing]** If a signal contains frequencies above the Nyquist frequency, those frequencies will be incorrectly represented as lower frequencies -- a phenomenon called **aliasing**. To avoid this, the sampling rate must be at least twice the highest frequency present in the signal (the **Nyquist-Shannon sampling theorem**).

---

<br>

## 16. Fast Fourier Transform (FFT)

The **Fast Fourier Transform** reduces the computational cost from $O(n^2)$ (DFT) to $O(n \log n)$ by exploiting:

1. **Symmetry** of trigonometric functions
2. **Periodicity** of complex exponentials

Key properties of the FFT:

- FFT splits the input into **even and odd indices** recursively
- FFT exploits **symmetry and periodicity** in the DFT formula
- FFT dramatically **cuts down repeated calculations**

```python
import numpy as np

# Using NumPy's built-in FFT
t = np.linspace(0, 1, 256, endpoint=False)
signal = np.sin(2 * np.pi * 10 * t) + 0.5 * np.sin(2 * np.pi * 20 * t)

# Compute FFT
F = np.fft.fft(signal)

# Compute frequencies
freqs = np.fft.fftfreq(len(t), d=t[1] - t[0])

# Inverse FFT to recover original signal
recovered = np.fft.ifft(F)
```

> **[Why $n = 2^p$?]** The Cooley-Tukey FFT algorithm works most efficiently when $n$ is a power of 2, as the recursive splitting into halves can be done $\log_2 n$ times evenly.

---

<br>

## 17. Power Spectrum

The **power spectrum** shows how the power of a signal is distributed across different frequencies.

**Discrete case** (power at discrete frequency $k\omega_0$):

$$P_k = |\tilde{C}_k|^2$$

**Continuous case** (power at continuous frequency $\omega$):

$$P(\omega) = |\hat{F}(\omega)|^2$$

```python
import numpy as np
import matplotlib.pyplot as plt

# Generate a signal with two frequency components
n = 256
dt = 0.001
t = np.arange(n) * dt
signal = np.sin(2 * np.pi * 50 * t) + 0.5 * np.sin(2 * np.pi * 120 * t)

# Compute FFT and power spectrum
F = np.fft.fft(signal)
freqs = np.fft.fftfreq(n, d=dt)
power = np.abs(F) ** 2

# Plot only positive frequencies
pos_mask = freqs >= 0
plt.plot(freqs[pos_mask], power[pos_mask])
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power |F(w)|^2")
plt.title("Power Spectrum")
plt.show()
```

---

<br>

## 18. Summary Table

| Topic | Key Formula | Notes |
|:---|:---|:---|
| **Fourier Series** ($2\pi$ period) | $f(x) = \frac{a_0}{2} + \sum_{n=1}^{\infty}(a_n\cos nx + b_n\sin nx)$ | Periodic functions on $[-\pi,\pi]$ |
| **Coefficients** ($2\pi$ period) | $a_n = \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\cos nx\,dx$, $b_n = \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\sin nx\,dx$ | Derived from orthogonality |
| **General Period** $2\ell$ | $f(x) = \frac{a_0}{2} + \sum a_n\cos\frac{n\pi}{\ell}x + \sum b_n\sin\frac{n\pi}{\ell}x$ | Change of variable $x' = \frac{\pi x}{\ell}$ |
| **Fourier Integral** | $f(x) = \int_0^{\infty}[a(\omega)\cos\omega x + b(\omega)\sin\omega x]\,d\omega$ | Non-periodic; $\ell \to \infty$ |
| **Fourier Transform** | $\hat{f}(\omega) = \int_{-\infty}^{\infty}f(\xi)e^{-i\omega\xi}\,d\xi$ | Complex exponential form |
| **Inverse FT** | $f(x) = \frac{1}{2\pi}\int_{-\infty}^{\infty}\hat{f}(\omega)e^{i\omega x}\,d\omega$ | Recovers original function |
| **Sinusoidal Fit** | $f(t) = A_0 + A_1\cos(\omega_0 t) + B_1\sin(\omega_0 t)$ | $C_1 = \sqrt{A_1^2+B_1^2}$, $\theta = \arctan(-B_1/A_1)$ |
| **Discrete Coefficients** | $A_0 = \bar{y}$, $A_j = \frac{2}{n}\sum y_i\cos(j\omega_0 t_i)$, $B_j = \frac{2}{n}\sum y_i\sin(j\omega_0 t_i)$ | Equispaced data |
| **Euler's Formula** | $e^{i\theta} = \cos\theta + i\sin\theta$ | Connects trig and exponential |
| **Complex Fourier Series** | $f(t) = \sum_{k=-\infty}^{\infty}\tilde{C}_k e^{ik\omega_0 t}$ | $\tilde{C}_k = \frac{1}{T}\int_{-T/2}^{T/2}f(t)e^{-ik\omega_0 t}\,dt$ |
| **DFT** | $\hat{f}_k = \sum_{j=0}^{n-1}f_j e^{-ik\omega_0 j}$ | $O(n^2)$ complexity |
| **FFT** | Same result as DFT | $O(n\log n)$ complexity |
| **Nyquist Frequency** | $f_{\max} = \frac{1}{2}f_s = \frac{1}{2\Delta t}$ | Maximum detectable frequency |
| **Power Spectrum** | $P_k = \lvert\tilde{C}_k\rvert^2$, $P(\omega) = \lvert\hat{F}(\omega)\rvert^2$ | Energy distribution over frequencies |
| **FT Linearity** | $\mathcal{F}\{\alpha f + \beta g\} = \alpha\mathcal{F}\{f\} + \beta\mathcal{F}\{g\}$ | Linear operator |
| **FT Derivative** | $\mathcal{F}\{f^{(n)}\} = (i\omega)^n\mathcal{F}\{f\}$ | Differentiation becomes multiplication |
| **FT Convolution** | $\mathcal{F}\{f*g\} = \hat{f}(\omega)\cdot\hat{g}(\omega)$ | Convolution becomes product |
| **FT $x$-shift** | $\mathcal{F}^{-1}\{e^{-ia\omega}\hat{f}(\omega)\} = f(x-a)$ | Spatial shift |
| **FT $\omega$-shift** | $\mathcal{F}^{-1}\{\hat{f}(\omega-a)\} = e^{iax}f(x)$ | Frequency modulation |
