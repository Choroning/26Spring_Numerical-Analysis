# HW9 Assignment Explanation — Fourier Series and Fourier Analysis

[Repository Guide](../../README.md) | [English](Assignment-Explanation.md) | [한국어](Assignment-Explanation.ko.md)

> **Last Updated:** 2026-11-16
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 16

> **Learning Objectives:**
> 1. Explain how the indices of the normalized DFT coefficients correspond to the entries of the frequency array
> 2. Complete a function that computes the frequencies and powers of the one-sided power spectrum
> 3. Interpret the frequency components and amplitudes of a signal from the peaks of its power spectrum
> 4. Explain how the sampling frequency and record length determine the frequency resolution and the Nyquist frequency

---

## Table of Contents

- [1. Assignment Overview](#1-assignment-overview)
- [2. Requirements](#2-requirements)
- [3. Problem: One-Sided Power Spectrum](#3-problem-one-sided-power-spectrum)
  - [3.1 Approach](#31-approach)
  - [3.2 Code Written](#32-code-written)
  - [3.3 Results and Verification](#33-results-and-verification)
  - [3.4 Results Used by Later Cells](#34-results-used-by-later-cells)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Assignment Overview

> **Try It First:** This document contains the answers and solutions. Solve the problems in the [assignment handout](Original/Assignment.md) yourself before reading on.

The assignment is to fill in the single blank of the lab notebook `Lab10-Fourier.ipynb` distributed by the instructor and run it. The notebook has six parts: time and frequency domains, Fourier series, the discrete Fourier transform (DFT), the power spectrum, signal filtering and noise removal, and vibration analysis, and the other cells were distributed with their code filled in. The function completed in the blank is used not only in the example of Section 4 but also in Exercise 4.1, the noisy signal of Section 5, and the vibration analysis of Section 6. The problem statement is in the [assignment handout](Original/Assignment.md), the distributed blank notebook is [Original/Lab10-Fourier.ipynb](Original/Lab10-Fourier.ipynb), and the notebook with the blank filled in and run is [Solution/Lab10-Fourier_Completed.ipynb](Solution/Lab10-Fourier_Completed.ipynb).

| Problem | Notebook section | Blanks to fill |
|:-----|:-----|:-----|
| Problem | 4.1 Definition of the power spectrum | `fP` and `Pyy` in the function `power_spectrum` |

*Table 1. The blank that makes up the assignment*

The discrete Fourier transform, the Nyquist frequency, and the power spectrum are covered in Sections 14, 15, and 17 of the [Chapter 16 lecture note](../../Lecture-Notes/Ch16_Fourier-Analysis/Concepts_Lecture.md). This document repeats only what the solution needs.

---

<br>

## 2. Requirements

The completed notebook must satisfy the following requirements. The rest of the cell stays as distributed, and only the blank is filled in.

1. `fP` must hold the frequencies of the one-sided spectrum, that is, the frequencies for $`k = 1, 2, \ldots, N/2`$, from $`\Delta f`$ to the Nyquist frequency.
2. `Pyy` must hold the powers $`|Y_k|^2`$ for the same values of $`k`$, with the same length as `fP` and with entries that refer to the same $`k`$ position by position.
3. The function takes as input `Y` (length $`N`$) and `freq` (the frequencies for indices $`1 \ldots N-1`$, length $`N-1`$) as returned by `compute_dft`. This format is followed by the example of Section 4, by Exercise 4.1, and by the cells of Sections 5 and 6 that build `Y` and `freq` themselves.

---

<br>

## 3. Problem: One-Sided Power Spectrum

### 3.1 Approach

The key is to align the indices of the two arrays. `compute_dft` builds `freq = np.arange(1, N) / T_total`, so the $`i`$th entry of `freq` is the frequency $`(i+1)/T_{\text{total}}`$ for $`k = i + 1`$. The $`k`$th entry of `Y`, in contrast, is $`Y_k`$ itself. The values $`k = 1, \ldots, n_2`$ ($`n_2`$ = `N // 2`) needed for the one-sided spectrum sit at the following positions in the two arrays.

| Value needed | Position in `freq` | Position in `Y` |
|:-----|:-----|:-----|
| $`k = 1`$ ($`\Delta f`$) | `freq[0]` | `Y[1]` |
| $`k = n_2`$ (the Nyquist frequency for even $`N`$) | `freq[n2 - 1]` | `Y[n2]` |
| All | `freq[:n2]` | `Y[1:n2+1]` |

*Table 2. Positions of the entries needed for the one-sided spectrum*

Hence `fP` is `freq[:n2]`, and `Pyy` is the squared absolute value of `Y[1:n2+1]`. For complex numbers `np.abs` returns the magnitude $`|Y_k| = \sqrt{\mathrm{Re}(Y_k)^2 + \mathrm{Im}(Y_k)^2}`$. The DC component $`Y_0`$ and the mirrored half $`Y_{N/2+1}, \ldots, Y_{N-1}`$ are excluded. The definitions of the DFT and the power spectrum are covered in Sections 14 and 17 of the [Chapter 16 lecture note](../../Lecture-Notes/Ch16_Fourier-Analysis/Concepts_Lecture.md).

### 3.2 Code Written

```python
    fP  = freq[:n2]            #frequencies from delta_f to Nyquist
    Pyy = np.abs(Y[1:n2+1])**2 #power at each positive frequency
```

### 3.3 Results and Verification

This is the result for the signal of Example 16.3, $`f(t) = 5 + \cos(2\pi \cdot 12.5\, t) + \sin(2\pi \cdot 18.75\, t)`$, sampled with $`N = 8`$ and $`\Delta t = 0.02`$ s.

```text
=== One-Sided Power Spectrum (Ex. 16.3) ===
   Freq [Hz]    Power |Yk|^2
------------------------------
      6.2500      0.00000000
     12.5000      0.25000000  <-- peak
     18.7500      0.25000000  <-- peak
     25.0000      0.00000000

Both signal frequencies are clearly identified as peaks.
```

The output matches the output in the distributed notebook and is confirmed as follows.

- **Frequencies:** $`T_{\text{total}} = 8 \times 0.02 = 0.16`$ s, so $`\Delta f = 6.25`$ Hz, and the $`n_2 = 4`$ frequencies are 6.25, 12.5, 18.75, and 25 Hz. The last one, 25 Hz, equals $`f_s / 2 = 50 / 2`$, the Nyquist frequency.
- **Powers:** The DFT output gives $`Y_2 = 0.5`$ and $`Y_3 = -0.5i`$, so $`|Y_2|^2 = |Y_3|^2 = 0.25`$ and the rest are 0. The two peaks match the two signal frequencies exactly.
- **Amplitude and power:** A sinusoid of amplitude $`A`$ splits into the positive-frequency coefficient $`Y_k`$ and the mirror $`Y_{N-k}`$, each of magnitude $`A/2`$, so under this notebook's definition $`P_k = |Y_k|^2`$ its power is $`A^2/4`$. Both components have amplitude 1, giving 0.25.
- **Parseval's relation:** The sum of $`|Y_k|^2`$ over all $`k`$ is $`25 + 4 \times 0.25 = 26`$, which equals the mean square of the samples, $`\frac{1}{N}\sum y_j^2 = 26`$. The one-sided sum 0.5 is half of the power 1 that remains after removing the DC component, and the other half lies in the mirror.

### 3.4 Results Used by Later Cells

Exercise 4.1 applies the same function to the 10 points of Example 16.1 ($`T = 1.5`$ s, $`\Delta t = 0.15`$ s).

```text
=== Power Spectrum for Example 16.1 ===
   Freq [Hz]    Power |Yk|^2
------------------------------
      0.6667        0.250047  <-- peak
      1.3333        0.000000
      2.0000        0.000000
      2.6667        0.000000
      3.3333        0.000000

Dominant frequency from spectrum : 0.6667 Hz
Expected (1/T = w0/(2*pi))       : 0.6667 Hz
Match? YES
```

The peak lies at $`1/T = 0.6667`$ Hz, and the power 0.250047 equals $`C_1^2/4`$ computed from the combined amplitude $`C_1 = 1.0001`$ of the first-order Fourier series in Section 2 of the notebook. This holds because sampling one period at equal spacing gives $`Y_1 = (A_1 - iB_1)/2`$, so $`|Y_1|^2 = (A_1^2 + B_1^2)/4`$. Since $`N = 10`$, $`n_2 = 5`$, and the last frequency 3.3333 Hz is $`f_s / 2 = (1/0.15)/2`$, the Nyquist frequency.

The vibration analysis of Section 6 applies the function to an acceleration signal with $`f_s = 500`$ Hz and $`N = 512`$ and finds the frequency with the largest power in the first half of the spectrum (up to 125 Hz).

```text
Dominant frequency: 25.4 Hz  ->  1523 RPM
```

The rotation component of the signal is 25 Hz, yet the result is 25.4 Hz. The frequency resolution is $`\Delta f = 500/512 = 0.977`$ Hz, so 25 Hz falls between the neighboring frequencies 24.41 Hz and 25.39 Hz. Repeating the computation gives powers of 0.250 and 0.557 at these two frequencies: the power of the 25 Hz component is split between the neighboring frequencies, and the larger one, 25.39 Hz, is selected.

> **Note:** The noisy-signal cell of Section 5 uses 64 samples on $`t \in [0, 2\pi)`$, so the unit of `freq` is cycles per unit time. The two peaks of the signal $`\cos(10t) + \sin(3t)`$ therefore appear in `fP` at $`3/(2\pi) = 0.477`$ and $`10/(2\pi) = 1.592`$. The values 3 and 10 mentioned by the axis label (cycles per 2pi) and the printed sentence are not these values but the DFT indices $`k`$.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| One-sided power spectrum | Only the frequencies and powers $`P_k = \vert Y_k \vert^2`$ for $`k = 1, \ldots, N/2`$ are kept, and the DC component and the mirrored half are excluded. |
| Index alignment | `freq` starts at $`k = 1`$, so `fP = freq[:n2]`, while `Y` starts at $`k = 0`$, so `Pyy` comes from `Y[1:n2+1]`. |
| Amplitude and power | A sinusoid of amplitude $`A`$ appears as $`A/2`$ in both $`Y_k`$ and $`Y_{N-k}`$, so its power under this definition is $`A^2/4`$. |
| Frequency resolution | $`\Delta f = 1/(N \Delta t)`$, and when a signal frequency is not an integer multiple of this spacing its power is split between neighboring frequencies and the peak location shifts. |
| Nyquist frequency | For even $`N`$ the last frequency of the one-sided spectrum is $`f_s/2`$, and components above it cannot be distinguished by the sampling. |

---

<br>

## Self-Check Questions

1. **Number of entries:** In Example 16.3 ($`N = 8`$, $`\Delta t = 0.02`$ s), how many entries does `fP` have, and what are their frequencies in Hz?

   <details>
   <summary>Answer</summary>

   There are $`n_2 = 8 // 2 = 4`$ entries, the multiples of $`\Delta f = 1/0.16 = 6.25`$ Hz: 6.25, 12.5, 18.75, and 25 Hz. The 25 Hz entry is the Nyquist frequency.

   </details>

2. **Off-by-one index:** How would the output of Example 16.3 change if `fP = freq[1:n2+1]` were written by mistake?

   <details>
   <summary>Answer</summary>

   The frequencies would shift by one position to 12.5, 18.75, 25, and 31.25 Hz, so the two peaks actually at 12.5 Hz and 18.75 Hz would be printed as if they were at 18.75 Hz and 25 Hz. The power values stay the same, and only the frequency labels are off by $`\Delta f`$.

   </details>

3. **Power of 0.25:** The two components of Example 16.3 have amplitude 1, so why is the power 0.25 rather than 1?

   <details>
   <summary>Answer</summary>

   A real sinusoid splits equally between the positive and negative frequencies, so in the normalized DFT both $`Y_k`$ and the mirror $`Y_{N-k}`$ have magnitude $`A/2 = 0.5`$. The notebook defines the one-sided power as $`|Y_k|^2`$ without doubling, which gives $`0.5^2 = 0.25`$.

   </details>

4. **25.4 Hz:** In the vibration analysis the rotation frequency is 25 Hz, yet the result is 25.4 Hz. Explain why.

   <details>
   <summary>Answer</summary>

   With $`N = 512`$ and $`f_s = 500`$ Hz, the frequency resolution is $`500/512 = 0.977`$ Hz, and the frequencies of the spectrum are its integer multiples. 25 Hz falls between 24.41 Hz and 25.39 Hz, its power is split between these two frequencies, and the one with the larger power, 25.39 Hz, is printed.

   </details>
