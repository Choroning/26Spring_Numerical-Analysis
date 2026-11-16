# HW9 Assignment Handout — Fourier Series and Fourier Analysis

[Repository Guide](../../../README.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **Last Updated:** 2026-11-16

This is the Week 11 homework of Numerical Analysis. The task is to fill in the blank of the lab notebook `Lab10-Fourier.ipynb` distributed by the instructor and run it. The notebook has six parts: time and frequency domains, Fourier series, the discrete Fourier transform (DFT), the power spectrum, signal filtering and noise removal, and an engineering application in vibration analysis. The other cells, including the other functions marked `COMPLETE` (`fourier_series_1`, `compute_dft`, and `lowpass_filter`) and the cells of Exercises 2.1, 4.1, and 5.1, were distributed with their code filled in, and the only blank is the function `power_spectrum` in Section 4. The first code cell of the notebook imports `numpy` as `np` and `matplotlib.pyplot` as `plt`, and imports `fft` and `ifft` from `scipy.fft`.

The code blocks below reproduce cells of the notebook, and the place where code is missing is the blank to fill. The distributed blank notebook is [Lab10-Fourier.ipynb](Lab10-Fourier.ipynb).

- **Deadline:** one week after distribution
- **Submission:** Submit the notebook online after filling in the blank and running it.

---

<br>

## Problem. One-Sided Power Spectrum (Notebook Section 4.1)

Section 3 of the notebook summarizes the sampling parameters of the discrete Fourier transform as follows.

| Parameter | Symbol | Formula | Meaning |
|:-----|:-----|:-----|:-----|
| Sampling interval | $`\Delta t`$ | given | Time between samples |
| Sampling frequency | $`f_s`$ | $`1/\Delta t`$ | Samples per second |
| Total record length | $`T_{\text{total}}`$ | $`N \cdot \Delta t`$ | Duration of the signal |
| Frequency resolution | $`\Delta f`$ | $`1/T_{\text{total}}`$ | Smallest distinguishable frequency |
| Nyquist frequency | $`f_{\text{Nyq}}`$ | $`f_s/2`$ | Maximum detectable frequency |

*Table 1. Sampling parameters of the discrete Fourier transform*

To detect a frequency $`f`$, one must sample at $`f_s \geq 2f`$, and frequencies above $`f_{\text{Nyq}}`$ appear as false lower-frequency signals (aliasing). `scipy.fft.fft(y)` returns $`N`$ complex coefficients, which after normalization are as follows.

$$
Y_k = \frac{1}{N}\sum_{j=0}^{N-1} y_j e^{-2\pi i k j / N}
$$

In code this is `Y = fft(y) / N`. $`Y[0]`$ is the mean value (DC component), $`Y[1 \ldots N/2]`$ are the positive-frequency components, and $`Y[N/2+1 \ldots N-1]`$ mirror the negative frequencies and are redundant for real signals. The notebook computes the normalized DFT with the following function.

```python
def compute_dft(y, dt):
    """
    Compute the normalized DFT of signal y.

    Parameters
    ----------
    y  : array  -- time-domain signal (N points)
    dt : float  -- sampling interval [s]

    Returns
    -------
    Y    : complex array  -- normalized DFT coefficients  Y = fft(y)/N
    freq : array          -- frequency array [Hz] for indices 1..N-1
    fs   : float          -- sampling frequency [Hz]
    """
    N       = len(y)
    fs      = 1.0 / dt
    T_total = N * dt

    # --- COMPLETE ---
    Y    =       fft(y) / N
    freq =       np.arange(1, N) / T_total
    # --- COMPLETE ---

    return Y, freq, fs
```

The power spectrum shows how the signal energy is distributed across frequencies. The one-sided spectrum (frequencies from $`\Delta f`$ to the Nyquist frequency) is defined as follows.

$$
P_k = |Y_k|^2 \quad \text{for } k = 1, 2, \ldots, N/2
$$

- A large $`P_k`$ means a strong sinusoidal component at frequency $`f_k`$.
- A small $`P_k`$ means no significant component there.
- Only $`k = 1 \ldots N/2`$ is used because $`Y[N/2+1 \ldots N-1]`$ mirrors the positive-frequency half.

| Feature | Interpretation |
|:-----|:-----|
| Sharp peak at $`f_k`$ | The signal has a strong periodic component at $`f_k`$ Hz |
| Broad, flat floor | Random noise (white noise spreads equally over all frequencies) |
| Peak at $`f = 0`$ | DC offset (nonzero mean) |

*Table 2. Interpreting the power spectrum*

Fill in the blanks of the function `power_spectrum` to compute the frequencies `fP` of the one-sided spectrum and the power `Pyy` at each frequency.

```python
def power_spectrum(Y, freq, fs):
    """
    Compute the one-sided power spectrum from normalized DFT coefficients Y.

    Parameters
    ----------
    Y    : complex array -- normalized DFT (Y = fft(y)/N), length N
    freq : array         -- frequency array for indices 1..N-1
    fs   : float         -- sampling frequency [Hz]

    Returns
    -------
    fP  : array -- frequency values for one-sided spectrum
    Pyy : array -- power at each frequency (|Y_k|^2)
    """
    N  = len(Y)
    n2 = N // 2

    # --- COMPLETE ---
    fP  =        #frequencies from delta_f to Nyquist
    Pyy =        #power at each positive frequency
    # --- COMPLETE ---

    return fP, Pyy
```

The notebook samples the signal of Example 16.3 as follows and computes its DFT with `Y_dft, freq_dft, _ = compute_dft(y_dft, dt_dft)`.

```python
# Example 16.3 — Known signal: f(t) = 5 + cos(2pi*12.5*t) + sin(2pi*18.75*t)
n_dft     = 8
dt_dft    = 0.02          # s
fs_dft    = 1.0 / dt_dft  # Hz
T_total   = n_dft * dt_dft
t_dft     = np.arange(n_dft) / fs_dft

freq1_ex  = 12.5          # Hz
freq2_ex  = 18.75         # Hz
func_ex   = lambda t: 5 + np.cos(2*np.pi*freq1_ex*t) + np.sin(2*np.pi*freq2_ex*t)
y_dft     = func_ex(t_dft)
```

The output of the DFT coefficients in the distributed notebook is as follows.

```text
=== DFT Coefficients Y[k] ===
k    freq [Hz]       Re(Y[k])       Im(Y[k])       |Y[k]|
----------------------------------------------------------
0       0.0000       5.000000      -0.000000     5.000000
1       6.2500       0.000000       0.000000     0.000000
2      12.5000       0.500000      -0.000000     0.500000
3      18.7500      -0.000000      -0.500000     0.500000
4      25.0000       0.000000      -0.000000     0.000000
5      31.2500      -0.000000       0.500000     0.500000
6      37.5000       0.500000       0.000000     0.500000
7      43.7500       0.000000      -0.000000     0.000000

Y[0].real = 5.0000  <- mean (expected 5)
Y[2].real = 0.5000  <- cosine at 12.50 Hz (expected 0.5)
Y[3].imag = -0.5000  <- sine   at 18.75 Hz (expected -0.5)
```

The cell after the blank uses the completed function to compute, print, and plot the power spectrum of this signal. The beginning of the cell is as follows.

```python
fP_dft, Pyy_dft = power_spectrum(Y_dft, freq_dft, fs_dft)

print('=== One-Sided Power Spectrum (Ex. 16.3) ===')
print(f'{"Freq [Hz]":>12}  {"Power |Yk|^2":>14}')
print('-' * 30)
for f_k, p_k in zip(fP_dft, Pyy_dft):
    tag = '  <-- peak' if p_k > 0.1 else ''
    print(f'{f_k:>12.4f}  {p_k:>14.8f}{tag}')
```

The distributed notebook shows the following output for this cell.

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

---

<br>

After solving the problem, see the [assignment explanation](../Assignment-Explanation.md) for the answer and solution and the [completed notebook](../Solution/Lab10-Fourier_Completed.ipynb) for the full code.
