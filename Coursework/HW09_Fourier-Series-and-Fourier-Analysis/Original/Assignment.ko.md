# HW9 과제 안내문 — 푸리에 급수와 푸리에 해석

[레포지토리 안내](../../../README.ko.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **최종 수정일:** 2026-11-24

수치해석 11주차 과제이다. 교수자가 배포한 실습 노트북 `Lab10-Fourier.ipynb`의 빈칸을 채워 실행한다. 노트북은 시간 영역과 주파수 영역, 푸리에 급수, 이산 푸리에 변환(DFT), 파워 스펙트럼, 신호 필터링과 잡음 제거, 공학 응용인 진동 분석의 여섯 부분으로 이루어진다. `COMPLETE` 표시가 있는 다른 함수(`fourier_series_1`, `compute_dft`, `lowpass_filter`)와 Exercise 2.1, 4.1, 5.1의 셀을 포함해 나머지 셀은 코드가 채워진 채로 배포되었으며, 빈칸은 4절의 함수 `power_spectrum` 한 곳이다. 노트북의 첫 코드 셀은 `numpy`를 `np`로, `matplotlib.pyplot`을 `plt`로 불러오고 `scipy.fft`에서 `fft`와 `ifft`를 불러온다.

아래 코드 블록은 노트북의 셀을 옮긴 것이며, 코드가 비어 있는 자리가 채울 곳이다. 배포된 빈 노트북은 [Lab10-Fourier.ipynb](Lab10-Fourier.ipynb)에 있다.

---

<br>

## 문제. 단측 파워 스펙트럼 (노트북 4.1절)

노트북 3절은 이산 푸리에 변환의 표본 추출 매개변수를 다음과 같이 정리한다.

| 매개변수 | 기호 | 공식 | 의미 |
|:-----|:-----|:-----|:-----|
| 표본 추출 간격 | $`\Delta t`$ | 주어짐 | 표본 사이의 시간 |
| 표본 추출 주파수 | $`f_s`$ | $`1/\Delta t`$ | 초당 표본 수 |
| 전체 기록 길이 | $`T_{\text{total}}`$ | $`N \cdot \Delta t`$ | 신호의 지속 시간 |
| 주파수 분해능 | $`\Delta f`$ | $`1/T_{\text{total}}`$ | 구별할 수 있는 가장 작은 주파수 |
| 나이퀴스트 주파수 | $`f_{\text{Nyq}}`$ | $`f_s/2`$ | 검출할 수 있는 최대 주파수 |

*표 1. 이산 푸리에 변환의 표본 추출 매개변수*

주파수 $`f`$를 검출하려면 $`f_s \geq 2f`$로 표본을 추출해야 하며, $`f_{\text{Nyq}}`$보다 높은 주파수는 더 낮은 주파수의 거짓 신호로 나타난다(앨리어싱). `scipy.fft.fft(y)`는 복소 계수 $`N`$개를 돌려주며, 정규화하면 다음과 같다.

$$
Y_k = \frac{1}{N}\sum_{j=0}^{N-1} y_j e^{-2\pi i k j / N}
$$

코드로는 `Y = fft(y) / N`이다. $`Y[0]`$은 평균값(직류 성분), $`Y[1 \ldots N/2]`$는 양의 주파수 성분이고, $`Y[N/2+1 \ldots N-1]`$은 음의 주파수의 거울상으로 실수 신호에서는 중복된다. 노트북은 정규화한 DFT를 다음 함수로 계산한다.

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

파워 스펙트럼은 신호의 에너지가 주파수에 따라 어떻게 분포하는지 보여 준다. 단측 스펙트럼($`\Delta f`$부터 나이퀴스트 주파수까지의 주파수)은 다음과 같이 정의한다.

$$
P_k = |Y_k|^2 \quad \text{for } k = 1, 2, \ldots, N/2
$$

- $`P_k`$가 크면 주파수 $`f_k`$에 강한 정현파 성분이 있다.
- $`P_k`$가 작으면 그 주파수에 의미 있는 성분이 없다.
- $`Y[N/2+1 \ldots N-1]`$은 양의 주파수 절반의 거울상이므로 $`k = 1 \ldots N/2`$만 쓴다.

| 특징 | 해석 |
|:-----|:-----|
| $`f_k`$의 날카로운 피크 | 신호에 $`f_k`$ Hz의 강한 주기 성분이 있다 |
| 넓고 평평한 바닥 | 무작위 잡음(백색 잡음은 모든 주파수에 고르게 퍼진다) |
| $`f = 0`$의 피크 | 직류 오프셋(0이 아닌 평균) |

*표 2. 파워 스펙트럼의 해석*

함수 `power_spectrum`의 빈칸을 채워, 단측 스펙트럼의 주파수 `fP`와 각 주파수의 파워 `Pyy`를 구하라.

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

노트북은 예제 16.3의 신호를 다음과 같이 표본 추출하고, `Y_dft, freq_dft, _ = compute_dft(y_dft, dt_dft)`로 DFT를 계산해 둔다.

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

배포 노트북에 실린 DFT 계수의 출력은 다음과 같다.

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

빈칸 다음 셀은 완성한 함수로 이 신호의 파워 스펙트럼을 구해 출력하고 그림으로 그린다. 셀의 앞부분은 다음과 같다.

```python
fP_dft, Pyy_dft = power_spectrum(Y_dft, freq_dft, fs_dft)

print('=== One-Sided Power Spectrum (Ex. 16.3) ===')
print(f'{"Freq [Hz]":>12}  {"Power |Yk|^2":>14}')
print('-' * 30)
for f_k, p_k in zip(fP_dft, Pyy_dft):
    tag = '  <-- peak' if p_k > 0.1 else ''
    print(f'{f_k:>12.4f}  {p_k:>14.8f}{tag}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

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

문제를 풀어 본 뒤 [과제 해설](../Assignment-Explanation.ko.md)에서 정답과 풀이를, [완성한 노트북](../Solution/Lab10-Fourier_Completed.ipynb)에서 작성한 코드 전체를 볼 수 있다.
