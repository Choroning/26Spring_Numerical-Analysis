# 13장 실습 — 대칭 고유값 문제와 응용

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lab.md) | [한국어](Concepts_Lab.ko.md)

> **최종 수정일:** 2026-10-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 13

> **학습 목표:**
> 1. `np.linalg.eigh`를 사용해 대칭 행렬의 특별한 성질을 활용할 수 있다
> 2. 스펙트럼 분해 $`A = V\Lambda V^T`$를 사용해 저계수(low-rank) 근사를 만들 수 있다
> 3. `scipy.linalg.eigh(K, M)`로 일반화 고유값 문제 $`K\mathbf{v} = \lambda M\mathbf{v}`$를 풀 수 있다
> 4. 질량-스프링 시스템의 고유진동수와 모드 형상(mode shape)을 계산할 수 있다
> 5. 선형 ODE 시스템 $`\dot{\mathbf{y}} = A\mathbf{y}`$의 안정성을 고유값으로 분석할 수 있다

---

<br>

## 목차

- [1. 대칭 행렬과 스펙트럼 분해](#1-대칭-행렬과-스펙트럼-분해)
  - [1.1 특별한 성질](#11-특별한-성질)
  - [1.2 `np.linalg.eigh`와 재구성](#12-nplinalgeigh와-재구성)
  - [1.3 저계수 근사](#13-저계수-근사)
- [2. 응용: 진동 해석](#2-응용-진동-해석)
  - [2.1 운동 방정식과 일반화 고유값 문제](#21-운동-방정식과-일반화-고유값-문제)
  - [2.2 `scipy.linalg.eigh(K, M)`로 고유진동수 계산](#22-scipylinalgeighk-m로-고유진동수-계산)
  - [2.3 모드 형상 시각화](#23-모드-형상-시각화)
  - [2.4 자유 진동 시간 응답](#24-자유-진동-시간-응답)
- [3. 응용: ODE 시스템과 안정성](#3-응용-ode-시스템과-안정성)
  - [3.1 고유쌍을 이용한 일반해](#31-고유쌍을-이용한-일반해)
  - [3.2 안정성 분류](#32-안정성-분류)
  - [3.3 예제: 문제 13.13](#33-예제-문제-1313)
- [4. 방법 비교와 공학적 맥락](#4-방법-비교와-공학적-맥락)
  - [4.1 방법 개요](#41-방법-개요)
  - [4.2 주요 공식](#42-주요-공식)
  - [4.3 공학적 맥락](#43-공학적-맥락)
  - [4.4 실전 팁](#44-실전-팁)
- [요약](#요약)
- [점검 문제](#점검-문제)

---

<br>

## 1. 대칭 행렬과 스펙트럼 분해

### 1.1 특별한 성질

$`A = A^T`$(대칭)이면 고유구조에 일반적인 경우에는 없는 여러 보장이 생긴다.

| 성질 | 설명 |
|:---------|:----------|
| 실수 고유값 | 모든 $`\lambda_i \in \mathbb{R}`$ |
| 직교 고유벡터 | $`i \neq j`$일 때 $`\mathbf{v}_i^T \mathbf{v}_j = 0`$ (중복 고유값이면 직교하도록 고를 수 있음) |
| 정규직교 고유기저 | $`V^T V = V V^T = I`$ |
| 스펙트럼 분해 | $`A = V \Lambda V^T = \displaystyle\sum_{i=1}^n \lambda_i\, \mathbf{v}_i\mathbf{v}_i^T`$ |

*표 1. 실수 대칭 행렬의 고유구조가 갖는 성질*

이러한 보장 덕분에 대칭 고유값 전용 알고리즘(야코비법, 삼중대각 형태에 대한 QR)은 일반 알고리즘보다 훨씬 효율적이다.

행렬이 대칭(복소수이면 에르미트)임을 알고 있을 때는 `eig` 대신 **`np.linalg.eigh`** 를 사용한다. 이 함수의 특징은 다음과 같다.

- 실수 고유값만 반환한다(반올림 오차에 의한 가짜 허수부가 없다).
- 고유값을 오름차순으로 정렬한다.
- `eig`보다 약 2배 빠르다.
- 정규직교 고유벡터 행렬을 반환한다.

---

<br>

### 1.2 `np.linalg.eigh`와 재구성

스펙트럼 분해 $`A = V \Lambda V^T`$를 이용하면 고유쌍만으로 대칭 행렬을 정확하게 재구성할 수 있다.

```python
import numpy as np

A = np.array([[20.0,  3.0,  2.0],
              [ 3.0,  9.0,  4.0],
              [ 2.0,  4.0, 12.0]])

# eigh: 정렬됨, 실수, 정규직교
eigvals_h, V = np.linalg.eigh(A)   # 오름차순
print('eigh eigenvalues (ascending):', np.round(eigvals_h, 8))

print('\nOrthonormality check  V^T V:')
print(np.round(V.T @ V, 10))

# 스펙트럼 분해  A = V Lambda V^T
Lambda = np.diag(eigvals_h)
A_reconstructed = V @ Lambda @ V.T
print(f'\nmax|A - V Lambda V^T| = {np.max(np.abs(A - A_reconstructed)):.2e}')
```

재구성 오차는 기계 정밀도 수준($`10^{-14}`$ 이하)이어야 한다. 이 항등식은 많은 응용의 기초가 된다. $`V`$와 $`\Lambda`$를 얻으면 $`A^k`$, $`e^A`$, $`A^{-1}`$ 같은 표현이 모두 $`V f(\Lambda) V^T`$, 즉 값싼 직교 곱 사이에 낀 대각 연산으로 줄어든다.

---

<br>

### 1.3 저계수 근사

고유값을 $`|\lambda_1| \ge |\lambda_2| \ge \cdots \ge |\lambda_n|`$ 순서로 놓으면, 스펙트럼 합 $`A = \sum_{i=1}^n \lambda_i \mathbf{v}_i \mathbf{v}_i^T`$를 절댓값이 가장 큰 $`k`$개 항으로 절단할 수 있다.

```math
A_k = \sum_{i=1}^{k} \lambda_i\, \mathbf{v}_i\mathbf{v}_i^T, \qquad
\|A - A_k\|_F = \sqrt{\sum_{i=k+1}^{n} \lambda_i^2}
```

이것은 Frobenius 노름에서의 **최적 랭크-$`k`$ 근사**이다(대칭 $`A`$에 대한 Eckart-Young 정리). 이 원리는 주성분 분석(PCA), 이미지 압축, 노이즈 제거의 기반이 된다.

```python
print('Low-rank approximation error ||A - A_k||_F:')
for k in [1, 2, 3]:
    idx = np.argsort(np.abs(eigvals_h))[::-1][:k]
    A_k = sum(eigvals_h[i] * np.outer(V[:, i], V[:, i]) for i in idx)
    err_F = np.linalg.norm(A - A_k, 'fro')
    err_theory = np.sqrt(np.sum(eigvals_h[np.argsort(np.abs(eigvals_h))[:-k]] ** 2))
    print(f'  k={k}: ||A - A_k||_F = {err_F:.6f}  (theory: {err_theory:.6f})')
```

$`k`$가 커질수록 크기가 작은 고유값만 버리게 되어 근사 오차가 줄어든다. 이 행렬에서는 $`k = 1`$일 때 14.520460, $`k = 2`$일 때 6.088518이며, $`k = n`$이면 재구성이 정확해진다.

---

<br>

## 2. 응용: 진동 해석

### 2.1 운동 방정식과 일반화 고유값 문제

세 개의 질량이 네 개의 스프링으로 연결되고 양쪽 끝이 벽에 고정된 시스템을 생각한다.

```text
Wall ─[k₁]─ m₁ ─[k₂]─ m₂ ─[k₃]─ m₃ ─[k₄]─ Wall
```

*그림 1. 양쪽 벽에 고정된 세 질량, 네 스프링 시스템*

뉴턴의 제2법칙에서 감쇠가 없는 자유 진동 방정식을 얻는다.

```math
M\ddot{\mathbf{x}} + K\mathbf{x} = \mathbf{0}
```

조화 운동 $`\mathbf{x}(t) = \mathbf{v}\sin(\omega t)`$를 가정하면 **일반화 고유값 문제**가 유도된다.

```math
K\mathbf{v} = \omega^2 M\mathbf{v}
```

여기서 $`\lambda = \omega^2`$는 **고유진동수**의 제곱이고 $`\mathbf{v}`$는 대응하는 **모드 형상(mode shape)** 이다. 강성 행렬과 질량 행렬은 다음과 같다.

```math
K = \begin{bmatrix} k_1+k_2 & -k_2 & 0 \\ -k_2 & k_2+k_3 & -k_3 \\ 0 & -k_3 & k_3+k_4 \end{bmatrix},
\qquad
M = \begin{bmatrix} m_1 & & \\ & m_2 & \\ & & m_3 \end{bmatrix}
```

$`K`$와 $`M`$은 모두 대칭 양의 정부호(SPD) 행렬이므로 `scipy.linalg.eigh(K, M)`이 알맞은 도구이다. 이 함수는 실수 고유값을 오름차순으로 반환하고, $`M`$-정규화된 모드 형상($`\phi_i^T M \phi_i = 1`$)을 제공한다.

---

<br>

### 2.2 `scipy.linalg.eigh(K, M)`로 고유진동수 계산

$`\omega^2`$로부터 각진동수 $`\omega = \sqrt{\lambda}`$ (rad/s), 선형 진동수 $`f = \omega / (2\pi)`$ (Hz), 주기 $`T = 1/f`$ (s)를 얻는다.

```python
import scipy as sc

# 시스템 매개변수
k1, k2, k3, k4 = 40.0, 30.0, 20.0, 10.0   # N/m
m1, m2, m3     =  2.0,  1.5,  1.0          # kg

K_vib = np.array([
    [k1+k2,   -k2,      0],
    [  -k2, k2+k3,    -k3],
    [    0,   -k3,  k3+k4]
])
M_vib = np.diag([m1, m2, m3])

# K v = omega^2 M v 풀이
omega_sq, modes = sc.linalg.eigh(K_vib, M_vib)   # 고유값 오름차순

omega_n = np.sqrt(omega_sq)
freq_hz = omega_n / (2 * np.pi)

print('=== Natural Frequencies ===')
print(f'{"Mode":>5}  {"omega^2":>14}  {"omega (rad/s)":>14}  {"f (Hz)":>10}  {"T (s)":>10}')
print('-' * 60)
for i in range(3):
    print(f'{i+1:>5}  {omega_sq[i]:>14.4f}  {omega_n[i]:>14.4f}  {freq_hz[i]:>10.4f}  {1/freq_hz[i]:>10.4f}')

# 일반화 고유값 방정식 검증
print('\nVerification  max|K v - omega^2 M v|:')
for i in range(3):
    v   = modes[:, i]
    res = K_vib @ v - omega_sq[i] * M_vib @ v
    print(f'  Mode {i+1}: {np.max(np.abs(res)):.2e}')
```

실행하면 고유진동수 $`\omega = 3.0078`$, $`5.6889`$, $`7.5447`$ rad/s(진동수 0.4787, 0.9054, 1.2008 Hz)를 얻는다.

> **참고:** 일반화 고유값 문제는 $`\tilde{K} = M^{-1/2} K M^{-1/2}`$와 같은 변환으로 표준 문제로 환원할 수 있다. `scipy.linalg.eigh(K, M)`은 같은 종류의 환원을 $`M`$의 촐레스키 분해로 내부에서 처리하고, $`M`$-정규직교인 모드 형상($`\Phi^T M \Phi = I`$)을 반환한다.

---

<br>

### 2.3 모드 형상 시각화

`modes`의 각 열은 해당 진동 모드에서 세 질량이 어떻게 변위되는지를 나타내는 벡터이다. 노드 위치에 대해 변위를 그리면 모드 구조가 드러난다.

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(13, 5))
x_pos    = [0, 1, 2, 3, 4]   # wall, m1, m2, m3, wall
x_masses = [1, 2, 3]

for i, ax in enumerate(axes):
    v      = modes[:, i]
    v_norm = v / np.max(np.abs(v))   # 최대 진폭 = 1로 스케일
    disp   = np.array([0.0, v_norm[0], v_norm[1], v_norm[2], 0.0])

    ax.axhline(0, color='gray', ls='--', lw=1)
    ax.plot(x_pos, disp, 'b-o', ms=13, lw=2,
            markerfacecolor='steelblue', markeredgecolor='navy')
    ax.plot([0, 4], [0, 0], 'ks', ms=14)   # 벽

    ax.set(xlim=[-0.5, 4.5], ylim=[-1.7, 1.7],
           title=f'Mode {i+1}\nω = {omega_n[i]:.2f} rad/s  ({freq_hz[i]:.3f} Hz)',
           xlabel='Node', ylabel='Normalized displacement')
    ax.set_xticks(x_pos, ['Wall', 'm₁', 'm₂', 'm₃', 'Wall'])
    ax.grid(True, alpha=0.3)

plt.suptitle('Three-Mass Spring System — Mode Shapes', fontsize=13)
plt.tight_layout()
plt.show()
```

가장 낮은 진동수의 모드에서는 일반적으로 모든 질량이 같은 방향으로 움직이고, 높은 모드일수록 체인을 따라 부호 변화(node)가 많아진다. 부호 변화의 개수는 모드 번호와 함께 증가하며, 이는 진동하는 현의 고조파와 비슷한 패턴이다.

---

<br>

### 2.4 자유 진동 시간 응답

초기 조건 $`\mathbf{x}(0) = \mathbf{x}_0`$와 $`\dot{\mathbf{x}}(0) = \mathbf{v}_0`$가 주어지면 시간 응답은 각 모드의 중첩이다.

```math
\mathbf{x}(t) = \sum_{i=1}^{n}\! \left[c_i \cos(\omega_i t) + s_i \sin(\omega_i t)\right] \boldsymbol\phi_i
```

모드 진폭은 초기 조건을 모드 형상에 $`M`$-사영하여 $`c_i = \boldsymbol\phi_i^T M \mathbf{x}_0`$, $`s_i = \boldsymbol\phi_i^T M \mathbf{v}_0 / \omega_i`$로 얻는다.

```python
# 초기 조건: m1을 50 mm 변위시킨 뒤 정지 상태에서 놓기
x0_vib  = np.array([0.05, 0.0, 0.0])
dx0_vib = np.zeros(3)

c_amp = modes.T @ M_vib @ x0_vib
s_amp = modes.T @ M_vib @ dx0_vib / (omega_n + 1e-30)

t = np.linspace(0, 3.0, 1000)
x_t = np.zeros((3, len(t)))
for i in range(3):
    x_t += np.outer(modes[:, i],
                    c_amp[i] * np.cos(omega_n[i] * t)
                    + s_amp[i] * np.sin(omega_n[i] * t))

fig, ax = plt.subplots(figsize=(10, 4))
for j in range(3):
    ax.plot(t, x_t[j] * 1000, lw=1.8, label=f'm_{j+1}')
ax.set(xlabel='Time (s)', ylabel='Displacement (mm)',
       title='Free Vibration: Initial Displacement of m₁ = 50 mm')
ax.legend(); ax.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()
```

응답은 세 고유진동수에서의 순수한 사인파의 합이다. 감쇠가 없으면 운동이 무한히 계속되지만, 실제 시스템에서는 감쇠로 에너지가 서서히 빠져나가며 대개 높은 모드가 더 빨리 감쇠한다.

---

<br>

## 3. 응용: ODE 시스템과 안정성

### 3.1 고유쌍을 이용한 일반해

선형 1계 시스템 $`\dot{\mathbf{y}} = A\mathbf{y}`$, $`\mathbf{y}(0) = \mathbf{y}_0`$에서 $`A`$가 일차독립인 고유벡터 $`n`$개를 가지면(대각화 가능하면) 일반해는 고유 모드들의 합이다.

```math
\mathbf{y}(t) = \sum_{i=1}^{n} c_i\, \mathbf{v}_i\, e^{\lambda_i t}, \qquad
\mathbf{c} = V^{-1}\mathbf{y}_0
```

각 모드는 고유벡터 방향 $`\mathbf{v}_i`$를 따라 독립적으로 변하며 $`e^{\lambda_i t}`$로 지수적으로 스케일된다. 계수 벡터 $`\mathbf{c}`$는 $`V\mathbf{c} = \mathbf{y}_0`$을 풀어 얻는다.

---

<br>

### 3.2 안정성 분류

모든 모드의 장기 거동은 고유값의 실수부로 결정된다.

| 조건 | 거동 |
|:----------|:---------|
| 모든 $`\text{Re}(\lambda_i) < 0`$ | 해가 **감쇠** → **점근 안정(asymptotically stable)** |
| 하나라도 $`\text{Re}(\lambda_i) > 0`$ | 해가 **증가** → **불안정(unstable)** |
| 모든 $`\text{Re}(\lambda_i) \leq 0`$이고 일부 $`= 0`$ | **경계 안정(marginally stable)** (유계이나 감쇠하지 않음) |

*표 2. 고유값의 실수부에 따른 안정성 분류*

세 번째 행은 3.1절처럼 $`A`$가 대각화 가능하다고 가정한 것이다. 실수부가 0인 고유값에 대응하는 일차독립 고유벡터가 부족하면 해가 $`t`$에 비례해 커질 수 있다.

$`\lambda_i`$의 허수부는 각진동수 $`\text{Im}(\lambda_i)`$의 진동을 만들고, 이 진동은 지수 포락선 $`e^{\text{Re}(\lambda_i) t}`$로 변조된다.

**강성비(stiffness ratio)** $`|\lambda_{\max} / \lambda_{\min}|`$(실수부 기준)는 감쇠 시간 척도가 얼마나 차이 나는지를 나타낸다. 강성비가 크면 시스템이 수치적으로 **강성(stiff)** 을 띤다. 느린 모드가 답을 지배하는데도 명시적 ODE 적분기는 가장 빠른 모드를 따라가기 위해 매우 작은 스텝을 써야 한다.

---

<br>

### 3.3 예제: 문제 13.13

다음 시스템을 생각한다.

```math
\begin{cases} \dot{y}_1 = -5 y_1 + 3 y_2 \\ \dot{y}_2 = 100 y_1 - 301 y_2 \end{cases},
\qquad \mathbf{y}(0) = [50,\; 100]^T
```

```python
A_ode = np.array([[ -5.0,    3.0],
                  [100.0, -301.0]])
y0    = np.array([50.0, 100.0])

lam_ode, V_ode = np.linalg.eig(A_ode)
print('Eigenvalues:')
for i, lam in enumerate(lam_ode):
    status = 'decays' if np.real(lam) < 0 else 'GROWS'
    print(f'  lambda_{i+1} = {lam:.6f}   -> mode {status}')
stable = np.all(np.real(lam_ode) < 0)
print(f'System is {"STABLE" if stable else "UNSTABLE"}')

# 모드 상수  c = V^{-1} y0
c_ode = np.linalg.solve(V_ode, y0)
print(f'\nModal constants: c1 = {c_ode[0]:.4f},  c2 = {c_ode[1]:.4f}')

# 시간 해
t = np.linspace(0, 1.5, 800)
y_t = np.zeros((2, len(t)))
for i in range(2):
    y_t += np.real(c_ode[i] * np.outer(V_ode[:, i], np.exp(lam_ode[i] * t)))

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(t, y_t[0], 'b-', lw=2, label='y₁(t)')
ax.plot(t, y_t[1], 'r-', lw=2, label='y₂(t)')
ax.set(xlabel='Time t', ylabel='y', title='ODE Solution  ẏ = A y')
ax.legend(); ax.grid(True, alpha=0.4)

ax = axes[1]
ax.plot(y_t[0], y_t[1], 'k-', lw=2)
ax.plot(y0[0], y0[1], 'go', ms=12, label='Start y(0)', zorder=5)
ax.plot(y_t[0, -1], y_t[1, -1], 'rs', ms=12, label='End', zorder=5)
ax.set(xlabel='y₁', ylabel='y₂', title='Phase Portrait')
ax.legend(); ax.grid(True, alpha=0.4)

plt.tight_layout()
plt.show()

print(f'\nTime constants:  tau_1 = 1/|lambda_1| = {1/abs(lam_ode[0]):.4f} s  (slow mode)')
print(f'                 tau_2 = 1/|lambda_2| = {1/abs(lam_ode[1]):.6f} s  (fast mode)')
print(f'Stiffness ratio  |lambda_2 / lambda_1| = {abs(lam_ode[1]/lam_ode[0]):.1f}')
```

이 행렬의 두 고유값은 모두 음의 실수이며 크게 떨어져 있다($`\lambda_1 \approx -4`$, $`\lambda_2 \approx -302`$). 따라서 시스템은 안정하지만 강성이다. 출력되는 시간 상수는 $`\tau_1 \approx 0.2506`$ s, $`\tau_2 \approx 0.0033`$ s이고 강성비는 약 75.7이다. 빠른 모드는 100분의 몇 초 안에 감쇠하고, 느린 모드가 눈에 보이는 과도 응답을 지배한다. 이 예는 23장에서 다루는 강성 ODE 방법의 전형적인 도입 사례이다.

---

<br>

## 4. 방법 비교와 공학적 맥락

### 4.1 방법 개요

| 방법 | 찾는 값 | 비용 | 비고 |
|:-------|:------|:-----|:------|
| `np.linalg.eig(A)` | 모든 $`n`$개의 고유쌍 | $`O(n^3)`$ | 일반 행렬; 복소값 가능 |
| `np.linalg.eigh(A)` | 모든 $`n`$개의 고유쌍 | $`O(n^3)`$ | 대칭 전용; 실수, 오름차순 정렬 |
| **거듭제곱법** (12장 실습) | 지배 $`\lambda`$ | 반복당 $`O(n^2)`$ | 속도 $`\approx \|\lambda_2/\lambda_1\|`$ |
| **역거듭제곱법과 편이 거듭제곱법** (12장 실습) | 최소 또는 $`\sigma`$에 가까운 $`\lambda`$ | LU 이후 반복당 $`O(n^2)`$ | 내부 고유값에 사용 |
| `sc.linalg.eigh(K, M)` | 일반화 $`K\mathbf{v} = \lambda M\mathbf{v}`$ | $`O(n^3)`$ | SPD 쌍; 고유진동수 |

*표 3. 고유값 계산 방법 개요*

### 4.2 주요 공식

| 개념 | 공식 |
|:--------|:--------|
| 스펙트럼 분해 | $`A = V \Lambda V^T`$ (대칭 $`A`$, $`V^T V = I`$) |
| 저계수 근사 | $`A_k = \sum_{i=1}^k \lambda_i\, \mathbf{v}_i \mathbf{v}_i^T`$ |
| Frobenius 오차 | $`\Vert A - A_k\Vert_F = \sqrt{\sum_{i=k+1}^n \lambda_i^2}`$ (고유값을 절댓값 내림차순으로 정렬) |
| 고유진동수 | $`\omega_i = \sqrt{\lambda_i}`$, $`K\mathbf{v} = \lambda M\mathbf{v}`$로부터 |
| ODE 해 | $`\mathbf{y}(t) = \sum_i c_i\, \mathbf{v}_i\, e^{\lambda_i t}, \quad \mathbf{c} = V^{-1}\mathbf{y}_0`$ |
| 안정성 기준 | 점근 안정 $`\Leftrightarrow`$ 모든 $`\text{Re}(\lambda_i) < 0`$ |

*표 4. 대칭 고유값 문제와 응용의 주요 공식*

### 4.3 공학적 맥락

| 응용 | 고유값의 의미 | 고유벡터의 의미 |
|:------------|:-------------------|:--------------------|
| 구조 진동 | $`\lambda = \omega^2`$ (고유진동수 제곱) | 모드 형상 |
| ODE 시스템 | 감쇠 / 증가율 $`e^{\lambda t}`$ | 상태 공간의 모드 방향 |
| 반복 해법 (12장) | 반복행렬 $`T`$(12장의 $`M`$)의 스펙트럼 반경 $`\rho(T)`$ | 수렴 방향 |
| 조건수 (11장) | $`\kappa_2 = \vert\lambda\vert_{\max} / \vert\lambda\vert_{\min}`$ (대칭 행렬, 2-노름) | 해당 없음 |

*표 5. 공학 응용에서 고유값과 고유벡터의 의미*

### 4.4 실전 팁

1. 대칭 행렬에는 항상 `eig`보다 `eigh`를 사용한다. 더 빠르고, 실수 값을 주며, 결과가 정렬되어 있다.
2. 대형 행렬에서 고유값 몇 개만 필요하면 전체 `eig` 대신 거듭제곱법이나 편이 거듭제곱법(12장 실습)을 사용한다.
3. 일반화 고유값 문제는 $`K\mathbf{v} - \lambda M\mathbf{v} \approx \mathbf{0}`$인지 확인하여 검증한다.
4. 강성비 $`|\lambda_{\max} / \lambda_{\min}|`$가 크면 강성 ODE라는 신호이며, 암시적 적분기가 필요하다.

---

<br>

## 요약

| 개념 | 핵심 요약 |
|:-----|:-----|
| 대칭 행렬의 성질 | 실수 대칭 행렬은 실수 고유값과 정규직교 고유벡터를 가지며 $`A = V \Lambda V^T`$로 분해된다. |
| `np.linalg.eigh` | 대칭 행렬 전용 함수로, 실수 고유값을 오름차순으로 정렬하고 정규직교 고유벡터를 반환한다. |
| 저계수 근사 | 절댓값이 큰 $`k`$개 항만 남긴 $`A_k = \sum_{i=1}^k \lambda_i \mathbf{v}_i \mathbf{v}_i^T`$가 Frobenius 노름에서 최적 랭크-$`k`$ 근사이며, 오차는 버린 고유값 제곱합의 제곱근이다. |
| 일반화 고유값 문제 | 진동계의 $`K\mathbf{v} = \omega^2 M\mathbf{v}`$는 `scipy.linalg.eigh(K, M)`으로 풀며, 고유진동수는 $`\omega_i = \sqrt{\lambda_i}`$이다. |
| 모드 중첩 | 자유 진동 응답은 각 모드의 사인파 합이며, 진폭은 초기 조건의 $`M`$-사영으로 정해진다. |
| ODE 해 | 대각화 가능한 $`A`$에서 $`\dot{\mathbf{y}} = A\mathbf{y}`$의 해는 $`\mathbf{y}(t) = \sum_i c_i \mathbf{v}_i e^{\lambda_i t}`$, $`\mathbf{c} = V^{-1}\mathbf{y}_0`$이다. |
| 안정성 | 모든 고유값의 실수부가 음수이면 점근 안정하고, 하나라도 양수이면 불안정하다. |
| 강성비 | 실수부 기준 $`\vert\lambda_{\max}/\lambda_{\min}\vert`$가 크면 강성 시스템이며, 명시적 적분기는 매우 작은 스텝이 필요하다. |

---

<br>

## 점검 문제

1. **저계수 근사 오차:** 1.2절의 행렬은 고유값이 21.729156, 13.182326, 6.088518이다. 랭크-1 근사와 랭크-2 근사의 Frobenius 오차를 구하라.

   <details>
   <summary>정답</summary>

   $`\Vert A - A_1\Vert_F = \sqrt{13.182326^2 + 6.088518^2} = \sqrt{173.774 + 37.070} = \sqrt{210.844} \approx 14.520`$이고, $`\Vert A - A_2\Vert_F = \sqrt{6.088518^2} = 6.088518`$이다. 1.3절의 코드 출력 14.520460, 6.088518과 일치한다.

   </details>

2. **일반화 고유값 문제:** 질량 $`m_1 = m_2 = 2`$ kg인 두 질량이 강성 $`k = 10`$ N/m인 스프링 세 개로 두 벽 사이에 연결되어 있다. $`K`$와 $`M`$을 쓰고 고유진동수와 모드 형상을 구하라.

   <details>
   <summary>정답</summary>

   $`K = \begin{bmatrix} 20 & -10 \\ -10 & 20 \end{bmatrix}`$, $`M = \begin{bmatrix} 2 & 0 \\ 0 & 2 \end{bmatrix}`$이다. $`\det(K - \omega^2 M) = (20 - 2\omega^2)^2 - 100 = 0`$에서 $`20 - 2\omega^2 = \pm 10`$이므로 $`\omega^2 = 5`$ 또는 $`15`$이다. 따라서 $`\omega_1 = \sqrt{5} \approx 2.236`$ rad/s($`f_1 \approx 0.356`$ Hz), $`\omega_2 = \sqrt{15} \approx 3.873`$ rad/s($`f_2 \approx 0.616`$ Hz)이다. 모드 형상은 $`\omega_1`$에서 $`[1, 1]^T`$(동위상), $`\omega_2`$에서 $`[1, -1]^T`$(역위상)이며, 13장 예제 13.4와 같은 결과이다.

   </details>

3. **안정성 분류:** 다음 행렬로 정의한 $`\dot{\mathbf{y}} = A\mathbf{y}`$의 안정성을 판정하라. (a) $`\begin{bmatrix} -2 & 1 \\ 0 & -5 \end{bmatrix}`$ (b) $`\begin{bmatrix} 0 & 1 \\ -4 & 0 \end{bmatrix}`$ (c) $`\begin{bmatrix} 1 & 2 \\ 2 & 1 \end{bmatrix}`$

   <details>
   <summary>정답</summary>

   (a)는 삼각행렬이므로 고유값이 대각 원소 $`-2`$, $`-5`$이고, 모두 음수이므로 점근 안정하다(강성비 2.5). (b)는 $`\lambda^2 + 4 = 0`$에서 $`\lambda = \pm 2i`$이므로 경계 안정하며, 각진동수 2(주기 $`\pi`$)로 감쇠 없이 진동한다. (c)는 $`(1 - \lambda)^2 - 4 = 0`$에서 고유값이 3과 $`-1`$이므로, 양수 고유값 3 때문에 불안정하다.

   </details>

4. **강성비와 시간 상수:** 3.3절의 시스템에서 시간 상수 $`\tau_i = 1/|\lambda_i|`$와 강성비를 구하고, 명시적 적분기에 어떤 어려움이 생기는지 설명하라.

   <details>
   <summary>정답</summary>

   대각합은 $`-306`$, 행렬식은 $`1505 - 300 = 1205`$이므로 $`\lambda = \frac{-306 \pm \sqrt{306^2 - 4 \cdot 1205}}{2}`$에서 $`\lambda_1 \approx -3.990`$, $`\lambda_2 \approx -302.010`$이다. 시간 상수는 $`\tau_1 \approx 0.2506`$ s, $`\tau_2 \approx 0.00331`$ s이고, 강성비는 $`302.010/3.990 \approx 75.7`$이다. 빠른 모드는 $`5\tau_2 \approx 0.017`$ s 안에 사라지지만, 명시적 적분기는 안정성을 위해 $`\tau_2`$ 규모의 작은 스텝을 계속 써야 한다. 따라서 느린 모드가 사라지는 $`5\tau_1 \approx 1.25`$ s까지 적분하려면 매우 많은 스텝이 필요하다.

   </details>
