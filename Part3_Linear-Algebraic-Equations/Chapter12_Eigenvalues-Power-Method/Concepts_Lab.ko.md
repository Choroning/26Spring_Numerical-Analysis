# 제12장 실습 — 고유값 계산을 위한 거듭제곱법

> **최종 수정일:** 2026-04-14
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 13 (Power Method)

> **선수 지식**: [프로그래밍언어] Python, NumPy, SciPy. [선형대수학] 고유값 / 고유벡터, LU 분해 (제10장).
>
> **학습 목표**:
> 1. `numpy.linalg.eig`로 고유쌍을 계산하고 대각합 / 행렬식 항등식을 검증할 수 있다
> 2. 행렬의 지배 고유값(dominant eigenvalue)을 찾기 위해 거듭제곱법(Power Method)을 구현할 수 있다
> 3. 역거듭제곱법(Inverse Power Method)을 적용해 가장 작은 고유값을 찾을 수 있다
> 4. 편이 거듭제곱법(Shifted Power Method)을 사용해 임의 목표값 근처의 고유값을 찾을 수 있다
> 5. 스펙트럼 간격 $|\lambda_2 / \lambda_1|$에 따른 수렴 속도를 분석할 수 있다

---

<br>

## 목차

- [1. 고유값과 고유벡터](#1-고유값과-고유벡터)
  - [1.1 특성 다항식과 주요 항등식](#11-특성-다항식과-주요-항등식)
  - [1.2 기하학적 해석](#12-기하학적-해석)
  - [1.3 `numpy.linalg.eig` 사용](#13-numpylinalgeig-사용)
- [2. 거듭제곱법](#2-거듭제곱법)
  - [2.1 알고리즘](#21-알고리즘)
  - [2.2 단계별 추적](#22-단계별-추적)
  - [2.3 파이썬 구현](#23-파이썬-구현)
  - [2.4 수렴 분석](#24-수렴-분석)
- [3. 역 / 편이 거듭제곱법](#3-역--편이-거듭제곱법)
  - [3.1 역거듭제곱법 — 최소 고유값](#31-역거듭제곱법--최소-고유값)
  - [3.2 편이 거듭제곱법 — 임의의 고유값](#32-편이-거듭제곱법--임의의-고유값)
  - [3.3 파이썬 구현](#33-파이썬-구현)
  - [3.4 서로 다른 편이로 모든 고유값 찾기](#34-서로-다른-편이로-모든-고유값-찾기)
- [요약](#요약)

---

<br>

## 1. 고유값과 고유벡터

### 1.1 특성 다항식과 주요 항등식

정방행렬 $A$의 **고유값(eigenvalue)** $\lambda$와 **고유벡터(eigenvector)** $\mathbf{v}$는 $\mathbf{v} \neq \mathbf{0}$인 조건에서 $A\mathbf{v} = \lambda\mathbf{v}$를 만족합니다. 이를 정리하면 $(A - \lambda I)\mathbf{v} = \mathbf{0}$이 되고, 비자명(nontrivial)한 해가 존재하려면 행렬이 특이(singular)해야 합니다:

$$\det(A - \lambda I) = 0 \quad \leftarrow \text{특성 방정식}$$

$n \times n$ 행렬의 경우 이 식은 $\lambda$에 대한 $n$차 다항식이며, 정확히 $n$개의 고유값을 제공합니다 (중복도 포함, 복소수일 수 있음). 검증에 유용한 대수적 항등식은 다음과 같습니다:

| 성질 | 공식 |
|:---------|:--------|
| 대각합 항등식 | $\text{tr}(A) = \sum_i a_{ii} = \sum_i \lambda_i$ |
| 행렬식 항등식 | $\det(A) = \prod_i \lambda_i$ |
| 계수(rank) | 0이 아닌 고유값의 개수와 같음 |

대칭 행렬의 경우, NumPy가 반환하는 고유벡터는 정규직교(orthonormal)입니다: $V^T V = I$.

---

<br>

### 1.2 기하학적 해석

대부분의 벡터 $\mathbf{u}$에 대해 상(image) $A\mathbf{u}$는 $\mathbf{u}$와 다른 방향을 가리킵니다. 고유벡터는 특별합니다: $A\mathbf{v}$는 $\mathbf{v}$와 정확히 평행하며 단지 $\lambda$만큼 스케일됩니다. 기하학적으로, 사상 $\mathbf{x} \mapsto A\mathbf{x}$는 단위원을 타원으로 변환하며, 이 타원의 축은 고유벡터와 정렬됩니다:

```python
import numpy as np
import matplotlib.pyplot as plt

A_geo = np.array([[3.0, 1.0], [1.0, 3.0]])
eigvals_geo, eigvecs_geo = np.linalg.eig(A_geo)

theta   = np.linspace(0, 2 * np.pi, 300)
circle  = np.array([np.cos(theta), np.sin(theta)])
ellipse = A_geo @ circle

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
colors = ['crimson', 'steelblue']

ax = axes[0]
ax.plot(circle[0], circle[1], 'b-', lw=1.5)
for i, (lam, c) in enumerate(zip(eigvals_geo, colors)):
    v = eigvecs_geo[:, i]
    ax.annotate('', xy=v, xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color=c, lw=2.5))
ax.set(xlim=[-1.5, 1.5], ylim=[-1.5, 1.5], aspect='equal',
       title='Before: unit circle + eigenvectors')

ax = axes[1]
ax.plot(ellipse[0], ellipse[1], 'b-', lw=1.5)
for i, (lam, c) in enumerate(zip(eigvals_geo, colors)):
    av = A_geo @ eigvecs_geo[:, i]
    ax.annotate('', xy=av, xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color=c, lw=2.5))
ax.set(xlim=[-5, 5], ylim=[-5, 5], aspect='equal',
       title='After A·: ellipse + scaled eigenvectors')

plt.tight_layout()
plt.show()
print(f'Eigenvalues: {eigvals_geo}')
```

고유벡터 $\mathbf{v}_1, \mathbf{v}_2$는 변환이 순수 스케일링으로만 작용하는 유일한 방향입니다. 다른 모든 방향은 회전됩니다.

---

<br>

### 1.3 `numpy.linalg.eig` 사용

`numpy.linalg.eig(A)`는 모든 $n$개의 고유쌍을 $O(n^3)$ 시간에 반환합니다. 반환된 고유벡터 행렬의 열은 (순서대로) 고유값 배열의 항목과 대응됩니다:

```python
# 강의노트 예제 (13장, 문제 13.1)
A = np.array([[20.0,  3.0,  2.0],
              [ 3.0,  9.0,  4.0],
              [ 2.0,  4.0, 12.0]])

eigvals, eigvecs = np.linalg.eig(A)

print('Eigenvalues:')
for i, lam in enumerate(eigvals):
    print(f'  lambda_{i+1} = {lam:.8f}')

print('\nEigenvectors (columns):')
print(np.round(eigvecs, 6))

# A v = lambda v 검증
print('\nVerification  max|A v - lambda v|:')
for i in range(len(eigvals)):
    res = A @ eigvecs[:, i] - eigvals[i] * eigvecs[:, i]
    print(f'  lambda_{i+1}: {np.max(np.abs(res)):.2e}')
```

고유쌍을 계산한 후에는 반드시 대각합과 행렬식 항등식을 검증(sanity check)해야 합니다:

```python
print('=== Key Properties ===')
print(f'trace(A)               = {np.trace(A):.6f}')
print(f'sum of eigenvalues     = {np.sum(eigvals):.6f}   <- must match')
print(f'det(A)                 = {np.linalg.det(A):.6f}')
print(f'product of eigenvalues = {np.prod(eigvals):.6f}   <- must match')

# 대칭 행렬의 경우 고유벡터는 정규직교
VTV = eigvecs.T @ eigvecs
print('\nV^T V (= I for symmetric A):')
print(np.round(VTV, 8))
```

`trace ≠ sum(lambda)`이거나 `det ≠ prod(lambda)`라면 무언가 잘못된 것입니다 — 행렬 입력이 틀렸거나 수치 문제가 발생한 것입니다.

---

<br>

## 2. 거듭제곱법

### 2.1 알고리즘

`np.linalg.eig`는 모든 $n$개의 고유쌍을 계산합니다 — $O(n^3)$ 연산. 대형 행렬에서는 가장 큰 $|\lambda|$를 가지는 **지배(dominant)** 고유값만 필요한 경우가 많습니다. **거듭제곱법(Power Method)**은 행렬-벡터 곱만 사용하여 반복당 $O(n^2)$로 이를 찾습니다.

**왜 작동하는가?** 시작 벡터를 고유기저로 전개하면:

$$\mathbf{x}^{(0)} = c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_n\mathbf{v}_n$$

$A$를 $k$번 곱한 후:

$$A^k\mathbf{x}^{(0)} = \lambda_1^k\!\left[c_1\mathbf{v}_1 + c_2\!\left(\frac{\lambda_2}{\lambda_1}\right)^{\!k}\mathbf{v}_2 + \cdots\right]$$

$|\lambda_1| > |\lambda_2| \ge \cdots \ge |\lambda_n|$이면 비율 $(\lambda_2 / \lambda_1)^k \to 0$이 되고 $A^k\mathbf{x}^{(0)}$는 $\mathbf{v}_1$에 정렬됩니다. 반복당 수렴 속도는 대략 $|\lambda_2 / \lambda_1|$입니다.

**알고리즘** (오버플로우 방지를 위해 $\infty$-노름으로 정규화). $\mathbf{x}^{(0)} = \mathbf{1}$로 시작하여, $k = 1, 2, \ldots$에 대해:

1. $\mathbf{y}^{(k)} = A\,\mathbf{x}^{(k-1)}$
2. $\lambda^{(k)} = y^{(k)}_{i^*}$, 여기서 $i^* = \arg\max_i |y_i^{(k)}|$ — 고유값 추정값
3. $\mathbf{x}^{(k)} = \mathbf{y}^{(k)} / \lambda^{(k)}$ — $\|\mathbf{x}^{(k)}\|_\infty = 1$이 되도록 정규화

---

<br>

### 2.2 단계별 추적

수동 추적(hand trace)은 반복을 구체적으로 보여줍니다. 위와 동일한 $A$와 $\mathbf{x}^{(0)} = \mathbf{1}$ 사용:

```python
x_k = np.ones(3)
print(f'{"Iter":>5}  {"lambda est":>12}  {"x1":>10}  {"x2":>10}  {"x3":>10}')
print('-' * 57)
for k in range(1, 8):
    y     = A @ x_k
    lam_k = y[np.argmax(np.abs(y))]
    x_k   = y / lam_k
    print(f'{k:>5}  {lam_k:>12.6f}  {x_k[0]:>10.6f}  {x_k[1]:>10.6f}  {x_k[2]:>10.6f}')

print(f'\nTrue dominant lambda (numpy): {np.max(np.abs(eigvals)):.6f}')
```

몇 번의 반복 안에 `lambda est`는 실제 지배 고유값 근처로 수렴하고 정규화된 벡터 `x_k`는 대응되는 고유벡터($\|\mathbf{x}\|_\infty = 1$)에 접근합니다.

---

<br>

### 2.3 파이썬 구현

상대 변화 중단 조건과 플롯용 이력(history)을 가진 재사용 가능한 함수로 감쌉니다:

```python
def power_method(A, x0=None, tol=1e-8, max_iter=500):
    """
    Find the dominant eigenvalue and eigenvector of A using the Power Method.

    Returns
    -------
    lam      : float        — dominant eigenvalue
    v        : (n,) ndarray — eigenvector (inf-norm normalized)
    lam_hist : list         — eigenvalue estimate history
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]
    x = np.ones(n) if x0 is None else np.array(x0, dtype=float)
    lam_hist = []
    lam_old  = np.inf

    for _ in range(max_iter):
        y   = A @ x
        lam = y[np.argmax(np.abs(y))]   # 절댓값이 가장 큰 원소
        x   = y / lam                    # 정규화

        lam_hist.append(float(lam))
        if abs(lam - lam_old) / abs(lam) < tol:
            break
        lam_old = lam

    return float(lam), x, lam_hist
```

> **참고:** 절댓값이 가장 큰 원소($\infty$-노름)로 정규화하면 오버플로우를 방지하고 `lam`이 고유값 추정값 자체가 됩니다 — 별도의 레일리 몫(Rayleigh quotient)이 필요 없습니다.

---

<br>

### 2.4 수렴 분석

거듭제곱법의 이론적 수렴 속도는 $|\lambda_2 / \lambda_1|$입니다: 이 비율이 1에 가까울수록 수렴이 느립니다.

```python
lam_pm, v_pm, hist_pm = power_method(A, tol=1e-12)
lam_true = float(eigvals[np.argmax(np.abs(eigvals))])

print(f'Converged in {len(hist_pm)} iterations')
print(f'Dominant eigenvalue (Power Method): {lam_pm:.12f}')
print(f'True value (numpy)                : {lam_true:.12f}')
print(f'Absolute error                    : {abs(lam_pm - lam_true):.2e}')

sorted_lam = np.sort(np.abs(eigvals))[::-1]
rate = sorted_lam[1] / sorted_lam[0]
print(f'\nTheoretical convergence rate |λ₂/λ₁| = {rate:.6f}')
print(f'  -> error halves every {-np.log(2)/np.log(rate):.1f} iterations')

# 수렴 그래프
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
iters = range(1, len(hist_pm) + 1)

ax1.plot(iters, hist_pm, 'b.-', ms=6)
ax1.axhline(lam_true, color='r', ls='--', lw=2, label=f'True λ = {lam_true:.4f}')
ax1.set(xlabel='Iteration', ylabel='λ estimate', title='Eigenvalue Convergence')
ax1.legend(); ax1.grid(True, alpha=0.4)

errors = np.abs(np.array(hist_pm) - lam_true)
ax2.semilogy(iters, errors, 'r.-', ms=6, label='Actual error')
ref = errors[0] * rate ** np.array(range(len(errors)))
ax2.semilogy(iters, ref, 'k--', lw=1.5, label='|λ₂/λ₁|^k reference')
ax2.set(xlabel='Iteration', ylabel='|λ^(k) - λ_true|',
        title='Absolute Error (log scale)')
ax2.legend(); ax2.grid(True, which='both', ls=':', alpha=0.5)

plt.tight_layout()
plt.show()
```

로그-오차 그래프는 기울기가 $\log|\lambda_2/\lambda_1|$인 직선이어야 합니다. 고유값들이 잘 분리되어 있으면 빠른 선형 수렴을 보이지만, 스펙트럼 간격이 작으면 (예: $|\lambda_2/\lambda_1| = 0.99$) 수천 번의 반복이 필요합니다.

---

<br>

## 3. 역 / 편이 거듭제곱법

### 3.1 역거듭제곱법 — 최소 고유값

$A\mathbf{v} = \lambda\mathbf{v}$이고 $A$가 가역이면 $A^{-1}\mathbf{v} = (1/\lambda)\mathbf{v}$입니다. $A$의 가장 작은 $|\lambda_i|$는 $A^{-1}$의 가장 큰 $|1/\lambda_i|$가 됩니다. 따라서 $A^{-1}$에 거듭제곱법을 적용하면 $1/\lambda_{\min}$으로 수렴하고, $\lambda_{\min} = 1 / \mu_{\max}$로 역산됩니다.

**구현:** $A^{-1}$을 명시적으로 만들지 마십시오. 각 반복은 한 번 계산된 $A$의 LU 분해를 사용해 $A\mathbf{y} = \mathbf{x}$를 풉니다. $O(n^3)$ 분해 후 반복당 비용은 $O(n^2)$입니다.

---

<br>

### 3.2 편이 거듭제곱법 — 임의의 고유값

사용자가 선택한 **편이(shift)** $\sigma$에 대해 $(A - \sigma I)$에 역거듭제곱법을 적용합니다. $(A - \sigma I)^{-1}$의 고유값은 $1/(\lambda_i - \sigma)$입니다. $A$의 고유값 중 $\sigma$에 가장 가까운 것이 가장 큰 $|1/(\lambda_i - \sigma)|$를 가지므로 반복은 그 값으로 수렴합니다:

$$\lambda = \sigma + \frac{1}{\mu}, \qquad \mu = (A - \sigma I)^{-1}\text{의 지배 고유값}$$

$\sigma = 0$이면 역거듭제곱법이 됩니다 (최소 $|\lambda|$ 탐색). $\sigma$를 임의 목표값 근처로 설정하면 그와 가장 가까운 내부 고유값을 선택합니다.

---

<br>

### 3.3 파이썬 구현

단일 함수가 역거듭제곱법($\sigma = 0$)과 일반 편이 변형 모두를 처리합니다:

```python
import scipy as sc

def shifted_power_method(A, sigma=0.0, x0=None, tol=1e-8, max_iter=500):
    """
    Find the eigenvalue of A nearest to `sigma` (Shifted Inverse Power Method).

    Algorithm:
      1. Build B = A - sigma * I
      2. Factor B = LU  (once)
      3. Power iteration on B^{-1}: solve B y = x, normalize
      4. True eigenvalue:  lambda = sigma + 1/mu   where mu = max(|y|)

    sigma = 0  ->  Inverse Power Method  (finds smallest |lambda|)
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]
    x = np.ones(n) if x0 is None else np.array(x0, dtype=float)

    B      = A - sigma * np.eye(n)
    lu_obj = sc.linalg.lu_factor(B)

    lam_hist = []
    lam_old  = np.inf

    for _ in range(max_iter):
        y   = sc.linalg.lu_solve(lu_obj, x)
        mu  = y[np.argmax(np.abs(y))]     # = 1 / (lambda - sigma)
        x   = y / mu
        lam = sigma + 1.0 / mu

        lam_hist.append(float(lam))
        if abs(lam - lam_old) / (abs(lam) + 1e-30) < tol:
            break
        lam_old = lam

    return float(lam), x, lam_hist
```

> **참고:** LU 분해는 반복문 바깥에서 한 번만 수행되므로, 각 반복은 $O(n^2)$ 삼각 풀이만 필요합니다. 반복문 안에서 `np.linalg.inv(B) @ x`를 사용하지 마십시오 — 더 느리고 수치적으로도 덜 정확합니다.

---

<br>

### 3.4 서로 다른 편이로 모든 고유값 찾기

관심 있는 각 고유값 근처의 편이를 선택하면, 편이 거듭제곱법은 `eig`를 호출하지 않고도 소형 행렬의 모든 고유값을 복원할 수 있습니다:

```python
lam_sorted = np.sort(np.real(eigvals))
print(f'All eigenvalues (numpy, ascending): {np.round(lam_sorted, 6)}\n')

tests = [
    ('Power Method   (dominant)', power_method,         dict(),            np.max(np.abs(eigvals))),
    ('sigma=0        (smallest)', shifted_power_method, dict(sigma=0.0),   lam_sorted[0]),
    ('sigma=10       (middle)',   shifted_power_method, dict(sigma=10.0),  lam_sorted[1]),
    ('sigma=24       (largest)',  shifted_power_method, dict(sigma=24.0),  lam_sorted[2]),
]

print(f'{"Method":<34}  {"Found λ":>14}  {"True λ":>12}  {"Error":>10}  {"Iters":>6}')
print('-' * 84)
for label, func, kwargs, lam_t in tests:
    lam_f, _, hist = func(A, tol=1e-12, **kwargs)
    print(f'{label:<34}  {lam_f:>14.8f}  {lam_t:>12.8f}  {abs(lam_f-lam_t):>10.2e}  {len(hist):>6}')
```

고유값 근처에 편이를 배치하면 수렴이 극적으로 빨라집니다: $|1/(\lambda_i - \sigma)|$가 다른 역수들을 압도하므로 유효 스펙트럼 간격이 매우 커집니다. 실제로는 대략적인 $\sigma$ 추정(예: Gershgorin 원반 사용)만으로도 몇 번의 반복으로 수렴합니다.

---

<br>

## 요약

| 방법 | 찾는 값 | 반복당 비용 | 수렴 속도 |
|:-------|:------|:------------------|:-----------------|
| `np.linalg.eig(A)` | 모든 $n$개의 고유쌍 | 총 $O(n^3)$ | — |
| **거듭제곱법** | 지배 $\lambda$ (최대 $\|\lambda\|$) | $O(n^2)$ | $\|\lambda_2 / \lambda_1\|$ |
| **역거듭제곱법** | 최소 $\|\lambda\|$ | LU 이후 $O(n^2)$ | $\|\lambda_{\min} / \lambda_{\text{second-smallest}}\|$ |
| **편이 거듭제곱법** | $\sigma$에 가장 가까운 $\lambda$ | LU 이후 $O(n^2)$ | 편이 품질에 의존 |

주요 공식:

1. **특성 방정식**: $\det(A - \lambda I) = 0$
2. **대각합 항등식**: $\text{tr}(A) = \sum_i \lambda_i$ (검증용)
3. **행렬식 항등식**: $\det(A) = \prod_i \lambda_i$
4. **거듭제곱 단계**: $\lambda^{(k)} = \max\|A\mathbf{x}^{(k-1)}\|_\infty, \quad \mathbf{x}^{(k)} = A\mathbf{x}^{(k-1)} / \lambda^{(k)}$
5. **편이 복원**: $\lambda = \sigma + 1/\mu$, 여기서 $\mu$는 $(A - \sigma I)^{-1}$의 지배 고유값

실전 팁:

- 항상 대각합 / 행렬식 항등식과 잔차 $\|A\mathbf{v} - \lambda\mathbf{v}\|$로 고유쌍을 검증하십시오
- $(A - \sigma I)$는 반복문 바깥에서 한 번만 인수분해하십시오 — 명시적으로 역행렬을 구하지 마십시오
- 스펙트럼 간격이 작으면 거듭제곱법 수렴이 느립니다; 편이로 유효 간격을 증폭하십시오
- 전체 스펙트럼을 계산하지 않고 내부 고유값을 탐색할 때 편이 거듭제곱법을 사용하십시오

---
