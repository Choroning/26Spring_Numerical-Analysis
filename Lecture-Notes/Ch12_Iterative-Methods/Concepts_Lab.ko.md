# 12장 실습 — 고유값 계산을 위한 거듭제곱법

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lab.md) | [한국어](Concepts_Lab.ko.md)

> **최종 수정일:** 2026-10-20
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 13

> **학습 목표:**
> 1. `numpy.linalg.eig`로 고유쌍을 계산하고 대각합 항등식과 행렬식 항등식으로 검증할 수 있다
> 2. 행렬의 지배 고유값(dominant eigenvalue)을 찾기 위해 거듭제곱법(power method)을 구현할 수 있다
> 3. 역거듭제곱법(inverse power method)을 적용해 가장 작은 고유값을 찾을 수 있다
> 4. 편이 거듭제곱법(shifted power method)을 사용해 임의의 목표값 근처의 고유값을 찾을 수 있다
> 5. 고유값 비 $`|\lambda_2 / \lambda_1|`$에 따른 수렴 속도를 분석할 수 있다

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
  - [2.3 Python 구현](#23-python-구현)
  - [2.4 수렴 분석](#24-수렴-분석)
- [3. 역거듭제곱법과 편이 거듭제곱법](#3-역거듭제곱법과-편이-거듭제곱법)
  - [3.1 최소 고유값을 찾는 역거듭제곱법](#31-최소-고유값을-찾는-역거듭제곱법)
  - [3.2 임의의 고유값을 찾는 편이 거듭제곱법](#32-임의의-고유값을-찾는-편이-거듭제곱법)
  - [3.3 Python 구현](#33-python-구현)
  - [3.4 서로 다른 편이로 모든 고유값 찾기](#34-서로-다른-편이로-모든-고유값-찾기)
- [4. 방법 비교와 실전 지침](#4-방법-비교와-실전-지침)
  - [4.1 방법 비교](#41-방법-비교)
  - [4.2 주요 공식](#42-주요-공식)
  - [4.3 실전 팁](#43-실전-팁)
- [요약](#요약)
- [점검 문제](#점검-문제)

---

<br>

## 1. 고유값과 고유벡터

### 1.1 특성 다항식과 주요 항등식

정방행렬 $`A`$의 **고유값(eigenvalue)** $`\lambda`$와 **고유벡터(eigenvector)** $`\mathbf{v}`$는 $`\mathbf{v} \neq \mathbf{0}`$인 조건에서 $`A\mathbf{v} = \lambda\mathbf{v}`$를 만족한다. 이를 정리하면 $`(A - \lambda I)\mathbf{v} = \mathbf{0}`$이 되고, 비자명(nontrivial)한 해가 존재하려면 이 행렬이 특이(singular)해야 한다.

```math
\det(A - \lambda I) = 0 \quad \leftarrow \text{특성 방정식}
```

$`n \times n`$ 행렬의 경우 이 식은 $`\lambda`$에 대한 $`n`$차 다항식이며, 중복도를 포함하여 정확히 $`n`$개의 고유값을 준다(복소수일 수 있다). 검증에 유용한 대수적 항등식은 다음과 같다.

| 성질 | 공식 |
|:---------|:--------|
| 대각합 항등식 | $`\text{tr}(A) = \sum_i a_{ii} = \sum_i \lambda_i`$ |
| 행렬식 항등식 | $`\det(A) = \prod_i \lambda_i`$ |
| 랭크(rank) | 대각화 가능한 행렬(예: 대칭 행렬)에서는 0이 아닌 고유값의 개수(중복도 포함)와 같음 |

*표 1. 고유값 검증에 쓰는 대수적 항등식*

대칭 행렬의 고유값이 서로 다르면 `eig`가 반환하는 단위 고유벡터는 정규직교(orthonormal)이다($`V^T V = I`$). 중복 고유값이 있어도 정규직교 고유벡터를 얻으려면 `eigh`를 사용한다(13장 실습).

---

<br>

### 1.2 기하학적 해석

대부분의 벡터 $`\mathbf{u}`$에 대해 상(image) $`A\mathbf{u}`$는 $`\mathbf{u}`$와 다른 방향을 가리킨다. 고유벡터는 특별해서, $`A\mathbf{v}`$가 $`\mathbf{v}`$와 정확히 평행하고 $`\lambda`$배로 스케일될 뿐이다. 기하학적으로, 대칭 행렬 $`A`$에 대한 사상 $`\mathbf{x} \mapsto A\mathbf{x}`$는 단위원을 축이 고유벡터 방향과 일치하는 타원으로 변환한다(일반 행렬에서는 타원의 축이 특이벡터 방향이다).

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

이 예에서 고유벡터 $`\mathbf{v}_1, \mathbf{v}_2`$는 변환이 순수한 스케일링으로만 작용하는 유일한 방향이며, 다른 모든 방향은 회전된다. 코드는 고유값 4와 2를 출력한다.

---

<br>

### 1.3 `numpy.linalg.eig` 사용

`numpy.linalg.eig(A)`는 $`n`$개의 고유쌍을 모두 $`O(n^3)`$ 시간에 반환한다. 반환된 고유벡터 행렬의 열은 순서대로 고유값 배열의 항목에 대응한다.

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

실행하면 고유값 $`21.72915618`$, $`13.18232626`$, $`6.08851756`$을 얻는다. 고유쌍을 계산한 뒤에는 대각합 항등식과 행렬식 항등식으로 결과를 검증(sanity check)하는 것이 좋다.

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

`trace ≠ sum(lambda)`이거나 `det ≠ prod(lambda)`이면 행렬 입력이 틀렸거나 수치 문제가 발생한 것이다. 이 행렬에서는 대각합 41과 행렬식 1744가 각각 고유값의 합, 곱과 일치한다.

---

<br>

## 2. 거듭제곱법

### 2.1 알고리즘

`np.linalg.eig`는 $`n`$개의 고유쌍을 모두 계산하며, 이는 $`O(n^3)`$ 연산이다. 대형 행렬에서는 $`|\lambda|`$가 가장 큰 **지배(dominant)** 고유값만 필요한 경우가 많다. **거듭제곱법(power method)** 은 행렬-벡터 곱만 사용하여 반복당 $`O(n^2)`$로 이 값을 찾는다.

**왜 작동하는가?** $`A`$가 선형독립인 고유벡터 $`n`$개를 가진다고 가정하고 시작 벡터를 고유기저로 전개하면 다음과 같다.

```math
\mathbf{x}^{(0)} = c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_n\mathbf{v}_n
```

$`A`$를 $`k`$번 곱하면 다음과 같다.

```math
A^k\mathbf{x}^{(0)} = \lambda_1^k\!\left[c_1\mathbf{v}_1 + c_2\!\left(\frac{\lambda_2}{\lambda_1}\right)^{\!k}\mathbf{v}_2 + \cdots\right]
```

$`|\lambda_1| > |\lambda_2| \ge \cdots \ge |\lambda_n|`$이고 $`c_1 \neq 0`$이면 비율 $`(\lambda_2 / \lambda_1)^k \to 0`$이 되어 $`A^k\mathbf{x}^{(0)}`$의 방향이 $`\mathbf{v}_1`$에 정렬된다. 반복당 수렴 비는 대략 $`|\lambda_2 / \lambda_1|`$이다.

**알고리즘** (오버플로를 막기 위해 $`\infty`$-노름으로 정규화). $`\mathbf{x}^{(0)} = \mathbf{1}`$에서 시작하여 $`k = 1, 2, \ldots`$에 대해 다음을 반복한다.

1. $`\mathbf{y}^{(k)} = A\,\mathbf{x}^{(k-1)}`$를 계산한다.
2. $`\lambda^{(k)} = y^{(k)}_{i^*}`$로 고유값 추정값을 정한다. 여기서 $`i^* = \arg\max_i |y_i^{(k)}|`$이다.
3. $`\mathbf{x}^{(k)} = \mathbf{y}^{(k)} / \lambda^{(k)}`$로 $`\|\mathbf{x}^{(k)}\|_\infty = 1`$이 되도록 정규화한다.

---

<br>

### 2.2 단계별 추적

손 추적(hand trace)으로 반복을 구체적으로 확인한다. 위와 같은 $`A`$와 $`\mathbf{x}^{(0)} = \mathbf{1}`$을 사용한다.

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

몇 번의 반복 안에 `lambda est`는 실제 지배 고유값 근처로 다가가고, 정규화된 벡터 `x_k`는 대응하는 고유벡터($`\|\mathbf{x}\|_\infty = 1`$)에 접근한다. 7번째 반복의 추정값은 21.838011로, 참값 21.729156과의 차이가 약 0.11이다.

---

<br>

### 2.3 Python 구현

상대 변화 중단 조건과 플롯용 이력(history)을 갖춘 재사용 가능한 함수로 감싼다.

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

> **참고:** 절댓값이 가장 큰 원소($`\infty`$-노름)로 정규화하면 오버플로를 막고 `lam`이 곧 고유값 추정값이 되므로 별도의 레일리 몫(Rayleigh quotient)이 필요 없다.

---

<br>

### 2.4 수렴 분석

거듭제곱법의 이론적 수렴 비는 $`|\lambda_2 / \lambda_1|`$이며, 이 비가 1에 가까울수록 수렴이 느리다.

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

로그 오차 그래프는 기울기가 $`\log|\lambda_2/\lambda_1|`$인 직선에 가까워야 한다. 이 행렬에서는 $`|\lambda_2/\lambda_1| = 0.606665`$이므로 오차가 약 1.4번 반복마다 절반으로 줄고, `tol=1e-12`에서 51번 만에 수렴한다. 고유값이 잘 분리되어 있으면 빠른 선형 수렴을 보이지만, 스펙트럼 간격이 작아 비가 1에 가까우면(예: $`|\lambda_2/\lambda_1| = 0.99`$) 수천 번의 반복이 필요하다.

---

<br>

## 3. 역거듭제곱법과 편이 거듭제곱법

### 3.1 최소 고유값을 찾는 역거듭제곱법

$`A\mathbf{v} = \lambda\mathbf{v}`$이고 $`A`$가 가역이면 $`A^{-1}\mathbf{v} = (1/\lambda)\mathbf{v}`$이다. $`A`$의 가장 작은 $`|\lambda_i|`$는 $`A^{-1}`$의 가장 큰 $`|1/\lambda_i|`$가 된다. 따라서 $`A^{-1}`$에 거듭제곱법을 적용하면 $`1/\lambda_{\min}`$으로 수렴하고, $`\lambda_{\min} = 1 / \mu_{\max}`$로 되돌린다.

**구현:** $`A^{-1}`$을 명시적으로 만들지 않는다. 각 반복은 한 번 계산해 둔 $`A`$의 LU 분해로 $`A\mathbf{y} = \mathbf{x}`$를 푼다. $`O(n^3)`$ 분해 후 반복당 비용은 $`O(n^2)`$이다.

---

<br>

### 3.2 임의의 고유값을 찾는 편이 거듭제곱법

사용자가 고른 **편이(shift)** $`\sigma`$에 대해 $`(A - \sigma I)`$에 역거듭제곱법을 적용하므로, 이 방법은 정확히는 편이 역거듭제곱법(shifted inverse power method)이다. $`(A - \sigma I)^{-1}`$의 고유값은 $`1/(\lambda_i - \sigma)`$이다. $`A`$의 고유값 중 $`\sigma`$에 가장 가까운 것이 가장 큰 $`|1/(\lambda_i - \sigma)|`$를 가지므로 반복은 그 값으로 수렴한다.

```math
\lambda = \sigma + \frac{1}{\mu}, \qquad \mu = (A - \sigma I)^{-1}\text{의 지배 고유값}
```

$`\sigma = 0`$이면 역거듭제곱법이 된다(최소 $`|\lambda|`$ 탐색). $`\sigma`$를 임의의 목표값 근처로 두면 그 값에 가장 가까운 고유값을 고르므로, 가장 크거나 가장 작은 고유값이 아닌 내부 고유값도 찾을 수 있다.

---

<br>

### 3.3 Python 구현

하나의 함수가 역거듭제곱법($`\sigma = 0`$)과 일반적인 편이 방법을 모두 처리한다.

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

> **참고:** LU 분해는 반복문 바깥에서 한 번만 수행하므로 각 반복에는 $`O(n^2)`$ 삼각 풀이만 필요하다. 반복문 안에서 `np.linalg.inv(B) @ x`를 쓰면 더 느리고 수치적으로도 덜 정확하므로 쓰지 않는다.

---

<br>

### 3.4 서로 다른 편이로 모든 고유값 찾기

관심 있는 각 고유값 근처에 편이를 두면, 편이 거듭제곱법은 `eig`를 호출하지 않고도 소형 행렬의 모든 고유값을 복원할 수 있다.

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

편이를 고유값 가까이에 둘수록 $`|1/(\lambda_i - \sigma)|`$가 다른 역수들을 압도하여 유효한 고유값 비가 작아지고 수렴이 빨라진다. 수렴 비는 $`|\lambda_j - \sigma| / |\lambda_k - \sigma|`$이며, $`\lambda_j`$는 $`\sigma`$에 가장 가까운 고유값, $`\lambda_k`$는 두 번째로 가까운 고유값이다. 위 출력에서 $`\sigma = 24`$는 이 비가 약 0.21이라 18번 만에 수렴하지만, $`\sigma = 10`$은 두 고유값 6.09와 13.18의 거의 가운데에 있어 비가 약 0.81이 되고 127번이 걸린다. 실제로는 게르슈고린 원판(Gershgorin disc) 등으로 대략적인 $`\sigma`$를 정하며, $`\sigma`$가 목표 고유값에 다른 고유값보다 충분히 가까우면 몇 번의 반복으로 수렴한다.

---

<br>

## 4. 방법 비교와 실전 지침

### 4.1 방법 비교

| 방법 | 찾는 값 | 반복당 비용 | 수렴 속도 |
|:-------|:------|:------------------|:-----------------|
| `np.linalg.eig(A)` | 모든 $`n`$개의 고유쌍 | 총 $`O(n^3)`$ | 해당 없음 (직접 계산) |
| **거듭제곱법** | 지배 $`\lambda`$ (최대 $`\vert\lambda\vert`$) | $`O(n^2)`$ | $`\vert\lambda_2 / \lambda_1\vert`$ |
| **역거듭제곱법** | 최소 $`\vert\lambda\vert`$ | LU 이후 $`O(n^2)`$ | $`\vert\lambda_{\min} / \lambda_{\text{second-smallest}}\vert`$ |
| **편이 거듭제곱법** | $`\sigma`$에 가장 가까운 $`\lambda`$ | LU 이후 $`O(n^2)`$ | 편이 품질에 의존 |

*표 2. 고유값 계산 방법의 비교*

### 4.2 주요 공식

1. **특성 방정식:** $`\det(A - \lambda I) = 0`$
2. **대각합 항등식:** $`\text{tr}(A) = \sum_i \lambda_i`$ (검증용)
3. **행렬식 항등식:** $`\det(A) = \prod_i \lambda_i`$
4. **거듭제곱 단계:** $`\lambda^{(k)} = y^{(k)}_{i^*}`$ ($`\mathbf{y}^{(k)} = A\mathbf{x}^{(k-1)}`$에서 절댓값이 가장 큰 성분, 부호 포함), $`\mathbf{x}^{(k)} = A\mathbf{x}^{(k-1)} / \lambda^{(k)}`$
5. **편이 복원:** $`\lambda = \sigma + 1/\mu`$, 여기서 $`\mu`$는 $`(A - \sigma I)^{-1}`$의 지배 고유값

### 4.3 실전 팁

- 고유쌍은 항상 대각합 항등식, 행렬식 항등식, 잔차 $`\|A\mathbf{v} - \lambda\mathbf{v}\|`$로 검증한다.
- $`(A - \sigma I)`$는 반복문 바깥에서 한 번만 인수분해하고, 역행렬을 명시적으로 구하지 않는다.
- 스펙트럼 간격이 작으면 거듭제곱법의 수렴이 느리므로, 편이로 유효 간격을 넓힌다.
- 전체 스펙트럼을 계산하지 않고 내부 고유값을 찾을 때는 편이 거듭제곱법을 사용한다.

---

<br>

## 요약

| 개념 | 핵심 요약 |
|:-----|:-----|
| 특성 방정식 | 고유값은 $`\det(A - \lambda I) = 0`$의 근이며, $`n \times n`$ 행렬은 중복도를 포함해 $`n`$개의 고유값을 가진다. |
| 대각합과 행렬식 항등식 | $`\text{tr}(A) = \sum_i \lambda_i`$와 $`\det(A) = \prod_i \lambda_i`$로 계산한 고유값을 검증한다. |
| `np.linalg.eig` | 모든 고유쌍을 $`O(n^3)`$에 계산하며, 고유벡터 행렬의 열이 고유값 배열의 순서와 대응한다. |
| 거듭제곱법 | $`A`$를 반복해서 곱하고 절댓값이 가장 큰 성분으로 정규화하여 지배 고유값과 고유벡터를 반복당 $`O(n^2)`$로 구한다. |
| 수렴 비 | 거듭제곱법의 오차는 반복마다 약 $`\vert\lambda_2/\lambda_1\vert`$배로 줄며, 이 비가 1에 가까울수록 수렴이 느리다. |
| 역거듭제곱법 | $`A^{-1}`$에 거듭제곱법을 적용하되 LU 분해로 $`A\mathbf{y} = \mathbf{x}`$를 풀어 최소 $`\vert\lambda\vert`$를 구한다. |
| 편이 거듭제곱법 | $`(A - \sigma I)^{-1}`$에 거듭제곱법을 적용하고 $`\lambda = \sigma + 1/\mu`$로 되돌려 $`\sigma`$에 가장 가까운 고유값을 구한다. |

---

<br>

## 점검 문제

1. **거듭제곱법 손 계산:** $`A = \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}`$에 $`\mathbf{x}^{(0)} = [1, 0]^T`$로 거듭제곱법을 세 번 반복하여 고유값 추정값을 구하라. 추정값은 어디로 수렴하며, 반복당 수렴 비는 얼마인가?

   <details>
   <summary>정답</summary>

   $`A\mathbf{x}^{(0)} = [2, 1]^T`$이므로 $`\lambda^{(1)} = 2`$, $`\mathbf{x}^{(1)} = [1, 0.5]^T`$이다. $`A\mathbf{x}^{(1)} = [2.5, 2]^T`$이므로 $`\lambda^{(2)} = 2.5`$, $`\mathbf{x}^{(2)} = [1, 0.8]^T`$이다. $`A\mathbf{x}^{(2)} = [2.8, 2.6]^T`$이므로 $`\lambda^{(3)} = 2.8`$이다. $`A`$의 고유값은 3과 1이므로 추정값은 3으로 수렴하고, 반복당 수렴 비는 $`\vert\lambda_2/\lambda_1\vert = 1/3`$이다.

   </details>

2. **대각합과 행렬식 항등식:** 1.3절의 행렬 $`A`$에서 두 고유값 21.729156과 13.182326을 알 때, 대각합 항등식으로 나머지 고유값을 구하고 행렬식 항등식으로 확인하라.

   <details>
   <summary>정답</summary>

   $`\text{tr}(A) = 20 + 9 + 12 = 41`$이므로 나머지 고유값은 $`41 - 21.729156 - 13.182326 = 6.088518`$이다. 세 고유값의 곱은 $`21.729156 \times 13.182326 \times 6.088518 \approx 1744.0`$이며 $`\det(A) = 1744`$와 일치한다.

   </details>

3. **편이 선택:** 같은 행렬에 $`\sigma = 15`$로 편이 거듭제곱법을 적용하면 어떤 고유값으로 수렴하는가? 수렴했을 때의 $`\mu`$와 반복당 수렴 비를 구하라.

   <details>
   <summary>정답</summary>

   각 고유값과 15의 거리는 $`\vert 21.729 - 15\vert = 6.729`$, $`\vert 13.182 - 15\vert = 1.818`$, $`\vert 6.089 - 15\vert = 8.911`$이므로 가장 가까운 13.182326으로 수렴한다. 이때 $`\mu = 1/(13.182326 - 15) \approx -0.5502`$이고, 수렴 비는 $`1.818/6.729 \approx 0.270`$이다.

   </details>

4. **시작 벡터의 조건:** 거듭제곱법의 시작 벡터에 지배 고유벡터 $`\mathbf{v}_1`$ 성분이 없으면($`c_1 = 0`$) 이론적으로 어떻게 되는가? 실제 부동 소수점 계산에서는 어떤가?

   <details>
   <summary>정답</summary>

   2.1절의 전개에서 $`c_1 = 0`$이면 $`\lambda_1^k`$ 항이 없으므로, 정확한 산술에서는 반복이 남은 성분 중 절댓값이 가장 큰 고유값(예: $`c_2 \neq 0`$이고 $`\vert\lambda_2\vert > \vert\lambda_3\vert`$이면 $`\lambda_2`$)으로 수렴한다. 부동 소수점 계산에서는 반올림 오차가 $`\mathbf{v}_1`$ 성분을 조금씩 만들어 내므로 대개 결국 $`\lambda_1`$로 수렴하지만, 그만큼 많은 반복이 필요할 수 있다.

   </details>

5. **역거듭제곱법의 비용:** $`n \times n`$ 밀집 행렬에서 역거듭제곱법을 $`m`$번 반복할 때, LU 분해를 한 번 해 두는 방법과 반복마다 행렬을 새로 분해하여 푸는 방법의 총비용을 비교하라.

   <details>
   <summary>정답</summary>

   LU 분해를 한 번 하면 $`O(n^3)`$ 분해 한 번과 반복당 $`O(n^2)`$ 삼각 풀이로 총 $`O(n^3 + mn^2)`$이다. 반복마다 분해하면 총 $`O(mn^3)`$이 들어, $`m`$이 클수록 차이가 커진다.

   </details>
