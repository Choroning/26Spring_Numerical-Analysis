# 13장 실습 — 고유값과 고유벡터

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lab.md) | [한국어](Concepts_Lab.ko.md)

> **최종 수정일:** 2026-11-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 13

> **학습 목표:**
> 1. `numpy.linalg.eig`로 고유쌍을 계산하고 대각합 항등식과 행렬식 항등식으로 검증할 수 있다
> 2. 거듭제곱법(power method)으로 지배 고유값(dominant eigenvalue)을 찾고, 고유값 비 $`|\lambda_2 / \lambda_1|`$에 따른 수렴 속도를 분석할 수 있다
> 3. 역거듭제곱법(inverse power method)과 편이 거듭제곱법(shifted power method)으로 가장 작은 고유값과 임의의 목표값 근처의 고유값을 찾을 수 있다
> 4. `np.linalg.eigh`로 대칭 행렬의 특별한 성질을 활용하고, 스펙트럼 분해 $`A = V\Lambda V^T`$로 저랭크(low-rank) 근사를 만들 수 있다
> 5. `scipy.linalg.eigh(K, M)`로 일반화 고유값 문제 $`K\mathbf{v} = \lambda M\mathbf{v}`$를 풀어 질량-스프링 시스템의 고유 주파수와 모드 형상(mode shape)을 계산할 수 있다
> 6. 선형 ODE 시스템 $`\dot{\mathbf{y}} = A\mathbf{y}`$의 안정성을 고유값으로 분석할 수 있다

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
- [5. 대칭 행렬과 스펙트럼 분해](#5-대칭-행렬과-스펙트럼-분해)
  - [5.1 특별한 성질](#51-특별한-성질)
  - [5.2 `np.linalg.eigh`와 재구성](#52-nplinalgeigh와-재구성)
  - [5.3 저랭크 근사](#53-저랭크-근사)
- [6. 응용: 진동 해석](#6-응용-진동-해석)
  - [6.1 운동 방정식과 일반화 고유값 문제](#61-운동-방정식과-일반화-고유값-문제)
  - [6.2 `scipy.linalg.eigh(K, M)`로 고유 주파수 계산](#62-scipylinalgeighk-m로-고유-주파수-계산)
  - [6.3 모드 형상 시각화](#63-모드-형상-시각화)
  - [6.4 자유 진동 시간 응답](#64-자유-진동-시간-응답)
- [7. 응용: ODE 시스템과 안정성](#7-응용-ode-시스템과-안정성)
  - [7.1 고유쌍을 이용한 일반해](#71-고유쌍을-이용한-일반해)
  - [7.2 안정성 분류](#72-안정성-분류)
  - [7.3 예제: 문제 13.13](#73-예제-문제-1313)
- [8. 방법 비교와 공학적 맥락](#8-방법-비교와-공학적-맥락)
  - [8.1 방법 개요](#81-방법-개요)
  - [8.2 주요 공식](#82-주요-공식)
  - [8.3 공학적 맥락](#83-공학적-맥락)
  - [8.4 실전 팁](#84-실전-팁)
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

대칭 행렬의 고유값이 서로 다르면 `eig`가 반환하는 단위 고유벡터는 정규직교(orthonormal)이다($`V^T V = I`$). 중복 고유값이 있어도 정규직교 고유벡터를 얻으려면 `eigh`를 사용한다(5.2절).

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

**작동하는 이유.** $`A`$가 선형독립인 고유벡터 $`n`$개를 가진다고 가정하고 시작 벡터를 고유기저로 전개하면 다음과 같다.

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

## 5. 대칭 행렬과 스펙트럼 분해

### 5.1 특별한 성질

$`A = A^T`$(대칭)이면 고유구조에 일반적인 경우에는 없는 여러 보장이 생긴다.

| 성질 | 설명 |
|:---------|:----------|
| 실수 고유값 | 모든 $`\lambda_i \in \mathbb{R}`$ |
| 직교 고유벡터 | $`i \neq j`$일 때 $`\mathbf{v}_i^T \mathbf{v}_j = 0`$ (중복 고유값이면 직교하도록 고를 수 있음) |
| 정규직교 고유기저 | $`V^T V = V V^T = I`$ |
| 스펙트럼 분해 | $`A = V \Lambda V^T = \displaystyle\sum_{i=1}^n \lambda_i\, \mathbf{v}_i\mathbf{v}_i^T`$ |

*표 3. 실수 대칭 행렬의 고유구조가 갖는 성질*

이러한 보장 덕분에 대칭 고유값 전용 알고리즘(야코비 고유값 방법, 삼중대각 형태에 대한 QR)은 일반 알고리즘보다 훨씬 효율적이다.

행렬이 대칭(복소수이면 에르미트)임을 알고 있을 때는 `eig` 대신 **`np.linalg.eigh`** 를 사용한다. 이 함수의 특징은 다음과 같다.

- 실수 고유값만 반환한다(반올림 오차에 의한 가짜 허수부가 없다).
- 고유값을 오름차순으로 정렬한다.
- `eig`보다 약 2배 빠르다.
- 정규직교 고유벡터 행렬을 반환한다.

---

<br>

### 5.2 `np.linalg.eigh`와 재구성

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

재구성 오차는 머신 엡실론에 가까운 수준($`10^{-14}`$ 이하)이어야 한다. 이 항등식은 많은 응용의 기초가 된다. $`V`$와 $`\Lambda`$를 얻으면 $`A^k`$, $`e^A`$, $`A^{-1}`$ 같은 표현이 모두 $`V f(\Lambda) V^T`$, 즉 값싼 직교 곱 사이에 낀 대각 연산으로 줄어든다.

---

<br>

### 5.3 저랭크 근사

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

## 6. 응용: 진동 해석

### 6.1 운동 방정식과 일반화 고유값 문제

세 개의 질량이 네 개의 스프링으로 연결되고 양쪽 끝이 벽에 고정된 시스템을 생각한다.

```text
벽 ─[k₁]─ m₁ ─[k₂]─ m₂ ─[k₃]─ m₃ ─[k₄]─ 벽
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

여기서 $`\lambda = \omega^2`$는 **고유 주파수**의 제곱이고 $`\mathbf{v}`$는 대응하는 **모드 형상(mode shape)** 이다. 강성 행렬과 질량 행렬은 다음과 같다.

```math
K = \begin{bmatrix} k_1+k_2 & -k_2 & 0 \\ -k_2 & k_2+k_3 & -k_3 \\ 0 & -k_3 & k_3+k_4 \end{bmatrix},
\qquad
M = \begin{bmatrix} m_1 & & \\ & m_2 & \\ & & m_3 \end{bmatrix}
```

$`K`$와 $`M`$은 모두 대칭 양정치(SPD) 행렬이므로 `scipy.linalg.eigh(K, M)`이 알맞은 도구이다. 이 함수는 실수 고유값을 오름차순으로 반환하고, $`M`$-정규화된 모드 형상($`\phi_i^T M \phi_i = 1`$)을 제공한다.

---

<br>

### 6.2 `scipy.linalg.eigh(K, M)`로 고유 주파수 계산

$`\omega^2`$로부터 각주파수 $`\omega = \sqrt{\lambda}`$ (rad/s), 주파수 $`f = \omega / (2\pi)`$ (Hz), 주기 $`T = 1/f`$ (s)를 얻는다.

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

실행하면 고유 주파수 $`\omega = 3.0078`$, $`5.6889`$, $`7.5447`$ rad/s(주파수 0.4787, 0.9054, 1.2008 Hz)를 얻는다.

> **참고:** 일반화 고유값 문제는 $`\tilde{K} = M^{-1/2} K M^{-1/2}`$와 같은 변환으로 표준 문제로 환원할 수 있다. `scipy.linalg.eigh(K, M)`은 같은 종류의 환원을 $`M`$의 촐레스키 분해로 내부에서 처리하고, $`M`$-정규직교인 모드 형상($`\Phi^T M \Phi = I`$)을 반환한다.

---

<br>

### 6.3 모드 형상 시각화

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

가장 낮은 주파수의 모드에서는 일반적으로 모든 질량이 같은 방향으로 움직이고, 높은 모드일수록 체인을 따라 부호 변화(node)가 많아진다. 부호 변화의 개수는 모드 번호와 함께 증가하며, 이는 진동하는 현의 고조파와 비슷한 패턴이다.

---

<br>

### 6.4 자유 진동 시간 응답

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

응답은 세 고유 주파수에서의 순수한 사인파의 합이다. 감쇠가 없으면 운동이 무한히 계속되지만, 실제 시스템에서는 감쇠로 에너지가 서서히 빠져나가며 대개 높은 모드가 더 빨리 감쇠한다.

---

<br>

## 7. 응용: ODE 시스템과 안정성

### 7.1 고유쌍을 이용한 일반해

선형 1계 시스템 $`\dot{\mathbf{y}} = A\mathbf{y}`$, $`\mathbf{y}(0) = \mathbf{y}_0`$에서 $`A`$가 선형독립인 고유벡터 $`n`$개를 가지면(대각화 가능하면) 일반해는 고유 모드들의 합이다.

```math
\mathbf{y}(t) = \sum_{i=1}^{n} c_i\, \mathbf{v}_i\, e^{\lambda_i t}, \qquad
\mathbf{c} = V^{-1}\mathbf{y}_0
```

각 모드는 고유벡터 방향 $`\mathbf{v}_i`$를 따라 독립적으로 변하며 $`e^{\lambda_i t}`$로 지수적으로 스케일된다. 계수 벡터 $`\mathbf{c}`$는 $`V\mathbf{c} = \mathbf{y}_0`$을 풀어 얻는다.

---

<br>

### 7.2 안정성 분류

모든 모드의 장기 거동은 고유값의 실수부로 결정된다.

| 조건 | 거동 |
|:----------|:---------|
| 모든 $`\text{Re}(\lambda_i) < 0`$ | 해가 **감쇠** → **점근 안정(asymptotically stable)** |
| 하나라도 $`\text{Re}(\lambda_i) > 0`$ | 해가 **증가** → **불안정(unstable)** |
| 모든 $`\text{Re}(\lambda_i) \leq 0`$이고 일부 $`= 0`$ | **경계 안정(marginally stable)** (유계이나 감쇠하지 않음) |

*표 4. 고유값의 실수부에 따른 안정성 분류*

세 번째 행은 7.1절처럼 $`A`$가 대각화 가능하다고 가정한 것이다. 실수부가 0인 고유값에 대응하는 선형독립 고유벡터가 부족하면 해가 $`t`$에 비례해 커질 수 있다.

$`\lambda_i`$의 허수부는 각주파수 $`\text{Im}(\lambda_i)`$의 진동을 만들고, 이 진동은 지수 포락선 $`e^{\text{Re}(\lambda_i) t}`$로 변조된다.

**강성비(stiffness ratio)** $`|\lambda_{\max} / \lambda_{\min}|`$(실수부 기준)는 감쇠 시간 척도가 얼마나 차이 나는지를 나타낸다. 강성비가 크면 시스템이 수치적으로 **강성(stiff)** 을 띤다. 느린 모드가 답을 지배하는데도 양적(explicit) ODE 적분기는 가장 빠른 모드를 따라가기 위해 매우 작은 스텝을 써야 한다.

---

<br>

### 7.3 예제: 문제 13.13

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

## 8. 방법 비교와 공학적 맥락

### 8.1 방법 개요

| 방법 | 찾는 값 | 비용 | 비고 |
|:-------|:------|:-----|:------|
| `np.linalg.eig(A)` | 모든 $`n`$개의 고유쌍 | $`O(n^3)`$ | 일반 행렬; 복소값 가능 |
| `np.linalg.eigh(A)` | 모든 $`n`$개의 고유쌍 | $`O(n^3)`$ | 대칭 전용; 실수, 오름차순 정렬 |
| **거듭제곱법** (2절) | 지배 $`\lambda`$ | 반복당 $`O(n^2)`$ | 속도 $`\approx \vert\lambda_2/\lambda_1\vert`$ |
| **역거듭제곱법과 편이 거듭제곱법** (3절) | 최소 또는 $`\sigma`$에 가까운 $`\lambda`$ | LU 이후 반복당 $`O(n^2)`$ | 내부 고유값에 사용 |
| `sc.linalg.eigh(K, M)` | 일반화 $`K\mathbf{v} = \lambda M\mathbf{v}`$ | $`O(n^3)`$ | SPD 쌍; 고유 주파수 |

*표 5. 고유값 계산 방법 개요*

### 8.2 주요 공식

| 개념 | 공식 |
|:--------|:--------|
| 스펙트럼 분해 | $`A = V \Lambda V^T`$ (대칭 $`A`$, $`V^T V = I`$) |
| 저랭크 근사 | $`A_k = \sum_{i=1}^k \lambda_i\, \mathbf{v}_i \mathbf{v}_i^T`$ |
| Frobenius 오차 | $`\Vert A - A_k\Vert_F = \sqrt{\sum_{i=k+1}^n \lambda_i^2}`$ (고유값을 절댓값 내림차순으로 정렬) |
| 고유 주파수 | $`\omega_i = \sqrt{\lambda_i}`$, $`K\mathbf{v} = \lambda M\mathbf{v}`$로부터 |
| ODE 해 | $`\mathbf{y}(t) = \sum_i c_i\, \mathbf{v}_i\, e^{\lambda_i t}, \quad \mathbf{c} = V^{-1}\mathbf{y}_0`$ |
| 안정성 기준 | 점근 안정 $`\Leftrightarrow`$ 모든 $`\text{Re}(\lambda_i) < 0`$ |

*표 6. 대칭 고유값 문제와 응용의 주요 공식*

### 8.3 공학적 맥락

| 응용 | 고유값의 의미 | 고유벡터의 의미 |
|:------------|:-------------------|:--------------------|
| 구조 진동 | $`\lambda = \omega^2`$ (고유 주파수 제곱) | 모드 형상 |
| ODE 시스템 | 감쇠 / 증가율 $`e^{\lambda t}`$ | 상태 공간의 모드 방향 |
| 반복 해법 (12장) | 반복행렬 $`M`$의 스펙트럼 반경 $`\rho(M)`$ | 수렴 방향 |
| 조건수 (11장) | $`\kappa_2 = \vert\lambda\vert_{\max} / \vert\lambda\vert_{\min}`$ (대칭 행렬, 2-노름) | 해당 없음 |

*표 7. 공학 응용에서 고유값과 고유벡터의 의미*

### 8.4 실전 팁

1. 대칭 행렬에는 항상 `eig`보다 `eigh`를 사용한다. 더 빠르고, 실수 값을 주며, 결과가 정렬되어 있다.
2. 대형 행렬에서 고유값 몇 개만 필요하면 전체 `eig` 대신 거듭제곱법(2절)이나 편이 거듭제곱법(3절)을 사용한다.
3. 일반화 고유값 문제는 $`K\mathbf{v} - \lambda M\mathbf{v} \approx \mathbf{0}`$인지 확인하여 검증한다.
4. 강성비 $`|\lambda_{\max} / \lambda_{\min}|`$가 크면 강성 ODE라는 신호이며, 음적(implicit) 적분기가 필요하다.

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
| 대칭 행렬의 성질 | 실수 대칭 행렬은 실수 고유값과 정규직교 고유벡터를 가지며 $`A = V \Lambda V^T`$로 분해된다. |
| `np.linalg.eigh` | 대칭 행렬 전용 함수로, 실수 고유값을 오름차순으로 정렬하고 정규직교 고유벡터를 반환한다. |
| 저랭크 근사 | 절댓값이 큰 $`k`$개 항만 남긴 $`A_k = \sum_{i=1}^k \lambda_i \mathbf{v}_i \mathbf{v}_i^T`$가 Frobenius 노름에서 최적 랭크-$`k`$ 근사이며, 오차는 버린 고유값 제곱합의 제곱근이다. |
| 일반화 고유값 문제 | 진동계의 $`K\mathbf{v} = \omega^2 M\mathbf{v}`$는 `scipy.linalg.eigh(K, M)`으로 풀며, 고유 주파수는 $`\omega_i = \sqrt{\lambda_i}`$이다. |
| 모드 중첩 | 자유 진동 응답은 각 모드의 사인파 합이며, 진폭은 초기 조건의 $`M`$-사영으로 정해진다. |
| ODE 해 | 대각화 가능한 $`A`$에서 $`\dot{\mathbf{y}} = A\mathbf{y}`$의 해는 $`\mathbf{y}(t) = \sum_i c_i \mathbf{v}_i e^{\lambda_i t}`$, $`\mathbf{c} = V^{-1}\mathbf{y}_0`$이다. |
| 안정성 | 모든 고유값의 실수부가 음수이면 점근 안정하고, 하나라도 양수이면 불안정하다. |
| 강성비 | 실수부 기준 $`\vert\lambda_{\max}/\lambda_{\min}\vert`$가 크면 강성 시스템이며, 양적 적분기는 매우 작은 스텝이 필요하다. |

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

6. **저랭크 근사 오차:** 5.2절의 행렬은 고유값이 21.729156, 13.182326, 6.088518이다. 랭크-1 근사와 랭크-2 근사의 Frobenius 오차를 구하라.

   <details>
   <summary>정답</summary>

   $`\Vert A - A_1\Vert_F = \sqrt{13.182326^2 + 6.088518^2} = \sqrt{173.774 + 37.070} = \sqrt{210.844} \approx 14.520`$이고, $`\Vert A - A_2\Vert_F = \sqrt{6.088518^2} = 6.088518`$이다. 5.3절의 코드 출력 14.520460, 6.088518과 일치한다.

   </details>

7. **일반화 고유값 문제:** 질량 $`m_1 = m_2 = 2`$ kg인 두 질량이 강성 $`k = 10`$ N/m인 스프링 세 개로 두 벽 사이에 연결되어 있다. $`K`$와 $`M`$을 쓰고 고유 주파수와 모드 형상을 구하라.

   <details>
   <summary>정답</summary>

   $`K = \begin{bmatrix} 20 & -10 \\ -10 & 20 \end{bmatrix}`$, $`M = \begin{bmatrix} 2 & 0 \\ 0 & 2 \end{bmatrix}`$이다. $`\det(K - \omega^2 M) = (20 - 2\omega^2)^2 - 100 = 0`$에서 $`20 - 2\omega^2 = \pm 10`$이므로 $`\omega^2 = 5`$ 또는 $`15`$이다. 따라서 $`\omega_1 = \sqrt{5} \approx 2.236`$ rad/s($`f_1 \approx 0.356`$ Hz), $`\omega_2 = \sqrt{15} \approx 3.873`$ rad/s($`f_2 \approx 0.616`$ Hz)이다. 모드 형상은 $`\omega_1`$에서 $`[1, 1]^T`$(동위상), $`\omega_2`$에서 $`[1, -1]^T`$(역위상)이며, 13장 예제 13.4와 같은 결과이다.

   </details>

8. **안정성 분류:** 다음 행렬로 정의한 $`\dot{\mathbf{y}} = A\mathbf{y}`$의 안정성을 판정하라. (a) $`\begin{bmatrix} -2 & 1 \\ 0 & -5 \end{bmatrix}`$ (b) $`\begin{bmatrix} 0 & 1 \\ -4 & 0 \end{bmatrix}`$ (c) $`\begin{bmatrix} 1 & 2 \\ 2 & 1 \end{bmatrix}`$

   <details>
   <summary>정답</summary>

   (a)는 삼각행렬이므로 고유값이 대각 원소 $`-2`$, $`-5`$이고, 모두 음수이므로 점근 안정하다(강성비 2.5). (b)는 $`\lambda^2 + 4 = 0`$에서 $`\lambda = \pm 2i`$이므로 경계 안정하며, 각주파수 2(주기 $`\pi`$)로 감쇠 없이 진동한다. (c)는 $`(1 - \lambda)^2 - 4 = 0`$에서 고유값이 3과 $`-1`$이므로, 양수 고유값 3 때문에 불안정하다.

   </details>

9. **강성비와 시간 상수:** 7.3절의 시스템에서 시간 상수 $`\tau_i = 1/|\lambda_i|`$와 강성비를 구하고, 양적 적분기에 어떤 어려움이 생기는지 설명하라.

   <details>
   <summary>정답</summary>

   대각합은 $`-306`$, 행렬식은 $`1505 - 300 = 1205`$이므로 $`\lambda = \frac{-306 \pm \sqrt{306^2 - 4 \cdot 1205}}{2}`$에서 $`\lambda_1 \approx -3.990`$, $`\lambda_2 \approx -302.010`$이다. 시간 상수는 $`\tau_1 \approx 0.2506`$ s, $`\tau_2 \approx 0.00331`$ s이고, 강성비는 $`302.010/3.990 \approx 75.7`$이다. 빠른 모드는 $`5\tau_2 \approx 0.017`$ s 안에 사라지지만, 양적 적분기는 안정성을 위해 $`\tau_2`$ 규모의 작은 스텝을 계속 써야 한다. 따라서 느린 모드가 사라지는 $`5\tau_1 \approx 1.25`$ s까지 적분하려면 매우 많은 스텝이 필요하다.

   </details>
