# HW5 과제 안내문 — 고유값과 고유벡터

[레포지토리 안내](../../../README.ko.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **최종 수정일:** 2026-11-17

수치해석 6주차 과제이다. 교수자가 배포한 실습 노트북 `Lab06-Eigen.ipynb`의 빈칸을 채워 실행한다. 노트북은 고유값과 고유벡터, 거듭제곱법, 역거듭제곱법과 편이 거듭제곱법, 대칭 행렬과 스펙트럼 분해, 진동 해석 응용, ODE 시스템과 안정성 응용의 여섯 부분으로 이루어지며, 빈칸은 `COMPLETE`라고 적힌 두 주석 줄 사이에 모두 5곳 있다. 노트북의 첫 코드 셀은 `numpy`를 `np`로, `scipy`를 `sc`로, `matplotlib.pyplot`을 `plt`로 불러온다.

아래 각 문제의 코드 블록은 노트북의 셀을 옮긴 것이며, 코드가 비어 있는 자리가 채울 곳이다. 배포된 빈 노트북은 [Lab06-Eigen.ipynb](Lab06-Eigen.ipynb)에 있다.

- **제출 기한:** 배포 후 1주
- **제출 방법:** 빈칸을 채워 실행한 노트북을 온라인으로 제출한다.

---

<br>

## 문제 1. 고유값과 고유벡터의 계산 (노트북 1.3절)

정사각 행렬 $`A`$의 고유값 $`\lambda`$와 고유벡터 $`\mathbf{v}`$는 다음을 만족한다.

$$
A\mathbf{v} = \lambda\mathbf{v}, \qquad \mathbf{v} \neq \mathbf{0}
$$

기하학적으로 $`A`$는 고유벡터의 크기만 바꾸고 방향은 돌리지 않는다. 식을 정리하면 $`(A - \lambda I)\mathbf{v} = \mathbf{0}`$이고, 자명하지 않은 해가 있으려면 $`A - \lambda I`$가 특이 행렬이어야 하므로 특성 방정식 $`\det(A - \lambda I) = 0`$을 얻는다. $`n \times n`$ 행렬에서 이는 $`\lambda`$에 관한 $`n`$차 다항식이며, 중복도를 세면 복소수를 포함해 정확히 $`n`$개의 고유값을 준다.

노트북 1.3절은 `numpy.linalg.eig`로 강의 노트의 예제(13장 문제 13.1) 행렬의 고유값과 고유벡터를 구한다. 다음 셀의 빈칸을 채워라.

```python
# Lecture-notes example  (Ch. 13, Problem 13.1)
A = np.array([[20.0,  3.0,  2.0],
              [ 3.0,  9.0,  4.0],
              [ 2.0,  4.0, 12.0]])

# --- COMPLETE ---
eigvals, eigvecs = 
# --- COMPLETE ---

print('Eigenvalues:')
for i, lam in enumerate(eigvals):
    print(f'  lambda_{i+1} = {lam:.8f}')

print('\nEigenvectors (columns):')
print(np.round(eigvecs, 6))

# Verify  A v = lambda v
print('\nVerification  max|A v - lambda v|:')
for i in range(len(eigvals)):
    res = A @ eigvecs[:, i] - eigvals[i] * eigvecs[:, i]
    print(f'  lambda_{i+1}: {np.max(np.abs(res)):.2e}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
Eigenvalues:
  lambda_1 = 21.72915618
  lambda_2 = 13.18232626
  lambda_3 = 6.08851756

Eigenvectors (columns):
[[ 0.89819   0.427012  0.104481]
 [ 0.30972  -0.44601  -0.839732]
 [ 0.311976 -0.786598  0.532855]]

Verification  max|A v - lambda v|:
  lambda_1: 3.55e-15
  lambda_2: 7.11e-15
  lambda_3: 1.33e-15
```

## 문제 2. 거듭제곱법 (노트북 2.1절)

`np.linalg.eig`는 $`n`$개의 고유쌍을 모두 계산하며 $`O(n^3)`$ 연산이 든다. 큰 행렬에서는 절댓값이 가장 큰 지배 고유값만 필요한 경우가 많고, **거듭제곱법(power method)** 은 반복마다 $`O(n^2)`$으로 이를 계산한다.

원리는 다음과 같다. $`\mathbf{x}^{(0)}`$을 고유기저로 $`\mathbf{x}^{(0)} = c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots`$와 같이 전개하면 다음이 성립한다.

```math
A^k\mathbf{x}^{(0)} = \lambda_1^k\!\left[c_1\mathbf{v}_1 + c_2\!\left(\frac{\lambda_2}{\lambda_1}\right)^{\!k}\mathbf{v}_2 + \cdots\right]
```

$`|\lambda_1| > |\lambda_2|`$이면 $`(\lambda_2/\lambda_1)^k \to 0`$이므로 $`A^k\mathbf{x}^{(0)}`$은 $`\mathbf{v}_1`$ 방향으로 정렬된다. 수렴 속도는 반복마다 대략 $`|\lambda_2/\lambda_1|`$이다.

알고리즘은 오버플로를 막기 위해 무한 노름으로 정규화한다. $`\mathbf{x}^{(0)} = \mathbf{1}`$에서 시작하여 $`k = 1, 2, \ldots`$에 대해 다음을 반복한다.

1. $`\mathbf{y}^{(k)} = A\,\mathbf{x}^{(k-1)}`$
2. $`\lambda^{(k)} = y^{(k)}_{i^*}`$, 여기서 $`i^* = \arg\max_i |y_i^{(k)}|`$ (고유값 추정)
3. $`\mathbf{x}^{(k)} = \mathbf{y}^{(k)} / \lambda^{(k)}`$ (정규화, $`\|\mathbf{x}^{(k)}\|_\infty = 1`$)

다음 함수의 빈칸을 채워 거듭제곱법을 완성하라.

```python
def power_method(A, x0=None, tol=1e-8, max_iter=500):
    """
    Find the dominant eigenvalue and eigenvector of A using the Power Method.

    Parameters
    ----------
    A        : (n, n) array
    x0       : (n,)   initial guess  (default: ones)
    tol      : float  relative change in lambda as stopping criterion
    max_iter : int

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
        # --- COMPLETE ---
        lam =                                         
        x   =                                         
        # --- COMPLETE ---

        lam_hist.append(float(lam))
        if abs(lam - lam_old) / abs(lam) < tol:
            break
        lam_old = lam

    return float(lam), x, lam_hist
```

다음 셀은 완성한 함수를 문제 1의 행렬 `A`에 적용하고 수렴을 그래프로 분석한다. 그래프를 그리는 부분은 생략했다.

```python
# Test power_method and analyze convergence
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
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
Converged in 51 iterations
Dominant eigenvalue (Power Method): 21.729156182199
True value (numpy)                : 21.729156182169
Absolute error                    : 3.02e-11

Theoretical convergence rate |λ₂/λ₁| = 0.606665
  -> error halves every 1.4 iterations
```

## 문제 3. 역거듭제곱법과 편이 거듭제곱법 (노트북 3.1, 3.2절)

$`A\mathbf{v} = \lambda\mathbf{v}`$이면 $`A^{-1}\mathbf{v} = (1/\lambda)\mathbf{v}`$이다. 따라서 $`A`$의 가장 작은 $`|\lambda_i|`$는 $`A^{-1}`$의 가장 큰 $`|1/\lambda_i|`$가 되고, $`A^{-1}`$에 거듭제곱법을 적용하면 $`1/\lambda_{\min}`$으로 수렴한다. 구현에서는 $`A^{-1}`$을 직접 만들지 않는다. 반복마다 $`A\mathbf{y} = \mathbf{x}`$를 푸는데, LU 분해를 한 번 계산해 두면 이 풀이는 $`O(n^2)`$만 든다.

사용자가 정한 **편이(shift)** $`\sigma`$에 대해 $`(A - \sigma I)`$에 역거듭제곱법을 적용하면 다음이 성립한다.

$$
\text{eigenvalues of }(A-\sigma I)^{-1}: \quad \frac{1}{\lambda_i - \sigma}
$$

$`\sigma`$에 가장 가까운 고유값이 가장 큰 $`|1/(\lambda_i - \sigma)|`$를 가지므로 방법은 그 고유값으로 수렴한다. $`\sigma = 0`$으로 두면 역거듭제곱법이 된다. 다음 함수의 빈칸을 채워라.

```python
def shifted_power_method(A, sigma=0.0, x0=None, tol=1e-8, max_iter=500):
    """
    Find the eigenvalue of A nearest to `sigma` (Shifted Inverse Power Method).

    Algorithm:
      1. Build B = A - sigma*I
      2. Factor B = LU  (once)
      3. Power iteration on B^{-1}: solve B y = x, normalize
      4. True eigenvalue:  lambda = sigma + 1/mu   where mu = max(|y|)

    sigma = 0  ->  Inverse Power Method  (finds smallest |lambda|)
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]
    x = np.ones(n) if x0 is None else np.array(x0, dtype=float)

    # --- COMPLETE ---
    B      = A - sigma * np.eye(n)                   # shifted matrix
    lu_obj = 
    # --- COMPLETE ---

    lam_hist = []
    lam_old  = np.inf

    for _ in range(max_iter):
        # --- COMPLETE ---
        y   = 
        mu  = y[np.argmax(np.abs(y))]                # largest element of y  (= 1 / (lambda - sigma))
        x   = y / mu
        lam = 
        # --- COMPLETE ---

        lam_hist.append(float(lam))
        if abs(lam - lam_old) / (abs(lam) + 1e-30) < tol:
            break
        lam_old = lam

    return float(lam), x, lam_hist
```

다음 셀은 서로 다른 편이로 세 고유값을 모두 찾는다.

```python
# Find all three eigenvalues using different shifts
lam_sorted = np.sort(np.real(eigvals))
print(f'All eigenvalues (numpy, ascending): {np.round(lam_sorted, 6)}')
print()

tests = [
    ('Power Method   (dominant)', power_method,           dict(),           np.max(np.abs(eigvals))),
    ('sigma=0        (smallest)', shifted_power_method,   dict(sigma=0.0),  lam_sorted[0]),
    ('sigma=10       (middle)',   shifted_power_method,   dict(sigma=10.0), lam_sorted[1]),
    ('sigma=24       (largest)', shifted_power_method,   dict(sigma=24.0), lam_sorted[2]),
]

print(f'{"Method":<34}  {"Found lambda":>14}  {"True lambda":>12}  {"Error":>10}  {"Iters":>6}')
print('-' * 84)
for label, func, kwargs, lam_t in tests:
    lam_f, _, hist = func(A, tol=1e-12, **kwargs)
    print(f'{label:<34}  {lam_f:>14.8f}  {lam_t:>12.8f}  {abs(lam_f-lam_t):>10.2e}  {len(hist):>6}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
All eigenvalues (numpy, ascending): [ 6.088518 13.182326 21.729156]

Method                                Found lambda   True lambda       Error   Iters
------------------------------------------------------------------------------------
Power Method   (dominant)              21.72915618   21.72915618    3.02e-11      51
sigma=0        (smallest)               6.08851756    6.08851756    2.68e-12      38
sigma=10       (middle)                13.18232626   13.18232626    5.06e-12     127
sigma=24       (largest)               21.72915618   21.72915618    1.34e-12      18
```

## 문제 4. 진동계의 고유 주파수 (노트북 5.1절)

양 끝이 벽에 고정된 스프링 네 개로 질량 세 개가 연결되어 있다.

```text
Wall ─[k₁]─ m₁ ─[k₂]─ m₂ ─[k₃]─ m₃ ─[k₄]─ Wall
```

뉴턴의 제2법칙에서 자유 진동 방정식 $`M\ddot{\mathbf{x}} + K\mathbf{x} = \mathbf{0}`$을 얻는다. 조화 운동 $`\mathbf{x}(t) = \mathbf{v}\sin(\omega t)`$를 가정하면 다음 **일반화 고유값 문제(generalized eigenvalue problem)** 가 된다.

$$
K\mathbf{v} = \omega^2 M\mathbf{v}
$$

여기서 $`\lambda = \omega^2`$은 고유 주파수의 제곱이고 $`\mathbf{v}`$는 모드 형상이다. 강성 행렬과 질량 행렬은 다음과 같다.

```math
K = \begin{bmatrix} k_1+k_2 & -k_2 & 0 \\ -k_2 & k_2+k_3 & -k_3 \\ 0 & -k_3 & k_3+k_4 \end{bmatrix},
\qquad
M = \begin{bmatrix} m_1 & & \\ & m_2 & \\ & & m_3 \end{bmatrix}
```

`scipy.linalg.eigh(K, M)`는 대칭 양의 정부호 행렬 쌍의 일반화 고유값 문제를 풀고 고유값을 오름차순으로 돌려준다. 다음 셀의 빈칸에 각진동수 `omega_n`(rad/s)과 주파수 `freq_hz`(Hz)를 채워라.

```python
# System parameters
k1, k2, k3, k4 = 40.0, 30.0, 20.0, 10.0   # N/m
m1, m2, m3     =  2.0,  1.5,  1.0          # kg

K_vib = np.array([
    [k1+k2,   -k2,      0],
    [  -k2, k2+k3,    -k3],
    [    0,   -k3,  k3+k4]
])
M_vib = np.diag([m1, m2, m3])

# Solve  K v = omega^2 M v
omega_sq, modes = sc.linalg.eigh(K_vib, M_vib)   # eigenvalues ascending

# --- COMPLETE ---
omega_n =                                        #  — rad/s
freq_hz =                                        #  — Hz
# --- COMPLETE ---

print('=== Natural Frequencies ===')
print(f'{"Mode":>5}  {"omega^2":>14}  {"omega (rad/s)":>14}  {"f (Hz)":>10}  {"T (s)":>10}')
print('-' * 60)
for i in range(3):
    print(f'{i+1:>5}  {omega_sq[i]:>14.4f}  {omega_n[i]:>14.4f}  {freq_hz[i]:>10.4f}  {1/freq_hz[i]:>10.4f}')

print('\nMode shapes (columns = eigenvectors, M-normalized):')
print(np.round(modes, 6))

# Verify generalized eigenvalue equation
print('\nVerification  max|K v - omega^2 M v|:')
for i in range(3):
    v   = modes[:, i]
    res = K_vib @ v - omega_sq[i] * M_vib @ v
    print(f'  Mode {i+1}: {np.max(np.abs(res)):.2e}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
=== Natural Frequencies ===
 Mode         omega^2   omega (rad/s)      f (Hz)       T (s)
------------------------------------------------------------
    1          9.0470          3.0078      0.4787      2.0889
    2         32.3636          5.6889      0.9054      1.1045
    3         56.9227          7.5447      1.2008      0.8328

Mode shapes (columns = eigenvectors, M-normalized):
[[-0.329371  0.484607 -0.395817]
 [-0.569877  0.085174  0.578492]
 [-0.543958 -0.720715 -0.429743]]

Verification  max|K v - omega^2 M v|:
  Mode 1: 4.44e-15
  Mode 2: 7.11e-15
  Mode 3: 1.42e-14
```

## 문제 5. ODE 시스템의 모드 상수 (노트북 6.1, 6.2절)

$`\dot{\mathbf{y}} = A\mathbf{y}`$, $`\mathbf{y}(0) = \mathbf{y}_0`$의 해는 다음과 같다.

```math
\mathbf{y}(t) = \sum_{i=1}^{n} c_i\, \mathbf{v}_i\, e^{\lambda_i t}, \qquad \mathbf{c} = V^{-1}\mathbf{y}_0
```

안정성은 고유값의 실수부로 판정한다.

| 조건 | 거동 |
|:-----|:-----|
| 모든 $`\text{Re}(\lambda_i) < 0`$ | 해가 감쇠하며 점근적으로 안정하다 |
| 어느 하나라도 $`\text{Re}(\lambda_i) > 0`$ | 해가 커지며 불안정하다 |
| 모든 $`\text{Re}(\lambda_i) \leq 0`$이고 일부가 $`= 0`$ | 한계 안정이다 |

*표 1. 고유값의 실수부에 따른 선형 ODE 시스템의 안정성*

예제(강의 문제 13.13)는 다음 시스템이다.

```math
\begin{cases} \dot{y}_1 = -5y_1 + 3y_2 \\ \dot{y}_2 = 100y_1 - 301y_2 \end{cases},
\qquad \mathbf{y}(0) = [50,\; 100]^T
```

다음 셀의 빈칸에 모드 상수 $`\mathbf{c} = V^{-1}\mathbf{y}_0`$을 구하는 코드를 채워라.

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

# --- COMPLETE: find modal constants c = V^{-1} y0 ---
c_ode = 
# --- COMPLETE ---

print(f'\nModal constants: c1 = {c_ode[0]:.4f},  c2 = {c_ode[1]:.4f}')

# Time solution
t = np.linspace(0, 1.5, 800)
y_t = np.zeros((2, len(t)))
for i in range(2):
    y_t += np.real(c_ode[i] * np.outer(V_ode[:, i], np.exp(lam_ode[i] * t)))

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(t, y_t[0], 'b-', lw=2, label='$y_1(t)$')
ax.plot(t, y_t[1], 'r-', lw=2, label='$y_2(t)$')
ax.set(xlabel='Time $t$', ylabel='$y$', title='ODE Solution $\\dot{\\mathbf{y}} = A\\mathbf{y}$')
ax.legend(); ax.grid(True, alpha=0.4)

ax = axes[1]
ax.plot(y_t[0], y_t[1], 'k-', lw=2)
ax.plot(y0[0], y0[1], 'go', ms=12, label='Start $\\mathbf{y}(0)$', zorder=5)
ax.plot(y_t[0,-1], y_t[1,-1], 'rs', ms=12, label='End', zorder=5)
ax.set(xlabel='$y_1$', ylabel='$y_2$', title='Phase Portrait')
ax.legend(); ax.grid(True, alpha=0.4)

plt.tight_layout()
plt.show()

print(f'\nTime constants:  tau_1 = 1/|lambda_1| = {1/abs(lam_ode[0]):.4f} s  (slow mode)')
print(f'                 tau_2 = 1/|lambda_2| = {1/abs(lam_ode[1]):.6f} s  (fast mode)')
print(f'Stiffness ratio  |lambda_2 / lambda_1| = {abs(lam_ode[1]/lam_ode[0]):.1f}')
```

배포 노트북에는 이 셀이 출력하는 글자가 다음과 같이 실려 있다. 두 출력 사이에는 해의 시간 변화와 위상 궤적 그래프가 있다.

```text
Eigenvalues:
  lambda_1 = -3.989933   -> mode decays
  lambda_2 = -302.010067   -> mode decays
System is STABLE

Modal constants: c1 = 53.6413,  c2 = 82.8879

Time constants:  tau_1 = 1/|lambda_1| = 0.2506 s  (slow mode)
                 tau_2 = 1/|lambda_2| = 0.003311 s  (fast mode)
Stiffness ratio  |lambda_2 / lambda_1| = 75.7
```

---

<br>

문제를 풀어 본 뒤 [과제 해설](../Assignment-Explanation.ko.md)에서 정답과 풀이를, [완성한 노트북](../Solution/Lab06-Eigen_Completed.ipynb)에서 작성한 코드 전체를 볼 수 있다.
