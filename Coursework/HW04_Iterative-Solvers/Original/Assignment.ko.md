# HW4 과제 안내문 — 선형계의 반복법

[레포지토리 안내](../../../README.ko.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **최종 수정일:** 2026-11-24

수치해석 5주차 과제이다. 교수자가 배포한 실습 노트북 `Lab05-IterativeSolver.ipynb`의 빈칸을 채워 실행한다. 노트북은 행렬 분할과 고정점 반복의 틀, 야코비 반복법, 가우스-자이델 방법, 수렴 분석, 과완화법(SOR), 2차원 정상상태 열전도 응용의 여섯 부분으로 이루어지며, 빈칸은 `COMPLETE`라고 적힌 두 주석 줄 사이에 모두 6곳 있다. 노트북의 첫 코드 셀은 `numpy`를 `np`로, `matplotlib.pyplot`을 `plt`로 불러오고 `time`을 불러온다.

아래 각 문제의 코드 블록은 노트북의 셀을 옮긴 것이며, 코드가 비어 있는 자리가 채울 곳이다. 빈칸 옆에 남은 주석은 배포된 노트북에 적혀 있던 설명이다. 배포된 빈 노트북은 [Lab05-IterativeSolver.ipynb](Lab05-IterativeSolver.ipynb)에 있다.

---

<br>

## 문제 1. 행렬 분할 (노트북 1.3절)

모든 정방 행렬은 세 부분으로 나눌 수 있다.

$$
A = D + L + U
$$

여기서 $`D`$는 $`A`$의 대각 부분, $`L`$은 대각과 그 위가 0인 엄밀 하삼각 부분, $`U`$는 대각과 그 아래가 0인 엄밀 상삼각 부분이다.

```math
A = \begin{bmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{bmatrix}
= \underbrace{\begin{bmatrix} a_{11} & & \\ & a_{22} & \\ & & a_{33} \end{bmatrix}}_{D}
+ \underbrace{\begin{bmatrix} 0 & & \\ a_{21} & 0 & \\ a_{31} & a_{32} & 0 \end{bmatrix}}_{L}
+ \underbrace{\begin{bmatrix} 0 & a_{12} & a_{13} \\ & 0 & a_{23} \\ & & 0 \end{bmatrix}}_{U}
```

고전적인 반복법은 모두 $`A\mathbf{x} = \mathbf{b}`$를 고정점 반복 $`\mathbf{x}^{(k+1)} = T\mathbf{x}^{(k)} + \mathbf{c}`$로 고쳐 쓴다. $`T`$는 반복행렬, $`\mathbf{c}`$는 상수 벡터이며, 반복은 스펙트럼 반경 $`\rho(T) = \max_i |\lambda_i(T)| < 1`$일 때 그리고 그때에만 수렴한다. 함수 `split_DLU`의 빈칸을 채워 $`L`$과 $`U`$를 구하라. 같은 셀의 뒷부분은 노트북의 예제 행렬로 함수를 시험한다.

```python
def split_DLU(A):
    """
    Split matrix A into diagonal D, strictly lower triangular L,
    and strictly upper triangular U such that A = D + L + U.

    Returns
    -------
    D, L, U : each (n, n) ndarray
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]

    # --- COMPLETE ---
    D = np.diag(np.diag(A))                          # diagonal matrix (use np.diag twice)
    L =                                              # (strictly lower)
    U =                                              # (strictly upper)
    # --- COMPLETE ---

    return D, L, U


# Test with the example from the lecture notes
A_ex = np.array([[ 3.0, -0.1, -0.2],
                 [ 0.1,  7.0, -0.3],
                 [ 0.3, -0.2, 10.0]])
b_ex = np.array([7.85, -19.3, 71.4])   # true solution: x = [3, -2.5, 7]

D, L, U = split_DLU(A_ex)
print('D =\n', D)
print('L =\n', L)
print('U =\n', U)
print('\nSanity check  D + L + U - A =\n', D + L + U - A_ex)
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
D =
 [[ 3.  0.  0.]
 [ 0.  7.  0.]
 [ 0.  0. 10.]]
L =
 [[ 0.   0.   0. ]
 [ 0.1  0.   0. ]
 [ 0.3 -0.2  0. ]]
U =
 [[ 0.  -0.1 -0.2]
 [ 0.   0.  -0.3]
 [ 0.   0.   0. ]]

Sanity check  D + L + U - A =
 [[0. 0. 0.]
 [0. 0. 0.]
 [0. 0. 0.]]
```

## 문제 2. 야코비 반복법 (노트북 2.2절)

$`i`$번째 식 $`\sum_j a_{ij} x_j = b_i`$에서 $`x_i`$를 분리하면 다음을 얻는다.

```math
x_i^{(k+1)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j \neq i} a_{ij}\, x_j^{(k)}\right), \quad i = 0, 1, \ldots, n-1
```

오른쪽의 값은 모두 이전 반복값 $`\mathbf{x}^{(k)}`$를 쓰므로 새 벡터 전체를 동시에 계산할 수 있다(병렬화가 쉽다). 행렬로 쓰면 $`A = D + (L + U)`$로 분할한 뒤 $`D\mathbf{x} = \mathbf{b} - (L+U)\mathbf{x}`$에서 다음을 얻는다.

$$
\mathbf{x}^{(k+1)} = \underbrace{-D^{-1}(L+U)}_{T_J}\mathbf{x}^{(k)} + D^{-1}\mathbf{b}
$$

수렴 판정에는 흔히 근사 상대 오차를 쓰며, $`\varepsilon_a`$가 사용자가 정한 허용 오차 $`\varepsilon_s`$보다 작아지면 반복을 멈춘다.

```math
\varepsilon_a = \frac{\|\mathbf{x}^{(k+1)} - \mathbf{x}^{(k)}\|_2}{\|\mathbf{x}^{(k+1)}\|_2} \times 100\%
```

함수 `jacobi`의 빈칸을 채워 한 번의 야코비 반복으로 `x_new`를 구하라.

```python
def jacobi(A, b, x0=None, tol=1e-6, max_iter=200):
    """
    Solve Ax = b using the Jacobi iteration.

    Parameters
    ----------
    A        : (n, n) array — coefficient matrix (must be diagonally dominant)
    b        : (n,)   array — right-hand side
    x0       : (n,)   array — initial guess (default: zeros)
    tol      : float        — relative error tolerance (stopping criterion)
    max_iter : int          — maximum number of iterations

    Returns
    -------
    x        : (n,) ndarray — approximate solution
    errors   : list of float — relative approximate error at each iteration
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)

    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)
    errors = []

    d = np.diag(A)                  # diagonal entries a_{ii}
    R = A - np.diag(d)              # off-diagonal part  (L + U)

    for _ in range(max_iter):
        # --- COMPLETE ---
        x_new =                                       
        # --- COMPLETE ---

        # Relative approximate error
        if np.linalg.norm(x_new) > 0:
            ea = np.linalg.norm(x_new - x) / np.linalg.norm(x_new) * 100
        else:
            ea = np.inf
        errors.append(ea)

        x = x_new
        if ea < tol:
            break

    return x, errors
```

이어지는 셀은 노트북의 예제로 함수를 시험하고 반복표를 출력한다.

```python
# Test Jacobi on the lecture-notes example  (true solution: x = [3, -2.5, 7])
x_true = np.array([3.0, -2.5, 7.0])

x_jac, err_jac = jacobi(A_ex, b_ex, tol=1e-6)
x_ref          = np.linalg.solve(A_ex, b_ex)

print(f'Jacobi converged in {len(err_jac)} iterations')
print(f'Solution : {np.round(x_jac, 8)}')
print(f'Reference: {np.round(x_ref, 8)}')
print(f'Max error: {np.max(np.abs(x_jac - x_ref)):.2e}')

# Print iteration table
print(f'\n{"Iter":>5}  {"x1":>12}  {"x2":>12}  {"x3":>12}  {"ea (%)": >12}')
print('-' * 60)
x_k = np.zeros(3)
d_A = np.diag(A_ex)
R_A = A_ex - np.diag(d_A)
for k in range(min(8, len(err_jac))):
    x_k = (b_ex - R_A @ x_k) / d_A
    print(f'{k+1:>5}  {x_k[0]:>12.6f}  {x_k[1]:>12.6f}  {x_k[2]:>12.6f}  {err_jac[k]:>12.6f}')
```

## 문제 3. 가우스-자이델 방법 (노트북 3.1절)

야코비 반복법은 $`\mathbf{x}^{(k+1)}`$의 모든 성분을 이전 값만으로 계산한다. 가우스-자이델 방법은 새로 계산한 성분을 곧바로 쓴다.

```math
x_i^{(k+1)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j=0}^{i-1} a_{ij}\, x_j^{(k+1)} - \sum_{j=i+1}^{n-1} a_{ij}\, x_j^{(k)}\right)
```

첫째 합은 이번 반복에서 이미 갱신된 값을, 둘째 합은 이전 값을 쓴다. 행렬로 쓰면 $`A = (D+L) + U`$로 분할하여 다음을 얻는다.

$$
\mathbf{x}^{(k+1)} = \underbrace{-(D+L)^{-1}U}_{T_{GS}}\mathbf{x}^{(k)} + (D+L)^{-1}\mathbf{b}
$$

| 특징 | 야코비 반복법 | 가우스-자이델 방법 |
|:-----|:-----|:-----|
| 쓰는 값 | 모두 이전 값 $`\mathbf{x}^{(k)}`$ | $`j < i`$는 새 값, $`j > i`$는 이전 값 |
| 수렴 속도 | 느리다 | 보통 약 2배 빠르다 |
| 병렬성 | 쉽다(성분이 서로 독립) | 순차적이다(각 성분이 앞 성분에 의존) |
| 메모리 | $`\mathbf{x}`$의 사본 두 개가 필요하다 | 제자리에서 갱신할 수 있다 |

*표 1. 야코비 반복법과 가우스-자이델 방법의 비교*

노트북은 $`\mathbf{x}^{(0)} = \mathbf{0}`$에서 예제의 첫 반복을 다음과 같이 보인다. 굵게 쓴 값이 이번 반복에서 갱신된 값이다.

```math
\begin{aligned}
x_1^{(1)} &= \frac{7.85 - (-0.1)(0) - (-0.2)(0)}{3} = 2.6167 \\
x_2^{(1)} &= \frac{-19.3 - (0.1)(\mathbf{2.6167}) - (-0.3)(0)}{7} = -2.7945 \\
x_3^{(1)} &= \frac{71.4 - (0.3)(\mathbf{2.6167}) - (-0.2)(\mathbf{-2.7945})}{10} = 7.0056
\end{aligned}
```

함수 `gauss_seidel`의 빈칸을 채워 각 성분을 가장 최근의 값으로 갱신하라.

```python
def gauss_seidel(A, b, x0=None, tol=1e-6, max_iter=200):
    """
    Solve Ax = b using the Gauss-Seidel iteration.

    Parameters
    ----------
    A, b, x0, tol, max_iter  — same as jacobi()

    Returns
    -------
    x      : (n,) ndarray — approximate solution
    errors : list of float — relative approximate error per iteration
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)
    errors = []

    for _ in range(max_iter):
        x_old = x.copy()

        # Update each component using the most recently available values
        for i in range(n):
            # --- COMPLETE ---
            sigma =                                  
            x[i]  =                                  
            # --- COMPLETE ---

        # Relative approximate error
        if np.linalg.norm(x) > 0:
            ea = np.linalg.norm(x - x_old) / np.linalg.norm(x) * 100
        else:
            ea = np.inf
        errors.append(ea)

        if ea < tol:
            break

    return x, errors
```

이어지는 셀은 같은 예제에서 두 방법의 반복 횟수를 비교하고 수렴 이력을 그린다.

```python
# Compare Jacobi vs Gauss-Seidel on the lecture-notes example
x_gs, err_gs = gauss_seidel(A_ex, b_ex, tol=1e-6)

print(f'Jacobi       converged in {len(err_jac):3d} iterations')
print(f'Gauss-Seidel converged in {len(err_gs):3d} iterations')
print(f'GS speedup (iterations): {len(err_jac)/len(err_gs):.1f}x')
print(f'\nGS solution: {np.round(x_gs, 8)}')
print(f'True value : {x_true}')

plt.figure(figsize=(8, 4))
plt.semilogy(range(1, len(err_jac)+1), err_jac, 'b.-', ms=8, label='Jacobi')
plt.semilogy(range(1, len(err_gs)+1),  err_gs,  'r.-', ms=8, label='Gauss-Seidel')
plt.axhline(1e-6, color='gray', ls='--', label='Tolerance $10^{-6}$')
plt.xlabel('Iteration $k$')
plt.ylabel('Relative error $\\varepsilon_a$ (%)')
plt.title('Jacobi vs Gauss-Seidel Convergence')
plt.legend()
plt.grid(True, which='both', ls=':', alpha=0.5)
plt.tight_layout()
plt.show()
```

## 문제 4. 대각 우위 검사 (노트북 4.1절)

모든 행 $`i`$에서 다음이 성립하면 행렬 $`A`$를 엄밀한 대각 우위(strictly diagonally dominant)라 한다.

$$
|a_{ii}| > \sum_{j \neq i} |a_{ij}|
$$

$`A`$가 엄밀한 대각 우위이면 야코비 반복법과 가우스-자이델 방법은 어떤 초기값 $`\mathbf{x}^{(0)}`$에서도 수렴한다. 이는 충분조건일 뿐이므로, 엄밀한 대각 우위가 아닌 행렬에서도 반복법이 수렴할 수 있다. 노트북은 문제 1의 예제 행렬을 행마다 다음과 같이 검사한다.

- **0행:** $`|3| > |-0.1| + |-0.2| = 0.3`$
- **1행:** $`|7| > |0.1| + |-0.3| = 0.4`$
- **2행:** $`|10| > |0.3| + |-0.2| = 0.5`$

함수 `is_diag_dominant`의 빈칸을 채워 각 행의 비대각 성분 절댓값의 합 `offdiag_sum`과 여유 `margins`를 구하라. 여유가 양수인 행은 대각 우위를 만족한다. 같은 셀의 뒷부분은 세 행렬로 함수를 시험한다.

```python
def is_diag_dominant(A, strict=True):
    """
    Check whether A is (strictly) diagonally dominant.

    Returns
    -------
    dominant : bool
    margins  : (n,) array — |a_ii| - sum_{j!=i}|a_ij| for each row
                             Positive margin means the row is dominant.
    """
    A = np.array(A, dtype=float)

    # --- COMPLETE ---
    diag_abs = np.abs(np.diag(A))                    # |a_ii|
    offdiag_sum =                                    # row sums of |A| minus the diagonal
    margins =
    # --- COMPLETE ---

    if strict:
        dominant = bool(np.all(margins > 0))
    else:
        dominant = bool(np.all(margins >= 0))
    return dominant, margins


# Test on well-conditioned and poorly-conditioned matrices
A_good = np.array([[ 4.0, -1.0,  0.5],
                   [ 0.2,  5.0, -1.0],
                   [-0.5,  0.3,  6.0]])

A_bad  = np.array([[ 1.0,  2.0, -3.0],
                   [ 1.5,  2.0,  0.5],
                   [-1.0,  1.0,  1.0]])

dom_good, mar_good = is_diag_dominant(A_good)
dom_bad,  mar_bad  = is_diag_dominant(A_bad)

print('A_good — strictly diagonally dominant:', dom_good)
print('  Row margins:', np.round(mar_good, 4))

print('\nA_bad  — strictly diagonally dominant:', dom_bad)
print('  Row margins:', np.round(mar_bad, 4))

# Check our lecture-notes example
dom_ex, mar_ex = is_diag_dominant(A_ex)
print(f'\nA_ex   — strictly diagonally dominant: {dom_ex},  margins: {np.round(mar_ex, 4)}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
A_good — strictly diagonally dominant: True
  Row margins: [2.5 3.8 5.2]

A_bad  — strictly diagonally dominant: False
  Row margins: [-4.  0. -1.]

A_ex   — strictly diagonally dominant: True,  margins: [2.7 6.6 9.5]
```

## 문제 5. 과완화법(SOR) (노트북 5.1절)

가우스-자이델 방법은 $`x_i^{(k)}`$에서 새 값 $`\tilde{x}_i^{(k+1)}`$로 한 걸음 나아간다. 과완화법(successive over-relaxation, SOR)은 이 갱신을 가중 인수 $`\omega > 1`$만큼 더 나아가게 하여 수렴을 가속한다.

```math
x_i^{(k+1)} = (1 - \omega)\,x_i^{(k)}
+ \frac{\omega}{a_{ii}}\left(b_i - \sum_{j < i} a_{ij}\, x_j^{(k+1)} - \sum_{j > i} a_{ij}\, x_j^{(k)}\right)
```

| $`\omega`$ | 이름 | 효과 |
|:-----|:-----|:-----|
| $`\omega = 1`$ | 가우스-자이델 방법 | 가속이 없다 |
| $`1 < \omega < 2`$ | 과완화(SOR) | 수렴이 빨라진다 |
| $`0 < \omega < 1`$ | 저완화 | 느려지지만 발산하는 경우를 안정시키는 데 쓸 수 있다 |
| $`\omega \geq 2`$ | - | 발산한다 |

*표 2. 가중 인수의 범위에 따른 완화의 종류*

최적의 $`\omega`$는 $`\rho(T_{SOR}(\omega))`$를 최소로 한다. 삼중대각 모형 문제($`n \times n`$ 격자의 5점 라플라스 스텐실)에서는 다음과 같다.

$$
\omega_{\text{opt}} = \frac{2}{1 + \sqrt{1 - \rho(T_J)^2}}
$$

함수 `sor`의 빈칸을 채워 가우스-자이델 값 `x_gs_i`와 이전 값을 섞는 SOR 갱신을 완성하라.

```python
def sor(A, b, omega=1.0, x0=None, tol=1e-6, max_iter=500):
    """
    Solve Ax = b using Successive Over-Relaxation (SOR).

    Parameters
    ----------
    A, b, x0, tol, max_iter — same as gauss_seidel()
    omega  : float — relaxation parameter (0 < omega < 2)

    Returns
    -------
    x      : (n,) ndarray
    errors : list of float
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)
    errors = []

    for _ in range(max_iter):
        x_old = x.copy()

        for i in range(n):
            # Gauss-Seidel update (x_gs_i)
            sigma  = A[i, :i] @ x[:i] + A[i, i+1:] @ x_old[i+1:]
            x_gs_i = (b[i] - sigma) / A[i, i]

            # --- COMPLETE: apply SOR blending ---
            x[i] =
            # --- COMPLETE ---

        if np.linalg.norm(x) > 0:
            ea = np.linalg.norm(x - x_old) / np.linalg.norm(x) * 100
        else:
            ea = np.inf
        errors.append(ea)

        if ea < tol:
            break

    return x, errors
```

이어지는 셀은 노트북 3.2절의 열 네트워크로 야코비 반복법, 가우스-자이델 방법, 여러 $`\omega`$의 SOR을 비교한다. 3.2절은 내부 절점 세 개의 시스템 $`K = \begin{bmatrix} 4 & -1 & 0 \\ -1 & 4 & -1 \\ 0 & -1 & 4 \end{bmatrix}`$, $`\mathbf{f} = [200, 0, 50]^T`$를 `K_heat`, `f_heat`로 세우고, 가우스-자이델 방법으로 허용 오차 $`10^{-8}`$까지 푼 오차 이력을 `err_heat_gs`에 담아 둔다.

```python
# Compare Jacobi, GS, and SOR (ω = 1.1, 1.3, 1.5) on the heat network
results = {}
x_jac_h, err_jac_h = jacobi(K_heat, f_heat, tol=1e-8)
results['Jacobi'] = err_jac_h
results['GS (ω=1.0)'] = err_heat_gs

for w in [1.1, 1.3, 1.5, 1.7]:
    _, err = sor(K_heat, f_heat, omega=w, tol=1e-8)
    results[f'SOR ω={w}'] = err

plt.figure(figsize=(9, 5))
styles = ['b--', 'r-', 'g-', 'm-', 'c-', 'orange']
for (label, errs), style in zip(results.items(), styles):
    plt.semilogy(range(1, len(errs)+1), errs, style, lw=1.8, label=f'{label} ({len(errs)} iters)')
plt.axhline(1e-8, color='gray', ls=':', label='Tolerance')
plt.xlabel('Iteration $k$')
plt.ylabel('Relative error $\\varepsilon_a$ (%)')
plt.title('Jacobi / Gauss-Seidel / SOR Comparison (Heat Network)')
plt.legend(fontsize=8)
plt.grid(True, which='both', ls=':', alpha=0.5)
plt.tight_layout()
plt.show()

print(f'\n{"Method":<14}  {"Iterations":>12}')
print('-' * 28)
for label, errs in results.items():
    print(f'{label:<14}  {len(errs):>12}')
```

## 문제 6. 2차원 열전도의 격자 SOR (노트북 6.2절)

정사각형 판 $`[0,1]\times[0,1]`$의 경계 조건은 아래쪽($`y=0`$)이 $`T = 0`$°C, 위쪽($`y=1`$)이 $`T = 100`$°C, 왼쪽($`x=0`$)과 오른쪽($`x=1`$)이 $`T = 0`$°C이다. 정상상태 온도는 2차원 라플라스 방정식을 만족한다.

$$
\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0
$$

간격 $`h = 1/(n+1)`$인 $`(n+2)\times(n+2)`$ 격자에서 유한 차분으로 이산화하면 다음 식을 얻는다.

$$
T_{i-1,j} + T_{i+1,j} + T_{i,j-1} + T_{i,j+1} - 4T_{i,j} = 0
$$

이 5점 스텐실은 미지수가 $`n^2`$개이고 0이 아닌 성분이 $`5n^2`$개뿐인 크고 희소한 대칭 시스템을 만들어 반복법에 알맞다. 노트북은 $`n^2 \times n^2`$ 행렬을 조립하지 않고 스텐실로 격자를 직접 훑는다. 함수 `solve_laplace_2d`의 빈칸을 채워 `method='sor'`일 때의 갱신을 완성하라.

```python
def solve_laplace_2d(n, T_top=100.0, T_bot=0.0, T_left=0.0, T_right=0.0,
                     method='gs', omega=1.0, tol=1e-6, max_iter=5000):
    """
    Solve the 2D Laplace equation on an (n+2) x (n+2) grid using
    a point iterative method (Gauss-Seidel or SOR).

    Interior nodes: indices 1..n, 1..n
    Returns T : (n+2, n+2) temperature array, iters : iteration count
    """
    T = np.zeros((n + 2, n + 2))

    # Boundary conditions
    T[0,  :]  = T_bot
    T[-1, :]  = T_top
    T[:, 0]   = T_left
    T[:, -1]  = T_right

    iters = 0
    for it in range(max_iter):
        T_old = T.copy()

        # Sweep interior nodes
        for i in range(1, n + 1):          # row index
            for j in range(1, n + 1):      # column index
                # Gauss-Seidel update for the 5-point stencil
                T_gs = 0.25 * (T[i-1, j] + T[i+1, j] + T[i, j-1] + T[i, j+1])

                if method == 'sor':
                    # --- COMPLETE: apply SOR blending ---
                    T[i, j] =
                    # --- COMPLETE ---
                else:
                    T[i, j] = T_gs

        iters += 1

        # Convergence check on interior nodes
        diff = np.max(np.abs(T[1:-1, 1:-1] - T_old[1:-1, 1:-1]))
        if diff < tol:
            break

    return T, iters
```

이어지는 셀은 $`n = 30`$에서 가우스-자이델 방법과 $`\omega = 1.7`$인 SOR을 비교하고 온도 분포를 그린다.

```python
# Solve and visualize for n = 30
n_grid = 30

T_gs_2d,  it_gs  = solve_laplace_2d(n_grid, method='gs',  tol=1e-6)
T_sor_2d, it_sor = solve_laplace_2d(n_grid, method='sor', omega=1.7, tol=1e-6)

print(f'2D Laplace (n={n_grid}, {n_grid**2} unknowns):')
print(f'  Gauss-Seidel : {it_gs:4d} iterations')
print(f'  SOR (ω=1.7)  : {it_sor:4d} iterations')
print(f'  SOR speedup  : {it_gs / it_sor:.1f}x')
print(f'  Max |GS - SOR|: {np.max(np.abs(T_gs_2d - T_sor_2d)):.2e}')

# Plot
x = y = np.linspace(0, 1, n_grid + 2)
X, Y = np.meshgrid(x, y)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, (T_plot, title) in zip(axes,
    [(T_gs_2d,  f'Gauss-Seidel ({it_gs} iters)'),
     (T_sor_2d, f'SOR ω=1.7 ({it_sor} iters)')]):
    cf = ax.contourf(X, Y, T_plot, 20, cmap='hot')
    fig.colorbar(cf, ax=ax, label='T (°C)')
    ax.contour(X, Y, T_plot, 10, colors='k', linewidths=0.5, alpha=0.4)
    ax.set(title=title, xlabel='x', ylabel='y')

plt.suptitle('2D Steady-State Heat Conduction: $\\nabla^2 T = 0$', fontsize=13)
plt.tight_layout()
plt.show()
```

---

<br>

문제를 풀어 본 뒤 [과제 해설](../Assignment-Explanation.ko.md)에서 정답과 풀이를, [완성한 노트북](../Solution/Lab05-IterativeSolver_Completed.ipynb)에서 작성한 코드 전체를 볼 수 있다.
