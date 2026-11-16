# HW7 과제 안내문 — 최적화 방법

[레포지토리 안내](../../../README.ko.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **최종 수정일:** 2026-11-16

수치해석 9주차 과제이다. 교수자가 배포한 실습 노트북 `Lab08-Optimization.ipynb`의 빈칸을 채워 실행한다. 노트북은 국소 최솟값과 전역 최솟값, 황금 분할 탐색, 포물선 보간, 최급 하강법, `scipy.optimize`와 공학 응용, 제약 최적화(선택)의 여섯 부분으로 이루어지며, 빈칸은 `COMPLETE`라고 적힌 두 주석 줄 사이에 모두 3곳 있다. 원기둥 설계 문제를 확장하는 5.4절의 Exercise 셀은 코드가 채워진 채 배포되었고, 코드가 비어 있는 두 Exercise 셀은 이 안내문 끝의 [노트북의 추가 연습 문제](#노트북의-추가-연습-문제)에 옮겼다. 노트북의 첫 코드 셀은 `numpy`를 `np`로, `matplotlib.pyplot`을 `plt`로 불러오고, `scipy.optimize`의 `minimize_scalar`와 `minimize`를 불러온다.

아래 각 문제의 코드 블록은 노트북의 셀을 옮긴 것이며, 코드가 비어 있는 자리가 채울 곳이다. 배포된 빈 노트북은 [Lab08-Optimization.ipynb](Lab08-Optimization.ipynb)에 있다.

- **제출 기한:** 배포 후 1주
- **제출 방법:** 빈칸을 채워 실행한 노트북을 온라인으로 제출한다.

---

<br>

## 노트북의 예제 함수 (노트북 1절)

**최적화(optimization)** 는 목적 함수를 최소화(또는 최대화)하는 결정 변수의 값을 찾는다. 매끄러운 $`f`$에서 $`x^*`$가 국소 최솟값이기 위한 **필요조건**은 1차원에서 $`f'(x^*) = 0`$, 다차원에서 $`\nabla f(x^*) = \mathbf{0}`$이다. 1차원에서 국소 **최솟값**이기 위한 **충분조건**은 $`f'(x^*) = 0`$이면서 $`f''(x^*) > 0`$인 것이다. 대부분의 수치 방법은 **국소** 최솟값을 찾으며, **전역** 최솟값을 찾기는 훨씬 어렵다.

노트북은 다음 함수를 예제로 계속 쓴다.

$$
f(x) = \frac{x^2}{10} - 2\sin(x)
$$

```python
# Running example: f(x) = x²/10 - 2sin(x)
f1d = lambda x: x**2 / 10.0 - 2.0 * np.sin(x)
```

노트북은 이 함수의 그래프에 $`[0, 4]`$의 전역 최솟값과 $`[5, 8]`$의 국소 최솟값을 표시하고, $`[0, 4]`$에서 시작한 수치 방법은 전역 최솟값을, $`[5, 8]`$에서 시작한 방법은 국소 최솟값만 찾는다고 설명한다.

## 문제 1. 황금 분할 탐색 (노트북 2.1~2.3절)

근 찾기의 이분법처럼, 황금 분할 탐색은 **단봉** 함수의 최솟값을 담은 구간 $`[x_l, x_u]`$를 줄여 간다. 단계마다 두 내부 점 $`x_1 > x_2`$를 대칭으로 둔다.

$$
d = (\varphi - 1)(x_u - x_l), \qquad x_1 = x_l + d, \quad x_2 = x_u - d
$$

여기서 $`\varphi = \frac{1+\sqrt{5}}{2} \approx 1.618`$은 **황금비**이다. 최소화에서는 $`f(x_1)`$과 $`f(x_2)`$를 비교해 구간을 다음과 같이 갱신한다.

| 조건 | 남기는 구간 | 근거 |
|:-----|:-----|:-----|
| $`f(x_1) < f(x_2)`$ | $`[x_2, x_u]`$ | 최솟값이 $`x_2`$의 오른쪽에 있다 |
| $`f(x_1) \geq f(x_2)`$ | $`[x_l, x_1]`$ | 최솟값이 $`x_1`$의 왼쪽에 있다 |

*표 1. 황금 분할 탐색의 구간 갱신 규칙*

갱신한 뒤에는 이전 내부 점 하나가 새 내부 점이 되므로 반복마다 함수를 **한 번만** 새로 계산하면 된다. 반복마다 구간은 $`(\varphi - 1) \approx 0.618`$배로 줄어드는 선형 수렴이며, 근사 상대 오차는 $`\varepsilon_a = (2 - \varphi) \left|\dfrac{\Delta x}{x_{\rm opt}}\right|`$, $`\Delta x = x_u - x_l`$이다. 다음 함수의 두 빈칸을 채워라.

```python
def goldmin(f, xl, xu, es=1e-7, maxit=100):
    """
    Find the minimum of a unimodal f on [xl, xu] using Golden-Section Search.

    Parameters
    ----------
    f       : callable — unimodal objective function
    xl, xu  : float   — initial bracket
    es      : float   — stopping criterion (approx. relative error)
    maxit   : int     — maximum iterations

    Returns
    -------
    xopt    : float  — minimizer estimate
    fopt    : float  — f(xopt)
    ea      : float  — final approximate relative error
    n_iter  : int    — iterations performed
    history : list of (iter, xopt, f(xopt), ea, interval_width)
    """
    phi = (1.0 + np.sqrt(5)) / 2.0

    d  = (phi - 1) * (xu - xl)
    x1 = xl + d;  f1 = f(x1)
    x2 = xu - d;  f2 = f(x2)

    history = []
    xopt = x1 if f1 < f2 else x2

    for i in range(maxit):
        xint = xu - xl   # current interval width

        if f1 < f2:
            # minimum is to the right of x2: new bracket = [x2, xu]
            xopt = x1
            # --- COMPLETE: update bracket and interior points ---
            xl = x2          # shift left bound to old x2
            x2 = x1          # old x1 becomes new left interior point
            f2 = f1
            x1 =              #   — new right interior point
            f1 = f(x1)
            # --- COMPLETE ---
        else:
            # minimum is to the left of x1: new bracket = [xl, x1]
            xopt = x2
            # --- COMPLETE: update bracket and interior points ---
            xu = x1          # shift right bound to old x1
            x1 = x2          # old x2 becomes new right interior point
            f1 = f2
            x2 =              #  — new left interior point
            f2 = f(x2)
            # --- COMPLETE ---

        ea = (2.0 - phi) * abs(xint / xopt) if xopt != 0 else (2.0 - phi) * abs(xint)
        history.append((i + 1, xopt, f(xopt), ea, xu - xl))

        if ea < es:
            break

    return xopt, f(xopt), ea, i + 1, history
```

다음 셀은 완성한 함수를 예제 함수에 적용한다.

```python
# Test goldmin on the running example
x_true = 1.4275517797083348   # known global minimum (from scipy)

xopt, fopt, ea, n_iter, history = goldmin(f1d, xl=0.0, xu=4.0, es=1e-9, maxit=100)

print(f'Minimum found : x* = {xopt:.10f}')
print(f'f(x*)         : {fopt:.10f}')
print(f'True minimum  : {x_true:.10f}')
print(f'|error|       : {abs(xopt - x_true):.3e}')
print(f'Iterations    : {n_iter}')

# Iteration table
print(f'\n{"Iter":>5}  {"xopt":>14}  {"f(xopt)":>14}  {"ea":>12}  {"interval":>10}')
print('-' * 62)
for row in history[::3]:   # print every 3rd row to keep table short
    print(f'{row[0]:>5}  {row[1]:>14.8f}  {row[2]:>14.8f}  {row[3]:>12.6e}  {row[4]:>10.6f}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
Minimum found : x* = 1.4275517674
f(x*)         : -1.7757256531
True minimum  : 1.4275517797
|error|       : 1.229e-08
Iterations    : 45

 Iter            xopt         f(xopt)            ea    interval
--------------------------------------------------------------
    1      1.52786405     -1.76472025  1.000000e+00    2.472136
    4      1.52786405     -1.76472025  2.360680e-01    0.583592
    7      1.44271910     -1.77547480  5.901699e-02    0.137767
   10      1.42261911     -1.77569914  1.412887e-02    0.032522
   13      1.42736407     -1.77572561  3.324285e-03    0.007678
   16      1.42736407     -1.77572561  7.847573e-04    0.001812
   19      1.42752750     -1.77572565  1.852349e-04    0.000428
   22      1.42756608     -1.77572565  4.372684e-05    0.000101
   25      1.42755134     -1.77572565  1.032261e-05    0.000024
   28      1.42755134     -1.77572565  2.436838e-06    0.000006
   31      1.42755216     -1.77572565  5.752591e-07    0.000001
   34      1.42755185     -1.77572565  1.358003e-07    0.000000
   37      1.42755177     -1.77572565  3.205810e-08    0.000000
   40      1.42755177     -1.77572565  7.567891e-09    0.000000
   43      1.42755177     -1.77572565  1.786537e-09    0.000000
```

## 문제 2. 포물선 보간 (노트북 3.1~3.3절)

포물선 보간은 구간의 끝점만 쓰는 대신 세 점 $`(x_1, f_1)`$, $`(x_2, f_2)`$, $`(x_3, f_3)`$을 지나는 **이차 다항식**을 맞추고, 그 포물선의 꼭짓점으로 바로 건너뛴다.

```math
x_4 = x_2 - \frac{1}{2}\cdot
\frac{(x_2 - x_1)^2[f_2 - f_3] - (x_2 - x_3)^2[f_2 - f_1]}
{(x_2 - x_1)[f_2 - f_3] - (x_2 - x_3)[f_2 - f_1]}
```

최솟값 근처에서 $`f`$가 포물선에 가까우면 황금 분할 탐색보다 훨씬 빠를 수 있다. $`x_4`$를 구한 뒤에는 현재 세 점 가운데 가장 나쁜 점을 바꾸고 반복하며, $`|x_4^{\rm new} - x_4^{\rm old}|`$가 충분히 작으면 멈춘다. 매끄러운 최솟값 근처에서는 차수 약 1.32의 초선형 수렴을 보이지만, 세 점이 거의 한 직선 위에 있으면(분모 $`\approx 0`$) 실패할 수 있다. 다음 함수의 빈칸을 채워라.

```python
def parabolic(f, x1, x2, x3, es=1e-7, maxit=50):
    """
    Minimize f using Parabolic Interpolation.

    Requires three initial points x1 < x2 < x3 with f(x2) < f(x1) and f(x2) < f(x3).

    Parameters
    ----------
    f       : callable
    x1, x2, x3 : float — initial three points  (x1 < x2 < x3)
    es      : float   — stopping criterion
    maxit   : int

    Returns
    -------
    x4, f(x4), ea, n_iter, history
    """
    history = []
    x4_old = x2

    for i in range(maxit):
        y1, y2, y3 = f(x1), f(x2), f(x3)

        num = (x2 - x1)**2 * (y2 - y3) - (x2 - x3)**2 * (y2 - y1)
        den = (x2 - x1) * (y2 - y3)    - (x2 - x3) * (y2 - y1)

        if abs(den) < 1e-30:
            break   # denominator near zero — cannot proceed

        # --- COMPLETE: parabolic interpolation formula ---
        x4 = 
        # --- COMPLETE ---

        ea = abs((x4 - x4_old) / x4) if x4 != 0 else abs(x4 - x4_old)
        history.append((i + 1, x4, f(x4), ea))

        if ea < es:
            break

        # Update: replace the worst of the three points
        pts = sorted([(f(p), p) for p in [x1, x2, x3, x4]], key=lambda t: t[0])
        x1, x2, x3 = sorted([pts[0][1], pts[1][1], pts[2][1]])

        x4_old = x4

    return x4, f(x4), ea, i + 1, history
```

다음 셀은 황금 분할 탐색과 포물선 보간을 비교한다. 그래프를 그리는 부분은 생략했다.

```python
# Compare Golden-Section vs Parabolic Interpolation
xopt_gs, fopt_gs, _, n_gs, hist_gs = goldmin(f1d, 0.0, 4.0, es=1e-10, maxit=100)
xopt_pi, fopt_pi, _, n_pi, hist_pi = parabolic(f1d, 0.0, 1.0, 4.0, es=1e-10, maxit=50)

print(f'Golden-Section:         {n_gs:3d} iterations,  final |error| = {abs(xopt_gs - x_true):.3e}')
print(f'Parabolic Interpolation:{n_pi:3d} iterations,  final |error| = {abs(xopt_pi - x_true):.3e}')
print('\nNote: Parabolic converges faster when f is smooth near the minimum.')
print('\nUse: curvature information of a function')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
Golden-Section:          49 iterations,  final |error| = 1.229e-08
Parabolic Interpolation:  9 iterations,  final |error| = 2.222e-09

Note: Parabolic converges faster when f is smooth near the minimum.

Use: curvature information of a function
```

## 문제 3. 최급 하강법 (노트북 4.1~4.4절)

$`\mathbf{x} \in \mathbb{R}^n`$인 함수 $`f(\mathbf{x})`$에서 **기울기** $`\nabla f`$는 가장 가파르게 증가하는 방향을 가리킨다. 따라서 **음의 기울기 방향**으로 움직이면 $`f`$가 가장 빠르게 줄어든다. 알고리즘은 $`\mathbf{x}^{(0)}`$에서 시작하여 $`k = 0, 1, 2, \ldots`$에 대해 다음을 반복한다.

1. 기울기 $`\mathbf{g}^{(k)} = \nabla f(\mathbf{x}^{(k)})`$를 계산한다.
2. 스텝 크기 $`\alpha^{(k)}`$를 고른다(선 탐색).
3. $`\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} - \alpha^{(k)} \mathbf{g}^{(k)}`$로 갱신한다.
4. $`\|\nabla f(\mathbf{x}^{(k)})\| < \varepsilon`$이거나 $`\|\mathbf{x}^{(k+1)} - \mathbf{x}^{(k)}\| < \delta`$이면 멈춘다.

스텝 크기가 너무 크면 $`f`$가 커지고(최솟값을 지나친다), 너무 작으면 수렴이 느리다. 백트래킹(Armijo) 규칙은 $`\alpha = 1`$에서 시작해 $`f(\mathbf{x} - \alpha\mathbf{g}) \geq f(\mathbf{x})`$인 동안 $`\alpha \leftarrow c\,\alpha`$ $`(c < 1)`$로 줄인다. 최급 하강법은 선형 수렴하며, 강볼록 함수에서는 $`\|\mathbf{x}^{(k)} - \mathbf{x}^*\| \leq \left(\frac{\kappa - 1}{\kappa + 1}\right)^k \|\mathbf{x}^{(0)} - \mathbf{x}^*\|`$이다. 여기서 $`\kappa = \lambda_{\max}/\lambda_{\min}`$은 헤세 행렬의 **조건수**이며, 조건수가 크면 지그재그로 느리게 수렴한다.

노트북은 2차원 예제 $`f(x, y) = 2 + x - y + 2x^2 + 2xy + y^2`$를 쓴다. 해석적 최솟값은 $`(-1, 1.5)`$이다. 다음 셀이 함수와 기울기를 정의한다. 그래프를 그리는 부분은 생략했다.

```python
# Running 2D example: f(x, y) = 2 + x - y + 2x² + 2xy + y²
# Analytical minimum: ∂f/∂x = 1+4x+2y=0, ∂f/∂y = -1+2x+2y=0 → x=-1, y=3/2

def f2d(x):
    return 2.0 + x[0] - x[1] + 2*x[0]**2 + 2*x[0]*x[1] + x[1]**2

def grad_f2d(x):
    return np.array([1.0 + 4*x[0] + 2*x[1],
                    -1.0 + 2*x[0] + 2*x[1]])

x2_true = np.array([-1.0, 1.5])   # exact minimum
```

다음 셀의 빈칸에 최급 하강법의 갱신 식을 채워라.

```python
def backtrack_linesearch(f, x, grad, c=0.5, alpha0=1.0, max_ls=50):
    """
    Backtracking line search: find alpha so that f(x - alpha*grad) < f(x).

    Starts from alpha0 and multiplies by c until descent condition is satisfied.
    """
    alpha = alpha0
    fx = f(x)
    for _ in range(max_ls):
        if f(x - alpha * grad) < fx:
            break
        alpha *= c
    return alpha


def steepest_descent(f, grad_f, x0, tol=1e-6, maxit=1000, record=True):
    """
    Minimize f using steepest descent with backtracking line search.

    Parameters
    ----------
    f      : callable — objective function, accepts 1-D array
    grad_f : callable — gradient function, returns 1-D array
    x0     : array-like — initial guess
    tol    : float    — stopping criterion: ||grad f|| < tol
    maxit  : int      — maximum iterations
    record : bool     — record iteration history

    Returns
    -------
    x       : array  — minimizer estimate
    f(x)    : float
    n_iter  : int
    history : list of (iter, x.copy(), f(x), ||grad||, alpha)
    """
    x = np.array(x0, dtype=float)
    history = []

    for i in range(maxit):
        g = grad_f(x)
        grad_norm = np.linalg.norm(g)

        if record:
            history.append((i, x.copy(), f(x), grad_norm, None))

        if grad_norm < tol:
            break

        # --- COMPLETE: gradient descent update ---
        alpha = backtrack_linesearch(f, x, g)   # backtracking line search
        x =                                    # — steepest descent step
        # --- COMPLETE ---

        if record and history:
            history[-1] = (history[-1][0], history[-1][1],
                           history[-1][2], history[-1][3], alpha)

    return x, f(x), i + 1, history
```

다음 셀은 완성한 함수를 2차원 예제에 적용한다. 그래프를 그리는 부분은 생략했다.

```python
# Test steepest descent on the 2D example
x0 = np.array([-2.0, 3.0])
xopt2d, fopt2d, n_sd, hist_sd = steepest_descent(f2d, grad_f2d, x0, tol=1e-8, maxit=200)

print(f'Starting point    : {x0}')
print(f'Minimizer found   : {np.round(xopt2d, 8)}')
print(f'True minimum      : {x2_true}')
print(f'f(minimizer)      : {fopt2d:.10f}')
print(f'f(true min)       : {f2d(x2_true):.10f}')
print(f'|error|           : {np.linalg.norm(xopt2d - x2_true):.3e}')
print(f'Iterations        : {n_sd}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다. 출력 뒤에는 하강 경로와 기울기 노름의 그래프가 있다.

```text
Starting point    : [-2.  3.]
Minimizer found   : [-1.00000001  1.50000001]
True minimum      : [-1.   1.5]
f(minimizer)      : 0.7500000000
f(true min)       : 0.7500000000
|error|           : 2.107e-08
Iterations        : 200
```

## 노트북의 추가 연습 문제

노트북에는 빈칸 외에 코드 셀이 비어 있는 Exercise가 두 개 있다.

**연습 A (노트북 2절 끝).** 황금 분할 탐색을 완성한 뒤 다음 연습이 이어진다. `goldmin`으로 구간 $`[5, 8]`$에서 $`f(x) = x^2/10 - 2\sin(x)`$의 **국소** 최솟값을 찾아라. 그 값은 앞에서 찾은 전역 최솟값과 비교해 어떠한가?

```python
# --- YOUR CODE HERE ---
# Apply goldmin to find the local minimum in [5, 8]
```

**연습 B (노트북 4절 끝).** 최급 하강법을 시험한 뒤 다음 연습이 이어진다. 백트래킹 대신 **고정된** 스텝 크기(예: $`\alpha = 0.01, 0.1, 0.5`$)로 최급 하강법을 실행하고, $`\alpha`$의 선택이 수렴에 어떤 영향을 주는지 관찰하라. 힌트: `backtrack_linesearch(...)`를 상수 `alpha`로 바꾼다.

```python
# --- YOUR CODE HERE ---
# Try fixed step sizes alpha = 0.01, 0.1, 0.5 on f2d
# Plot convergence curves for each
```

---

<br>

문제를 풀어 본 뒤 [과제 해설](../Assignment-Explanation.ko.md)에서 정답과 풀이를, [완성한 노트북](../Solution/Lab08-Optimization_Completed.ipynb)에서 작성한 코드 전체를 볼 수 있다.
