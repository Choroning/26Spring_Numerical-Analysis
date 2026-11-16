# HW6 과제 안내문 — 근 찾기 방법

[레포지토리 안내](../../../README.ko.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **최종 수정일:** 2026-11-16

수치해석 7주차 과제이다. 교수자가 배포한 실습 노트북 `Lab07-RootSolving.ipynb`의 빈칸을 채워 실행한다. 노트북은 동기와 그래프 방법, 증분 탐색, 이분법, 가위치법, 뉴턴-랩슨법, 할선법, `scipy.optimize`와 방법 비교의 일곱 부분으로 이루어지며, 빈칸은 `COMPLETE`라고 적힌 두 주석 줄 사이에 모두 5곳 있다. 이분법과 뉴턴-랩슨법을 적금 문제에 적용하는 두 Exercise 셀은 코드가 채워진 채 배포되었다. 노트북의 첫 코드 셀은 `numpy`를 `np`로, `matplotlib.pyplot`을 `plt`로 불러오고, `scipy.optimize`의 `brentq`, `root_scalar`와 `functools`의 `partial`을 불러온다.

아래 각 문제의 코드 블록은 노트북의 셀을 옮긴 것이며, 코드가 비어 있는 자리가 채울 곳이다. 배포된 빈 노트북은 [Lab07-RootSolving.ipynb](Lab07-RootSolving.ipynb)에 있다.

- **제출 기한:** 배포 후 1주
- **제출 방법:** 빈칸을 채워 실행한 노트북을 온라인으로 제출한다.

---

<br>

## 노트북의 예제 문제 (노트북 1절)

함수 $`f`$의 **근(root)** 은 $`f(x^*) = 0`$을 만족하는 값 $`x^*`$이다. 노트북은 두 공학 문제를 예제로 쓴다.

- **적금:** Shin은 매달 $`p`$달러를 연이율 $`r`$로 적립하며, $`n`$번 납입한 뒤의 총액은 $`A = \frac{12p}{r}\left[\left(1+\frac{r}{12}\right)^n - 1\right]`$이다. 매달 60달러씩 $`n = 24`$개월 납입한 뒤 $`A = 1500`$달러가 되는 $`r`$을 찾는 것이 목표이며, $`f(r) = \frac{12p}{r}\left[(1+r/12)^n - 1\right] - A = 0`$으로 정의한다.
- **번지 점프 선수:** 질량 $`m`$인 번지 점프 선수의 시각 $`t`$에서의 낙하 속도는 $`v(m) = \sqrt{mg/c_d} \tanh\left(\sqrt{c_d g/m}\, t\right)`$이다. 항력 계수가 $`c_d = 0.25`$ kg/m일 때 $`t = 4`$ s에서 $`v = 36`$ m/s가 되는 질량 $`m`$을 찾는 것이 목표이며, $`f(m) = \sqrt{mg/c_d} \tanh\left(\sqrt{c_d g/m}\, t\right) - v = 0`$으로 정의한다.

다음 셀이 두 함수를 정의한다. 그래프를 그리는 부분은 생략했다.

```python
# Problem parameters
cd, g_acc, v_target, t_fall = 0.25, 9.81, 36.0, 4.0
deposit, duration, target_amount = 60.0, 24, 1500.0

# Define the two root-finding functions
def f_bungee(m):
    """f(m) = v(m) - v_target  (bungee jumper mass problem)"""
    return np.sqrt(m * g_acc / cd) * np.tanh(np.sqrt(cd * g_acc / m) * t_fall) - v_target

def f_savings(r):
    """f(r) = A(r) - target  (savings account interest rate)"""
    ratio = 1.0 + r / 12.0
    return 12.0 * deposit / r * (ratio**duration - 1.0) - target_amount
```

그래프에서 번지 점프 선수의 근은 $`[50, 250]`$ 안에 있으며, $`f(50) < 0`$이고 $`f(250) > 0`$이다. $`f`$가 연속이고 $`f(x_l) \cdot f(x_u) < 0`$이면 $`[x_l, x_u]`$에 근이 적어도 하나 있다(중간값 정리).

## 문제 1. 증분 탐색 (노트북 2절)

근 찾기 방법을 적용하기 전에 근을 감싸는 구간을 찾아야 한다. 증분 탐색은 $`[x_{\min}, x_{\max}]`$를 $`n_s - 1`$개의 같은 소구간으로 나누고 각 소구간에서 부호 변화를 확인한다.

1. $`[x_{\min}, x_{\max}]`$에 같은 간격의 점 $`n_s`$개 $`x_0, x_1, \ldots, x_{n_s-1}`$을 만든다.
2. 각 쌍 $`(x_k, x_{k+1})`$에서 $`f(x_k) \cdot f(x_{k+1}) < 0`$이면 그 구간을 기록한다.
3. 찾은 구간을 모두 돌려준다.

$`n_s`$가 너무 작으면 한 구간에 근이 짝수 개 들어가 놓칠 수 있다. 다음 함수의 빈칸에 부호 변화를 판정하는 조건을 채워라.

```python
def incsearch(func, xmin, xmax, ns=50):
    """
    Find all sign-change brackets of func in [xmin, xmax].

    Parameters
    ----------
    func : callable
    xmin, xmax : float  — search interval
    ns   : int          — number of points (ns-1 subintervals)

    Returns
    -------
    nb   : int          — number of brackets found
    xb   : list of (xl, xu) tuples
    """
    x = np.linspace(xmin, xmax, ns)
    nb = 0
    xb = []

    for k in range(ns - 1):
        # --- COMPLETE: detect sign change between x[k] and x[k+1] ---
        if :
            nb += 1
            xb.append((x[k], x[k+1]))
        # --- COMPLETE ---

    return nb, xb
```

다음 셀은 근이 여러 개인 함수로 완성한 함수를 시험한다. 그래프를 그리는 부분은 생략했다.

```python
# Test incsearch on a function with multiple roots
f_multi = lambda x: np.sin(10 * x) + np.cos(3 * x)
xmin, xmax = 3.0, 6.0

for ns in [20, 50, 100]:
    nb, xb = incsearch(f_multi, xmin, xmax, ns=ns)
    print(f'ns={ns:3d}  ->  {nb} bracket(s) found')

# Best result
nb, xb = incsearch(f_multi, xmin, xmax, ns=100)
print(f'\nBrackets (ns=100):')
for xl_b, xu_b in xb:
    print(f'  [{xl_b:.4f}, {xu_b:.4f}]')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
ns= 20  ->  7 bracket(s) found
ns= 50  ->  5 bracket(s) found
ns=100  ->  9 bracket(s) found

Brackets (ns=100):
  [3.2424, 3.2727]
  [3.3636, 3.3939]
  [3.7273, 3.7576]
  [4.2121, 4.2424]
  [4.2424, 4.2727]
  [4.6970, 4.7273]
  [5.1515, 5.1818]
  [5.1818, 5.2121]
  [5.6667, 5.6970]
```

## 문제 2. 이분법 (노트북 3.1절)

$`f(x_l) \cdot f(x_u) < 0`$인 구간 $`[x_l, x_u]`$가 주어지면 이분법은 다음을 반복한다.

1. 중점 $`x_r = \dfrac{x_l + x_u}{2}`$를 계산한다.
2. 근사 상대 오차 $`\varepsilon_a = \left|\dfrac{x_r^{\rm new} - x_r^{\rm old}}{x_r^{\rm new}}\right| \times 100\%`$를 계산한다.
3. 구간을 갱신한다. $`f(x_l) \cdot f(x_r) < 0`$이면 근이 $`[x_l, x_r]`$에 있으므로 $`x_u = x_r`$로 두고, 그렇지 않으면 근이 $`[x_r, x_u]`$에 있으므로 $`x_l = x_r`$로 둔다.
4. $`\varepsilon_a < \varepsilon_s`$가 되거나 최대 반복 횟수에 이를 때까지 반복한다.

구간의 폭은 단계마다 절반이 되므로 $`n`$번 반복한 뒤의 오차는 $`(x_u - x_l)/2^n`$ 이하이다. 이는 비율 0.5의 선형 수렴이며, 절대 허용 오차 $`\varepsilon_{\rm abs}`$를 얻는 데 필요한 최소 반복 횟수는 $`n \geq \left\lceil \log_2\!\left(\frac{x_u - x_l}{\varepsilon_{\rm abs}}\right) \right\rceil`$이다. 처음 구간이 올바르면 항상 수렴하지만 수렴이 느리다. 다음 함수의 빈칸을 채워라.

```python
def bisect(func, xl, xu, es=1e-7, maxit=50):
    """
    Find a root of func in [xl, xu] using the bisection method.

    Parameters
    ----------
    func : callable — continuous function
    xl, xu : float  — initial bracket  (must satisfy f(xl)*f(xu) < 0)
    es     : float  — stopping criterion (approx. relative error)
    maxit  : int    — maximum iterations

    Returns
    -------
    xr      : float — root estimate
    fxr     : float — f(xr)
    ea      : float — final approximate relative error
    n_iter  : int   — iterations performed
    history : list of (iter, xr, f(xr), ea)
    """
    if func(xl) * func(xu) > 0:
        raise ValueError('f(xl) and f(xu) must have opposite signs (bracket check failed)')

    xr_old = xl       # previous midpoint (initialized to xl)
    history = []

    for i in range(maxit):
        # --- COMPLETE ---
        xr =                                # — midpoint
        ea = abs((xr - xr_old) / xr) if xr != 0 else abs(xr - xr_old)  # approx. relative error

        if func(xl) * func(xr) < 0:     #  — sign test: root in left half?
            xu = xr
        else:
            xl =                             #      — root in right half
        # --- COMPLETE ---

        history.append((i + 1, xr, func(xr), ea))
        xr_old = xr

        if ea < es:
            break

    return xr, func(xr), ea, i + 1, history
```

다음 셀은 완성한 함수를 번지 점프 선수 문제에 적용한다.

```python
# Apply bisection to the bungee jumper problem
m_true = 142.7376  # known true root (for error analysis)

xr, fxr, ea, n_iter, history = bisect(f_bungee, 50, 250, es=1e-9, maxit=60)

print(f'Root found : m = {xr:.8f} kg')
print(f'f(root)    : {fxr:.3e}')
print(f'Rel. error : {ea:.3e}')
print(f'Iterations : {n_iter}')

# Print iteration table
print(f'\n{"Iter":>5}  {"xr":>14}  {"f(xr)":>14}  {"ea":>12}')
print('-' * 52)
for row in history:
    print(f'{row[0]:>5}  {row[1]:>14.8f}  {row[2]:>14.6e}  {row[3]:>12.6e}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
Root found : m = 142.73763308 kg
f(root)    : -5.203e-10
Rel. error : 6.525e-10
Iterations : 31

 Iter              xr           f(xr)            ea
----------------------------------------------------
    1    150.00000000    1.420430e-01  6.666667e-01
    2    100.00000000   -1.197380e+00  5.000000e-01
    3    125.00000000   -4.086015e-01  2.000000e-01
    4    137.50000000   -1.108041e-01  9.090909e-02
    5    143.75000000    2.057713e-02  4.347826e-02
    6    140.62500000   -4.379968e-02  2.222222e-02
    7    142.18750000   -1.129258e-02  1.098901e-02
    8    142.96875000    4.720768e-03  5.464481e-03
    9    142.57812500   -3.266139e-03  2.739726e-03
   10    142.77343750    7.322383e-04  1.367989e-03
   11    142.67578125   -1.265717e-03  6.844627e-04
   12    142.72460938   -2.664315e-04  3.421143e-04
   13    142.74902344    2.329804e-04  1.710279e-04
   14    142.73681641   -1.670630e-05  8.552125e-05
   15    142.74291992    1.081418e-04  4.275880e-05
   16    142.73986816    4.571898e-05  2.137986e-05
   17    142.73834229    1.450664e-05  1.069004e-05
   18    142.73757935   -1.099755e-06  5.345050e-06
   19    142.73796082    6.703461e-06  2.672518e-06
   20    142.73777008    2.801858e-06  1.336261e-06
   21    142.73767471    8.510526e-07  6.681308e-07
   22    142.73762703   -1.243508e-07  3.340655e-07
   23    142.73765087    3.633510e-07  1.670327e-07
   24    142.73763895    1.195001e-07  8.351637e-08
   25    142.73763299   -2.425352e-09  4.175818e-08
   26    142.73763597    5.853737e-08  2.087909e-08
   27    142.73763448    2.805601e-08  1.043955e-08
   28    142.73763373    1.281533e-08  5.219773e-09
   29    142.73763336    5.194991e-09  2.609887e-09
   30    142.73763318    1.384819e-09  1.304943e-09
   31    142.73763308   -5.202665e-10  6.524716e-10
```

## 문제 3. 가위치법 (노트북 4.1절)

가위치법은 구간을 **중점**에서 나누는 대신 $`(x_l, f(x_l))`$과 $`(x_u, f(x_u))`$를 잇는 **직선**의 x 절편을 찾는다.

$$
x_r = x_u - f(x_u)\,\frac{x_l - x_u}{f(x_l) - f(x_u)}
$$

이는 선형 보간과 같아서, $`[x_l, x_u]`$에서 $`f`$가 거의 직선이면 중점보다 좋은 추정을 준다. 구간의 갱신 규칙(부호 판정)은 이분법과 같다. 그림으로 나타내면 다음과 같다.

```text
  f(xl) *
          \                 <- True f(x)
           \    ___---
  0 --------\--/--------  x-axis
              \/ <- xr (line intercept)
          *  f(xu)
```

$`f`$가 크게 휘어 있으면 구간의 한 끝점이 여러 반복 동안 고정되는 일이 많으며, 이 때문에 가위치법이 이분법보다 느려지는 경우가 있다. 다음 함수의 빈칸을 채워라.

```python
def false_position(func, xl, xu, es=1e-7, maxit=50):
    """
    Find a root of func in [xl, xu] using the False Position (Regula Falsi) method.

    Parameters
    ----------
    func : callable
    xl, xu : float  — initial bracket
    es     : float  — stopping criterion (approx. relative error)
    maxit  : int    — maximum iterations

    Returns
    -------
    xr, f(xr), ea, n_iter, history
    """
    if func(xl) * func(xu) > 0:
        raise ValueError('f(xl) and f(xu) must have opposite signs')

    xr_old = xl
    history = []

    for i in range(maxit):
        fl, fu = func(xl), func(xu)

        # --- COMPLETE: linear interpolation (secant line x-intercept) ---
        xr = 
        # --- COMPLETE ---

        fr = func(xr)
        ea = abs((xr - xr_old) / xr) if xr != 0 else abs(xr - xr_old)

        history.append((i + 1, xr, fr, ea))
        xr_old = xr

        if ea < es:
            break

        if fl * fr < 0:
            xu = xr
        else:
            xl = xr

    return xr, func(xr), ea, i + 1, history
```

다음 셀은 번지 점프 선수 문제에서 이분법과 가위치법을 비교한다. 그래프를 그리는 부분은 생략했다.

```python
# Compare bisection vs. false position on the bungee jumper problem
xr_bs, _, _, _, hist_bs = bisect(f_bungee, 50, 250, es=1e-10, maxit=60)
xr_fp, _, _, _, hist_fp = false_position(f_bungee, 50, 250, es=1e-10, maxit=60)

err_bs = [abs(row[1] - m_true) for row in hist_bs]
err_fp = [abs(row[1] - m_true) for row in hist_fp]

print(f'Bisection:      {len(hist_bs):2d} iterations,  final |error| = {err_bs[-1]:.3e}')
print(f'False Position: {len(hist_fp):2d} iterations,  final |error| = {err_fp[-1]:.3e}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
Bisection:      34 iterations,  final |error| = 3.312e-05
False Position: 42 iterations,  final |error| = 3.313e-05
```

## 문제 4. 뉴턴-랩슨법 (노트북 5.1절)

뉴턴-랩슨법은 도함수 $`f'(x)`$로 현재 추정값에서의 접선을 만들고, 그 x 절편을 다음 추정값으로 삼는다.

$$
x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}
$$

$`x_i`$에서 접선 $`y = f(x_i) + f'(x_i)(x - x_i)`$를 그리고 $`y = 0`$으로 놓아 $`x`$를 풀면 $`x_{i+1}`$이 된다. 근 $`x^*`$ 근처에서 오차 $`e_i = x_i - x^*`$는 $`e_{i+1} \approx -\frac{f''(x^*)}{2f'(x^*)} e_i^2`$을 만족하는 이차 수렴이며, 맞는 소수 자릿수가 반복마다 대략 두 배가 된다. 다만 $`f'(x)`$를 해석적으로 알아야 하고, $`f'(x_i) \approx 0`$이거나 시작값이 나쁘면 발산할 수 있다.

| 실패 유형 | 원인 | 대책 |
|:-----|:-----|:-----|
| 발산 | 시작값이 근에서 너무 멀다 | 먼저 구간법을 쓴다 |
| 진동 | 근 근처에 극값이 있다 | 작은 스텝을 쓰거나 이분법으로 바꾼다 |
| 도함수 0 | $`f'(x_i) \approx 0`$ | 작은 섭동을 더하거나 방법을 바꾼다 |

*표 1. 뉴턴-랩슨법이 실패하는 경우*

다음 함수의 빈칸을 채워라.

```python
def newtraph(f, fp, x0, es=1e-7, maxit=30):
    """
    Find a root of f using the Newton-Raphson method.

    Parameters
    ----------
    f    : callable — function
    fp   : callable — derivative of f
    x0   : float   — initial guess
    es   : float   — stopping criterion (approx. relative error)
    maxit: int     — maximum iterations

    Returns
    -------
    x1, f(x1), ea, n_iter, history
    """
    history = []
    xi = x0

    for i in range(maxit):
        fxi  = f(xi)
        fpxi = fp(xi)

        if abs(fpxi) < 1e-30:
            raise ZeroDivisionError('Derivative near zero — Newton-Raphson fails')

        # --- COMPLETE ---
        x1 =                                # — Newton-Raphson update
        ea = abs((x1 - xi) / x1) if x1 != 0 else abs(x1 - xi)
        # --- COMPLETE ---

        history.append((i + 1, x1, f(x1), ea))

        if ea < es:
            break

        xi = x1

    return x1, f(x1), ea, i + 1, history
```

앞 셀은 시험용 함수를 `f_nr = lambda x: x - np.exp(-x)`(근 약 0.5671), `fp_nr = lambda x: 1.0 + np.exp(-x)`, `x_true = 0.5671432904097838`로 정의해 둔다. 다음 셀은 완성한 함수를 이 함수에 적용한다.

```python
# Test newtraph on f(x) = x - e^{-x}, root = 0.5671...
xr_nr, fxr_nr, ea_nr, n_nr, hist_nr = newtraph(f_nr, fp_nr, x0=0.0, es=1e-12, maxit=20)

print(f'Root found : {xr_nr:.15f}')
print(f'True root  : {x_true:.15f}')
print(f'|error|    : {abs(xr_nr - x_true):.3e}')
print(f'Iterations : {n_nr}')

print(f'\n{"Iter":>5}  {"xi":>18}  {"f(xi)":>14}  {"ea":>14}')
print('-' * 58)
for row in hist_nr:
    print(f'{row[0]:>5}  {row[1]:>18.12f}  {row[2]:>14.6e}  {row[3]:>14.6e}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
Root found : 0.567143290409784
True root  : 0.567143290409784
|error|    : 1.110e-16
Iterations : 5

 Iter                  xi           f(xi)              ea
----------------------------------------------------------
    1      0.500000000000   -1.065307e-01    1.000000e+00
    2      0.566311003197   -1.304510e-03    1.170929e-01
    3      0.567143165035   -1.964805e-07    1.467287e-03
    4      0.567143290410   -4.440892e-15    2.210639e-07
    5      0.567143290410    1.110223e-16    5.089684e-15
```

## 문제 5. 할선법 (노트북 6.1절)

$`f'(x)`$를 해석적으로 구하기 어렵거나 계산 비용이 클 때가 있다. 할선법은 앞선 두 반복값으로 도함수를 근사한다.

$$
f'(x_i) \approx \frac{f(x_i) - f(x_{i-1})}{x_i - x_{i-1}}
$$

이를 뉴턴-랩슨 공식에 대입하면 다음을 얻는다.

```math
x_{i+1} = x_i - f(x_i)\,\frac{x_{i-1} - x_i}{f(x_{i-1}) - f(x_i)}
```

할선법은 황금비 $`\phi = \frac{1+\sqrt{5}}{2} \approx 1.618`$을 차수로 하여 수렴한다. 초선형 수렴이지만 뉴턴-랩슨법의 이차 수렴보다는 느리다.

| 방법 | 차수 | 단계마다 계산하는 함수 | 실질 속도 |
|:-----|:-----|:-----|:-----|
| 뉴턴-랩슨법 | 2 (이차) | $`f`$와 $`f'`$ | 빠름 |
| 할선법 | 1.618 | $`f`$만 | 중간 |
| 이분법 | 1 (선형) | $`f`$만 | 느림 |

*표 2. 뉴턴-랩슨법, 할선법, 이분법의 수렴 비교*

다음 함수의 빈칸을 채워라.

```python
def secant(f, x0, x1, es=1e-7, maxit=30):
    """
    Find a root of f using the Secant method.

    Parameters
    ----------
    f     : callable — function (no derivative needed)
    x0, x1: float   — two initial guesses (do NOT need to bracket the root)
    es    : float   — stopping criterion
    maxit : int     — maximum iterations

    Returns
    -------
    x2, f(x2), ea, n_iter, history
    """
    history = []

    for i in range(maxit):
        f0, f1 = f(x0), f(x1)

        if abs(f1 - f0) < 1e-30:
            raise ZeroDivisionError('f(x1) ≈ f(x0) — secant method fails (zero denominator)')

        # --- COMPLETE: secant update formula ---
        x2 = 
        # --- COMPLETE ---

        ea = abs((x2 - x1) / x2) if x2 != 0 else abs(x2 - x1)
        history.append((i + 1, x2, f(x2), ea))

        if ea < es:
            break

        x0, x1 = x1, x2

    return x2, f(x2), ea, i + 1, history
```

다음 셀은 같은 시험용 함수에서 뉴턴-랩슨법과 할선법을 비교한다. 그래프를 그리는 부분은 생략했다.

```python
# Compare Newton-Raphson vs Secant on f(x) = x - e^{-x}
xr_sc, fxr_sc, ea_sc, n_sc, hist_sc = secant(f_nr, x0=0.0, x1=1.0, es=1e-12, maxit=30)

print('=== f(x) = x - e^{-x},  true root = 0.5671432904097838 ===')
print(f'Newton-Raphson:  {n_nr} iters,  |error| = {abs(xr_nr - x_true):.3e}')
print(f'Secant:          {n_sc} iters,  |error| = {abs(xr_sc - x_true):.3e}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
=== f(x) = x - e^{-x},  true root = 0.5671432904097838 ===
Newton-Raphson:  5 iters,  |error| = 1.110e-16
Secant:          6 iters,  |error| = 0.000e+00
```

---

<br>

문제를 풀어 본 뒤 [과제 해설](../Assignment-Explanation.ko.md)에서 정답과 풀이를, [완성한 노트북](../Solution/Lab07-RootSolving_Completed.ipynb)에서 작성한 코드 전체를 볼 수 있다.
