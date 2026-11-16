# HW1 과제 안내문 — 수치 오차의 원인

[레포지토리 안내](../../../README.ko.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **최종 수정일:** 2026-11-16

수치해석 2주차 과제이다. 교수자가 배포한 실습 노트북 `Lab02-SourceOfErrors.ipynb`의 빈칸을 채워 실행한다. 노트북은 오차의 정의와 지표, 절단 오차와 테일러 급수, 반올림 오차와 부동 소수점 연산, 유한 차분과 오차 분석의 네 부분으로 이루어지며, 빈칸은 `COMPLETE` 또는 `Complete`라고 적힌 두 주석 줄 사이에 모두 5곳 있다. 노트북의 첫 코드 셀은 `numpy`를 `np`로, `matplotlib.pyplot`을 `plt`로 불러오고 `math`를 불러온다.

아래 각 문제의 코드 블록은 노트북의 셀을 옮긴 것이며, 코드가 비어 있는 자리가 채울 곳이다. 배포된 빈 노트북은 [Lab02-SourceOfErrors.ipynb](Lab02-SourceOfErrors.ipynb)에 있다.

- **제출 기한:** 배포 후 1주
- **제출 방법:** 빈칸을 채워 실행한 노트북을 온라인으로 제출한다.

---

<br>

## 문제 1. 지수 함수의 매클로린 부분합 (노트북 2.2절)

지수 함수는 다음과 같이 간결한 매클로린 전개를 갖는다.

$$
e^x = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \cdots = \sum_{k=0}^{\infty} \frac{x^k}{k!}
$$

항을 하나씩 더해 $`e^{0.5}`$를 근사하고 오차가 줄어드는 모습을 관찰한다. 앞 셀은 참값을 다음과 같이 정한다.

```python
# True value
x = 0.5
f_true = np.exp(x)
print(f'True value e^{x} = {f_true:.15f}')
```

다음 셀의 빈칸을 채워 남은 부분합 `f_approx[3]`, `f_approx[4]`, `f_approx[5]`를 구하라.

```python
# Build approximations term by term  (k=0,1,2,...,5 → 6 terms)
# f_approx[k] = sum_{i=0}^{k} x^i / i!

# --- COMPLETE: fill in the missing partial sums ---
f_approx = np.zeros(6)
f_approx[0] = 1                               # k=0: just the constant term
f_approx[1] = f_approx[0] + x                # k=1: + x
f_approx[2] = f_approx[1] + x**2 / math.factorial(2)   # k=2: + x^2/2!
f_approx[3] =                                # k=3: + x^3/3!
f_approx[4] =                                # k=4: + x^4/4!
f_approx[5] =                                # k=5: + x^5/5!
# --- COMPLETE ---

print(f'{"Terms":>6}  {"Approximation":>20}  {"True Rel Err (%)":>18}')
print('-' * 50)
for k in range(6):
    et = abs(f_true - f_approx[k]) / f_true * 100
    print(f'{k+1:>6}  {f_approx[k]:>20.15f}  {et:>18.10f}')
```

## 문제 2. 중지 기준을 쓴 급수 계산 (노트북 2.3절)

항의 개수를 미리 정하는 대신, 실용적인 알고리즘은 $`\varepsilon_a < \varepsilon_s`$가 될 때까지 항을 계속 더한다. 노트북 1절은 근사 백분율 상대 오차를 다음과 같이 정의한다.

```math
\varepsilon_a = \frac{| x_{\text{new}} - x_{\text{prev}} |}{|x_{\text{new}}|} \times 100\%
```

참값을 모르는 실제 계산에서는 스카보로 기준(Scarborough criterion)이 결과의 유효 숫자 $`n`$자리를 보장하려면 $`\varepsilon_a`$가 얼마나 작아야 하는지 알려 준다.

```math
\varepsilon_a < \varepsilon_s = (0.5 \times 10^{2-n})\%
```

노트북은 이 기준을 `scarborough = lambda n: 0.5 * 10**(2 - n)`으로 정의해 둔다. 함수 `exp_series`의 빈칸을 채워, 반복마다 부분합에 다음 항을 더하게 하라.

```python
def exp_series(x, es=1e-4, maxit=50):
    """
    Approximate e^x using the Maclaurin series, stopping when
    approximate relative error < es (%) or after maxit iterations.

    Parameters
    ----------
    x     : evaluation point
    es    : stopping criterion in % (default 1e-4 %)
    maxit : maximum number of terms (default 50)

    Returns
    -------
    fx   : approximated value
    ea   : final approximate relative error (%)
    k    : number of terms used
    """
    sol = 1.0       # k=0 term
    ea  = 100.0

    print(f'{"k":>4}  {"Approximation":>20}  {"ea (%)":>15}')
    print('-' * 45)

    for k in range(1, maxit + 1):
        sol_old = sol
        # - - - Complete - - -
        sol = sol + 
        # - - - Complete - - -
        if sol != 0:
            ea = abs(sol - sol_old) / abs(sol) * 100
        print(f'{k:>4}  {sol:>20.15f}  {ea:>15.10f}')
        if ea < es:
            break

    return sol, ea, k
```

다음 셀은 유효 숫자 4자리를 목표로 함수를 실행한다.

```python
# Target: 4 significant figures → es = 0.005 %
es = scarborough(4)
print(f'Stopping criterion: {es} %\n')
fx, ea, k = exp_series(0.5, es=es)
print(f'\nResult : {fx:.10f}')
print(f'True   : {np.exp(0.5):.10f}')
print(f'Terms  : {k}')
```

노트북은 이 셀의 기대 출력을 다음과 같이 제시한다.

```text
Stopping criterion: 0.005 %

   k         Approximation           ea (%)
---------------------------------------------
   1     1.500000000000000    33.3333333333
   2     1.625000000000000     7.6923076923
   3     1.645833333333333     1.2658227848
   4     1.648437500000000     0.1579778831
   5     1.648697916666667     0.0157952930
   6     1.648719618055555     0.0013162571

Result : 1.6487196181
True   : 1.6487212707
Terms  : 6
```

## 문제 3. 이차 방정식의 근과 재앙적 소거 (노트북 3.5절)

거의 같은 두 수를 빼면 유효 숫자가 상쇄되어 결과가 반올림 오차에 지배된다. 이를 **재앙적 소거(catastrophic cancellation, subtractive cancellation)** 라 한다. 예를 들어 근의 공식 $`x = \dfrac{-b \pm \sqrt{b^2 - 4ac}}{2a}`$에서 $`b^2 \gg 4ac`$이면, 한 근은 거의 같은 두 큰 수의 뺄셈을 거친다. 노트북은 먼저 방정식 $`x^2 - 100000.0001x + 1 = 0`$으로 표준 공식과 안정한 형태를 비교한다.

```python
# Equation: x^2 - 100000.0001 x + 1 = 0
# Roots: x1 ≈ 100000.0001, x2 ≈ 0.000010000
a, b, c = 1, -100000.0001, 1

discriminant = math.sqrt(b**2 - 4*a*c)

# Standard quadratic formula
x1_std = (-b + discriminant) / (2*a)
x2_std = (-b - discriminant) / (2*a)    # cancellation occurs here!

# Numerically stable form: use x2 = c / (a * x1)
x2_stable = c / (a * x1_std)

# True roots (high precision)
x1_true = 100000.00000999990
x2_true = 0.000010000000100000

print(f'Standard   x1 = {x1_std:.10f}   error = {abs(x1_std - x1_true):.2e}')
print(f'Standard   x2 = {x2_std:.10f}   error = {abs(x2_std - x2_true):.2e}')
print(f'Stable     x2 = {x2_stable:.10f}   error = {abs(x2_stable - x2_true):.2e}')
```

```text
Standard   x1 = 100000.0000900000   error = 8.00e-05
Standard   x2 = 0.0000100000   error = 3.99e-12
Stable     x2 = 0.0000100000   error = 1.09e-13
```

이어지는 연습 셀은 $`|b|`$가 커질 때 표준 공식으로 구한 $`x_2`$의 오차가 어떻게 변하는지 조사한다. 아래의 각 $`b`$ 값에 대해 두 근과 $`x_2`$의 상대 오차를 구하도록 빈칸 `x1_s`와 `x2_stab`을 채워라.

```python
# --- Exercise ---
# Investigate how the error in x2_std changes as |b| grows.
# For each value of b below, compute both roots and the relative error in x2.

print(f'{"b":>12}  {"x2 standard":>16}  {"x2 stable":>16}  {"rel err std (%)":>17}  {"rel err stable (%)":>20}')
print('-' * 90)

for b_val in [1e2, 1e4, 1e6, 1e8, 1e10]:
    a_, b_, c_ = 1, -b_val, 1
    disc = math.sqrt(b_**2 - 4*a_*c_)
    # --- COMPLETE: compute x2_std and x2_stable ---
    x1_s = 
    x2_s = (-b_ - disc) / (2*a_)  # standard formula
    x2_stab =   # stable formula: c / (a * x1)
    # --- COMPLETE ---
    x2_ref = 1.0 / x1_s             # reference (stable)
    err_std  = abs(x2_s    - x2_ref) / abs(x2_ref) * 100
    err_stab = abs(x2_stab - x2_ref) / abs(x2_ref) * 100
    print(f'{b_val:>12.0e}  {x2_s:>16.10f}  {x2_stab:>16.10f}  {err_std:>17.6e}  {err_stab:>20.6e}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
           b       x2 standard         x2 stable    rel err std (%)    rel err stable (%)
------------------------------------------------------------------------------------------
       1e+02      0.0100010002      0.0100010002       1.221123e-11          0.000000e+00
       1e+04      0.0001000000      0.0001000000       1.117663e-07          0.000000e+00
       1e+06      0.0000010000      0.0000010000       7.614492e-04          0.000000e+00
       1e+08      0.0000000075      0.0000000100       2.549419e+01          0.000000e+00
       1e+10      0.0000000000      0.0000000001       1.000000e+02          0.000000e+00
```

## 문제 4. 세 가지 유한 차분 공식 (노트북 4.3절)

테일러 전개 $`f(x + h) = f(x) + f'(x)h + \frac{f''(x)}{2!}h^2 + \cdots`$에서 $`f'(x)`$를 분리하면 다음 근사식을 얻는다.

| 근사 | 공식 | 오차의 차수 |
|:-----|:-----|:-----|
| 전진 차분 | $`f'(x) \approx \dfrac{f(x+h)-f(x)}{h}`$ | $`O(h)`$ |
| 후진 차분 | $`f'(x) \approx \dfrac{f(x)-f(x-h)}{h}`$ | $`O(h)`$ |
| 중심 차분 | $`f'(x) \approx \dfrac{f(x+h)-f(x-h)}{2h}`$ | $`O(h^2)`$ |

*표 1. 테일러 전개에서 얻는 세 가지 유한 차분 공식*

예제 함수는 $`f(x) = -0.1x^4 - 0.15x^3 - 0.5x^2 - 0.25x + 1.2`$이고, 해석적 도함수는 $`f'(x) = -0.4x^3 - 0.45x^2 - x - 0.25`$이다. 앞 셀은 $`x = 0.5`$, $`h = 0.5`$에서 필요한 함숫값을 구해 둔다.

```python
f     = lambda x: -0.1*x**4 - 0.15*x**3 - 0.5*x**2 - 0.25*x + 1.2
dfdx  = lambda x: -0.4*x**3 - 0.45*x**2 -       x - 0.25

x  = 0.5
h  = 0.5

fm = f(x - h)   # f(x-h)
fi = f(x)       # f(x)
fp = f(x + h)   # f(x+h)

print(f'f(x-h) = {fm:8.4f}')
print(f'f(x)   = {fi:8.4f}')
print(f'f(x+h) = {fp:8.4f}')
print(f"True f'(x) = {dfdx(x):8.4f}")
```

```text
f(x-h) =   1.2000
f(x)   =   0.9250
f(x+h) =   0.2000
True f'(x) =  -0.9125
```

다음 셀의 빈칸에 세 가지 유한 차분 공식을 채워라.

```python
# --- COMPLETE: fill in the three finite difference formulas ---
dfdx_forward  =                 
dfdx_backward =                 
dfdx_central  =                 
# --- COMPLETE ---

true_val = dfdx(x)
et_fwd = abs(true_val - dfdx_forward)  / abs(true_val) * 100
et_bwd = abs(true_val - dfdx_backward) / abs(true_val) * 100
et_cen = abs(true_val - dfdx_central)  / abs(true_val) * 100

print(f'{"Method":>12}  {"Estimate":>10}  {"True Rel Err (%)":>18}')
print('-' * 45)
print(f'{"Forward":>12}  {dfdx_forward:>10.5f}  {et_fwd:>18.4f}')
print(f'{"Backward":>12}  {dfdx_backward:>10.5f}  {et_bwd:>18.4f}')
print(f'{"Central":>12}  {dfdx_central:>10.5f}  {et_cen:>18.4f}')
print(f'{"True":>12}  {true_val:>10.5f}')
```

## 문제 5. 스텝 크기에 따른 전진 차분과 중심 차분의 오차 (노트북 4.5절)

노트북 4.4절은 스텝 크기 $`h`$를 줄일 때 중심 차분의 절단 오차는 $`h^2`$에 비례해 줄지만 반올림 오차는 $`\varepsilon_{\text{mach}} / h`$에 비례해 커지므로, 전체 오차를 최소로 하는 최적 스텝 크기 $`h^* \approx \varepsilon_{\text{mach}}^{1/3}`$이 있다고 설명하고, 다음 셀로 중심 차분의 스텝 크기 분석을 수행한다.

```python
def step_size_analysis(func, dfunc, x, n_steps=15):
    """
    Compute central-difference derivative for step sizes h = 10^0, 10^-1, ..., 10^{-(n-1)}
    and record the absolute error against the analytical derivative.
    """
    true_val = dfunc(x)
    H = np.zeros(n_steps)
    E = np.zeros(n_steps)

    h = 1.0
    for i in range(n_steps):
        H[i] = h
        approx = (func(x + h) - func(x - h)) / (2 * h)
        E[i]   = abs(true_val - approx)
        h /= 10

    return H, E

H, E = step_size_analysis(f, dfdx, x=0.5)

print(f'{"Step size":>14}  {"Approx deriv":>18}  {"Abs error":>14}')
print('-' * 52)
for i in range(len(H)):
    approx = (f(0.5 + H[i]) - f(0.5 - H[i])) / (2 * H[i])
    print(f'{H[i]:>14.1e}  {approx:>18.14f}  {E[i]:>14.4e}')
```

연습 문제는 다음과 같다. 전진 차분과 중심 차분에 대해 이 스텝 크기 분석을 반복하고 두 결과를 한 그림에 그려, 중심 차분은 $`O(h^2)`$로, 전진 차분은 $`O(h)`$로 수렴함을 눈으로 확인하라. 다음 셀의 빈칸 `E_fwd[i]`와 `E_cen[i]`를 채운다.

```python
n_steps = 15
H  = np.zeros(n_steps)
E_fwd = np.zeros(n_steps)
E_cen = np.zeros(n_steps)

true_val = dfdx(0.5)
h = 1.0
for i in range(n_steps):
    H[i] = h
    # - - - - Complete - - - 
    E_fwd[i] =  
    E_cen[i] =  
    # - - - - Complete - - - 
    h /= 10

plt.figure(figsize=(8, 5))
plt.loglog(H, E_fwd, 'rs-', markersize=6, label='Forward difference')
plt.loglog(H, E_cen, 'bo-', markersize=6, label='Central difference')

h_ref = np.logspace(0, -8, 80)
plt.loglog(h_ref, 0.5  * h_ref,    'r--', alpha=0.5, label=r'$O(h)$')
plt.loglog(h_ref, 0.05 * h_ref**2, 'b--', alpha=0.5, label=r'$O(h^2)$')

plt.xlabel('Step size $h$')
plt.ylabel('Absolute error')
plt.title('Forward vs. central difference: convergence rate')
plt.legend()
plt.grid(True, which='both', ls=':', alpha=0.5)
plt.tight_layout()
```

---

<br>

문제를 풀어 본 뒤 [과제 해설](../Assignment-Explanation.ko.md)에서 정답과 풀이를, [완성한 노트북](../Solution/Lab02-SourceOfErrors_Completed.ipynb)에서 작성한 코드 전체를 볼 수 있다.
