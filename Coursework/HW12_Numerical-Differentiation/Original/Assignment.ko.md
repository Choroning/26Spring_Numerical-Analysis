# HW12 과제 안내문 — 수치 미분

[레포지토리 안내](../../../README.ko.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **최종 수정일:** 2026-11-16

수치해석 14주차 과제이다. 교수자가 배포한 실습 노트북 `Lab13-diff.ipynb`의 빈칸을 채워 실행한다. 노트북은 기준 함수와 테일러 전개 복습, 전진, 후진, 중심 차분, 고차 스텐실과 리처드슨 외삽법, 2계도함수와 혼합 편도함수, 스텝 크기 선택과 복소 스텝 방법, 다변수 미분(야코비 행렬과 헤세 행렬), FEM 도함수 복원, 공학 응용의 여덟 부분으로 이루어진다. 나머지 코드 셀은 7절의 `COMPLETE` 주석이 붙은 셀과 Exercise 셀까지 모두 코드가 채워진 채 배포되었으며, 빈칸은 2.1절의 `COMPLETE`라고 적힌 두 주석 줄 사이의 한 곳뿐이다. 노트북의 첫 코드 셀은 `numpy`를 `np`로, `matplotlib.pyplot`을 `plt`로 불러오고, `scipy.linalg`의 `solve`를 불러오며 경고 출력을 끈다.

아래 코드 블록은 노트북의 셀을 옮긴 것이며, 코드가 비어 있는 자리가 채울 곳이다. 배포된 빈 노트북은 [Lab13-diff.ipynb](Lab13-diff.ipynb)에 있다.

- **제출 기한:** 배포 후 1주
- **제출 방법:** 빈칸을 채워 실행한 노트북을 온라인으로 제출한다.

---

<br>

## 문제. 중심 차분 함수 (노트북 2.1절)

노트북 1절은 교재(Chapra와 Canale)와 같은 다항식을 기준 함수로 쓴다.

$$
f(x) = 0.1x^5 - 0.2x^4 + 3.0x^3 + 2.0x^2 - x + 0.5
$$

해석적 도함수는 다음과 같다.

$$
f'(x) = 0.5x^4 - 0.8x^3 + 9.0x^2 + 4.0x - 1
$$

$$
f''(x) = 2x^3 - 2.4x^2 + 18x + 4
$$

$$
f'''(x) = 6x^2 - 4.8x + 18
$$

모든 유한 차분 공식은 테일러 정리에서 나온다.

$$
f(x+h) = f(x) + h f'(x) + \frac{h^2}{2!} f''(x) + \frac{h^3}{3!} f'''(x) + \cdots
$$

$$
f(x-h) = f(x) - h f'(x) + \frac{h^2}{2!} f''(x) - \frac{h^3}{3!} f'''(x) + \cdots
$$

이 전개식을 정리하면 절단 오차를 조절하면서 $`f'(x)`$, $`f''(x)`$ 등을 분리할 수 있다. 다음 셀은 기준 함수와 도함수를 정하고 $`x_0 = 0.5`$에서의 값을 출력한 뒤 그래프를 그린다.

```python
# Benchmark function and its analytic derivatives
f    = lambda x: 0.1*x**5 - 0.2*x**4 + 3.0*x**3 + 2.0*x**2 - x + 0.5
df   = lambda x: 0.5*x**4 - 0.8*x**3 + 9.0*x**2 + 4.0*x - 1.0
d2f  = lambda x: 2.0*x**3 - 2.4*x**2 + 18.0*x + 4.0
d3f  = lambda x: 6.0*x**2 - 4.8*x + 18.0

# Evaluation point and step size
x0 = 0.5
print(f'Benchmark at x = {x0}:')
print(f'  f({x0})    = {f(x0):.8f}')
print(f'  f\'({x0})   = {df(x0):.8f}')
print(f'  f\'\'({x0})  = {d2f(x0):.8f}')
print(f'  f\'\'\'({x0}) = {d3f(x0):.8f}')

# Visualize f and its derivatives
x_plot = np.linspace(-1, 2, 400)
fig, axes = plt.subplots(1, 2, figsize=(13, 4))

ax = axes[0]
ax.plot(x_plot, f(x_plot),   'b-',  lw=2, label='$f(x)$')
ax.plot(x_plot, df(x_plot),  'r-',  lw=2, label="$f'(x)$")
ax.plot(x_plot, d2f(x_plot), 'g-',  lw=2, label="$f''(x)$")
ax.axvline(x0, color='k', ls=':', lw=1.5, label=f'$x_0={x0}$')
ax.set(xlabel='x', ylabel='value', title='Benchmark Function and Derivatives')
ax.legend(); ax.grid(True, alpha=0.3)

# Taylor expansion visualization at x0
ax = axes[1]
h_demo = 0.4
x_expand = np.linspace(x0 - h_demo - 0.1, x0 + h_demo + 0.1, 200)
T1 = f(x0) + df(x0)*(x_expand - x0)
T2 = T1 + d2f(x0)/2*(x_expand - x0)**2
T3 = T2 + d3f(x0)/6*(x_expand - x0)**3
ax.plot(x_expand, f(x_expand),  'b-',  lw=2.5, label='$f(x)$')
ax.plot(x_expand, T1,           'r--', lw=1.5, label='Taylor $O(h)$')
ax.plot(x_expand, T2,           'g--', lw=1.5, label='Taylor $O(h^2)$')
ax.plot(x_expand, T3,           'm--', lw=1.5, label='Taylor $O(h^3)$')
ax.scatter([x0], [f(x0)], c='k', s=80, zorder=6)
ax.set(xlabel='x', ylabel='value',
       title=f'Taylor Approximation around $x_0={x0}$')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

plt.tight_layout(); plt.show()
```

```text
Benchmark at x = 0.5:
  f(0.5)    = 0.86562500
  f'(0.5)   = 3.18125000
  f''(0.5)  = 12.65000000
  f'''(0.5) = 17.10000000
```

노트북 2.1절은 테일러 전개 $`f(x+h) = f(x) + hf'(x) + \frac{h^2}{2}f''(\xi)`$를 $`f'(x)`$에 대해 풀어 다음 세 공식을 얻는다.

| 이름 | 공식 | 절단 오차 | 사용하는 값 |
|:-----|:-----|:-----|:-----|
| 전진 차분 | $`f'(x) \approx \dfrac{f(x+h) - f(x)}{h}`$ | $`O(h)`$, 1차 | $`f(x)`$, $`f(x+h)`$ |
| 후진 차분 | $`f'(x) \approx \dfrac{f(x) - f(x-h)}{h}`$ | $`O(h)`$, 1차 | $`f(x-h)`$, $`f(x)`$ |
| 중심 차분 | $`f'(x) \approx \dfrac{f(x+h) - f(x-h)}{2h}`$ | $`O(h^2)`$, 2차 | $`f(x-h)`$, $`f(x+h)`$ |

*표 1. 1계도함수의 세 가지 기본 유한 차분 공식*

중심 차분이 2차인 이유는 다음과 같다. 두 테일러 전개를 빼면 $`h^2`$ 항이 상쇄된다.

$$
f(x+h) - f(x-h) = 2h f'(x) + \frac{2h^3}{6} f'''(x) + \cdots
$$

$$
\Rightarrow f'(x) = \frac{f(x+h)-f(x-h)}{2h} - \frac{h^2}{6}f'''(x) + O(h^4)
$$

주도 오차는 $`O(h^2)`$이며, 전진 차분과 후진 차분의 $`O(h)`$와 대조된다. 다음 셀의 빈칸을 채워 함수 `central_diff`를 완성하라. 같은 셀의 `forward_diff`와 `backward_diff`는 채워진 채 배포되었다.

```python
def forward_diff(func, x, h):
    """
    First-order forward difference approximation of f'(x).
    Truncation error: O(h).

    Parameters
    ----------
    func : callable — function f
    x    : float   — point of differentiation
    h    : float   — step size

    Returns
    -------
    float : approximation of f'(x)
    """
    fdiff = (func(x + h) - func(x)) / h
    return fdiff 


def backward_diff(func, x, h):
    """
    First-order backward difference approximation of f'(x).
    Truncation error: O(h).
    """
    bdiff = (func(x) - func(x - h)) / h
    return bdiff


def central_diff(func, x, h):
    """
    Second-order central difference approximation of f'(x).
    Truncation error: O(h^2).
    """
    # - - - - - COMPLETE - - - - - 
    cdiff = 
    # - - - - - COMPLETE - - - - - 
    return cdiff
    


# Compare at x0 = 0.5
exact = df(x0)
print(f'Analytic f\'({x0}) = {exact:.10f}\n')
print(f'{"h":>8}  {"Forward":>14}  {"err_F":>10}  {"Backward":>14}  {"err_B":>10}  {"Central":>14}  {"err_C":>10}')
print('-' * 95)
for h in [0.5, 0.25, 0.1, 0.05, 0.01, 0.001]:
    fwd  = forward_diff(f, x0, h)
    bwd  = backward_diff(f, x0, h)
    cnt  = central_diff(f, x0, h)
    print(f'{h:>8.4f}  {fwd:>14.8f}  {abs(fwd-exact):>10.2e}  '
          f'{bwd:>14.8f}  {abs(bwd-exact):>10.2e}  '
          f'{cnt:>14.8f}  {abs(cnt-exact):>10.2e}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
Analytic f'(0.5) = 3.1812500000

       h         Forward       err_F        Backward       err_B         Central       err_C
-----------------------------------------------------------------------------------------------
  0.5000      7.06875000    3.89e+00      0.73125000    2.45e+00      3.90000000    7.19e-01
  0.2500      4.94179687    1.76e+00      1.77773437    1.40e+00      3.35976563    1.79e-01
  0.1000      3.84231000    6.61e-01      2.57721000    6.04e-01      3.20976000    2.85e-02
  0.0500      3.50463188    3.23e-01      2.87211937    3.09e-01      3.18837563    7.13e-03
  0.0100      3.24478505    6.35e-02      3.11828495    6.30e-02      3.18153500    2.85e-04
  0.0010      3.18757785    6.33e-03      3.17492785    6.32e-03      3.18125285    2.85e-06
```

---

<br>

문제를 풀어 본 뒤 [과제 해설](../Assignment-Explanation.ko.md)에서 정답과 풀이를, [완성한 노트북](../Solution/Lab13-diff_Completed.ipynb)에서 작성한 코드 전체를 볼 수 있다.
