# HW11 과제 안내문 — 수치 적분

[레포지토리 안내](../../../README.ko.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **최종 수정일:** 2026-11-24

수치해석 13주차 과제이다. 교수자가 배포한 실습 노트북 `Lab12-Int.ipynb`의 빈칸을 채워 실행한다. 노트북은 기준 함수 소개, 사다리꼴 법칙, 심프슨 법칙, 리처드슨 외삽법과 롬버그 적분, 가우스-르장드르 구적법, 가우스-로바토 구적법, `scipy.integrate`, 공학 응용의 여덟 부분으로 이루어진다. 나머지 코드 셀은 Exercise 셀까지 모두 코드가 채워진 채 배포되었으며, 빈칸은 `COMPLETE`라고 적힌 두 주석 줄 사이의 한 곳뿐이다. 노트북의 첫 코드 셀은 `numpy`를 `np`로, `matplotlib.pyplot`을 `plt`로 불러오고, `scipy`의 `integrate`와 `numpy.polynomial.legendre`의 `leggauss`를 불러온다.

아래 코드 블록은 노트북의 셀을 옮긴 것이며, 코드가 비어 있는 자리가 채울 곳이다. 배포된 빈 노트북은 [Lab12-Int.ipynb](Lab12-Int.ipynb)에 있다.

---

<br>

## 문제. 복합 사다리꼴 법칙 (노트북 2.2절)

노트북 1절은 정확한 적분값을 아는 하나의 **기준 함수**를 모든 방법의 비교에 쓴다.

```math
f(x) = 0.2 + 25x - 200x^2 + 675x^3 - 900x^4 + 400x^5, \quad x \in [0,\, 0.8]
```

```math
I_{\text{exact}} = \int_0^{0.8} f(x)\,dx = 1.640533\ldots
```

참 백분율 상대 오차는 다음과 같이 정의한다.

```math
\varepsilon_t = \left|\frac{I_{\text{exact}} - I_{\text{approx}}}{I_{\text{exact}}}\right| \times 100\,\%
```

다음 셀은 기준 함수와 정확한 적분값, 오차 함수를 정하고 함수의 그래프를 그린다.

```python
# Benchmark function and exact integral
f   = lambda x: 0.2 + 25*x - 200*x**2 + 675*x**3 - 900*x**4 + 400*x**5
a, b  = 0.0, 0.8
# Exact integral, computed in closed form:
#   F(x) = 0.2x + 12.5 x^2 - (200/3) x^3 + (675/4) x^4 - 180 x^5 + (200/3) x^6
#   F(0.8) - F(0) = 49216/30000 = 3076/1875 = 1.6405333333...
I_exact = 49216 / 30000     # = 1.6405333333... (full double precision)

def pct_err(I_approx):
    """Percent true relative error."""
    return abs(I_exact - I_approx) / abs(I_exact) * 100

# Visualise
x_plot = np.linspace(a, b, 300)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(x_plot, f(x_plot), 'b-', lw=2)
ax.fill_between(x_plot, f(x_plot), alpha=0.15, color='blue')
ax.set(xlabel='x', ylabel='f(x)', title='Benchmark Function  (shaded area = exact integral)')
ax.text(0.5, 1.5, f'$I_{{\\rm exact}} = {I_exact:.10f}$', fontsize=13)
ax.grid(True, alpha=0.3)
plt.tight_layout(); plt.show()
```

노트북 2.1절의 단일 구간 사다리꼴 법칙은 $`(a, f(a))`$와 $`(b, f(b))`$를 지나는 일차 다항식으로 $`f`$를 근사한다.

```math
I \approx (b-a)\,\frac{f(a)+f(b)}{2}
```

이 공식은 1차 이하의 다항식에서 정확하며, 절단 오차는 다음과 같다.

$$
E_t = -\frac{1}{12}f''(\xi)(b-a)^3, \quad \xi \in [a,b]
$$

2.2절의 복합 사다리꼴 법칙은 $`[a,b]`$를 폭 $`h = (b-a)/n`$인 같은 크기의 소구간 $`n`$개로 나눈다.

$$
I \approx \frac{h}{2}\left[f(x_0) + 2\sum_{i=1}^{n-1}f(x_i) + f(x_n)\right]
$$

전체 오차는 $`h^2 \propto n^{-2}`$에 비례하므로, $`n`$을 두 배로 늘리면 오차가 약 4분의 1로 줄어든다. 다음 셀의 빈칸을 채워 함수 `composite_trap`을 완성하라.

```python
# ---- Single-segment trapezoidal rule ----
I_trap1 = (b - a) * (f(a) + f(b)) / 2.0
print(f'Single-segment Trapezoidal: I = {I_trap1:.6f},  error = {pct_err(I_trap1):.2f}%')

# ---- Composite trapezoidal rule ----
def composite_trap(func, a, b, n):
    """
    Composite trapezoidal rule.

    Parameters
    ----------
    func : callable
    a, b : integration limits
    n    : number of sub-intervals

    Returns
    -------
    I : float
    """
    x = np.linspace(a, b, n + 1)
    y = func(x)
    h = (b - a) / n
    # - - - - - - COMPLETE - - - - - - - 
    I = h/2 * 
    # - - - - - - COMPLETE - - - - - - - 
    return I

# Test with increasing n
print('\nComposite Trapezoidal Rule:')
print(f'{"n":>6}  {"I":>12}  {"error %":>10}')
for n in [1, 2, 4, 8, 16, 32, 64]:
    I = composite_trap(f, a, b, n)
    print(f'{n:>6}  {I:>12.8f}  {pct_err(I):>10.4f}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
Single-segment Trapezoidal: I = 0.172800,  error = 89.47%

Composite Trapezoidal Rule:
     n             I     error %
     1    0.17280000     89.4668
     2    1.06880000     34.8505
     4    1.48480000      9.4928
     8    1.60080000      2.4220
    16    1.63055000      0.6085
    32    1.63803438      0.1523
    64    1.63990840      0.0381
```

---

<br>

문제를 풀어 본 뒤 [과제 해설](../Assignment-Explanation.ko.md)에서 정답과 풀이를, [완성한 노트북](../Solution/Lab12-Int_Completed.ipynb)에서 작성한 코드 전체를 볼 수 있다.
