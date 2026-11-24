# HW10 과제 안내문 — 데이터 보간과 스플라인

[레포지토리 안내](../../../README.ko.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **최종 수정일:** 2026-11-24

수치해석 12주차 과제이다. 교수자가 배포한 실습 노트북 `Lab11-Spline.ipynb`의 빈칸을 채워 실행한다. 노트북은 보간의 필요성, 뉴턴의 분할 차분, 라그랑주 보간, 룽게 현상, 이차 스플라인, 삼차 스플라인, `scipy.interpolate.CubicSpline`의 일곱 부분으로 이루어진다. Exercise 2.1과 3.1의 셀을 포함해 나머지 셀은 코드가 채워진 채로 배포되었으며, 빈칸은 3절의 함수 `lagrange_basis`에 `COMPLETE`라고 적힌 두 주석 줄 사이의 한 곳이다. 노트북의 첫 코드 셀은 `numpy`를 `np`로, `matplotlib.pyplot`을 `plt`로 불러오고 `scipy.interpolate`에서 `CubicSpline`과 `interp1d`를 불러온다.

아래 코드 블록은 노트북의 셀을 옮긴 것이며, 코드가 비어 있는 자리가 채울 곳이다. 배포된 빈 노트북은 [Lab11-Spline.ipynb](Lab11-Spline.ipynb)에 있다.

---

<br>

## 문제. 라그랑주 기저 함수 (노트북 3.1절)

노트북 1.2절의 셀은 데이터셋 A를 다음과 같이 정의한다.

```python
# Dataset A: Air density vs. temperature (from textbook Table 17.1)
T_air  = np.array([300., 400., 500.])        # temperature [K]
rho_air = np.array([0.616, 0.525, 0.457])   # density [kg/m^3]
```

노트북 2절은 이 데이터로 뉴턴의 분할 차분표 `b_air = newtdd_table(T_air, rho_air)`를 만들고, `newtdd_eval`로 뉴턴 보간 다항식을 계산해 둔다. 3절의 라그랑주 형태는 보간 다항식을 기저 함수 $`L_i(x)`$로 직접 쓴다.

$$
p_n(x) = \sum_{i=0}^{n} y_i L_i(x)
$$

각 기저 함수는 $`L_i(x_j) = \delta_{ij}`$, 곧 자기 노드에서 1이고 다른 모든 노드에서 0이며, 다음과 같다.

```math
L_i(x) = \prod_{\substack{j=0 \\ j \neq i}}^{n} \frac{x - x_j}{x_i - x_j}
```

노트북은 두 방법을 다음과 같이 비교하고, 두 방법 모두 $`n+1`$개의 점을 지나는 $`n`$차 이하의 유일한 다항식을 준다고 덧붙인다.

| 특징 | 뉴턴 | 라그랑주 |
|:-----|:-----|:-----|
| 데이터 점을 추가할 때 | 쉽다(새 열 하나) | 모두 다시 계산해야 한다 |
| 개념의 중심 | 분할 차분 | 기저 함수 |
| 결과 | 같은 다항식 | 같은 다항식 |

*표 1. 뉴턴 보간과 라그랑주 보간의 비교*

함수 `lagrange_basis`의 빈칸을 채워 $`i`$번째 라그랑주 기저 함수의 값을 계산하라. 같은 셀의 `lagrange_interp`는 이 함수로 보간값을 구한다.

```python
def lagrange_basis(x, i, xx):
    """
    Evaluate the i-th Lagrange basis polynomial at xx.

    Parameters
    ----------
    x  : array of interpolation nodes
    i  : int — which basis polynomial
    xx : scalar or array — evaluation point(s)

    Returns
    -------
    Li : value(s) of L_i(xx)
    """
    x  = np.array(x, dtype=float)
    xx = np.atleast_1d(np.array(xx, dtype=float))
    n  = len(x)
    Li = np.ones(xx.shape)

    for j in range(n):
        if j != i:
            # - - - - - - - COMPLETE - - - - - - - 
            Li =                # Lagrange basis
            # - - - - - - - COMPLETE - - - - - - - 
    return Li


def lagrange_interp(x, y, xx):
    """
    Evaluate the Lagrange interpolating polynomial at xx.

    Parameters
    ----------
    x, y : array — interpolation nodes and values
    xx   : scalar or array — evaluation point(s)

    Returns
    -------
    yy : interpolated value(s)
    """
    x  = np.array(x, dtype=float)
    y  = np.array(y, dtype=float)
    xx = np.atleast_1d(np.array(xx, dtype=float))
    n  = len(x)
    yy = np.zeros(xx.shape)

    for i in range(n):
        yy = yy + y[i] * lagrange_basis(x, i, xx)
    return yy if yy.size > 1 else yy[0]
```

다음 셀은 데이터셋 A의 세 기저 함수를 그리고, 라그랑주 보간과 뉴턴 보간의 곡선을 비교한 뒤 350 K에서 두 값을 출력한다. 셀의 마지막 부분은 다음과 같다.

```python
rho_L = lagrange_interp(T_air, rho_air, 350.)
rho_N = newtdd_eval(T_air, b_air, 350.)
print(f'Lagrange at 350 K = {rho_L:.8f}')
print(f'Newton   at 350 K = {rho_N:.8f}')
print(f'Difference        = {abs(rho_L - rho_N):.2e}  (should be ~machine precision)')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
Lagrange at 350 K = 0.56762500
Newton   at 350 K = 0.56762500
Difference        = 0.00e+00  (should be ~machine precision)
```

---

<br>

문제를 풀어 본 뒤 [과제 해설](../Assignment-Explanation.ko.md)에서 정답과 풀이를, [완성한 노트북](../Solution/Lab11-Spline_Completed.ipynb)에서 작성한 코드 전체를 볼 수 있다.
