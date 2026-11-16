# HW8 과제 안내문 — 회귀 분석

[레포지토리 안내](../../../README.ko.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **최종 수정일:** 2026-11-16

수치해석 10주차 과제이다. 교수자가 배포한 실습 노트북 `Lab09-Regression.ipynb`의 빈칸을 채워 실행한다. 노트북은 회귀의 필요성, 단순 선형 회귀, 선형화할 수 있는 비선형 회귀, 일반 선형 최소제곱, `scipy.optimize.curve_fit`을 이용한 비선형 회귀, 모델 비교와 선택의 여섯 부분으로 이루어지며, 빈칸은 `COMPLETE`라고 적힌 두 주석 줄 사이에 모두 5곳 있다. 노트북의 첫 코드 셀은 `numpy`를 `np`로, `matplotlib.pyplot`을 `plt`로, `scipy.stats`를 `stats`로 불러오고 `scipy.optimize`에서 `curve_fit`을 불러온다.

아래 각 문제의 코드 블록은 노트북의 셀을 옮긴 것이며, 코드가 비어 있는 자리가 채울 곳이다. 배포된 빈 노트북은 [Lab09-Regression.ipynb](Lab09-Regression.ipynb)에 있다.

- **제출 기한:** 배포 후 1주
- **제출 방법:** 빈칸을 채워 실행한 노트북을 온라인으로 제출한다.

노트북 1.2절의 셀은 실습 전체에서 쓰는 데이터셋을 다음과 같이 정의하고 산점도로 그린다.

```python
# Dataset 1: Beam deflection (clearly linear)
load       = np.array([0, 10, 20, 30, 40, 50, 60, 70], dtype=float)   # kN
deflection = np.array([0.0, 2.1, 4.3, 6.2, 8.5, 10.4, 12.6, 14.8])   # mm

# Dataset 2: Stress-response data (power law behavior)
stress   = np.linspace(1, 8, 8) * 10.0                                # MPa
response = np.array([25, 70, 380, 550, 610, 1220, 830, 1450], dtype=float)

# Dataset 3: Borosilicate glass refractive index (Cauchy equation)
wavelength = np.array([0.6563, 0.6439, 0.5890, 0.5338,
                       0.5086, 0.4861, 0.4340, 0.3988])  # micrometers
ridx       = np.array([1.50883, 1.50917, 1.51124, 1.51386,
                       1.51534, 1.51690, 1.52136, 1.52546])
```

---

<br>

## 문제 1. 단순 선형 회귀의 계수 (노트북 2.1절)

$`S_r = \sum(y_i - a_0 - a_1 x_i)^2`$을 최소화하여 $`\hat{y} = a_0 + a_1 x`$를 적합한다. 편도함수를 0으로 놓으면 다음을 얻는다.

$$
a_1 = \frac{n\sum x_i y_i - \sum x_i \sum y_i}{n\sum x_i^2 - \left(\sum x_i\right)^2} = \frac{\sum x_i y_i - n\bar{x}\bar{y}}{\sum x_i^2 - n\bar{x}^2}, \qquad a_0 = \bar{y} - a_1 \bar{x}
$$

함수 `linreg`의 빈칸을 채워 기울기 `a1`과 절편 `a0`을 구하라.

```python
def linreg(xi, yi):
    """
    Simple linear regression: fit y = a0 + a1*x via least squares.

    Returns
    -------
    a0 : intercept
    a1 : slope
    """
    xi = np.array(xi, dtype=float)
    yi = np.array(yi, dtype=float)
    n  = len(xi)

    # --- COMPLETE ---
    a1 =                  
    a0 =                  
    # --- COMPLETE ---

    return a0, a1
```

다음 셀은 이 함수를 보 처짐 데이터에 적용하고, 적합 직선과 잔차 그림을 그린다. 계산과 출력 부분은 다음과 같다.

```python
a0, a1 = linreg(load, deflection)
print(f'Slope     a1 = {a1:.6f} mm/kN')
print(f'Intercept a0 = {a0:.6f} mm')
print(f'Stiffness  k = 1/a1 = {1/a1:.4f} kN/mm')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
Slope     a1 = 0.210357 mm/kN
Intercept a0 = 0.000000 mm
Stiffness  k = 1/a1 = 4.7538 kN/mm
```

## 문제 2. 적합도 통계량 (노트북 2.2절)

노트북은 적합도 통계량을 다음과 같이 정리한다.

| 양 | 공식 | 의미 |
|:-----|:-----|:-----|
| $`S_T`$ | $`\sum(y_i - \bar{y})^2`$ | 총 제곱합 |
| $`S_r`$ (SSE) | $`\sum(y_i - \hat{y}_i)^2`$ | 잔차(설명되지 않는) 제곱합 |
| SSR | $`S_T - S_r`$ | 설명된 제곱합 |
| $`R^2`$ | $`1 - S_r/S_T`$ | 설명된 분산의 비율 |
| $`s_e`$ | $`\sqrt{S_r/(n-2)}`$ | 추정의 표준 오차 |

*표 1. 적합도 통계량*

함수 `goodness_of_fit`의 빈칸을 채워 `SST`, `SSE`, `SSR`, `R2`를 구하라.

```python
def goodness_of_fit(yi, yi_hat, n_params=2):
    """
    Compute goodness-of-fit statistics.

    Returns dict: SST, SSE, SSR, R2, R2adj, se
    """
    n     = len(yi)
    y_bar = np.mean(yi)

    # --- COMPLETE ---
    SST =             # total sum of squares
    SSE =             # residual sum of squares
    SSR =             # regression sum of squares
    R2  =             # coefficient of determination
    # --- COMPLETE ---

    se    = np.sqrt(SSE / (n - n_params))
    R2adj = 1 - (SSE / (n - n_params)) / (SST / (n - 1))
    return dict(SST=SST, SSE=SSE, SSR=SSR, R2=R2, R2adj=R2adj, se=se)
```

다음 셀은 문제 1의 직선에 이 함수를 적용한다.

```python
yhat_lin = a0 + a1 * load
gof_lin  = goodness_of_fit(deflection, yhat_lin, n_params=2)

print('=== Goodness of Fit: Beam Deflection (Linear) ===')
print(f'SST    = {gof_lin["SST"]:.6f}   (total variation)')
print(f'SSE    = {gof_lin["SSE"]:.6f}   (residual/unexplained)')
print(f'SSR    = {gof_lin["SSR"]:.6f}   (explained by model)')
print(f'R2     = {gof_lin["R2"]:.8f}')
print(f'Adj R2 = {gof_lin["R2adj"]:.8f}')
print(f'se     = {gof_lin["se"]:.6f} mm')
print(f'SST = SSR + SSE? {np.isclose(gof_lin["SST"], gof_lin["SSR"] + gof_lin["SSE"])}')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
=== Goodness of Fit: Beam Deflection (Linear) ===
SST    = 185.898750   (total variation)
SSE    = 0.048214   (residual/unexplained)
SSR    = 185.850536   (explained by model)
R2     = 0.99974064
Adj R2 = 0.99969742
se     = 0.089642 mm
SST = SSR + SSE? True
```

## 문제 3. 거듭제곱 모델의 선형화 (노트북 3.1절)

거듭제곱 모델 $`y = \alpha x^\beta`$의 양변에 $`\log_{10}`$을 취하면 선형식을 얻는다.

$$
\underbrace{\log_{10} y}_{Y} = \underbrace{\log_{10} \alpha}_{A_0} + \beta \underbrace{\log_{10} x}_{X} \quad\Rightarrow\quad Y = A_0 + \beta X
$$

매개변수는 $`\beta = A_1`$, $`\alpha = 10^{A_0}`$으로 되찾는다. 노트북은 로그 선형화가 원래의 $`y`$ 공간이 아니라 변환된(로그) 공간에서 SSE를 최소화하므로, 잡음이 크면 5절의 `curve_fit`으로 정확한 비선형 최소제곱을 구하라고 주의를 준다. 함수 `power_fit`의 빈칸을 채워 자료를 로그 공간으로 변환하고 거듭제곱 모델의 매개변수를 되찾아라.

```python
def power_fit(xi, yi):
    """
    Fit y = alpha * x^beta using log10 linearization.

    Returns
    -------
    alpha, beta : power law parameters
    A0, A1      : intercept and slope in log10 space
    """
    xi = np.array(xi, dtype=float)
    yi = np.array(yi, dtype=float)

    # --- COMPLETE: transform to log10 space ---
    X =              
    Y =              
    # --- COMPLETE ---

    A0, A1 = linreg(X, Y)      # fit in log space

    # --- COMPLETE: recover power law parameters ---
    alpha =          
    beta  =          
    # --- COMPLETE ---

    return alpha, beta, A0, A1
```

다음 셀은 응력과 응답 데이터에 이 함수를 적용하고, 앞 셀에서 같은 데이터에 직선을 적합한 결과(`gof_bad`)와 함께 그림으로 비교한 뒤 두 적합의 $`R^2`$과 $`s_e`$를 출력한다. 셀의 앞부분은 다음과 같다.

```python
alpha_pw, beta_pw, A0, A1 = power_fit(stress, response)
print(f'Power law fit:  y = {alpha_pw:.4f} * x^{beta_pw:.4f}')
print(f'Exponent beta ~ {beta_pw:.2f}  (response ~ stress squared)')

yhat_pw = alpha_pw * stress**beta_pw
gof_pw  = goodness_of_fit(response, yhat_pw, n_params=2)
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
Power law fit:  y = 0.2741 * x^1.9842
Exponent beta ~ 1.98  (response ~ stress squared)

Linear:    R2=0.88049,  se=189.79
Power law: R2=0.80882,  se=240.04
```

## 문제 4. 다항 회귀의 정규 방정식 (노트북 4.1절)

매개변수 $`\mathbf{c}`$에 대해 선형인 모델은 모두 $`\hat{y} = \mathbf{X}\mathbf{c}`$로 쓸 수 있다. 여기서 설계 행렬(design matrix) $`\mathbf{X}`$($`n\times(m+1)`$)의 행은 $`[z_0(x_i), z_1(x_i), \ldots, z_m(x_i)]`$이다. 최소제곱 해는 다음 정규 방정식을 만족한다.

$$
\mathbf{X}^T\mathbf{X}\mathbf{c} = \mathbf{X}^T\mathbf{y} \quad\Rightarrow\quad \mathbf{c} = (\mathbf{X}^T\mathbf{X})^{-1}\mathbf{X}^T\mathbf{y}
$$

| 기저 함수 $`z_j(x)`$ | 모델 |
|:-----|:-----|
| $`1, x`$ | 단순 선형 |
| $`1, x, x^2, \ldots, x^m`$ | $`m`$차 다항식 |
| $`1, 1/\lambda^2, 1/\lambda^4`$ | 코시 방정식 |
| $`1, x_1, x_2`$ | 다중 선형 회귀 |

*표 2. 기저 함수와 모델*

$`n > m+1`$이면 연립방정식은 과결정(overdetermined)이므로 정규 방정식이나 수치적으로 더 안정한 `np.linalg.lstsq`를 쓴다. 함수 `poly_lstsq`의 빈칸을 채워 설계 행렬을 만들고 정규 방정식을 풀어라.

```python
def poly_lstsq(xi, yi, m):
    """
    Polynomial regression: fit y = c0 + c1*x + ... + cm*x^m.

    Parameters
    ----------
    xi, yi : array-like
    m      : int  — polynomial degree

    Returns
    -------
    c : array  — coefficients [c0, c1, ..., cm]
    X : array  — design matrix  (n x (m+1))
    """
    xi = np.array(xi, dtype=float)
    yi = np.array(yi, dtype=float)
    n  = len(xi)

    # --- COMPLETE: build design matrix ---
    X = np.zeros((n, m + 1))
    for j in range(m + 1):
        X[:, j] =            # (j-th column = x to the power j)
    # --- COMPLETE ---

    # --- COMPLETE: solve normal equations X^T X c = X^T y ---
    c =                      
    # --- COMPLETE ---

    return c, X
```

다음 셀은 속도 데이터에 1, 2, 3차 다항식을 적합하여 그린다.

```python
# Polynomial regression on velocity data (quadratic model expected)
x_vel = np.array([0., 1., 2., 3., 4., 5.], dtype=float)  # time [s]
y_vel = np.array([2.1, 7.7, 13.6, 27.2, 40.9, 61.1])     # velocity [m/s]

fig, ax = plt.subplots(figsize=(8, 4))
ax.scatter(x_vel, y_vel, c='k', marker='s', s=80, zorder=5, label='Data')

colors = ['blue', 'red', 'green']
for m, col in zip([1, 2, 3], colors):
    c, X = poly_lstsq(x_vel, y_vel, m)
    x_plt = np.linspace(-0.2, 5.5, 300)
    y_plt = sum(c[j] * x_plt**j for j in range(m + 1))
    yhat  = X @ c
    gof_p = goodness_of_fit(y_vel, yhat, n_params=m + 1)
    ax.plot(x_plt, y_plt, color=col, lw=2,
            label=f'Deg {m}: $R^2={gof_p["R2"]:.5f}$, $s_e={gof_p["se"]:.3f}$')

ax.set(xlabel='Time [s]', ylabel='Velocity [m/s]',
       title='Polynomial Regression: Degrees 1, 2, 3')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
plt.tight_layout(); plt.show()

c2, _ = poly_lstsq(x_vel, y_vel, 2)
print(f'Quadratic: y = {c2[0]:.4f} + {c2[1]:.4f}t + {c2[2]:.4f}t²')
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
Quadratic: y = 2.4786 + 2.3593t + 1.8607t²
```

## 문제 5. 앙투안 방정식의 비선형 회귀 (노트북 5.1절)

미카엘리스-멘텐 식 $`y = V_{\max} x / (K_m + x)`$나 앙투안 방정식 $`\log_{10} P = A - B/(C + T)`$처럼 선형화할 수 없는 모델이 있다. `curve_fit`은 감쇠 뉴턴 방법인 레벤버그-마쿼트 알고리즘으로 $`\sum(y_i - f(x_i;\mathbf{c}))^2`$을 직접 최소화한다. 노트북은 사용법을 다음과 같이 제시한다.

```python
popt, pcov = curve_fit(f, xdata, ydata, p0=initial_guess)
# popt : optimal parameters
# pcov : covariance matrix  ->  perr = np.sqrt(np.diag(pcov))  gives std errors
```

다음 셀의 빈칸을 채워 물의 증기압 데이터에 앙투안 방정식을 적합하라.

```python
# Antoine equation: log10(P) = A - B/(C+T)  -> P = 10^(A - B/(C+T))
T_ant = np.array([40, 50, 60, 70, 80, 90, 100], dtype=float)    # degC
P_ant = np.array([55.3, 92.5, 149.4, 233.7, 355.1, 525.8, 760.0])  # mmHg

def antoine_model(T, A, B, C):
    return 10.0**(A - B / (C + T))

# --- COMPLETE: fit the Antoine equation ---
p0_ant = [8.0, 1500.0, 200.0]
popt_ant, pcov_ant =     
# --- COMPLETE ---

A_opt, B_opt, C_opt = popt_ant
perr_ant = np.sqrt(np.diag(pcov_ant))
gof_ant  = goodness_of_fit(P_ant, antoine_model(T_ant, *popt_ant), n_params=3)

print('=== Antoine Equation Fit (water vapor pressure) ===')
print(f'A = {A_opt:.6f} +/- {perr_ant[0]:.6f}')
print(f'B = {B_opt:.4f} +/- {perr_ant[1]:.4f}')
print(f'C = {C_opt:.4f} +/- {perr_ant[2]:.4f}')
print(f'R2 = {gof_ant["R2"]:.8f}')
print(f'Literature: A=8.07131, B=1730.63, C=233.426')

T_fit = np.linspace(35, 105, 300)
fig, ax = plt.subplots(figsize=(8, 4))
ax.scatter(T_ant, P_ant, c='k', marker='s', s=80, zorder=5, label='Data')
ax.plot(T_fit, antoine_model(T_fit, *popt_ant), 'r-', lw=2,
        label=f'Antoine fit $R^2={gof_ant["R2"]:.6f}$')
ax.set(xlabel='Temperature [°C]', ylabel='Vapor pressure [mmHg]',
       title='Antoine Equation: Vapor Pressure of Water')
ax.legend(); ax.grid(True, alpha=0.3)
plt.tight_layout(); plt.show()
```

배포 노트북에는 이 셀의 출력이 다음과 같이 실려 있다.

```text
=== Antoine Equation Fit (water vapor pressure) ===
A = 7.997574 +/- 0.007401
B = 1687.0509 +/- 4.6041
C = 229.7113 +/- 0.4237
R2 = 0.99999999
Literature: A=8.07131, B=1730.63, C=233.426
```

---

<br>

문제를 풀어 본 뒤 [과제 해설](../Assignment-Explanation.ko.md)에서 정답과 풀이를, [완성한 노트북](../Solution/Lab09-Regression_Completed.ipynb)에서 작성한 코드 전체를 볼 수 있다.
