# 제15장 강의 — 비선형 회귀(Nonlinear Regression)

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 15

> **선수 지식**: [통계학] 선형 회귀 (제14장).
>
> **학습 목표**:
> 1. 다항식 및 다중 선형 회귀를 적용할 수 있다
> 2. 비선형 회귀를 위한 가우스-뉴턴 방법을 구현할 수 있다
> 3. 모델 선택 기준을 비교할 수 있다

---

<br>

## 목차

- [1. 복습: 표준 오차와 결정 계수](#1-복습-표준-오차와-결정-계수)
  - [1.1 표본 표준 편차(SD)](#11-표본-표준-편차sd)
  - [1.2 평균의 표준 오차(SE)](#12-평균의-표준-오차se)
  - [1.3 회귀의 표준 오차](#13-회귀의-표준-오차)
  - [1.4 결정 계수(R-제곱)](#14-결정-계수r-제곱)
- [2. 비선형 관계의 선형화(Linearization)](#2-비선형-관계의-선형화linearization)
  - [2.1 지수 모델(Exponential Model)](#21-지수-모델exponential-model)
  - [2.2 거듭제곱 모델(Power Model)](#22-거듭제곱-모델power-model)
  - [2.3 포화 성장률 모델(Saturation-Growth-Rate Model)](#23-포화-성장률-모델saturation-growth-rate-model)
- [3. 다항 회귀(Polynomial Regression)](#3-다항-회귀polynomial-regression)
  - [3.1 문제 정의](#31-문제-정의)
  - [3.2 정규 방정식 유도](#32-정규-방정식-유도)
  - [3.3 행렬 형태의 정규 방정식](#33-행렬-형태의-정규-방정식)
  - [3.4 추정의 표준 오차](#34-추정의-표준-오차)
- [4. 다중 선형 회귀(Multiple Linear Regression)](#4-다중-선형-회귀multiple-linear-regression)
  - [4.1 문제 정의](#41-문제-정의)
  - [4.2 정규 방정식](#42-정규-방정식)
  - [4.3 추정의 표준 오차](#43-추정의-표준-오차)
  - [4.4 다변수 거듭제곱 모델로의 확장](#44-다변수-거듭제곱-모델로의-확장)
- [5. 일반 선형 최소제곱 회귀(General Linear Least Square Regression)](#5-일반-선형-최소제곱-회귀general-linear-least-square-regression)
  - [5.1 일반 모델 정의](#51-일반-모델-정의)
  - [5.2 행렬 정의](#52-행렬-정의)
  - [5.3 정규 방정식 유도(행렬 미적분)](#53-정규-방정식-유도행렬-미적분)
  - [5.4 해 계산](#54-해-계산)
  - [5.5 적합도 통계량](#55-적합도-통계량)
  - [5.6 매개변수 추정치의 공분산 행렬](#56-매개변수-추정치의-공분산-행렬)
- [6. 비선형 회귀(가우스-뉴턴, Gauss-Newton)](#6-비선형-회귀가우스-뉴턴-gauss-newton)
  - [6.1 동기: 앙투안 방정식(Antoine Equation)](#61-동기-앙투안-방정식antoine-equation)
  - [6.2 scipy.optimize를 이용한 Python 구현](#62-scipyoptimize를-이용한-python-구현)
- [7. 요약 표](#7-요약-표)

---

<br>

## 1. 복습: 표준 오차와 결정 계수

### 1.1 표본 표준 편차(SD)

**표본 표준 편차**(SD)는 개별 데이터 점이 데이터 세트의 평균으로부터 얼마나 변동(퍼짐)하는지를 측정한다.

$$
s_y = \sqrt{\frac{1}{n-1} \sum_{j=1}^{n} (y_j - \bar{y})^2}
$$

여기서:
- $n$ = 데이터 점의 수
- $y_j$ = 개별 관측값
- $\bar{y}$ = 표본 평균

**예제:** 키 = [170, 172, 168, 174, 176] cm

- 표본 평균 $\bar{y} = 172$ cm
- $s_y = 3.16$ cm

> **해석:** 평균적으로, 개별 키 측정값은 평균과 약 3.16 cm 차이가 난다.

### 1.2 평균의 표준 오차(SE)

**표준 오차**(SE)는 표본 평균이 참 모집단 평균을 얼마나 정밀하게 추정하는지를 측정한다.

$$
SE = \frac{s_y}{\sqrt{n}}
$$

**예제 (계속):**
- $SE = \frac{3.16}{\sqrt{5}} = 1.41$ cm

> **해석:** 평균 키는 172 cm로 추정되며, 표준 오차는 $\pm 1.41$ cm이다.

### 1.3 회귀의 표준 오차

회귀의 **표준 오차**(SE)는 관측값이 선형 회귀선으로부터 얼마나 벗어나는지를 측정한다.

$$
s_e = \sqrt{\frac{1}{n-2} \sum_{j=1}^{n} (y_j - \hat{y}_j)^2} = \sqrt{\frac{SSE}{n-2}}
$$

여기서:
- $\hat{y}_j$ = 회귀선에서의 예측값
- $n - 2$ = 자유도(두 매개변수 $\beta_0, \beta_1$이 추정됨)

**예제:** 키 = [170, 172, 168, 174, 176] cm, $x = [1, 2, 3, 4, 5]$

- $SE = \sqrt{\frac{SSE}{n-2}} = 2.61$ cm

> **해석:** 관측된 $y$ 값은 회귀선으로부터 약 2.61 cm 벗어난다.

### 1.4 결정 계수(R-제곱)

**결정 계수** $R^2$는 회귀 모델이 종속 변수의 변동성을 얼마나 잘 설명하는지를 측정한다.

총 제곱합은 다음과 같이 분해된다:

$$
SST = SSR + SSE
$$

양변을 $SST$로 나누면:

$$
1 = \frac{SSR}{SST} + \frac{SSE}{SST}
$$

따라서:

$$
R^2 = \frac{SSR}{SST} = 1 - \frac{SSE}{SST}
$$

여기서:
- $SST = \sum (y_i - \bar{y})^2$ — 총 제곱합(Total Sum of Squares)
- $SSR = \sum (\hat{y}_i - \bar{y})^2$ — 회귀 제곱합(설명된 변동)
- $SSE = \sum (y_i - \hat{y}_i)^2$ — 오차 제곱합(설명되지 않은 변동)

| $R^2$ 값 | 의미 |
|---|---|
| $R^2 = 1$ | 완벽한 적합 |
| $R^2 = 0$ | 적합 없음 |
| $0 < R^2 < 1$ | 부분적 적합 |

**예제 (ex 14.6):** $R^2 = 88\%$

> **해석:** 두 변수 간의 관계는 데이터 변동의 88%(약 90%)를 설명한다.

---

<br>

## 2. 비선형 관계의 선형화(Linearization)

핵심 아이디어는 **비선형 함수를 직선 형태로 변환** 하여 선형 회귀 기법을 적용할 수 있게 하는 것이다.

### 2.1 지수 모델(Exponential Model)

**원래 모델:**

$$
y = \alpha_1 e^{\beta_1 x}
$$

**선형화:** 양변에 자연로그를 취한다:

$$
\ln y = \ln(\alpha_1 e^{\beta_1 x}) = \ln \alpha_1 + \ln e^{\beta_1 x} = \ln \alpha_1 + \beta_1 x
$$

이것은 선형 형태를 준다: $\ln y = \underbrace{\ln \alpha_1}_{\text{절편}} + \underbrace{\beta_1}_{\text{기울기}} \cdot x$

> **그래프 효과:** $y$ vs $x$ 그래프는 지수 곡선을 보인다. $\ln y$ vs $x$ 그래프는 기울기 $\beta_1$, 절편 $\ln \alpha_1$인 직선이 된다.

### 2.2 거듭제곱 모델(Power Model)

**원래 모델:**

$$
y = \alpha_2 x^{\beta_2}
$$

**선형화:** 양변에 자연로그를 취한다:

$$
\ln y = \ln(\alpha_2 x^{\beta_2}) = \ln \alpha_2 + \ln x^{\beta_2} = \ln \alpha_2 + \beta_2 \ln x
$$

이것은 선형 형태를 준다: $\ln y = \underbrace{\ln \alpha_2}_{\text{절편}} + \underbrace{\beta_2}_{\text{기울기}} \cdot \ln x$

> **그래프 효과:** $y$ vs $x$ 그래프는 거듭제곱 곡선을 보인다. $\ln y$ vs $\ln x$ 그래프는 기울기 $\beta_2$, 절편 $\ln \alpha_2$인 직선이 된다.

### 2.3 포화 성장률 모델(Saturation-Growth-Rate Model)

**원래 모델:**

$$
y = \frac{\alpha_3 x}{\beta_3 + x}
$$

**선형화:** 양변의 역수를 취한다:

$$
\frac{1}{y} = \frac{\beta_3 + x}{\alpha_3 x} = \frac{\beta_3}{\alpha_3} \cdot \frac{1}{x} + \frac{1}{\alpha_3}
$$

이것은 선형 형태를 준다: $\frac{1}{y} = \underbrace{\frac{\beta_3}{\alpha_3}}_{\text{기울기}} \cdot \frac{1}{x} + \underbrace{\frac{1}{\alpha_3}}_{\text{절편}}$

> **그래프 효과:** $y$ vs $x$ 그래프는 포화 곡선을 보인다. $\frac{1}{y}$ vs $\frac{1}{x}$ 그래프는 기울기 $\frac{\beta_3}{\alpha_3}$, 절편 $\frac{1}{\alpha_3}$인 직선이 된다.

| 모델 | 원래 형태 | 변환된 변수 | 기울기 | 절편 |
|---|---|---|---|---|
| 지수 모델 | $y = \alpha_1 e^{\beta_1 x}$ | $\ln y$ vs $x$ | $\beta_1$ | $\ln \alpha_1$ |
| 거듭제곱 모델 | $y = \alpha_2 x^{\beta_2}$ | $\ln y$ vs $\ln x$ | $\beta_2$ | $\ln \alpha_2$ |
| 포화 성장 모델 | $y = \frac{\alpha_3 x}{\beta_3 + x}$ | $\frac{1}{y}$ vs $\frac{1}{x}$ | $\frac{\beta_3}{\alpha_3}$ | $\frac{1}{\alpha_3}$ |

---

<br>

## 3. 다항 회귀(Polynomial Regression)

### 3.1 문제 정의

데이터가 곡선 형태(예: 포물선)를 보일 때, 직선만으로는 충분하지 않다. 선형 회귀의 아이디어를 **다항 회귀** 로 확장한다.

**목표:**
- 선형 회귀 아이디어를 다항 회귀로 확장
- 신뢰 구간 추정

데이터 집합 $\{x_i, y_i\}_{i=1}^{n}$이 주어졌을 때, 다음을 만족하는 $\beta_0, \beta_1, \beta_2$를 구한다:

$$
\hat{y} = \beta_0 + \beta_1 x + \beta_2 x^2
$$

목적은 오차 제곱합을 최소화하는 것이다:

$$
\min_{\beta_0, \beta_1, \beta_2} SSE
$$

여기서:

$$
SSE = \sum_{i=1}^{n} (y_i - \hat{y}_i)^2 = \sum_{i=1}^{n} (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)^2
$$

### 3.2 정규 방정식 유도

$\beta_0, \beta_1, \beta_2$를 결정하기 위해 편도함수를 취하고 0으로 놓는다:

$$
\frac{\partial SSE}{\partial \beta_0} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)(-1) = 0
$$

$$
\frac{\partial SSE}{\partial \beta_1} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)(-x_i) = 0
$$

$$
\frac{\partial SSE}{\partial \beta_2} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_i - \beta_2 x_i^2)(-x_i^2) = 0
$$

정리하면, **정규 방정식** 은 다음과 같다:

$$
\sum y_i = \beta_0 \cdot n + \beta_1 \sum x_i + \beta_2 \sum x_i^2
$$

$$
\sum y_i x_i = \beta_0 \sum x_i + \beta_1 \sum x_i^2 + \beta_2 \sum x_i^3
$$

$$
\sum y_i x_i^2 = \beta_0 \sum x_i^2 + \beta_1 \sum x_i^3 + \beta_2 \sum x_i^4
$$

### 3.3 행렬 형태의 정규 방정식

정규 방정식은 다음과 같이 간결하게 표현할 수 있다:

$$
\begin{bmatrix} \sum y_i \\ \sum y_i x_i \\ \sum y_i x_i^2 \end{bmatrix}
=
\begin{bmatrix} n & \sum x_i & \sum x_i^2 \\ \sum x_i & \sum x_i^2 & \sum x_i^3 \\ \sum x_i^2 & \sum x_i^3 & \sum x_i^4 \end{bmatrix}
\begin{bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \end{bmatrix}
$$

> **핵심 관찰:** 이것은 Part 3의 방법(예: 가우스 소거법, LU 분해)으로 풀 수 있는 $3 \times 3$ 선형 시스템이다.

### 3.4 추정의 표준 오차

3개의 매개변수($\beta_0, \beta_1, \beta_2$)를 갖는 다항 회귀의 경우:

$$
e_i = y_i - \hat{y}_i
$$

$$
SSE = \sum e_i^2
$$

$$
s_e = \sqrt{\frac{SSE}{n - 3}}
$$

> **자유도:** 분모가 $n - 3$인 이유는 세 개의 매개변수($\beta_0, \beta_1, \beta_2$)가 $\hat{y}$를 구성하는 데 사용되었기 때문이다. 일반적으로 $m$차 다항식의 자유도는 $n - (m+1)$이다.

---

<br>

## 4. 다중 선형 회귀(Multiple Linear Regression)

다중 선형 회귀는 1차원에서 2차원(및 그 이상)으로 개념을 확장한다: **두 개 이상의 독립 변수** 를 사용하여 데이터를 적합한다.

### 4.1 문제 정의

모델은 데이터에 평면(또는 초평면)을 적합한다:

$$
\hat{y} = \beta_0 + \beta_1 x_1 + \beta_2 x_2
$$

목적은:

$$
SSE = \sum e_i^2 = \sum (y_i - \hat{y}_i)^2 = \sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})^2
$$

### 4.2 정규 방정식

$\min_{\beta_0, \beta_1, \beta_2} SSE$를 만족하는 $\beta_0, \beta_1, \beta_2$를 구한다:

$$
\frac{\partial SSE}{\partial \beta_0} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})(-1) = 0
$$

$$
\frac{\partial SSE}{\partial \beta_1} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})(-x_{1i}) = 0
$$

$$
\frac{\partial SSE}{\partial \beta_2} = 0 \implies 2\sum (y_i - \beta_0 - \beta_1 x_{1i} - \beta_2 x_{2i})(-x_{2i}) = 0
$$

행렬 형태:

$$
\begin{bmatrix} \sum y_i \\ \sum y_i x_{1i} \\ \sum y_i x_{2i} \end{bmatrix}
=
\begin{bmatrix} n & \sum x_{1i} & \sum x_{2i} \\ \sum x_{1i} & \sum x_{1i}^2 & \sum x_{1i} x_{2i} \\ \sum x_{2i} & \sum x_{1i} x_{2i} & \sum x_{2i}^2 \end{bmatrix}
\begin{bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \end{bmatrix}
$$

### 4.3 추정의 표준 오차

두 개의 독립 변수(3개의 매개변수)를 갖는 다중 선형 회귀의 경우:

$$
s_e = \sqrt{\frac{SSE}{n - 3}}
$$

### 4.4 다변수 거듭제곱 모델로의 확장

다중 선형 회귀는 **다변수 거듭제곱 모델** 로 확장될 수 있다:

$$
y = \beta_0 x_1^{\beta_1} x_2^{\beta_2} \cdots x_m^{\beta_m} + \varepsilon
$$

자연로그를 취하면:

$$
\ln y = \ln \beta_0 + \beta_1 \ln x_1 + \beta_2 \ln x_2 + \cdots + \beta_m \ln x_m
$$

> **핵심 통찰:** 로그 변환 후, 이것은 변환된 변수 $\ln y$, $\ln x_1$, $\ln x_2$, ..., $\ln x_m$에 대한 표준 다중 선형 회귀 문제가 된다.

---

<br>

## 5. 일반 선형 최소제곱 회귀(General Linear Least Square Regression)

### 5.1 일반 모델 정의

일반 선형 최소제곱 모델은 이전의 모든 경우를 통합한다. 모델은 매개변수 $\beta_j$에 대해 "선형"이지만, 기저 함수 $f_j$는 임의적(예: 삼각 함수)일 수 있다:

**예제 — 푸리에 유형 모델:**

$$
y = \beta_0 + \beta_1 \cos(\omega x) + \beta_2 \sin(\omega x)
$$

**일반 형태:**

$$
y = \beta_0 f_0(x_1, x_2, \ldots, x_p) + \beta_1 f_1(x_1, x_2, \ldots, x_p) + \cdots + \beta_m f_m(x_1, x_2, \ldots, x_p)
$$

> **중요:** 여기서 "선형"은 $x$가 아니라 **매개변수** $\beta_j$에 대해 선형이라는 의미이다. 함수 $f_j$는 $x$에 대해 비선형일 수 있다.

### 5.2 행렬 정의

$n$개의 데이터 점 $\{y_i, x_{1i}, x_{2i}, \ldots, x_{pi}\}_{i=1}^{n}$을 수집한다.

각 데이터 점 $i = 1, 2, \ldots, n$에 대해:

$$
y_i = \beta_0 f_0(x_{1i}, x_{2i}, \ldots, x_{pi}) + \beta_1 f_1(x_{1i}, x_{2i}, \ldots, x_{pi}) + \cdots + \beta_m f_m(x_{1i}, x_{2i}, \ldots, x_{pi})
$$

행렬 형태:

$$
\underset{(n \times 1)}{\mathbf{y}}
=
\underset{(n \times (m+1))}{\mathbf{X}}
\underset{((m+1) \times 1)}{\boldsymbol{\beta}}
+ \underset{(n \times 1)}{\mathbf{e}}
$$

여기서:

$$
\mathbf{y} = \begin{bmatrix} y_1 \\ y_2 \\ \vdots \\ y_n \end{bmatrix}, \quad
\mathbf{X} = \begin{bmatrix} f_{01} & f_{11} & f_{21} & \cdots & f_{m1} \\ f_{02} & f_{12} & f_{22} & \cdots & f_{m2} \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ f_{0n} & f_{1n} & f_{2n} & \cdots & f_{mn} \end{bmatrix}, \quad
\boldsymbol{\beta} = \begin{bmatrix} \beta_0 \\ \beta_1 \\ \beta_2 \\ \vdots \\ \beta_m \end{bmatrix}, \quad
\mathbf{e} = \begin{bmatrix} e_1 \\ e_2 \\ \vdots \\ e_n \end{bmatrix}
$$

여기서 $f_{ji}$는 데이터 점 $i$에서 평가된 기저 함수 $f_j$의 값이고, $\mathbf{e}$는 잔차(오차) 벡터이다.

> **표기법:** $\mathbf{X}$는 $n$개의 행(데이터 점)과 $m+1$개의 열(매개변수 수, $\beta_0$ 포함)을 갖는다.

### 5.3 정규 방정식 유도(행렬 미적분)

행렬 형태의 오차 제곱합:

$$
SSE = \sum e_i^2 = \mathbf{e}^T \mathbf{e} = (\mathbf{y} - \hat{\mathbf{y}})^T (\mathbf{y} - \hat{\mathbf{y}})
$$

최소화하기 위해 $\boldsymbol{\beta}$에 대한 도함수를 취하고 0으로 놓는다:

$$
\frac{\partial SSE}{\partial \boldsymbol{\beta}} = \frac{\partial}{\partial \boldsymbol{\beta}} (\mathbf{e}^T \mathbf{e})
= \frac{\partial}{\partial \boldsymbol{\beta}} (\mathbf{y} - \mathbf{X}\boldsymbol{\beta})^T (\mathbf{y} - \mathbf{X}\boldsymbol{\beta})
= (-2\mathbf{X})^T (\mathbf{y} - \mathbf{X}\boldsymbol{\beta}) = \mathbf{0}
$$

이로부터 **정규 방정식** 을 얻는다:

$$
\mathbf{X}^T \mathbf{X} \boldsymbol{\beta} = \mathbf{X}^T \mathbf{y}
$$

$\boldsymbol{\beta}$에 대해 풀면:

$$
\boxed{\boldsymbol{\beta} = (\mathbf{X}^T \mathbf{X})^{-1} \mathbf{X}^T \mathbf{y}}
$$

> **핵심 공식:** 이것은 이 장의 가장 중요한 결과이다. 임의의 일반 선형 모델에 적용 가능한, 오차 제곱합을 최소화하는 매개변수 벡터 $\boldsymbol{\beta}$에 대한 닫힌 형태의 해를 제공한다.

### 5.4 해 계산

$\boldsymbol{\beta}$가 계산되면:

| 양 | 수식 |
|---|---|
| 잔차 벡터 | $\mathbf{e} = \mathbf{y} - \mathbf{X}\boldsymbol{\beta}$ |
| 오차 제곱합 | $SSE = \mathbf{e}^T \mathbf{e}$ |
| 예측값 | $\hat{\mathbf{y}} = \mathbf{X}\boldsymbol{\beta}$ |

### 5.5 적합도 통계량

$$
SST = \sum (y_i - \bar{y})^2
$$

$$
SSR = \sum (\hat{y}_i - \bar{y})^2
$$

$$
R^2 = \frac{SSR}{SST}
$$

$$
SSE = SST - SSR
$$

**추정의 표준 오차** ($m+1$개의 매개변수를 갖는 일반적인 경우):

$$
s_e = \sqrt{\frac{SSE}{n - (m+1)}}
$$

> **자유도:** 분모 $n - (m+1)$은 모델 적합에 사용된 $m+1$개의 매개변수를 고려한다.

### 5.6 매개변수 추정치의 공분산 행렬

추정된 매개변수의 공분산 행렬:

$$
\text{cov}(\hat{\boldsymbol{\beta}}) = \sigma_e^2 (\mathbf{X}^T \mathbf{X})^{-1}
$$

여기서 $\sigma_e^2$는 $s_e^2 = \frac{SSE}{n-(m+1)}$로 추정된다.

> **실용적 활용:** $\text{cov}(\hat{\boldsymbol{\beta}})$의 대각 원소는 각 개별 매개변수 추정치의 분산을 제공하며, 그 제곱근은 매개변수의 표준 오차를 제공한다. 이를 사용하여 신뢰 구간을 구성할 수 있다.

---

<br>

## 6. 비선형 회귀(가우스-뉴턴, Gauss-Newton)

### 6.1 동기: 앙투안 방정식(Antoine Equation)

일부 모델은 매개변수에 대해 **본질적으로 비선형** 이며 변환으로 선형화할 수 없다. 예를 들어, 증기압에 대한 **앙투안 방정식(Antoine Equation)**:

$$
\log_{10} P_v = A - \frac{B}{C + T}
$$

여기서:
- $P_v$ = 증기압
- $T$ = 온도
- $A, B, C$ = 조정 가능한 매개변수

$C$가 비선형 표현(분모의 $C + T$) 안에 나타나므로, 이 방정식은 일반 선형 모델로 표현할 수 **없다**. 이를 적합하려면 비선형 회귀(최적화)가 필요하다.

### 6.2 scipy.optimize를 이용한 Python 구현

```python
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# Sample vapor pressure data (e.g., for water) in mmHg at different temperatures (deg C)
T_data = np.array([40, 50, 60, 70, 80])    # Temperature in Celsius
P_data = np.array([55.3, 92.5, 149.4, 233.7, 355.1])  # Vapor pressure in mmHg

# Convert pressure to log base 10
logP_data = np.log10(P_data)

# Antoine equation model: log10(P) = A - B / (C + T)
def antoine_model(params, T):
    A, B, C = params
    return A - B / (C + T)

# Objective function: sum of squared errors between predicted and actual logP
def objective(params):
    predictions = antoine_model(params, T_data)
    return np.sum((logP_data - predictions) ** 2)

# Initial guess for A, B, C
initial_guess = [8, 1500, 200]

# Minimize the objective function
result = minimize(objective, initial_guess, method='L-BFGS-B')

# Extract optimized parameters
A_opt, B_opt, C_opt = result.x
```

> **방법:** `scipy.optimize.minimize` 함수의 `method='L-BFGS-B'`는 유사 뉴턴(Quasi-Newton) 최적화 알고리즘을 사용한다. 초기 추정값이 필요하며, 알고리즘이 $SSE$를 최소화하도록 매개변수를 반복적으로 개선한다.

> **초기 추정값 민감도:** 비선형 회귀 방법은 초기 추정값에 의존한다. 부적절한 초기 추정값은 전역 최솟값이 아닌 지역 최솟값으로의 수렴을 유발할 수 있다.

---

<br>

## 7. 요약 표

| 주제 | 모델 | 매개변수 | 정규 방정식 / 방법 | $s_e$ |
|---|---|---|---|---|
| **선형 회귀** | $\hat{y} = \beta_0 + \beta_1 x$ | 2 | 닫힌 형태 (제14장) | $\sqrt{\frac{SSE}{n-2}}$ |
| **다항 회귀** | $\hat{y} = \beta_0 + \beta_1 x + \beta_2 x^2$ | 3 | $3 \times 3$ 정규 방정식 | $\sqrt{\frac{SSE}{n-3}}$ |
| **다중 선형 회귀** | $\hat{y} = \beta_0 + \beta_1 x_1 + \beta_2 x_2$ | 3 | $3 \times 3$ 정규 방정식 | $\sqrt{\frac{SSE}{n-3}}$ |
| **일반 선형 최소제곱** | $\hat{y} = \sum \beta_j f_j(\mathbf{x})$ | $m+1$ | $\boldsymbol{\beta} = (\mathbf{X}^T\mathbf{X})^{-1}\mathbf{X}^T\mathbf{y}$ | $\sqrt{\frac{SSE}{n-(m+1)}}$ |
| **선형화 모델** | 지수 / 거듭제곱 / 포화 | 2 | 변환 후 선형 회귀 | 변환에 따라 다름 |
| **비선형 회귀** | 예: 앙투안 방정식 | 다양 | 반복 최적화 (예: L-BFGS-B) | $\sqrt{\frac{SSE}{n-p}}$ |
