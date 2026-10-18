# 14장 — 선형 회귀(Linear Regression)

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **최종 수정일:** 2026-10-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 14

> **학습 목표:**
> 1. 평균, 중위수, 표준 편차, 변동 계수, MAD 같은 기술 통계량을 계산할 수 있다
> 2. 정규 분포의 성질과 표준화를 설명할 수 있다
> 3. 최소제곱 기준으로 정규 방정식을 유도하고 단순 선형 회귀의 기울기와 절편을 계산할 수 있다
> 4. 결정 계수($`r^2`$)와 추정의 표준 오차를 사용하여 회귀 품질을 평가할 수 있다

---

<br>

## 목차

- [1. 개요](#1-개요)
- [2. 통계 복습](#2-통계-복습)
  - [2.1 기술 통계(Descriptive Statistics)](#21-기술-통계descriptive-statistics)
  - [2.2 표준 편차와 분산](#22-표준-편차와-분산)
  - [2.3 변동 계수(Coefficient of Variation)](#23-변동-계수coefficient-of-variation)
  - [2.4 중위 절대 편차(Median Absolute Deviation, MAD)](#24-중위-절대-편차median-absolute-deviation-mad)
- [3. 정규 분포(Normal Distribution)](#3-정규-분포normal-distribution)
  - [3.1 중심극한정리(Central Limit Theorem)](#31-중심극한정리central-limit-theorem)
  - [3.2 표준 정규 분포(Standard Normal Distribution)](#32-표준-정규-분포standard-normal-distribution)
  - [3.3 대칭 성질](#33-대칭-성질)
  - [3.4 일반 정규 분포(General Normal Distribution)](#34-일반-정규-분포general-normal-distribution)
  - [3.5 표준화(Standardization)](#35-표준화standardization)
  - [3.6 정규 분포의 CDF와 PDF](#36-정규-분포의-cdf와-pdf)
- [4. 직선 최소제곱 회귀(Straight Line Least Square Regression)](#4-직선-최소제곱-회귀straight-line-least-square-regression)
  - [4.1 문제 설정](#41-문제-설정)
  - [4.2 대안적 오차 기준과 최소제곱법의 이유](#42-대안적-오차-기준과-최소제곱법의-이유)
  - [4.3 정규 방정식 유도](#43-정규-방정식-유도)
  - [4.4 계수 풀이](#44-계수-풀이)
- [5. 적합도 정량화](#5-적합도-정량화)
  - [5.1 SST, SSE, SSR 분해](#51-sst-sse-ssr-분해)
  - [5.2 SST = SSE + SSR 증명](#52-sst--sse--ssr-증명)
  - [5.3 결정 계수(Coefficient of Determination, $`r^2`$)](#53-결정-계수coefficient-of-determination-r2)
  - [5.4 추정의 표준 오차(Standard Error of the Estimate)](#54-추정의-표준-오차standard-error-of-the-estimate)
- [요약](#요약)
- [점검 문제](#점검-문제)

<br>

---

<br>

## 1. 개요

선형 회귀(Linear Regression)는 **곡선 맞춤(Curve Fitting)**(교재 Part IV, 14~18장)의 근본적인 기법이다. 목표는 **최소제곱(Least-Squares)** 기준을 사용하여 데이터 점들을 관통하는 **최적 직선(Best-Fit Straight Line)** 을 유도하는 것이다.

이 장은 기본 통계(평균, 표준 편차, 분포)를 복습한 뒤, 최적 직선의 기울기와 절편을 계산하고, 추정의 표준 오차를 계산하여 잔차 오차를 분석한다.

<br>

---

<br>

## 2. 통계 복습

### 2.1 기술 통계(Descriptive Statistics)

**통계학(Statistics)** 은 의미 있는 정보를 추론하기 위해 데이터를 다루는 학문이다.

**표본 추출(Sampling)** 은 전체 집단에 대한 통찰을 얻기 위해 큰 모집단에서 작은 그룹을 선택하는 과정이다.

#### 표본 평균(Sample Average, Mean)

$`n`$개의 관측값 $`y_1, y_2, \ldots, y_n`$이 있을 때, **표본 평균**은 다음과 같다:

```math
\bar{y} = \frac{1}{n} \sum_{j=1}^{n} y_j
```

#### 중위수(Median)

**중위수** $`\tilde{y}`$는 데이터를 정렬했을 때 가운데 값이다. 데이터가 짝수 개이면 가운데 두 값의 평균이다.

- 홀수 개: 데이터 [3, 5, 7]의 중위수는 5이다.
- 짝수 개: 데이터 [2, 4, 6, 8]의 중위수는 (4 + 6) / 2 = 5이다.

#### 최빈값(Mode)

**최빈값**은 데이터에서 가장 자주 나타나는 값이다.

- 데이터 [2, 4, 4, 4, 5]의 최빈값은 4이다.

#### 이상값(Outlier)

**이상값**은 나머지 데이터와 현저히 다른 데이터 점이다.

- 데이터 [2, 4, 4, 4, 100]에서 100은 이상값이다.

<br>

### 2.2 표준 편차와 분산

#### 모표준편차(Population Standard Deviation)

**표준 편차**는 데이터 점들이 평균으로부터 얼마나 퍼져 있는지를 측정한다:

```math
\sigma = \sqrt{\frac{1}{n} \sum_{j=1}^{n} (y_j - \mu)^2}
```

여기서 $`\mu`$는 전체 모집단의 **참 평균(True Mean)** 이고, $`n`$은 모집단의 총 수이다.

> **[확률및통계]** 모표준편차는 분모에 $`n`$을 사용하는데, 이는 모집단의 모든 구성원이 포함되어 자유도 손실이 없기 때문이다.

#### 표본 표준 편차(Sample Standard Deviation)

```math
s_y = \sqrt{\frac{1}{n-1} \sum_{j=1}^{n} (y_j - \bar{y})^2}
```

여기서 $`\bar{y}`$는 **표본 평균**이고 $`n`$은 표본 내 데이터 점의 총 수이다.

제곱근 내 분자에 해당하는 양을 **총 보정 제곱합(Total Corrected Sum of Squares, SST)** 이라 한다:

```math
SST = \sum_{j=1}^{n} (y_j - \bar{y})^2
```

> **[확률및통계]** 왜 $`n`$이 아니라 $`n - 1`$인가? 참 평균 $`\mu`$ 대신 표본 평균 $`\bar{y}`$를 사용하면 **변동성을 과소추정** 한다. $`\bar{y}`$는 이미 자유도 1을 사용한 확률 변수이므로, 모분산의 **비편향(Unbiased)** 추정치를 얻기 위해 $`n - 1`$로 나눈다(베셀 보정, Bessel's Correction).

<br>

### 2.3 변동 계수(Coefficient of Variation)

**변동 계수(c.v.)** 는 표준 편차 대 평균의 비율로, 백분율로 표현된다:

```math
c.v. = \frac{s_y}{\bar{y}} \cdot 100\%
```

> **[확률및통계]** 변동 계수는 산포도를 정규화하여, **스케일이 다른 데이터 세트도 비교할 수 있게** 한다. 낮은 c.v.는 데이터가 평균 주위에 밀집되어 있음을, 높은 c.v.는 평균에 비해 상당한 변동이 있음을 의미한다.

<br>

### 2.4 중위 절대 편차(Median Absolute Deviation, MAD)

**중위 절대 편차(MAD)** 는 값들이 중위수로부터 얼마나 벗어나는지를 측정한다:

```math
MAD = \text{median}(|y_i - \tilde{y}|)
```

**예제:**

- 데이터 [2, 4, 6, 8, 100]의 중위수는 6이다.
- 절대 편차는 |2-6| = 4, |4-6| = 2, |6-6| = 0, |8-6| = 2, |100-6| = 94이다.
- 편차를 정렬하면 [0, 2, 2, 4, 94]이므로 MAD = 2이다.

> **[확률및통계]** MAD는 **이상값에 강건(Robust)** 하다. 위 예제에서 큰 이상값 100은 94라는 편차를 만들지만, MAD는 편차의 중위수를 취하므로 이 극단값의 크기에 영향을 받지 않는다. 따라서 MAD는 이상값이 있을 때 표준 편차보다 더 신뢰할 수 있는 산포 측정치이다.

<br>

---

<br>

## 3. 정규 분포(Normal Distribution)

### 3.1 중심극한정리(Central Limit Theorem)

각각 유한한 평균과 분산을 가지는 많은 독립 확률 변수를 더하면, 많은 수의 i.i.d. 확률 변수의 합은 개별 확률 변수의 분포가 무엇이든 **근사적으로 정규 분포를 따른다**. 정확히는 표준화한 합의 분포가 표준 정규 분포로 수렴한다.

> **[확률및통계]** 중심극한정리(CLT)는 정규 분포가 널리 나타나는 근본적인 이유이다. 많은 작은 독립적 효과의 합인 측정 오차는 일반적으로 정규 분포를 따른다.

<br>

### 3.2 표준 정규 분포(Standard Normal Distribution)

**표준 정규 분포**는 중심극한정리 덕분에 매우 널리 사용된다.

**PDF (확률 밀도 함수, Probability Density Function):**

```math
\varphi(z) = \frac{1}{\sqrt{2\pi}} e^{-\frac{z^2}{2}}, \quad \text{for } -\infty < z < \infty
```

**CDF (누적 분포 함수, Cumulative Distribution Function):**

```math
\Phi(z) = \int_{-\infty}^{z} \varphi(t) \, dt, \quad \text{for } -\infty < z < \infty
```

이를 $`Z \sim \mathcal{N}(0, 1)`$로 표기하며, **평균 0, 분산 1**을 의미한다.

<br>

### 3.3 대칭 성질

**1. PDF의 대칭성:**

```math
\varphi(z) = \varphi(-z)
```

종 모양 곡선은 $`z = 0`$에 대해 대칭이다.

**2. 꼬리 면적의 대칭성:**

```math
P(Z \le -2) = P(Z \ge 2)
```

이는 다음과 동치이다:

```math
\Phi(-z) = 1 - \Phi(z)
```

**3. $`Z`$와 $`-Z`$의 대칭성:**

$`Z \sim \mathcal{N}(0, 1)`$이면, $`-Z \sim \mathcal{N}(0, 1)`$이기도 하다.

```math
P(-Z \le z) = P(Z \ge -z) = 1 - \Phi(-z) = \Phi(z)
```

<br>

### 3.4 일반 정규 분포(General Normal Distribution)

$`Z \sim \mathcal{N}(0, 1)`$이면:

```math
X = \mu + \sigma Z
```

는 **평균** $`\mu`$, **분산** $`\sigma^2`$인 정규 분포를 따르며, $`X \sim \mathcal{N}(\mu, \sigma^2)`$로 표기한다.

**평균과 분산의 증명:**

```math
E(\mu + \sigma Z) = E(\mu) + E(\sigma Z) = \mu
```

```math
Var(\mu + \sigma Z) = Var(\sigma Z) = \sigma^2 Var(Z) = \sigma^2
```

> **[확률및통계]** 기댓값의 선형성($`E[aX + b] = aE[X] + b`$)과 분산의 스케일링 성질($`Var(aX + b) = a^2 Var(X)`$)이 여기서 사용되는 핵심 성질이다. 상수를 더하면 평균은 이동하지만 분산은 변하지 않는다.

<br>

### 3.5 표준화(Standardization)

$`X \sim \mathcal{N}(\mu, \sigma^2)`$일 때, $`X`$의 **표준화**는 다음과 같다:

```math
Z = \frac{X - \mu}{\sigma} \sim \mathcal{N}(0, 1)
```

이 변환은 임의의 정규 확률 변수를 표준 정규로 변환하여, 표준 정규 분포표를 사용할 수 있게 한다.

<br>

### 3.6 정규 분포의 CDF와 PDF

$`X \sim \mathcal{N}(\mu, \sigma^2)`$라 하자.

**$`X`$의 CDF:**

```math
F(x) = P(X \le x) = P\left(\frac{X - \mu}{\sigma} \le \frac{x - \mu}{\sigma}\right) = \Phi\left(\frac{x - \mu}{\sigma}\right)
```

**$`X`$의 PDF:**

```math
f(x) = F'(x) = \varphi\left(\frac{x - \mu}{\sigma}\right) \cdot \frac{1}{\sigma} = \frac{1}{\sigma\sqrt{2\pi}} e^{-\frac{(x - \mu)^2}{2\sigma^2}}
```

> **[확률및통계]** PDF에서의 $`\frac{1}{\sigma}`$ 인자는 CDF를 미분할 때 연쇄 법칙(Chain Rule)에서 나온다. 이는 $`\sigma`$에 관계없이 PDF 아래 전체 면적이 1로 유지되도록 보장한다.

<br>

---

<br>

## 4. 직선 최소제곱 회귀(Straight Line Least Square Regression)

### 4.1 문제 설정

데이터 점 $`(x_1, y_1), (x_2, y_2), \ldots, (x_n, y_n)`$이 주어졌을 때, **선형 모델**을 구하고자 한다:

```math
\hat{y} = a_0 + a_1 x
```

여기서:
- $`a_0`$ = **y절편(y-intercept)**
- $`a_1`$ = **기울기(slope)**

각 관측값은 다음과 같이 표현할 수 있다:

```math
y_i = \hat{y}_i + e_i = a_0 + a_1 x_i + e_i
```

여기서 $`e_i = y_i - \hat{y}_i`$는 **잔차(Residual, Error)** 로, 관측값과 예측값 사이의 차이이다.

<br>

### 4.2 대안적 오차 기준과 최소제곱법의 이유

오차의 어떤 측도를 최소화하여 데이터를 관통하는 최적의 직선을 구하고자 한다. 여러 기준을 고려할 수 있다:

**기준 1 (오차의 합 최소화):**

```math
\min_{a_0, a_1} \sum e_i = \min_{a_0, a_1} \sum (y_i - a_0 - a_1 x_i)
```

문제: 양의 오차와 음의 오차가 **서로 상쇄** 되어 오도된 결과를 준다.

**기준 2 (절대 오차의 합 최소화):**

```math
\min_{a_0, a_1} \sum |e_i| = \min_{a_0, a_1} \sum |y_i - a_0 - a_1 x_i|
```

문제: 절대값 함수는 0에서 **미분 불가능** 하고, 해가 **유일하지 않을 수 있다**.

**기준 3 (미니맥스, Minimax, 최대 오차 최소화):**

```math
\min_{a_0, a_1} \max_i |y_i - a_0 - a_1 x_i|
```

문제: **미분 불가능** 하고, 단일 최대 오차에 집중하므로 **이상값에 민감** 하다.

**기준 4 (오차 제곱합 최소화, 최소제곱법, Least Squares):**

```math
\min_{a_0, a_1} \sum e_i^2 = \min_{a_0, a_1} \sum (y_i - a_0 - a_1 x_i)^2
```

이것이 **최소제곱(Least Squares)** 접근법이며, **오차 제곱합(Sum of Squared Errors, SSE)** 을 산출한다. 미분 가능하고, $`x_i`$가 모두 같지 않으면 주어진 데이터에 대해 **유일한** 직선을 산출한다.

> **[확률및통계]** 최소제곱 기준이 선호되는 이유: (1) 제곱은 부호 문제를 제거하고, (2) 목적 함수가 매끄럽고 미분 가능하며, (3) 유일한 닫힌 형태의 해를 가지고, (4) 오차가 정규 분포를 따른다는 가정 하에 최대 우도 추정(Maximum Likelihood Estimate)에 해당한다.

<br>

### 4.3 정규 방정식 유도

최소화할 목적 함수는 다음과 같다:

```math
SSE = \sum_{i=1}^{n} (y_i - \hat{y}_i)^2 = \sum_{i=1}^{n} (y_i - a_0 - a_1 x_i)^2
```

최솟값을 찾기 위해 SSE를 각 미지 계수에 대해 미분하고 0으로 놓는다:

```math
\frac{\partial SSE}{\partial a_0} = -2 \sum_{i=1}^{n} (y_i - a_0 - a_1 x_i) = 0
```

```math
\frac{\partial SSE}{\partial a_1} = -2 \sum_{i=1}^{n} (y_i - a_0 - a_1 x_i) x_i = 0
```

$`a_0`$와 $`a_1`$은 상수($`i`$에 의존하지 않음)이므로 합에서 인수 분해할 수 있다:

- $`\sum a_0 = a_0 \cdot n`$
- $`\sum a_0 x_i = a_0 \sum x_i`$
- $`\sum a_1 x_i = a_1 \sum x_i`$
- $`\sum a_1 x_i^2 = a_1 \sum x_i^2`$

편도함수를 0으로 놓으면 **정규 방정식(Normal Equations)** 을 얻는다:

```math
\sum y_i = n \cdot a_0 + a_1 \sum x_i
```

```math
\sum y_i x_i = a_0 \sum x_i + a_1 \sum x_i^2
```

**행렬 형태:**

```math
\begin{bmatrix} \sum y_i \\ \sum y_i x_i \end{bmatrix} = \begin{bmatrix} n & \sum x_i \\ \sum x_i & \sum x_i^2 \end{bmatrix} \begin{bmatrix} a_0 \\ a_1 \end{bmatrix} \quad (\ast)
```

> **[선형대수학]** 이것은 $`2 \times 2`$ 선형 시스템이다. 계수 행렬은 대칭이고 양의 정부호(Positive Definite)이므로($`x_i`$ 값이 모두 같지 않은 한), 유일한 해가 보장된다.

<br>

### 4.4 계수 풀이

#### 크래머 법칙(Cramer's Rule) 사용

계수 행렬의 행렬식:

```math
D = n \sum x_i^2 - \left(\sum x_i\right)^2
```

절편:

```math
a_0 = \frac{1}{D} \left( \sum y_i \cdot \sum x_i^2 - \sum x_i \cdot \sum y_i x_i \right)
```

기울기:

```math
a_1 = \frac{1}{D} \left( n \sum y_i x_i - \sum x_i \cdot \sum y_i \right)
```

#### 표본 평균을 사용한 대안적 형태

$`(\ast)`$의 첫 번째 행으로부터:

```math
\sum y_i = n \cdot a_0 + \left(\sum x_i\right) a_1
```

양변을 $`n`$으로 나누면:

```math
\frac{\sum y_i}{n} = a_0 + \frac{\sum x_i}{n} \cdot a_1
```

$`\bar{y} = \frac{\sum y_i}{n}`$과 $`\bar{x} = \frac{\sum x_i}{n}`$이 표본 평균이므로:

```math
\boxed{a_0 = \bar{y} - \bar{x} \cdot a_1} \quad (\ast\ast)
```

> **[확률및통계]** 이 결과는 우아한 기하학적 해석을 가진다. 회귀선은 데이터의 무게 중심인 **점 $`(\bar{x}, \bar{y})`$를 항상 통과** 한다.

<br>

---

<br>

## 5. 적합도 정량화

### 5.1 SST, SSE, SSR 분해

세 가지 핵심 양이 변동성의 서로 다른 측면을 측정한다:

| 양 | 전체 명칭 | 수식 | 해석 |
|----------|-----------|---------|----------------|
| **SST** | 총 제곱합(Sum of Squared Total, 총 보정 제곱합) | $`\displaystyle \sum_{i=1}^{n}(y_i - \bar{y})^2`$ | 평균으로부터 데이터 점의 총 편차 |
| **SSE** | 오차 제곱합(Sum of Squared Errors, 잔차 제곱합) | $`\displaystyle \sum_{i=1}^{n}(y_i - \hat{y}_i)^2`$ | 설명되지 않은 변동(관측값과 예측값 사이의 오차) |
| **SSR** | 회귀 제곱합(Sum of Squared Regression, 설명된 제곱합) | $`\displaystyle \sum_{i=1}^{n}(\hat{y}_i - \bar{y})^2`$ | 설명된 변동($`\hat{y}`$가 데이터를 얼마나 잘 맞추는지) |

*표 1. 변동성을 측정하는 세 제곱합 SST, SSE, SSR*

기본적인 분해:

```math
\boxed{SST = SSE + SSR}
```

```math
\text{총 변동} = \text{설명되지 않은 변동} + \text{설명된 변동}
```

> **[확률및통계]** $`SSR > SSE`$이면 $`r^2 > 0.5`$이므로(5.3절), 회귀 모델이 설명하는 변동이 설명하지 못하는 변동보다 크다. 이 비교는 적합의 정도를 나타낼 뿐이며, 모델의 통계적 유의성은 별도의 검정으로 판단한다.

<br>

### 5.2 SST = SSE + SSR 증명

**단계 1:** $`\hat{y}_i`$를 더하고 빼서 SST를 표현한다.

```math
SST = \sum (y_i - \bar{y})^2 = \sum (y_i - \hat{y}_i + \hat{y}_i - \bar{y})^2
```

**단계 2:** 제곱을 전개한다.

```math
= \sum \left[(y_i - \hat{y}_i)^2 + 2(y_i - \hat{y}_i)(\hat{y}_i - \bar{y}) + (\hat{y}_i - \bar{y})^2\right]
```

**단계 3:** 교차항이 소멸한다. 다음을 보여야 한다:

```math
\sum (y_i - \hat{y}_i)(\hat{y}_i - \bar{y}) = 0
```

**교차항이 0인 증명:**

$`e_i = y_i - \hat{y}_i`$로 놓으면:

```math
\sum e_i (\hat{y}_i - \bar{y}) = \sum e_i (a_0 + a_1 x_i - \bar{y})
```

```math
= (a_0 - \bar{y}) \sum e_i + a_1 \sum e_i x_i
```

정규 방정식으로부터 두 가지 성질이 필요하다:

**성질 1:** $`\sum e_i = 0`$

```math
\sum e_i = \sum (y_i - \hat{y}_i) = \sum (y_i - a_0 - a_1 x_i)
```

```math
= \sum y_i - a_0 n - a_1 \sum x_i = n\left(\bar{y} - a_0 - a_1 \bar{x}\right) = 0 \quad \text{(} (\ast\ast)\text{에 의해)}
```

**성질 2:** $`\sum e_i x_i = 0`$

```math
\sum e_i x_i = \sum (y_i - \hat{y}_i) x_i = \sum y_i x_i - \sum (a_0 + a_1 x_i) x_i
```

```math
= \sum y_i x_i - a_0 \sum x_i - a_1 \sum x_i^2 = 0 \quad \text{(두 번째 정규 방정식 } (\ast)\text{에 의해)}
```

따라서:

```math
(a_0 - \bar{y}) \underbrace{\sum e_i}_{= 0} + a_1 \underbrace{\sum e_i x_i}_{= 0} = 0
```

**결론:**

```math
SST = \sum (y_i - \hat{y}_i)^2 + \sum (\hat{y}_i - \bar{y})^2 = SSE + SSR
```

<br>

### 5.3 결정 계수(Coefficient of Determination, $`r^2`$)

**결정 계수**는 회귀선이 데이터를 얼마나 잘 적합하는지를 정량화한다:

```math
r^2 = \frac{SSR}{SST} = 1 - \frac{SSE}{SST}
```

| $`r^2`$ 값 | 해석 |
|-------------|----------------|
| $`r^2 = 1`$ | 완벽한 적합; 모든 점이 직선 위에 있다 ($`SSE = 0`$) |
| $`r^2 = 0`$ | 모델이 변동을 전혀 설명하지 못한다 ($`SSR = 0`$) |
| $`r^2 \approx 1`$ | 강한 선형 관계 |
| $`r^2 \approx 0`$ | 약하거나 선형 관계 없음 |

*표 2. 결정 계수 값의 해석*

> **[확률및통계]** 단순 선형 회귀에서 값 $`r = \sqrt{r^2}`$(적절한 부호 포함)는 **피어슨 상관 계수(Pearson Correlation Coefficient)** 이다. 이는 $`x`$와 $`y`$ 사이의 선형 관계의 강도와 방향을 측정한다. $`r`$의 부호는 기울기 $`a_1`$의 부호와 일치한다.

<br>

### 5.4 추정의 표준 오차(Standard Error of the Estimate)

**추정의 표준 오차**는 회귀선 주위로 데이터 점들이 퍼진 정도를 정량화한다:

```math
s_{y/x} = \sqrt{\frac{SSE}{n - 2}}
```

> **[확률및통계]** $`n`$이 아닌 $`n - 2`$로 나누는 이유는 데이터로부터 **두 개의 매개변수**($`a_0`$와 $`a_1`$)를 추정하여 자유도 2를 소모했기 때문이다. 이는 하나의 매개변수(평균)를 추정한 표본 표준 편차에서 $`n - 1`$로 나누는 것과 유사하다.

표준 오차 $`s_{y/x}`$는 잔차의 전형적인 크기를 정량화한다. 더 작은 $`s_{y/x}`$는 회귀선이 데이터에 더 밀착하여 적합됨을 나타낸다.

<br>

---

<br>

## 요약

| 개념 | 핵심 요약 |
|------|-----------|
| 선형 모델 | 최적 직선은 $`\hat{y} = a_0 + a_1 x`$로 나타낸다. |
| SSE (목적 함수) | 최소제곱법이 최소화하는 양은 $`\displaystyle \sum_{i=1}^{n}(y_i - a_0 - a_1 x_i)^2`$이다. |
| 정규 방정식 (행렬) | 계수는 연립방정식 $`\begin{bmatrix} n & \sum x_i \\\\ \sum x_i & \sum x_i^2 \end{bmatrix} \begin{bmatrix} a_0 \\\\ a_1 \end{bmatrix} = \begin{bmatrix} \sum y_i \\\\ \sum y_i x_i \end{bmatrix}`$을 풀어 구한다. |
| 기울기 ($`a_1`$) | 변화율을 나타내는 기울기는 $`\displaystyle \frac{n\sum x_i y_i - \sum x_i \sum y_i}{n\sum x_i^2 - (\sum x_i)^2}`$이다. |
| 절편 ($`a_0`$) | $`x = 0`$일 때의 y값인 절편은 $`\bar{y} - a_1 \bar{x}`$이다. |
| SST 분해 | $`SST = SSE + SSR`$이며, 총 변동은 설명되지 않은 변동과 설명된 변동의 합이다. |
| 결정 계수 | 적합도를 0에서 1 사이의 값으로 나타내는 지표는 $`r^2 = SSR / SST = 1 - SSE / SST`$이다. |
| 표준 오차 | 잔차의 산포는 $`s_{y/x} = \sqrt{SSE / (n - 2)}`$로 측정한다. |
| 표본 표준 편차 | 평균 주위의 데이터 산포는 $`s_y = \sqrt{SST / (n - 1)}`$로 측정한다. |
| 변동 계수 | 정규화된 산포는 $`c.v. = (s_y / \bar{y}) \cdot 100\%`$로 나타낸다. |
| MAD | 강건한 산포 측정치는 $`\text{median}(\vert y_i - \tilde{y} \vert)`$이다. |
| 표준 정규 PDF | 종 모양 곡선의 밀도는 $`\varphi(z) = \frac{1}{\sqrt{2\pi}} e^{-z^2/2}`$이다. |
| 일반 정규 PDF | 평균 $`\mu`$, 분산 $`\sigma^2`$인 정규 밀도는 $`f(x) = \frac{1}{\sigma\sqrt{2\pi}} e^{-(x-\mu)^2/(2\sigma^2)}`$이다. |
| 표준화 | $`Z = (X - \mu) / \sigma`$는 정규 확률 변수를 $`\mathcal{N}(0,1)`$로 변환한다. |

<br>

---

<br>

## 점검 문제

1. **최소제곱 직선:** 데이터 (1, 2), (2, 4), (3, 5), (4, 4), (5, 5)에 최소제곱 직선 $`\hat{y} = a_0 + a_1 x`$를 적합하여 $`a_0`$와 $`a_1`$을 구하라.

   <details>
   <summary>정답</summary>

   $`n = 5`$, $`\sum x_i = 15`$, $`\sum y_i = 20`$, $`\sum x_i^2 = 55`$, $`\sum x_i y_i = 66`$이다. $`D = 5 \cdot 55 - 15^2 = 50`$이므로 $`a_1 = (5 \cdot 66 - 15 \cdot 20) / 50 = 30 / 50 = 0.6`$이다. $`\bar{x} = 3`$, $`\bar{y} = 4`$이므로 $`a_0 = \bar{y} - a_1 \bar{x} = 4 - 0.6 \cdot 3 = 2.2`$이고, 회귀선은 $`\hat{y} = 2.2 + 0.6x`$이다.

   </details>

2. **적합도:** 1번 문제의 회귀선에 대해 SST, SSE, SSR, $`r^2`$, $`s_{y/x}`$를 구하라.

   <details>
   <summary>정답</summary>

   예측값은 2.8, 3.4, 4.0, 4.6, 5.2이고 잔차는 -0.8, 0.6, 1.0, -0.6, -0.2이다. $`SSE = 0.64 + 0.36 + 1.00 + 0.36 + 0.04 = 2.4`$이다. $`\bar{y} = 4`$이므로 $`SST = 4 + 0 + 1 + 0 + 1 = 6`$이고, $`SSR = SST - SSE = 3.6`$이다. 따라서 $`r^2 = 1 - 2.4 / 6 = 0.6`$이며, $`s_{y/x} = \sqrt{2.4 / (5 - 2)} = \sqrt{0.8} \approx 0.894`$이다.

   </details>

3. **잔차의 성질:** 최소제곱 회귀선의 잔차 합 $`\sum e_i`$가 0인 이유를 정규 방정식으로 설명하고, 이 성질이 회귀선과 점 $`(\bar{x}, \bar{y})`$의 관계와 어떻게 연결되는지 쓰라.

   <details>
   <summary>정답</summary>

   첫 번째 정규 방정식 $`\sum y_i = n a_0 + a_1 \sum x_i`$는 $`\partial SSE / \partial a_0 = -2 \sum (y_i - a_0 - a_1 x_i) = 0`$에서 나오며, 이 조건이 곧 $`\sum e_i = 0`$이다. 이 식을 $`n`$으로 나누면 $`\bar{y} = a_0 + a_1 \bar{x}`$이므로 회귀선은 항상 점 $`(\bar{x}, \bar{y})`$를 지난다. 1번 문제에서도 $`2.2 + 0.6 \cdot 3 = 4 = \bar{y}`$이다.

   </details>

4. **강건한 산포 측정:** 데이터 [1, 3, 3, 6, 50]의 평균, 중위수, MAD를 구하고, 이상값 50이 평균과 MAD에 미치는 영향을 비교하라.

   <details>
   <summary>정답</summary>

   평균은 63 / 5 = 12.6이고 중위수는 3이다. 중위수로부터의 절대 편차는 2, 0, 0, 3, 47이며, 정렬하면 [0, 0, 2, 3, 47]이므로 MAD = 2이다. 50을 10으로 바꾸면 평균은 4.6으로 크게 바뀌지만 MAD는 그대로 2이다. 평균은 이상값의 크기에 끌려가지만, MAD는 편차의 중위수만 보므로 이상값의 크기에 영향을 받지 않는다.

   </details>

5. **표준화:** $`X \sim \mathcal{N}(10, 4)`$일 때 $`P(X \le 12)`$와 $`P(X \le 8)`$을 표준 정규 CDF $`\Phi`$로 나타내고 값을 구하라. $`\Phi(1) \approx 0.8413`$이다.

   <details>
   <summary>정답</summary>

   $`\sigma = 2`$이므로 $`P(X \le 12) = \Phi((12 - 10) / 2) = \Phi(1) \approx 0.8413`$이다. 대칭 성질 $`\Phi(-z) = 1 - \Phi(z)`$에 의해 $`P(X \le 8) = \Phi(-1) = 1 - \Phi(1) \approx 0.1587`$이다.

   </details>
