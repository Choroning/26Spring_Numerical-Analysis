# 제14장 강의 — 선형 회귀(Linear Regression)

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 14

> **선수 지식**: [통계학] 기초 통계. [미적분학] 최소자승 개념 (제1-13장).
>
> **학습 목표**:
> 1. 단순 선형 회귀를 유도하고 적용할 수 있다
> 2. R-제곱과 표준 오차를 사용하여 회귀 품질을 평가할 수 있다
> 3. 비선형 모델의 선형화를 적용할 수 있다

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
  - [5.3 결정 계수(Coefficient of Determination, $r^2$)](#53-결정-계수coefficient-of-determination-r2)
  - [5.4 추정의 표준 오차(Standard Error of the Estimate)](#54-추정의-표준-오차standard-error-of-the-estimate)
- [6. 요약 표](#6-요약-표)

<br>

---

<br>

## 1. 개요

선형 회귀(Linear Regression)는 **곡선 맞춤(Curve Fitting)**(수치 해석 Part IV)의 근본적인 기법이다. 목표는 **최소제곱(Least-Squares)** 기준을 사용하여 데이터 점들을 관통하는 **최적 직선(Best-Fit Straight Line)** 을 유도하는 것이다.

**학습 목표:**

- 기본 통계 복습 (평균, 표준 편차, 분포)
- 최적 직선의 기울기와 절편 계산
- 추정의 표준 오차 계산 및 잔차 오차 분석

<br>

---

<br>

## 2. 통계 복습

### 2.1 기술 통계(Descriptive Statistics)

**통계학(Statistics)** 은 의미 있는 정보를 추론하기 위해 데이터를 다루는 학문이다.

**표본 추출(Sampling)** 은 전체 집단에 대한 통찰을 얻기 위해 큰 모집단에서 작은 그룹을 선택하는 과정이다.

#### 표본 평균(Sample Average, Mean)

$n$개의 관측값 $y_1, y_2, \ldots, y_n$이 있을 때, **표본 평균** 은 다음과 같다:

$$\bar{y} = \frac{1}{n} \sum_{j=1}^{n} y_j$$

#### 중위수(Median)

**중위수** $\tilde{y}$는 데이터를 정렬했을 때 가운데 값이다.

- 홀수 개: 데이터 = [3, 5, 7] --> 중위수 = 5
- 짝수 개: 데이터 = [2, 4, 6, 8] --> 중위수 = (4 + 6) / 2 = 5

#### 최빈값(Mode)

**최빈값** 은 데이터에서 가장 자주 나타나는 값이다.

- 데이터 = [2, 4, 4, 4, 5] --> 최빈값 = 4

#### 이상값(Outlier)

**이상값** 은 나머지 데이터와 현저히 다른 데이터 점이다.

- 데이터 = [2, 4, 4, 4, 100] --> 이상값 = 100

<br>

### 2.2 표준 편차와 분산

#### 모표준편차(Population Standard Deviation)

**표준 편차** 는 데이터 점들이 평균으로부터 얼마나 퍼져 있는지를 측정한다:

$$\sigma = \sqrt{\frac{1}{n} \sum_{j=1}^{n} (y_j - \mu)^2}$$

여기서 $\mu$는 전체 모집단의 **참 평균(True Mean)** 이고, $n$은 모집단의 총 수이다.

> **[통계]** 모표준편차는 분모에 $n$을 사용하는데, 이는 모집단의 모든 구성원이 포함되어 자유도 손실이 없기 때문이다.

#### 표본 표준 편차(Sample Standard Deviation)

$$s_y = \sqrt{\frac{1}{n-1} \sum_{j=1}^{n} (y_j - \bar{y})^2}$$

여기서 $\bar{y}$는 **표본 평균** 이고 $n$은 표본 내 데이터 점의 총 수이다.

제곱근 내 분자에 해당하는 양을 **총 보정 제곱합(Total Corrected Sum of Squares, SST)** 이라 한다:

$$SST = \sum_{j=1}^{n} (y_j - \bar{y})^2$$

> **[통계]** 왜 $n$이 아니라 $n - 1$인가? 참 평균 $\mu$ 대신 표본 평균 $\bar{y}$를 사용하면 **변동성을 과소추정** 한다. $\bar{y}$는 이미 자유도 1을 사용한 확률 변수이므로, **비편향(Unbiased)** 추정치를 얻기 위해 $n - 1$로 나눈다(베셀 보정, Bessel's Correction).

<br>

### 2.3 변동 계수(Coefficient of Variation)

**변동 계수(c.v.)** 는 표준 편차 대 평균의 비율로, 백분율로 표현된다:

$$c.v. = \frac{s_y}{\bar{y}} \cdot 100\%$$

> **[통계]** 변동 계수는 산포도를 정규화하여, **스케일이 다른 데이터 세트도 비교할 수 있게** 한다. 낮은 c.v.는 데이터가 평균 주위에 밀집되어 있음을, 높은 c.v.는 평균에 비해 상당한 변동이 있음을 의미한다.

<br>

### 2.4 중위 절대 편차(Median Absolute Deviation, MAD)

**중위 절대 편차(MAD)** 는 값들이 중위수로부터 얼마나 벗어나는지를 측정한다:

$$MAD = \text{median}(|y_i - \tilde{y}|)$$

**예제:**

- 데이터: [2, 4, 6, 8, 100] --> 중위수 = 6
- 절대 편차: |2-6|=4, |4-6|=2, |6-6|=0, |8-6|=2, |100-6|=94
- 정렬된 편차: [0, 2, 2, 4, 94] --> MAD = 2

> **[통계]** MAD는 **이상값에 강건(Robust)** 하다. 위 예제에서 큰 이상값 100은 94라는 편차를 만들지만, MAD는 편차의 중위수를 취하므로 이 극단값을 완전히 무시한다. 따라서 MAD는 이상값이 있을 때 표준 편차보다 더 신뢰할 수 있는 산포 측정치이다.

<br>

---

<br>

## 3. 정규 분포(Normal Distribution)

### 3.1 중심극한정리(Central Limit Theorem)

각각 유한한 평균과 분산을 가지는 많은 독립 확률 변수를 더하면, 많은 수의 i.i.d. 확률 변수의 합은 개별 확률 변수의 분포가 무엇이든 **정규 분포에 수렴** 한다.

> **[통계]** 중심극한정리(CLT)는 정규 분포가 널리 나타나는 근본적인 이유이다. 많은 작은 독립적 효과의 합인 측정 오차는 일반적으로 정규 분포를 따른다.

<br>

### 3.2 표준 정규 분포(Standard Normal Distribution)

**표준 정규 분포** 는 중심극한정리 덕분에 매우 널리 사용된다.

**PDF (확률 밀도 함수, Probability Density Function):**

$$\varphi(z) = \frac{1}{\sqrt{2\pi}} e^{-\frac{z^2}{2}}, \quad \text{for } -\infty < z < \infty$$

**CDF (누적 분포 함수, Cumulative Distribution Function):**

$$\Phi(z) = \int_{-\infty}^{z} \varphi(t) \, dt, \quad \text{for } -\infty < z < \infty$$

이를 $Z \sim \mathcal{N}(0, 1)$로 표기하며, **평균 0, 분산 1** 을 의미한다.

<br>

### 3.3 대칭 성질

**1. PDF의 대칭성:**

$$\varphi(z) = \varphi(-z)$$

종 모양 곡선은 $z = 0$에 대해 대칭이다.

**2. 꼬리 면적의 대칭성:**

$$P(Z \le -2) = P(Z \ge 2)$$

이는 다음과 동치이다:

$$\Phi(-z) = 1 - \Phi(z)$$

**3. $Z$와 $-Z$의 대칭성:**

$Z \sim \mathcal{N}(0, 1)$이면, $-Z \sim \mathcal{N}(0, 1)$이기도 하다.

$$P(-Z \le z) = P(Z \ge -z) = 1 - \Phi(-z) = \Phi(z)$$

<br>

### 3.4 일반 정규 분포(General Normal Distribution)

$Z \sim \mathcal{N}(0, 1)$이면:

$$X = \mu + \sigma Z$$

는 **평균** $\mu$, **분산** $\sigma^2$인 정규 분포를 따르며, $X \sim \mathcal{N}(\mu, \sigma^2)$로 표기한다.

**평균과 분산의 증명:**

$$E(\mu + \sigma Z) = E(\mu) + E(\sigma Z) = \mu$$

$$Var(\mu + \sigma Z) = Var(\sigma Z) = \sigma^2 Var(Z) = \sigma^2$$

> **[통계]** 기댓값의 선형성($E[aX + b] = aE[X] + b$)과 분산의 스케일링 성질($Var(aX + b) = a^2 Var(X)$)이 여기서 사용되는 핵심 성질이다. 상수를 더하면 평균은 이동하지만 분산은 변하지 않는다.

<br>

### 3.5 표준화(Standardization)

$X \sim \mathcal{N}(\mu, \sigma^2)$일 때, $X$의 **표준화** 는 다음과 같다:

$$Z = \frac{X - \mu}{\sigma} \sim \mathcal{N}(0, 1)$$

이 변환은 임의의 정규 확률 변수를 표준 정규로 변환하여, 표준 정규 분포표를 사용할 수 있게 한다.

<br>

### 3.6 정규 분포의 CDF와 PDF

$X \sim \mathcal{N}(\mu, \sigma^2)$라 하자.

**$X$의 CDF:**

$$F(x) = P(X \le x) = P\left(\frac{X - \mu}{\sigma} \le \frac{x - \mu}{\sigma}\right) = \Phi\left(\frac{x - \mu}{\sigma}\right)$$

**$X$의 PDF:**

$$f(x) = F'(x) = \varphi\left(\frac{x - \mu}{\sigma}\right) \cdot \frac{1}{\sigma} = \frac{1}{\sigma\sqrt{2\pi}} e^{-\frac{(x - \mu)^2}{2\sigma^2}}$$

> **[통계]** PDF에서의 $\frac{1}{\sigma}$ 인자는 CDF를 미분할 때 연쇄 법칙(Chain Rule)에서 나온다. 이는 $\sigma$에 관계없이 PDF 아래 전체 면적이 1로 유지되도록 보장한다.

<br>

---

<br>

## 4. 직선 최소제곱 회귀(Straight Line Least Square Regression)

### 4.1 문제 설정

데이터 점 $(x_1, y_1), (x_2, y_2), \ldots, (x_n, y_n)$이 주어졌을 때, **선형 모델** 을 구하고자 한다:

$$\hat{y} = a_0 + a_1 x$$

여기서:
- $a_0$ = **y절편(y-intercept)**
- $a_1$ = **기울기(slope)**

각 관측값은 다음과 같이 표현할 수 있다:

$$y_i = \hat{y}_i + e_i = a_0 + a_1 x_i + e_i$$

여기서 $e_i = y_i - \hat{y}_i$는 **잔차(Residual, Error)** — 관측값과 예측값 사이의 차이이다.

<br>

### 4.2 대안적 오차 기준과 최소제곱법의 이유

오차의 어떤 측도를 최소화하여 데이터를 관통하는 "최적" 직선을 구하고자 한다. 여러 기준을 고려할 수 있다:

**기준 1 — 오차의 합 최소화:**

$$\min_{a_0, a_1} \sum e_i = \min_{a_0, a_1} \sum (y_i - a_0 - a_1 x_i)$$

문제: 양의 오차와 음의 오차가 **서로 상쇄** 되어 오도된 결과를 준다.

**기준 2 — 절대 오차의 합 최소화:**

$$\min_{a_0, a_1} \sum |e_i| = \min_{a_0, a_1} \sum |y_i - a_0 - a_1 x_i|$$

문제: 절대값 함수는 0에서 **미분 불가능** 하고, 해가 **유일하지 않다**.

**기준 3 — 미니맥스(Minimax, 최대 오차 최소화):**

$$\min_{a_0, a_1} \max_i |y_i - a_0 - a_1 x_i|$$

문제: **미분 불가능** 하고, 단일 최대 오차에 집중하므로 **이상값에 민감** 하다.

**기준 4 — 오차 제곱합 최소화(최소제곱법, Least Squares):**

$$\min_{a_0, a_1} \sum e_i^2 = \min_{a_0, a_1} \sum (y_i - a_0 - a_1 x_i)^2$$

이것이 **최소제곱(Least Squares)** 접근법이며, **오차 제곱합(Sum of Squared Errors, SSE)** 을 산출한다. 미분 가능하고 주어진 데이터에 대해 **유일한** 직선을 산출한다.

> **[통계]** 최소제곱 기준이 선호되는 이유: (1) 제곱은 부호 문제를 제거하고, (2) 목적 함수가 매끄럽고 미분 가능하며, (3) 유일한 닫힌 형태의 해를 가지고, (4) 오차가 정규 분포를 따른다는 가정 하에 최대 우도 추정(Maximum Likelihood Estimate)에 해당한다.

<br>

### 4.3 정규 방정식 유도

최소화할 목적 함수는 다음과 같다:

$$SSE = \sum_{i=1}^{n} (y_i - \hat{y}_i)^2 = \sum_{i=1}^{n} (y_i - a_0 - a_1 x_i)^2$$

최솟값을 찾기 위해 SSE를 각 미지 계수에 대해 미분하고 0으로 놓는다:

$$\frac{\partial SSE}{\partial a_0} = -2 \sum_{i=1}^{n} (y_i - a_0 - a_1 x_i) = 0$$

$$\frac{\partial SSE}{\partial a_1} = -2 \sum_{i=1}^{n} (y_i - a_0 - a_1 x_i) x_i = 0$$

$a_0$와 $a_1$은 상수($i$에 의존하지 않음)이므로 합에서 인수 분해할 수 있다:

- $\sum a_0 = a_0 \cdot n$
- $\sum a_0 x_i = a_0 \sum x_i$
- $\sum a_1 x_i = a_1 \sum x_i$
- $\sum a_1 x_i^2 = a_1 \sum x_i^2$

편도함수를 0으로 놓으면 **정규 방정식(Normal Equations)** 을 얻는다:

$$\sum y_i = n \cdot a_0 + a_1 \sum x_i$$

$$\sum y_i x_i = a_0 \sum x_i + a_1 \sum x_i^2$$

**행렬 형태:**

$$\begin{bmatrix} \sum y_i \\ \sum y_i x_i \end{bmatrix} = \begin{bmatrix} n & \sum x_i \\ \sum x_i & \sum x_i^2 \end{bmatrix} \begin{bmatrix} a_0 \\ a_1 \end{bmatrix} \quad (\ast)$$

> **[선형대수]** 이것은 $2 \times 2$ 선형 시스템이다. 계수 행렬은 대칭이고 양의 정부호(Positive Definite)이므로($x_i$ 값이 모두 같지 않은 한), 유일한 해가 보장된다.

<br>

### 4.4 계수 풀이

#### 크래머 법칙(Cramer's Rule) 사용

계수 행렬의 행렬식:

$$D = n \sum x_i^2 - \left(\sum x_i\right)^2$$

절편:

$$a_0 = \frac{1}{D} \left( \sum y_i \cdot \sum x_i^2 - \sum x_i \cdot \sum y_i x_i \right)$$

기울기:

$$a_1 = \frac{1}{D} \left( n \sum y_i x_i - \sum x_i \cdot \sum y_i \right)$$

#### 표본 평균을 사용한 대안적 형태

$(\ast)$의 첫 번째 행으로부터:

$$\sum y_i = n \cdot a_0 + \left(\sum x_i\right) a_1$$

양변을 $n$으로 나누면:

$$\frac{\sum y_i}{n} = a_0 + \frac{\sum x_i}{n} \cdot a_1$$

$\bar{y} = \frac{\sum y_i}{n}$과 $\bar{x} = \frac{\sum x_i}{n}$이 표본 평균이므로:

$$\boxed{a_0 = \bar{y} - \bar{x} \cdot a_1} \quad (\ast\ast)$$

> **[통계]** 이 결과는 우아한 기하학적 해석을 가진다: 회귀선은 데이터의 무게 중심인 **점 $(\bar{x}, \bar{y})$를 항상 통과** 한다.

<br>

---

<br>

## 5. 적합도 정량화

### 5.1 SST, SSE, SSR 분해

세 가지 핵심 양이 변동성의 서로 다른 측면을 측정한다:

| 양 | 전체 명칭 | 수식 | 해석 |
|----------|-----------|---------|----------------|
| **SST** | 총 제곱합(Sum of Squared Total, 총 보정 제곱합) | $\displaystyle \sum_{i=1}^{n}(y_i - \bar{y})^2$ | 평균으로부터 데이터 점의 총 편차 |
| **SSE** | 오차 제곱합(Sum of Squared Errors, 잔차 제곱합) | $\displaystyle \sum_{i=1}^{n}(y_i - \hat{y}_i)^2$ | 설명되지 않은 변동(관측값과 예측값 사이의 오차) |
| **SSR** | 회귀 제곱합(Sum of Squared Regression, 설명된 제곱합) | $\displaystyle \sum_{i=1}^{n}(\hat{y}_i - \bar{y})^2$ | 설명된 변동($\hat{y}$가 데이터를 얼마나 잘 맞추는지) |

기본적인 분해:

$$\boxed{SST = SSE + SSR}$$

$$\text{총 변동} = \text{설명되지 않은 변동} + \text{설명된 변동}$$

> **[통계]** $SSR > SSE$이면, 회귀 모델이 설명하지 못하는 변동보다 더 많은 변동을 설명하는 것이므로, 모델이 유의미하다.

<br>

### 5.2 SST = SSE + SSR 증명

**단계 1:** $\hat{y}_i$를 더하고 빼서 SST를 표현:

$$SST = \sum (y_i - \bar{y})^2 = \sum (y_i - \hat{y}_i + \hat{y}_i - \bar{y})^2$$

**단계 2:** 제곱을 전개:

$$= \sum \left[(y_i - \hat{y}_i)^2 + 2(y_i - \hat{y}_i)(\hat{y}_i - \bar{y}) + (\hat{y}_i - \bar{y})^2\right]$$

**단계 3:** 교차항이 소멸한다. 다음을 보여야 한다:

$$\sum (y_i - \hat{y}_i)(\hat{y}_i - \bar{y}) = 0$$

**교차항이 0인 증명:**

$e_i = y_i - \hat{y}_i$로 놓으면:

$$\sum e_i (\hat{y}_i - \bar{y}) = \sum e_i (a_0 + a_1 x_i - \bar{y})$$

$$= (a_0 - \bar{y}) \sum e_i + a_1 \sum e_i x_i$$

정규 방정식으로부터 두 가지 성질이 필요하다:

**성질 1:** $\sum e_i = 0$

$$\sum e_i = \sum (y_i - \hat{y}_i) = \sum (y_i - a_0 - a_1 x_i)$$

$$= \sum y_i - a_0 n - a_1 \sum x_i = n\left(\bar{y} - a_0 - a_1 \bar{x}\right) = 0 \quad \text{(} (\ast\ast)\text{에 의해)}$$

**성질 2:** $\sum e_i x_i = 0$

$$\sum e_i x_i = \sum (y_i - \hat{y}_i) x_i = \sum y_i x_i - \sum (a_0 + a_1 x_i) x_i$$

$$= \sum y_i x_i - a_0 \sum x_i - a_1 \sum x_i^2 = 0 \quad \text{(두 번째 정규 방정식 } (\ast)\text{에 의해)}$$

따라서:

$$(a_0 - \bar{y}) \underbrace{\sum e_i}_{= 0} + a_1 \underbrace{\sum e_i x_i}_{= 0} = 0$$

**결론:**

$$SST = \sum (y_i - \hat{y}_i)^2 + \sum (\hat{y}_i - \bar{y})^2 = SSE + SSR$$

<br>

### 5.3 결정 계수(Coefficient of Determination, $r^2$)

**결정 계수** 는 회귀선이 데이터를 얼마나 잘 적합하는지를 정량화한다:

$$r^2 = \frac{SSR}{SST} = 1 - \frac{SSE}{SST}$$

| $r^2$ 값 | 해석 |
|-------------|----------------|
| $r^2 = 1$ | 완벽한 적합; 모든 점이 직선 위에 있다 ($SSE = 0$) |
| $r^2 = 0$ | 모델이 변동을 전혀 설명하지 못한다 ($SSR = 0$) |
| $r^2 \approx 1$ | 강한 선형 관계 |
| $r^2 \approx 0$ | 약하거나 선형 관계 없음 |

> **[통계]** 값 $r = \sqrt{r^2}$(적절한 부호 포함)는 **피어슨 상관 계수(Pearson Correlation Coefficient)** 이다. 이는 $x$와 $y$ 사이의 선형 관계의 강도와 방향을 측정한다. $r$의 부호는 기울기 $a_1$의 부호와 일치한다.

<br>

### 5.4 추정의 표준 오차(Standard Error of the Estimate)

**추정의 표준 오차** 는 회귀선 주위로 데이터 점들이 퍼진 정도를 정량화한다:

$$s_{y/x} = \sqrt{\frac{SSE}{n - 2}}$$

> **[통계]** $n$이 아닌 $n - 2$로 나누는 이유는 데이터로부터 **두 개의 매개변수**($a_0$와 $a_1$)를 추정하여 자유도 2를 소모했기 때문이다. 이는 하나의 매개변수(평균)를 추정한 표본 표준 편차에서 $n - 1$로 나누는 것과 유사하다.

표준 오차 $s_{y/x}$는 잔차의 전형적인 크기를 정량화한다. 더 작은 $s_{y/x}$는 회귀선이 데이터에 더 밀착하여 적합됨을 나타낸다.

<br>

---

<br>

## 6. 요약 표

| 개념 | 수식 | 목적 |
|---------|---------|---------|
| 선형 모델 | $\hat{y} = a_0 + a_1 x$ | 최적 직선 |
| SSE (목적 함수) | $\displaystyle \sum_{i=1}^{n}(y_i - a_0 - a_1 x_i)^2$ | 최소화할 양 |
| 정규 방정식 (행렬) | $\begin{bmatrix} n & \sum x_i \\\\ \sum x_i & \sum x_i^2 \end{bmatrix} \begin{bmatrix} a_0 \\\\ a_1 \end{bmatrix} = \begin{bmatrix} \sum y_i \\\\ \sum y_i x_i \end{bmatrix}$ | 계수를 위한 연립방정식 |
| 기울기 ($a_1$) | $\displaystyle \frac{n\sum x_i y_i - \sum x_i \sum y_i}{n\sum x_i^2 - (\sum x_i)^2}$ | 변화율 |
| 절편 ($a_0$) | $\bar{y} - a_1 \bar{x}$ | $x = 0$일 때의 y값 |
| SST 분해 | $SST = SSE + SSR$ | 총 변동 = 설명되지 않은 변동 + 설명된 변동 |
| 결정 계수 | $r^2 = SSR / SST = 1 - SSE / SST$ | 적합도 (0에서 1) |
| 표준 오차 | $s_{y/x} = \sqrt{SSE / (n - 2)}$ | 잔차의 산포 |
| 표본 표준 편차 | $s_y = \sqrt{SST / (n - 1)}$ | 평균 주위의 데이터 산포 |
| 변동 계수 | $c.v. = (s_y / \bar{y}) \cdot 100\%$ | 정규화된 산포 |
| MAD | $\text{median}(\|y_i - \tilde{y}\|)$ | 강건한 산포 측정치 |
| 표준 정규 PDF | $\varphi(z) = \frac{1}{\sqrt{2\pi}} e^{-z^2/2}$ | 종 모양 곡선 밀도 |
| 일반 정규 PDF | $f(x) = \frac{1}{\sigma\sqrt{2\pi}} e^{-(x-\mu)^2/(2\sigma^2)}$ | 평균 $\mu$, 분산 $\sigma^2$인 정규 밀도 |
| 표준화 | $Z = (X - \mu) / \sigma$ | $\mathcal{N}(0,1)$로 변환 |
