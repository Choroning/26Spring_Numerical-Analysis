# 7장 — 최적화

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **최종 수정일:** 2026-10-20
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 7

> **학습 목표:**
> 1. 1차원 최적화에 황금 분할 탐색을 적용할 수 있다
> 2. 포물선 보간으로 1차원 최적점을 추정할 수 있다
> 3. 국소 최솟값, 전역 최솟값, 볼록성의 정의를 설명할 수 있다
> 4. 최급 하강법 알고리즘을 구현할 수 있다

---

<br>

## 목차

- [1. 최적화 소개](#1-최적화-소개)
  - [1.1 최적화란?](#11-최적화란)
  - [1.2 동기 부여 예제: 알루미늄 캔 설계](#12-동기-부여-예제-알루미늄-캔-설계)
  - [1.3 수학적 표현](#13-수학적-표현)
  - [1.4 최소화 vs. 최대화](#14-최소화-vs-최대화)
  - [1.5 실행 가능 영역과 제약 조건](#15-실행-가능-영역과-제약-조건)
- [2. 핵심 위상 개념](#2-핵심-위상-개념)
  - [2.1 내부점(Interior Points)](#21-내부점interior-points)
  - [2.2 경계점(Boundary Points)](#22-경계점boundary-points)
  - [2.3 실행 가능 방향(Feasible Directions)](#23-실행-가능-방향feasible-directions)
- [3. 국소 최적성과 전역 최적성](#3-국소-최적성과-전역-최적성)
  - [3.1 전역 최솟값(Global Minimum)](#31-전역-최솟값global-minimum)
  - [3.2 하한과 상한(Infimum and Supremum)](#32-하한과-상한infimum-and-supremum)
  - [3.3 국소 최솟값(Local Minimum)](#33-국소-최솟값local-minimum)
  - [3.4 순 국소 최솟값(Strict Local Minimum)](#34-순-국소-최솟값strict-local-minimum)
- [4. 최적해의 존재](#4-최적해의-존재)
  - [4.1 바이어슈트라스 정리(Weierstrass Theorem)](#41-바이어슈트라스-정리weierstrass-theorem)
  - [4.2 강제 함수(Coercive Functions)](#42-강제-함수coercive-functions)
  - [4.3 강제 함수에 대한 증명 개요](#43-강제-함수에-대한-증명-개요)
- [5. 1차원 최적화](#5-1차원-최적화)
- [6. 황금 분할 탐색(Golden-Section Search)](#6-황금-분할-탐색golden-section-search)
  - [6.1 황금비(Golden Ratio)](#61-황금비golden-ratio)
  - [6.2 알고리즘: 구간 축소](#62-알고리즘-구간-축소)
  - [6.3 함수 평가의 재사용](#63-함수-평가의-재사용)
  - [6.4 오차 분석](#64-오차-분석)
- [7. 포물선 보간(Parabolic Interpolation)](#7-포물선-보간parabolic-interpolation)
  - [7.1 포물선 구성](#71-포물선-구성)
  - [7.2 최적점 공식](#72-최적점-공식)
  - [7.3 다음 반복을 위한 점 갱신](#73-다음-반복을-위한-점-갱신)
- [8. 등고선(Level Sets)과 기울기(Gradient)](#8-등고선level-sets과-기울기gradient)
  - [8.1 등고선 정의](#81-등고선-정의)
  - [8.2 기울기는 등고선에 직교한다](#82-기울기는-등고선에-직교한다)
- [9. 볼록 문제(Convex Problems)](#9-볼록-문제convex-problems)
  - [9.1 볼록 결합(Convex Combination)](#91-볼록-결합convex-combination)
  - [9.2 볼록 집합(Convex Set)](#92-볼록-집합convex-set)
  - [9.3 초평면(Hyperplane)과 다면체 집합(Polyhedral Set)](#93-초평면hyperplane과-다면체-집합polyhedral-set)
  - [9.4 꼭짓점(Extreme Point)](#94-꼭짓점extreme-point)
  - [9.5 볼록 함수(Convex Function)](#95-볼록-함수convex-function)
- [10. 최급 하강법(Steepest Gradient Descent Method)](#10-최급-하강법steepest-gradient-descent-method)
  - [10.1 방향 도함수와 코시-슈바르츠(Cauchy-Schwarz)](#101-방향-도함수와-코시-슈바르츠cauchy-schwarz)
  - [10.2 최급 하강 갱신 규칙](#102-최급-하강-갱신-규칙)
  - [10.3 수렴과 직교성 성질](#103-수렴과-직교성-성질)
  - [10.4 예제](#104-예제)
- [요약](#요약)
- [점검 문제](#점검-문제)

---

<br>

## 1. 최적화 소개

### 1.1 최적화란?

**최적화(Optimization)** 는 모든 **실행 가능한 해** 중에서 **최선**을 찾는 방법론이다. **최적해(optimal solution)** 는 최선의 목적 함수 값을 산출하는 실행 가능한 해이다.

함수 $`f(x)`$는 여러 유형의 극값점을 가질 수 있다:

- **전역 최댓값(Global maximum)**: 전체 정의역에서 $`f`$가 가장 큰 값을 갖는 점이다.
- **전역 최솟값(Global minimum)**: 전체 정의역에서 $`f`$가 가장 작은 값을 갖는 점이다.
- **국소 최댓값(Local maximum)**: 인근의 어떤 점에서보다 $`f`$의 값이 작지 않은 점이다.
- **국소 최솟값(Local minimum)**: 인근의 어떤 점에서보다 $`f`$의 값이 크지 않은 점이다.

### 1.2 동기 부여 예제: 알루미늄 캔 설계

한 회사가 지름 $`d`$와 높이 $`h`$를 가진 원기둥 형태의 알루미늄 캔을 설계해야 한다.

**공식:**

```math
\text{Volume} = \pi \left(\frac{d}{2}\right)^2 h = \frac{\pi d^2 h}{4}
```

```math
\text{Surface Area} = \frac{\pi d^2}{2} + \pi d h
```

부피와 겉넓이 모두 지름 $`d`$와 높이 $`h`$의 함수이다.

**목표:** 다음 조건을 만족하면서 **가능한 한 적은 알루미늄**을 사용하는 캔, 즉 겉넓이가 최소인 캔을 설계한다:

1. 높이는 지름보다 50% 이상 커야 한다: $`h \geq \frac{3}{2}d \iff \frac{3}{2}d - h \leq 0`$
2. 높이는 지름의 두 배를 넘을 수 없다: $`h \leq 2d \iff -2d + h \leq 0`$
3. 부피는 330 mL 이상이어야 한다: $`\frac{\pi d^2 h}{4} \geq 330`$
4. $`h > 0, \; d > 0`$이어야 하며, 이 조건은 이미 보장된다. 조건 1과 2에서 $`d \geq 0`$, $`h \geq 0`$이 나오고, 조건 3이 $`d = 0`$이나 $`h = 0`$을 배제하기 때문이다.

**최적화 문제 요약:**

```math
\text{minimize} \quad \frac{\pi d^2}{2} + \pi d h
```

```math
\text{subject to} \quad \frac{3}{2}d - h \leq 0
```

```math
-2d + h \leq 0
```

```math
\frac{\pi d^2 h}{4} \geq 330
```

### 1.3 수학적 표현

$`\mathbb{X} \subseteq \mathbb{R}^n`$을 $`n`$차원 실수 벡터의 집합이라 하자. 즉 $`\mathbb{X} \ni x = [x_1, x_2, \ldots, x_n]^T`$이다.

$`f : \mathbb{X} \to \mathbb{R}`$이라 하자.

최적화 문제의 일반적인 형태는 다음과 같다:

```math
\min_{x \in \mathbb{X}} f(x)
```

이는 다음과 같이 쓸 수도 있다:

```math
\text{minimize} \quad f(x) \quad \text{subject to} \quad x \in \mathbb{X}
```

여기서 각 기호의 뜻은 다음과 같다:

- $`\mathbb{X}`$는 **실행 가능(허용) 영역(feasible/admissible region)** 이다.
- $`f`$는 **목적 함수(objective function)** 이다.
- $`x`$는 **실행 가능한 해(feasible solution)** 또는 **결정 변수(decision variable)** 이다.

### 1.4 최소화 vs. 최대화

최소화 문제 $`\min_{x \in \mathbb{X}} f(x)`$는 동등한 최대화 문제로 변환할 수 있다:

```math
\max_{x \in \mathbb{X}} (-f(x))
```

그 역도 성립한다. 따라서 일반성을 잃지 않고 항상 최소화 문제만 고려할 수 있다.

### 1.5 실행 가능 영역과 제약 조건

실행 가능 영역 $`\mathbb{X}`$는 **제약 조건(constraints)** 의 집합으로 기술된다:

- **등식 제약 조건(Equality constraints):** $`h_i(x) = 0`$ 꼴이다. 여기서 $`h_i : \mathbb{R}^n \to \mathbb{R}`$이고 $`i \in \mathcal{E}_i`$(인덱스 집합)이다.
- **부등식 제약 조건(Inequality constraints):** $`g_j(x) \leq 0`$ 꼴이다. 여기서 $`g_j : \mathbb{R}^n \to \mathbb{R}`$이고 $`j \in \mathcal{E}_j`$(인덱스 집합)이다.

**예시:** 두 개의 부등식 제약 조건으로 정의된 다음 실행 가능 집합 $`\mathbb{X}`$를 생각한다:

```math
\mathbb{X} = \{x \in \mathbb{R}^2 : x_1^2 + x_2^2 \leq 1, \; (x_1 - 1)^2 + x_2^2 \leq 1\}
```

여기서 $`g_1(x_1, x_2) = x_1^2 + x_2^2 - 1 \leq 0`$이고 $`g_2(x_1, x_2) = (x_1 - 1)^2 + x_2^2 - 1 \leq 0`$이다. 실행 가능 영역은 두 원판의 교집합(겹치는 부분)이다.

---

<br>

## 2. 핵심 위상 개념

### 2.1 내부점(Interior Points)

**정의.** 점 $`\hat{x} \in \mathbb{X}`$가 $`\mathbb{X}`$의 **내부점(interior point)** 이란 $`\varepsilon > 0`$이 존재하여 $`\hat{x}`$를 중심으로 한 $`\varepsilon`$-공이 $`\mathbb{X}`$에 포함되는 것을 말한다:

```math
B(\hat{x}, \varepsilon) \subset \mathbb{X}
```

**정의.** $`\text{int}(\mathbb{X})`$는 $`\mathbb{X}`$의 모든 내부점의 집합이다.

### 2.2 경계점(Boundary Points)

**정의.** 점 $`\bar{x}`$가 다음을 만족하면 $`\mathbb{X}`$의 **경계점(boundary point)** 이라 한다:

```math
\bar{x} \in \mathbb{X} \quad \text{and} \quad \bar{x} \notin \text{int}(\mathbb{X})
```

다시 말해, $`\bar{x}`$는 집합에 속하지만 그 주위의 어떤 $`\varepsilon`$-공도 집합 안에 완전히 포함되지 않는다. 기하학적으로 경계점은 집합의 가장자리에 놓이고, 내부점은 모든 방향으로 여유를 두고 집합 안쪽에 놓인다.

### 2.3 실행 가능 방향(Feasible Directions)

**정의.** 벡터 $`d \in \mathbb{R}^n`$이 $`\bar{x} \in \mathbb{X}`$에서 집합 $`\mathbb{X} \subseteq \mathbb{R}^n`$에 대한 **실행 가능 방향(feasible direction)** 이란 $`\delta > 0`$이 존재하여 다음이 성립하는 것이다:

```math
\bar{x} + \alpha d \in \mathbb{X} \quad \text{for any } 0 \leq \alpha \leq \delta
```

> **참고:** **기하학적 의미.** 실행 가능 방향은 현재 점에서 그 방향으로 양의 작은 스텝만큼 이동해도 실행 가능 집합 안에 머무를 수 있는 방향이다. 내부점에서는 모든 방향이 실행 가능하다. 경계점에서는 바깥을 향하는 방향은 실행 가능하지 않다. 안쪽을 향하는 방향은 실행 가능하고, 경계를 따르는 방향은 경계가 평평한 부분(예: 다면체 집합의 면)에서만 실행 가능하다.

---

<br>

## 3. 국소 최적성과 전역 최적성

다음 최소화 문제를 생각한다:

```math
\text{minimize} \quad f(x) \quad \text{subject to} \quad x \in \mathbb{X}
```

여기서 $`f : \mathbb{R}^n \to \mathbb{R}`$이다.

### 3.1 전역 최솟값(Global Minimum)

**정의 (전역 최솟값).** 점 $`x^* \in \mathbb{X}`$가 다음을 만족하면 문제 $`\min_{x \in \mathbb{X}} f(x)`$의 **전역 최솟값(global minimum)** (전역 최소화자)이라 하고 $`x^* = \arg\min_{x \in \mathbb{X}} f(x)`$로 쓴다:

```math
f(x^*) \leq f(x) \quad \forall \; x \in \mathbb{X}
```

전역 최소화자 $`x^*`$가 다음을 만족하면 **순 전역 최소화자(strict global minimizer)** 라 한다:

```math
f(x^*) < f(x) \quad \forall \; x \in \mathbb{X} \setminus \{x^*\}
```

### 3.2 하한과 상한(Infimum and Supremum)

모든 함수가 전역 최소화자를 갖는 것은 아니다. $`x \in \mathbb{R}`$에서 $`f(x) = e^x`$를 생각하면 다음이 성립한다:

- $`f(x) \geq 0`$이므로 $`f`$는 아래로 유계이다.
- $`\lim_{x \to -\infty} f(x) = 0`$이고 0이 최대 하계이지만, $`f(x) = 0`$인 $`x`$는 존재하지 않는다.

**정의 (하한, 상한).** 함수 $`f : \mathbb{X} \to \mathbb{R}`$에 대해 다음과 같이 정의한다:

- $`f`$의 **최대 하계(greatest lower bound)** 를 **하한(infimum)** 이라 한다: $`\inf_{x \in \mathbb{X}} f(x)`$
- $`f`$의 **최소 상계(least upper bound)** 를 **상한(supremum)** 이라 한다: $`\sup_{x \in \mathbb{X}} f(x)`$

**예시:**

- $`\sup_{x \in \mathbb{R}} e^x = +\infty`$
- $`\inf_{x \in (0, +\infty)} \ln(x) = -\infty`$

> **참고:** **최솟값과 하한.** 함수의 최솟값(minimum)은 실제로 달성되는 하한이다. $`f(x^*) = \inf f`$인 $`x^*`$가 존재하면 $`\min f = \inf f`$이다. 하한이 달성되지 않으면, 하한은 잘 정의되더라도 최솟값은 존재하지 않는다.

### 3.3 국소 최솟값(Local Minimum)

**정의 (국소 최솟값).** $`x^* \in \mathbb{X}`$가 문제 $`\min_{x \in \mathbb{X}} f(x)`$의 **국소 최솟값(local minimum)** (국소 최소화자)이란 $`\varepsilon > 0`$이 존재하여 다음이 성립하는 것이다:

```math
f(x) \geq f(x^*) \quad \forall \; x \in \mathbb{X} \text{ with } \|x - x^*\| \leq \varepsilon
```

### 3.4 순 국소 최솟값(Strict Local Minimum)

국소 최소화자 $`x^*`$가 **순 국소 최소화자(strict local minimizer)** 이란 $`\varepsilon > 0`$이 존재하여 다음이 성립하는 것이다:

```math
f(x) > f(x^*) \quad \forall \; x \in \mathbb{X} \setminus \{x^*\} \text{ with } \|x - x^*\| \leq \varepsilon
```

> **핵심:** 전역 최소화자는 **전체** 실행 가능 집합에서 가장 작은 함수 값을 가지는 반면, 국소 최소화자는 어떤 이웃 안에서만 최적이면 된다. 모든 전역 최소화자는 국소 최소화자이지만, 그 역이 항상 성립하는 것은 아니다.

---

<br>

## 4. 최적해의 존재

전역 최소화자가 존재하는지 확인하는 것은 일반적으로 매우 어려운 문제이다. 그러나 일부 상황에서는 전역 최소화자의 존재를 보장할 수 있다.

### 4.1 바이어슈트라스 정리(Weierstrass Theorem)

**정리 (바이어슈트라스).** $`\mathbb{X} \subset \mathbb{R}^n`$이 공집합이 아닌 **콤팩트 집합**(닫히고 유계인 집합)이고 $`f : \mathbb{X} \to \mathbb{R}`$이 **연속 함수**일 때, 문제 $`\min_{x \in \mathbb{X}} f(x)`$와 $`\max_{x \in \mathbb{X}} f(x)`$는 **전역 최적해**를 가진다.

> **참고:** **콤팩트 집합.** $`\mathbb{R}^n`$에서의 콤팩트 집합은 닫힘(집합의 점들로 이루어진 수렴 수열의 극한이 모두 집합에 속함)과 유계(유한 반지름의 공 안에 들어감) 조건을 모두 만족하는 집합이다. 바이어슈트라스 정리는 연속 함수가 공집합이 아닌 콤팩트 집합에서 극값을 달성함을 보장한다.

### 4.2 강제 함수(Coercive Functions)

**정의 (강제 함수).** 함수 $`f : \mathbb{R}^n \to \mathbb{R}`$이 다음을 만족하면 **강제(coercive)** 라 한다:

```math
\lim_{\|x\| \to \infty} f(x) = +\infty
```

$`x`$가 원점에서 임의로 멀어지면 함수 값 $`f(x)`$가 무한대로 증가한다. 예를 들어 $`f(x) = x^2`$은 강제 함수이고, $`f(x) = \sin(x)`$는 강제 함수가 아니다.

**정리.** 모든 연속 강제 함수 $`f : \mathbb{R}^n \to \mathbb{R}`$은 $`\mathbb{R}^n`$에서 **전역 최소화자**를 가진다.

### 4.3 강제 함수에 대한 증명 개요

1. 임의의 점 $`\hat{x} \in \mathbb{R}^n`$을 하나 고정한다. $`f`$가 강제이므로, $`\|x\| > C`$이면 $`f(x) > f(\hat{x})`$가 되는 $`C > 0`$이 존재한다.
2. $`\mathbb{R}^n = \mathbb{X}_1 \cup \mathbb{X}_2`$를 서로소인 두 집합의 합집합으로 나타낸다. 여기서 $`\mathbb{X}_1 = \{x : \|x\| \leq C\}`$이고 $`\mathbb{X}_2 = \{x : \|x\| > C\}`$이다.
3. 그러면 $`\hat{x} \in \mathbb{X}_1`$이고, 모든 $`x \in \mathbb{X}_2`$에 대해 $`f(\hat{x}) < f(x)`$이다.
4. $`\mathbb{X}_1`$은 콤팩트 집합이므로, 바이어슈트라스 정리에 의해 $`f(x^*) = \min_{x \in \mathbb{X}_1} f(x)`$인 $`x^* \in \mathbb{X}_1`$이 존재한다.
5. 따라서 $`f(\hat{x}) \geq f(x^*)`$이고, 이로부터 모든 $`x \in \mathbb{R}^n`$에 대해 $`f(x) \geq f(x^*)`$가 된다.
6. 이는 $`x^*`$가 전역 최소화자임을 의미한다.

---

<br>

## 5. 1차원 최적화

1차원 최적화의 목표는 $`x \in \mathbb{R}`$인 함수 $`f(x)`$의 최솟값(또는 최댓값)을 찾는 것이다.

**목표:**

- 1차원 최적화와 다차원 최적화를 이해한다.
- 전역 최적과 국소 최적을 구별한다.
- 최대화 문제를 최소화 문제로 변환한다.
- **황금 분할 탐색(Golden-section search)** 과 **포물선 보간(Parabolic interpolation)** 방법을 학습한다.

---

<br>

## 6. 황금 분할 탐색(Golden-Section Search)

황금 분할 탐색은 닫힌 구간에서 단봉 함수(unimodal function, 구간 안에 극소점이 하나뿐인 함수)의 최솟값을 찾기 위한 효율적인 구간법이다. 최솟값을 포함하는 구간을 체계적으로 좁힌다.

### 6.1 황금비(Golden Ratio)

**정의 (황금비).** 선분을 두 부분 $`l_1`$과 $`l_2`$ ($`l_1 > l_2 > 0`$)로 나누되 다음을 만족하도록 한다:

```math
\frac{l_1 + l_2}{l_1} = \frac{l_1}{l_2} =: \phi \quad \text{(황금비)}
```

**전체 선분** 대 **큰 부분**의 비가 **큰 부분** 대 **작은 부분**의 비와 같다.

**유도:** 정의 관계식으로부터 다음을 얻는다:

```math
l_1 l_2 + l_2^2 = l_1^2
```

$`l_2^2`$으로 나누면 다음과 같다:

```math
\frac{l_1}{l_2} + 1 = \left(\frac{l_1}{l_2}\right)^2
```

```math
\phi^2 - \phi - 1 = 0
```

```math
\phi = \frac{1 \pm \sqrt{5}}{2}
```

$`\phi > 0`$을 취하면 다음을 얻는다:

```math
\phi = \frac{1 + \sqrt{5}}{2} \approx 1.61803
```

$`\phi^2 - \phi - 1 = 0`$의 양변을 $`\phi`$로 나누면 유용한 역수 관계가 나온다:

```math
\frac{1}{\phi} = \phi - 1 = \frac{\sqrt{5} - 1}{2} \approx 0.61803
```

### 6.2 알고리즘: 구간 축소

최솟값을 둘러싸는 구간 $`[x_l, x_u]`$가 주어지면, 다음을 정의한다:

```math
d = (\phi - 1)(x_u - x_l)
```

두 개의 내부점을 다음과 같이 선택한다:

```math
x_1^{(1)} = x_l + d = x_l(2 - \phi) + x_u(\phi - 1)
```

```math
x_2^{(1)} = x_u - d = x_u(2 - \phi) + x_l(\phi - 1)
```

$`x_2 < x_1`$이며, 두 내부점의 배치는 그림 1과 같다.

```text
x_l             x_2       x_1             x_u
 |---------------|---------|---------------|
 |<---------- d ---------->|
                 |<---------- d ---------->|
 |<------------ L = x_u - x_l ------------>|
```

*그림 1. 황금 분할 탐색의 두 내부점과 거리 d의 배치*

**결정 규칙:**

- $`f(x_2) > f(x_1)`$이면 최솟값은 $`[x_2, x_u]`$에 있으므로 $`x_l \leftarrow x_2^{(1)}`$로 설정한다.
- $`f(x_2) < f(x_1)`$이면 최솟값은 $`[x_l, x_1]`$에 있으므로 $`x_u \leftarrow x_1^{(1)}`$로 설정한다.

### 6.3 함수 평가의 재사용

이 방법의 장점은 각 반복 후 이전 단계의 두 내부점 중 하나가 새로운 내부점 중 하나와 일치한다는 것이다. 구체적으로, $`f(x_2) > f(x_1)`$이고 $`x_l \leftarrow x_2^{(1)}`$로 설정할 때 다음이 성립한다:

```math
x_2^{(2)} = x_u - d_{\text{new}} = x_1^{(1)}
```

이는 $`\phi^2 = \phi + 1`$ 성질을 사용하여 대수적으로 확인할 수 있다. 따라서 $`x_2^{(2)}`$에서의 함수 값은 $`f(x_1^{(1)})`$과 같으므로 **다시 계산할 필요가 없다**.

> **참고:** **효율성.** 황금 분할 탐색의 각 반복은 두 번이 아닌 **한 번의 새로운 함수 평가**만 필요하므로 매우 효율적이다. 구간은 매 단계 $`\frac{1}{\phi} \approx 0.61803`$의 일정한 비율로 축소된다.

### 6.4 오차 분석

한 번 반복한 뒤에는 두 내부점 가운데 함수 값이 더 작은 점을 최적점의 추정값으로 쓴다. $`f(x_2) > f(x_1)`$이면 추정값은 $`x_1`$이고, 참 최적점 $`x^*`$는 새 구간 $`[x_2, x_u]`$ 안에 있으며 다음 반복의 세 점은 $`x_2, x_1, x_u`$이다.

**경우 1:** $`x^*`$가 $`x_2`$와 $`x_1`$ 사이에 있으면($`x_2 \leq x^* \leq x_1`$), 추정값 $`x_1`$에서 $`x^*`$까지의 최대 거리는 다음과 같다:

```math
\Delta x_a = x_1 - x_2 = (2\phi - 3)(x_u - x_l) \approx 0.2361 \cdot (x_u - x_l)
```

**경우 2:** $`x^*`$가 $`x_1`$과 $`x_u`$ 사이에 있으면($`x_1 \leq x^* \leq x_u`$), 최대 거리는 다음과 같다:

```math
\Delta x_b = x_u - x_1 = (2 - \phi)(x_u - x_l) \approx 0.3820 \cdot (x_u - x_l)
```

$`f(x_2) < f(x_1)`$이면 추정값은 $`x_2`$이고 다음 반복의 세 점은 $`x_l, x_2, x_1`$이며, 대칭적으로 최대 거리는 $`x_2 - x_l = (2 - \phi)(x_u - x_l)`$ 또는 $`x_1 - x_2`$이다. 어느 경우든 최악의 오차는 $`(2 - \phi)(x_u - x_l)`$이므로, 이를 추정값의 크기로 나누어 **상대 오차(relative error)** 를 다음과 같이 정의한다:

```math
\varepsilon_a = (2 - \phi) \left| \frac{x_u - x_l}{x_{\text{opt}}} \right|
```

여기서 $`x_{\text{opt}}`$는 최적점의 현재 최선 추정값이다.

---

<br>

## 7. 포물선 보간(Parabolic Interpolation)

포물선 보간은 세 점을 지나는 이차식(포물선)을 맞추고, 포물선의 최적점을 다음 추정값으로 사용한다. 이 접근법은 매끄러운 함수에 대해 황금 분할 탐색보다 빠르게 수렴할 수 있다.

### 7.1 포물선 구성

세 점 $`(x_1, y_1)`$, $`(x_2, y_2)`$, $`(x_3, y_3)`$ ($`y_i = f(x_i)`$)이 주어지면, 다음 포물선 방정식을 구성한다:

```math
y = ax^2 + bx + c
```

이로부터 다음 선형 시스템을 얻는다:

```math
\begin{bmatrix} y_1 \\ y_2 \\ y_3 \end{bmatrix} = \begin{bmatrix} x_1^2 & x_1 & 1 \\ x_2^2 & x_2 & 1 \\ x_3^2 & x_3 & 1 \end{bmatrix} \begin{bmatrix} a \\ b \\ c \end{bmatrix}
```

**크래머 법칙(Cramer's rule)** 을 사용하여 $`D = \begin{vmatrix} x_1^2 & x_1 & 1 \\ x_2^2 & x_2 & 1 \\ x_3^2 & x_3 & 1 \end{vmatrix}`$이라 하면 다음을 얻는다:

```math
a = \frac{1}{D}\begin{vmatrix} y_1 & x_1 & 1 \\ y_2 & x_2 & 1 \\ y_3 & x_3 & 1 \end{vmatrix}, \quad b = \frac{1}{D}\begin{vmatrix} x_1^2 & y_1 & 1 \\ x_2^2 & y_2 & 1 \\ x_3^2 & y_3 & 1 \end{vmatrix}, \quad c = \frac{1}{D}\begin{vmatrix} x_1^2 & x_1 & y_1 \\ x_2^2 & x_2 & y_2 \\ x_3^2 & x_3 & y_3 \end{vmatrix}
```

### 7.2 최적점 공식

포물선의 최적점(꼭짓점)은 $`\frac{dy}{dx} = 0`$으로 놓아 구한다:

```math
2ax_4 + b = 0 \implies x_4 = -\frac{b}{2a}
```

크래머 법칙으로 구한 $`a`$, $`b`$를 대입하여 정리하면 다음과 같다:

```math
x_4 = x_2 - \frac{1}{2} \cdot \frac{(x_2 - x_1)^2(y_2 - y_3) - (x_2 - x_3)^2(y_2 - y_1)}{(x_2 - x_1)(y_2 - y_3) - (x_2 - x_3)(y_2 - y_1)}
```

> **참고:** **최대화** 문제의 경우 포물선은 아래로 열려야 하며($`a < 0`$), 꼭짓점이 근사 최댓값을 준다. **최소화** 문제의 경우 포물선은 위로 열려야 하며($`a > 0`$), 꼭짓점이 근사 최솟값을 준다.

### 7.3 다음 반복을 위한 점 갱신

$`x_4`$를 계산한 후, 다음 반복을 위한 새로운 세 점 $`(x_1, x_2, x_3)`$을 결정해야 한다. 결정은 $`x_4`$가 $`x_2`$에 대해 어디에 위치하는지에 따라 달라진다:

**$`x_4`$가 $`x_2`$와 $`x_3`$ 사이인 경우 (즉, $`x_2 \leq x_4 \leq x_3`$):**

| 조건 | 동작 |
|-----------|--------|
| $`f(x_4) < f(x_2)`$ (경우 4) | $`x_2 \Rightarrow x_1`$, $`x_4 \Rightarrow x_2`$, $`x_3`$ 유지 |
| $`f(x_4) \geq f(x_2)`$ (경우 3) | $`x_1`$ 유지, $`x_2`$ 유지, $`x_4 \Rightarrow x_3`$ |

*표 1. 새 점이 가운데 점의 오른쪽에 있을 때의 점 갱신*

**$`x_4`$가 $`x_1`$과 $`x_2`$ 사이인 경우 (즉, $`x_1 \leq x_4 \leq x_2`$):**

| 조건 | 동작 |
|-----------|--------|
| $`f(x_4) < f(x_2)`$ (경우 2) | $`x_2 \Rightarrow x_3`$, $`x_4 \Rightarrow x_2`$, $`x_1`$ 유지 |
| $`f(x_4) \geq f(x_2)`$ (경우 1) | $`x_4 \Rightarrow x_1`$, $`x_2`$ 유지, $`x_3`$ 유지 |

*표 2. 새 점이 가운데 점의 왼쪽에 있을 때의 점 갱신*

> **참고:** **수렴 속도.** 포물선 보간은 **초선형 수렴(superlinear convergence)** 을 보이며(차수 약 1.324), 이는 황금 분할 탐색의 선형 수렴보다 빠르다. 그러나 세 점이 최솟값을 잘 둘러싸지 않거나 맞춘 포물선이 퇴화하면 실패할 수 있다.

---

<br>

## 8. 등고선(Level Sets)과 기울기(Gradient)

최적화에서 등고선은 목적 함수를 분석하고 시각화하는 데 중요한 역할을 한다.

### 8.1 등고선 정의

**정의 (등고선).** 함수 $`f : \mathbb{R}^n \to \mathbb{R}`$과 상수 $`c`$에 대해 다음과 같이 정의한다:

- 집합 $`\{x \in \mathbb{R}^n : f(x) = c\}`$를 수준 $`c`$에서의 $`f`$의 **등고선(level set)** 이라 한다.
- 집합 $`\{x \in \mathbb{R}^n : f(x) \leq c\}`$를 $`c`$에서의 $`f`$의 **하위 등고선(lower level set)** 이라 한다.
- 집합 $`\{x \in \mathbb{R}^n : f(x) \geq c\}`$를 $`c`$에서의 $`f`$의 **상위 등고선(upper level set)** 이라 한다.

**예시:** $`f(x_1, x_2) = x_1^2 + x_2^2`$에 대해, $`c = 1`$에서의 등고선은 단위원 $`x_1^2 + x_2^2 = 1`$이다. 하위 등고선 $`f(x_1, x_2) \leq 1`$은 단위 원판이고, 상위 등고선 $`f(x_1, x_2) \geq 1`$은 단위 원판의 외부(경계인 단위원 포함)이다.

### 8.2 기울기는 등고선에 직교한다

등고선 $`S = \{x : f(x) = c\}`$ 위에 놓인 매개변수 곡선 $`\gamma = \{x(t) : t \in (a, b)\} \subset S`$를 생각한다. 여기서 $`x(t) : (a, b) \to S`$는 미분 가능한 함수이다.

$`t \in (a, b)`$에 대해 $`f(x(t)) = c`$이므로, $`t`$에 대해 미분하면 다음을 얻는다:

```math
\frac{df}{dt} = \nabla f^T \frac{dx}{dt} = 0
```

여기서 $`\nabla f`$와 $`\frac{dx}{dt}`$는 다음과 같다:

```math
\nabla f = \begin{bmatrix} \frac{\partial f}{\partial x_1} \\ \frac{\partial f}{\partial x_2} \end{bmatrix}, \quad \frac{dx}{dt} = \begin{bmatrix} \frac{dx_1}{dt} \\ \frac{dx_2}{dt} \end{bmatrix}
```

$`\frac{dx}{dt}`$는 곡선 $`x(t)`$의 접선 벡터이고 $`\nabla f^T \frac{dx}{dt} = 0`$이므로 다음이 성립한다:

```math
\nabla f \perp \frac{dx}{dt}
```

**기울기 $`\nabla f`$는 등고선에 직교한다.**

**예시:** $`f(x_1, x_2) = x_1^2 + x_2^2`$에 대해 기울기는 다음과 같다:

```math
\nabla f = \begin{bmatrix} 2x_1 \\ 2x_2 \end{bmatrix} = 2\begin{bmatrix} x_1 \\ x_2 \end{bmatrix}
```

기울기는 방사형으로 바깥을 향하며 원형 등고선에 수직이므로 직교성을 확인할 수 있다.

> **참고:** **기울기의 기하학적 의미.** 기울기가 0이 아닌 점에서 기울기는 가장 가파른 상승 방향을 가리킨다. 그런 점마다 등고선에 수직이며, $`f`$의 값이 더 큰 쪽, 즉 오르막을 가리킨다. 이것이 기울기 기반 최적화 방법의 기초이다.

---

<br>

## 9. 볼록 문제(Convex Problems)

### 9.1 볼록 결합(Convex Combination)

**정의 (볼록 결합).** $`\mathbb{R}^n`$의 점 $`x_1, x_2, \ldots, x_m`$과 $`\alpha_i \geq 0`$, $`\sum_{i=1}^m \alpha_i = 1`$인 실수가 주어질 때, 다음 점

```math
\sum_{i=1}^m \alpha_i x_i
```

을 이 점들의 **볼록 결합(convex combination)** 이라 한다.

### 9.2 볼록 집합(Convex Set)

**정의 (볼록 집합).** 집합 $`\mathbb{X} \subseteq \mathbb{R}^n`$이 **볼록(convex)** 이란 다음 관계가

```math
\alpha x + (1 - \alpha)y \in \mathbb{X}
```

모든 $`x, y \in \mathbb{X}`$와 모든 $`\alpha \in (0, 1)`$에 대해 성립하는 것이다.

기하학적으로, 집합 안의 임의의 두 점을 잇는 선분이 완전히 집합 안에 놓이면 그 집합은 볼록하다.

**예시:** $`\mathbb{X} = \{x \in \mathbb{R}^n : A_1 x = b_1, \; A_2 x \leq b_2\}`$가 볼록임을 보인다.

모든 $`x, y \in \mathbb{X}`$와 모든 $`\alpha \in (0, 1)`$에 대해 $`z = \alpha x + (1 - \alpha)y`$라 하면 다음이 성립한다:

```math
A_1 z = \alpha A_1 x + (1 - \alpha) A_1 y = \alpha b_1 + (1 - \alpha) b_1 = b_1
```

```math
A_2 z = \alpha A_2 x + (1 - \alpha) A_2 y \leq \alpha b_2 + (1 - \alpha) b_2 = b_2
```

따라서 $`z \in \mathbb{X}`$이며, $`\mathbb{X}`$는 볼록 집합이다.

### 9.3 초평면(Hyperplane)과 다면체 집합(Polyhedral Set)

**정의 (초평면).** 초평면 $`\mathbb{H} \subset \mathbb{R}^n`$은 다음 형태의 집합이다:

```math
\mathbb{H} = \{x \in \mathbb{R}^n : c^T x = b\}
```

여기서 $`c \in \mathbb{R}^n \setminus \{0\}`$이고 $`b \in \mathbb{R}`$이다.

초평면은 다음 두 **반공간(half-spaces)** 의 교집합이다:

```math
\mathbb{H}_+ = \{x \in \mathbb{R}^n : c^T x \geq b\}, \quad \mathbb{H}_- = \{x \in \mathbb{R}^n : c^T x \leq b\}
```

**$`\mathbb{H}`$가 볼록임의 증명:** $`x, y \in \mathbb{H}`$이면 $`c^T x = b`$이고 $`c^T y = b`$이다. $`z = \alpha x + (1 - \alpha) y`$에 대해 다음이 성립한다:

```math
c^T z = \alpha c^T x + (1 - \alpha) c^T y = \alpha b + (1 - \alpha)b = b
```

따라서 $`z \in \mathbb{H}`$이고, $`\mathbb{H}`$는 볼록하다.

**정의 (다면체 집합).** 유한 개의 선형 등식이나 부등식으로 정의된 집합을 **다면체 집합(polyhedral set)** (또는 다면체)이라 한다.

### 9.4 꼭짓점(Extreme Point)

**정의 (꼭짓점).** 볼록 집합 $`\mathbb{X}`$가 주어질 때, 점 $`x \in \mathbb{X}`$가 $`\mathbb{X}`$ 안의 두 **서로 다른** 점의 볼록 결합으로 표현될 수 없으면 그 점을 $`\mathbb{X}`$의 **꼭짓점(extreme point)** 이라 한다. 즉, $`x = \alpha x' + (1 - \alpha) x''`$인 서로 다른 점 $`x', x'' \in \mathbb{X}`$와 $`\alpha \in (0, 1)`$이 존재하지 않는다.

**예시:**

- 원(원판)의 경우 경계 위의 모든 점이 꼭짓점이다.
- 직사각형의 경우 네 **모퉁이(vertex)** 만이 꼭짓점이다. 변 위의 점(모퉁이 제외)은 두 끝점의 볼록 결합으로 쓸 수 있으므로 꼭짓점이 아니다.

실행 가능 영역이 닫히고 유계인 볼록 집합(예: 유계 다면체)이고 목적 함수가 선형이면, 최적해 가운데 적어도 하나는 꼭짓점에서 나타난다. 따라서 최적값을 찾을 때는 꼭짓점만 검사하면 된다.

### 9.5 볼록 함수(Convex Function)

**정의 (볼록 함수).** 볼록 집합 $`\mathbb{X} \subseteq \mathbb{R}^n`$이 주어질 때, 함수 $`f : \mathbb{X} \to \mathbb{R}`$이 다음을 만족하면 **볼록(convex)** 이라 한다:

```math
f(\alpha x + (1 - \alpha)y) \leq \alpha f(x) + (1 - \alpha) f(y) \quad \forall \; x, y \in \mathbb{X}, \; \alpha \in (0, 1)
```

$`f`$가 다음을 만족하면 **순볼록(strictly convex)** 이라 한다:

```math
f(\alpha x + (1 - \alpha)y) < \alpha f(x) + (1 - \alpha) f(y) \quad \forall \; x, y \in \mathbb{X}, \; x \neq y, \; \alpha \in (0, 1)
```

기하학적으로, 함수 그래프 위의 임의의 두 점을 잇는 선분이 항상 그래프와 같거나 그보다 높은 곳에 있으면 그 함수는 볼록하다.

**오목 함수(Concave functions):** $`-f`$가 볼록이면 함수 $`f`$를 **오목(concave)** 이라 한다.

**예시:**

- $`x > 0`$에서 $`f(x) = \frac{1}{x}`$는 **볼록**이다.
- $`x > 0`$에서 $`f(x) = \ln(x)`$는 **오목**이다.

**$`(0, \infty)`$에서 $`f(x) = \frac{1}{x}`$가 볼록임의 증명:** $`z = \alpha x_1 + (1 - \alpha)x_2`$일 때 $`\frac{1}{z} \leq \frac{\alpha}{x_1} + \frac{(1 - \alpha)}{x_2}`$를 보이면 된다.

양변에 $`x_1 x_2 z > 0`$을 곱하면 이는 $`x_1 x_2 \leq z(\alpha x_2 + (1 - \alpha)x_1)`$을 보이는 것과 동치이며, 정리하면 $`0 \leq (1 - \alpha)\alpha(x_1^2 - 2x_1 x_2 + x_2^2) = (1 - \alpha)\alpha(x_1 - x_2)^2`$이 되어 항상 참이다.

**볼록 함수의 성질:**

- 볼록 함수는 내부 $`\text{int}(\mathbb{X})`$에서 **연속**이다.
- 볼록 함수는 $`\text{int}(\mathbb{X})`$에서 **반드시 미분 가능한 것은 아니다**.
  - 예시: $`\mathbb{X} = \mathbb{R}`$에서 $`f(x) = |x|`$는 볼록이지만 $`x = 0`$에서 미분할 수 없다.

> **참고:** **볼록성의 중요성.** 최적화에서 볼록성의 중요성은 아무리 강조해도 지나치지 않다. 볼록 집합 위의 볼록 함수에서는 모든 국소 최소화자가 전역 최소화자이기도 하다. 이 성질 때문에 볼록 최적화 문제는 근본적으로 풀기 쉽다.

---

<br>

## 10. 최급 하강법(Steepest Gradient Descent Method)

### 10.1 방향 도함수와 코시-슈바르츠(Cauchy-Schwarz)

다음 비제약 최소화 문제를 생각한다:

```math
\min_{x \in \mathbb{R}^n} f(x)
```

여기서 $`f`$는 연속 미분 가능 함수이다.

$`x^{(0)} \in \mathbb{R}^n`$과 방향 $`d \in \mathbb{R}^n`$이 주어질 때, $`x^{(0)}`$에서 방향 $`d`$로의 $`f`$의 **방향 도함수(directional derivative)** 는 다음과 같다:

```math
\nabla f\big|_{x=x^{(0)}}^T d
```

> **[미적분학]** 방향 도함수는 점 $`x^{(0)}`$에서 방향 $`d`$로 움직일 때의 순간 변화율 $`\lim_{t \to 0} \frac{f(x^{(0)} + t d) - f(x^{(0)})}{t}`$으로 정의되며, $`f`$가 미분 가능하면 연쇄 법칙에 의해 기울기와 방향의 내적 $`\nabla f^T d`$와 같다. 따라서 가장 빠르게 감소하는 방향을 찾는 문제는 길이가 같은 $`d`$ 가운데 이 내적을 가장 작게 만드는 것을 찾는 문제가 되며, 아래에서는 코시-슈바르츠 부등식으로 그 방향을 정한다.

$`\|d\| = 1`$이면 이 값은 점 $`x^{(0)}`$에서 방향 $`d`$로 움직일 때 $`f`$의 증가율이다.

**코시-슈바르츠 부등식(Cauchy-Schwarz inequality)** 에 의해 다음이 성립한다:

```math
\nabla f^T d \leq \|\nabla f\| \cdot \|d\|
```

$`\alpha > 0`$에 대해 $`d = \alpha \nabla f`$를 선택하면 다음을 얻는다:

```math
\alpha \nabla f^T \nabla f = \alpha \|\nabla f\| \cdot \|\nabla f\|
```

이는 코시-슈바르츠 부등식에서 등호가 성립하는 경우이므로, $`\alpha \nabla f`$는 $`f`$의 **최대 증가율 방향**이다.

마찬가지로 **최대 감소율 방향**으로는 $`d = -\alpha \nabla f`$를 선택한다:

```math
-\nabla f^T d \leq \|\nabla f\| \cdot \|d\|
```

```math
\nabla f^T d \geq -\|\nabla f\| \cdot \|d\|
```

$`d = -\alpha \nabla f`$를 선택하면 $`-\alpha \nabla f^T \nabla f = -\alpha \|\nabla f\| \cdot \|\nabla f\|`$이므로 하한에 도달한다.

따라서 $`-\nabla f`$는 **최소화 방법에서 택할 최선의 방향**이다.

### 10.2 최급 하강 갱신 규칙

아이디어는 다음과 같이 반복적으로 갱신하는 것이다:

```math
x^{(k+1)} = x^{(k)} - \alpha_k \nabla f\big|_{x=x^{(k)}}
```

여기서 $`\alpha_k \geq 0`$은 다음을 만족하도록 선택한다:

```math
f(x^{(k+1)}) < f(x^{(k)}), \quad k \geq 0
```

스텝 크기 $`\alpha_k`$는 일반적으로 **직선 탐색(line search)** 으로 구한다. 즉, $`\alpha`$에 대해 $`\phi_k(\alpha) = f(x^{(k)} - \alpha \nabla f_k)`$를 최소화한다.

아래 구현은 직선 탐색에 6절의 황금 분할 탐색을 쓴다. 먼저 황금 분할 탐색을 다음과 같이 구현한다:

```python
import numpy as np

def golden_section_search(f, xl, xu, tol=1e-8, max_iter=200):
    """Minimize a unimodal function f on [xl, xu] by golden-section search."""
    phi = (1 + np.sqrt(5)) / 2
    d = (phi - 1) * (xu - xl)
    x1, x2 = xl + d, xu - d
    f1, f2 = f(x1), f(x2)
    for _ in range(max_iter):
        if xu - xl < tol:
            break
        if f2 > f1:                      # minimum in [x2, xu]
            xl, x2, f2 = x2, x1, f1      # old x1 becomes the new x2
            x1 = xl + (phi - 1) * (xu - xl)
            f1 = f(x1)
        else:                            # minimum in [xl, x1]
            xu, x1, f1 = x1, x2, f2      # old x2 becomes the new x1
            x2 = xu - (phi - 1) * (xu - xl)
            f2 = f(x2)
    return (xl + xu) / 2
```

최급 하강법은 이 함수로 매 반복의 스텝 크기를 구간 $`[0, 1]`$에서 찾는다:

```python
import numpy as np

def steepest_descent(f, grad_f, x0, tol=1e-8, max_iter=10000):
    """
    Steepest descent method with exact line search (golden-section).

    Parameters:
        f      : objective function
        grad_f : gradient function
        x0     : initial point (numpy array)
        tol    : tolerance for convergence
        max_iter: maximum number of iterations

    Returns:
        x      : approximate minimizer
        k      : number of iterations
    """
    x = np.array(x0, dtype=float)

    for k in range(max_iter):
        g = grad_f(x)
        if np.linalg.norm(g) < tol:
            break

        # Line search: minimize f(x - alpha * g)
        # Using golden-section search on alpha in [0, 1]
        phi = lambda a: f(x - a * g)
        alpha = golden_section_search(phi, 0, 1, tol=1e-10)

        x = x - alpha * g

    return x, k
```

직선 탐색 구간을 $`[0, 1]`$로 제한했으므로, 최적 스텝 크기가 1보다 큰 반복에서는 정확한 직선 탐색이 되지 않는다.

### 10.3 수렴과 직교성 성질

**정리.** $`\{x^{(k)} : k \geq 0\}`$을 정확한 직선 탐색을 쓰는 최급 하강법으로 생성한 수열이라 할 때, $`x^{(k)} \to x^*`$이면 다음이 성립한다:

```math
\nabla f(x^*) = 0
```

**증명 개요:** $`\phi_k(\alpha) = f(x^{(k)} - \alpha \nabla f_k)`$이고 $`\nabla f_k := \nabla f\big|_{x=x^{(k)}}`$라 하자. $`\phi_k(\alpha)`$를 최소화하는 $`\alpha`$에서는 다음이 성립한다:

```math
\frac{d\phi_k}{d\alpha} = \frac{df(x^{(k+1)})}{dx^{(k+1)}} \cdot \frac{dx^{(k+1)}}{d\alpha} = -\nabla f_{k+1}^T \nabla f_k = 0
```

이로부터 **직교성 성질(orthogonality property)** 이 도출된다:

```math
\nabla f_{k+1}^T \nabla f_k = 0, \quad k \geq 0
```

최급 하강법으로 생성된 **연속한 두 점**에서의 $`f`$의 기울기는 **서로 직교한다**.

$`f`$의 기울기가 연속이므로 $`k \to \infty`$로 극한을 취하면 $`\nabla f\big|_{x^*}^T \nabla f\big|_{x^*} = \|\nabla f(x^*)\|^2 = 0`$이다.

> **참고:** **지그재그 현상.** 직교성 성질은 최급 하강법 특유의 지그재그 패턴을 설명한다. 각 스텝이 이전 스텝에 수직이기 때문이다. 이 때문에 등고선이 길쭉한(악조건) 경우에는 수렴이 느려질 수 있는데, 방법이 골짜기를 가로질러 반복적으로 지나쳤다가 되돌아오기 때문이다.

### 10.4 예제

```math
f(x, y) = 3x^2 - 10x - 4xy + 2y^2 - 5y + 18
```

시작점은 $`(x^{(0)}, y^{(0)}) = (-2, 3)`$이고 허용 오차는 $`10^{-8}`$이다.

기울기는 다음과 같다:

```math
\nabla f = \begin{bmatrix} 6x - 10 - 4y \\ -4x + 4y - 5 \end{bmatrix}
```

최급 하강법은 최솟값을 향해 지그재그로 이동하는 반복값 수열을 생성하며, 연속한 기울기들은 서로 직교한다.

$`\nabla f = 0`$을 풀면 $`x = 7.5`$, $`y = 8.75`$이고, 이 점에서 $`f = -41.375`$이다. 헤세 행렬 $`\begin{bmatrix} 6 & -4 \\ -4 & 4 \end{bmatrix}`$는 첫 대각 성분이 양수이고 행렬식이 $`8 > 0`$이므로 양정치이며, 따라서 이 점이 유일한 최소점이다. 10.2절의 코드도 이 점으로 수렴한다.

---

<br>

## 요약

| 개념 | 핵심 요약 |
|-------|-------------|
| 최적화(Optimization) | 실행 가능한 해 중에서 최선(최솟값 또는 최댓값)을 찾는 방법론이다. |
| 실행 가능 영역(Feasible Region) | 등식 제약 조건($`h_i(x)=0`$)과 부등식 제약 조건($`g_j(x) \leq 0`$)으로 정의된다. |
| 최소화 vs. 최대화 | $`\min f(x) \iff \max(-f(x))`$이므로 모든 최대화 문제는 최소화 문제로 바꿔 쓸 수 있다. |
| 내부점과 경계점(Interior / Boundary) | 내부점은 어떤 $`\varepsilon`$-공과 함께 집합에 포함되는 점이고, 경계점은 집합에 속하지만 내부점이 아닌 점이다. |
| 전역 최솟값(Global Minimum) | 모든 $`x \in \mathbb{X}`$에 대해 $`f(x^*) \leq f(x)`$인 점 $`x^*`$이다. |
| 국소 최솟값(Local Minimum) | 어떤 이웃 안의 모든 $`x`$에 대해 $`f(x^*) \leq f(x)`$인 점 $`x^*`$이다. |
| 바이어슈트라스 정리(Weierstrass Theorem) | 공집합이 아닌 콤팩트 집합 $`\mathbb{X}`$ 위의 연속 함수 $`f`$는 전역 최적해를 가진다. |
| 강제 함수(Coercive Function) | $`\Vert x \Vert \to \infty`$일 때 $`f(x) \to \infty`$인 함수이며, 연속이면 $`\mathbb{R}^n`$에서 전역 최솟값을 가진다. |
| 황금 분할 탐색(Golden-Section Search) | 구간 기반 1차원 방법으로, 매 단계 구간을 $`\frac{1}{\phi} \approx 0.618`$배로 줄이고 반복마다 새 함수 평가를 한 번만 한다. |
| 포물선 보간(Parabolic Interpolation) | 세 점을 지나는 이차식의 꼭짓점 $`x_4 = -b/(2a)`$를 다음 추정값으로 쓰며 초선형 수렴한다. |
| 등고선(Level Set) | $`\{x : f(x) = c\}`$인 집합이며, 기울기는 등고선에 직교한다. |
| 볼록 집합(Convex Set) | 집합 안의 임의의 두 점을 잇는 선분이 집합 안에 머무는 집합이다. |
| 볼록 함수(Convex Function) | $`f(\alpha x + (1-\alpha)y) \leq \alpha f(x) + (1-\alpha)f(y)`$를 만족하며, 볼록 집합 위에서는 국소 최솟값이 곧 전역 최솟값이다. |
| 초평면(Hyperplane) | $`\{x : c^T x = b\}`$ 꼴의 집합이며 항상 볼록하다. |
| 꼭짓점(Extreme Point) | 집합 안의 서로 다른 두 점의 볼록 결합으로 표현할 수 없는 점이다. |
| 최급 하강법(Steepest Descent) | $`x^{(k+1)} = x^{(k)} - \alpha_k \nabla f_k`$로 갱신하며, $`-\nabla f`$가 가장 가파른 감소 방향이다. |
| 직교성 성질(Orthogonality Property) | 정확한 직선 탐색을 쓰면 연속한 최급 하강 기울기가 서로 수직이다: $`\nabla f_{k+1}^T \nabla f_k = 0`$. |

---

<br>

## 점검 문제

1. **황금 분할 탐색의 내부점:** 구간 $`[0, 4]`$에서 황금 분할 탐색을 시작할 때 두 내부점 $`x_1`$, $`x_2`$를 구하고, 한 번 반복한 뒤의 구간 길이를 구하라.

   <details>
   <summary>정답</summary>

   $`d = (\phi - 1)(4 - 0) \approx 0.61803 \times 4 = 2.4721`$이므로 $`x_1 = 0 + 2.4721 = 2.4721`$, $`x_2 = 4 - 2.4721 = 1.5279`$이다. 새 구간은 함수 값에 따라 $`[1.5279, 4]`$ 또는 $`[0, 2.4721]`$이며, 어느 쪽이든 길이는 $`d \approx 2.4721`$로 원래 길이의 약 0.618배이다.

   </details>

2. **포물선 보간:** $`f(x) = x^2 - 4x`$에 세 점 $`x_1 = 0`$, $`x_2 = 1`$, $`x_3 = 3`$을 써서 포물선 보간을 한 번 적용하라. $`x_4`$와 다음 반복의 세 점을 구하라.

   <details>
   <summary>정답</summary>

   $`y_1 = 0`$, $`y_2 = -3`$, $`y_3 = -3`$이다. 공식에 대입하면 분자는 $`(1)^2(0) - (-2)^2(-3) = 12`$, 분모는 $`(1)(0) - (-2)(-3) = -6`$이므로 $`x_4 = 1 - \frac{1}{2} \cdot \frac{12}{-6} = 2`$이다. $`f`$가 이차 함수이므로 한 번에 정확한 최소점 $`x = 2`$를 얻는다. $`x_4`$가 $`x_2`$와 $`x_3`$ 사이에 있고 $`f(x_4) = -4 < f(x_2) = -3`$이므로, 표 1의 경우 4에 따라 다음 세 점은 $`(1, 2, 3)`$이다.

   </details>

3. **국소 최소와 전역 최소:** 구간 $`[-3, 3]`$에서 $`f(x) = x^3 - 3x`$의 국소 최소화자와 전역 최소화자를 구하고, 전역 최소화자가 존재하는 이유를 설명하라.

   <details>
   <summary>정답</summary>

   $`f'(x) = 3x^2 - 3 = 0`$에서 $`x = \pm 1`$이고, $`f(1) = -2`$, $`f(-1) = 2`$이므로 내부의 국소 최소화자는 $`x = 1`$이다. 끝점에서는 $`f(-3) = -18`$, $`f(3) = 18`$이다. 따라서 전역 최소화자는 경계점 $`x = -3`$이고 $`f(-3) = -18`$이며, $`x = 1`$은 국소 최소화자이지만 전역 최소화자가 아니다. $`[-3, 3]`$은 공집합이 아닌 콤팩트 집합이고 $`f`$는 연속이므로 바이어슈트라스 정리에 의해 전역 최소화자가 존재한다.

   </details>

4. **최급 하강법 한 단계:** $`f(x, y) = x^2 + 4y^2`$에서 시작점 $`(2, 1)`$로 정확한 직선 탐색을 쓰는 최급 하강법을 한 번 적용하라. 스텝 크기, 새 점, 두 기울기의 내적을 구하라.

   <details>
   <summary>정답</summary>

   $`\nabla f = (2x, 8y)`$이므로 $`\nabla f(2, 1) = (4, 8)`$이다. $`\phi_0(\alpha) = (2 - 4\alpha)^2 + 4(1 - 8\alpha)^2`$에서 $`\phi_0'(\alpha) = -80 + 544\alpha = 0`$이므로 $`\alpha = 5/34 \approx 0.147`$이다. 새 점은 $`(2 - 20/34, \; 1 - 40/34) = (24/17, -3/17) \approx (1.412, -0.176)`$이고, 함수 값은 8에서 $`612/289 \approx 2.118`$로 줄어든다. 새 기울기는 $`(48/17, -24/17)`$이고 $`(4, 8) \cdot (48/17, -24/17) = (192 - 192)/17 = 0`$이므로 직교성 성질이 성립한다.

   </details>
