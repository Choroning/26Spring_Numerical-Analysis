# 13장 — 고유값 방법: 대칭 행렬

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **최종 수정일:** 2026-10-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 13

> **학습 목표:**
> 1. 특성 방정식으로 고유값과 고유벡터를 구할 수 있다
> 2. 고유값의 위치로 선형 ODE 시스템의 안정성과 진동 여부를 판정할 수 있다
> 3. 결합 진동계와 질량-스프링 시스템을 고유값 문제로 세우고 고유진동수와 모드를 구할 수 있다
> 4. 거듭제곱법과 역거듭제곱법으로 가장 큰 고유값과 가장 작은 고유값을 구할 수 있다

---

<br>

## 목차

- [1. 고유값과 고유벡터(Eigenvalues and Eigenvectors)](#1-고유값과-고유벡터eigenvalues-and-eigenvectors)
  - [1.1 정의](#11-정의)
  - [1.2 기하학적 해석](#12-기하학적-해석)
  - [1.3 예제: 투영 행렬(Projection Matrix)](#13-예제-투영-행렬projection-matrix)
  - [1.4 예제: 반사 행렬(Reflection Matrix)](#14-예제-반사-행렬reflection-matrix)
- [2. 고유값과 고유벡터 계산](#2-고유값과-고유벡터-계산)
  - [2.1 특성 방정식(Characteristic Equation)](#21-특성-방정식characteristic-equation)
  - [2.2 2x2 경우](#22-2x2-경우)
  - [2.3 예제: 고유값과 고유벡터 구하기](#23-예제-고유값과-고유벡터-구하기)
- [3. 응용: 미분방정식](#3-응용-미분방정식)
  - [3.1 1차 ODE: 스칼라 경우](#31-1차-ode-스칼라-경우)
  - [3.2 2차를 1차 시스템으로 변환](#32-2차를-1차-시스템으로-변환)
  - [3.3 행렬 지수함수와 안정성](#33-행렬-지수함수와-안정성)
  - [3.4 고유값 위치와 시스템 거동](#34-고유값-위치와-시스템-거동)
- [4. 예제 13.2: 1차 시스템의 안정성](#4-예제-132-1차-시스템의-안정성)
- [5. 고유값과 상미분방정식](#5-고유값과-상미분방정식)
  - [5.1 순수 진동: 단일 2차 ODE](#51-순수-진동-단일-2차-ode)
  - [5.2 진동 시스템의 고유벡터](#52-진동-시스템의-고유벡터)
- [6. 결합 2차 ODE 시스템](#6-결합-2차-ode-시스템)
  - [6.1 일반 형태](#61-일반-형태)
  - [6.2 주기적 해의 제안](#62-주기적-해의-제안)
  - [6.3 2차 시스템의 특성 방정식](#63-2차-시스템의-특성-방정식)
- [7. 예제 13.3: 결합 시스템의 거동](#7-예제-133-결합-시스템의-거동)
  - [7.1 고유값 계산](#71-고유값-계산)
  - [7.2 고유벡터 계산](#72-고유벡터-계산)
  - [7.3 진동수와 일반해](#73-진동수와-일반해)
  - [7.4 초기 조건 적용](#74-초기-조건-적용)
- [8. 물리적 설정: 질량-스프링 시스템](#8-물리적-설정-질량-스프링-시스템)
  - [8.1 두 질량, 세 스프링 시스템](#81-두-질량-세-스프링-시스템)
  - [8.2 질량-스프링 시스템의 고유값 공식](#82-질량-스프링-시스템의-고유값-공식)
- [9. 예제 13.4: 질량-스프링 시스템](#9-예제-134-질량-스프링-시스템)
- [10. 거듭제곱법(Power Method)](#10-거듭제곱법power-method)
  - [10.1 세 질량, 네 스프링 시스템 (예제 13.5)](#101-세-질량-네-스프링-시스템-예제-135)
  - [10.2 알고리즘](#102-알고리즘)
  - [10.3 반복 과정](#103-반복-과정)
  - [10.4 수렴과 최소 고유값 찾기](#104-수렴과-최소-고유값-찾기)
- [11. Python 구현](#11-python-구현)
  - [11.1 NumPy를 이용한 고유값 계산](#111-numpy를-이용한-고유값-계산)
  - [11.2 거듭제곱법 구현](#112-거듭제곱법-구현)
  - [11.3 질량-스프링 시스템 (예제 13.4)](#113-질량-스프링-시스템-예제-134)
  - [11.4 최소 고유값을 위한 역거듭제곱법](#114-최소-고유값을-위한-역거듭제곱법)
- [요약](#요약)
- [점검 문제](#점검-문제)

---

<br>

## 1. 고유값과 고유벡터(Eigenvalues and Eigenvectors)

### 1.1 정의

**정방행렬(square matrix)** 을 다룬다. $`n \times n`$ 행렬 $`[A]`$가 주어졌을 때 다음 방정식을 생각한다.

```math
[A]\{x\} = \lambda \{x\}
```

여기서 $`\{x\}`$는 $`n \times 1`$ 벡터이고 $`\lambda \in \mathbb{R}`$ 또는 $`\lambda \in \mathbb{C}`$이다.

임의의 $`\lambda`$ 값에 대해 **자명해(trivial solution)** $`\{x\} = \{0\}`$는 항상 존재한다.

$`[A]\{x\} = \lambda \{x\}`$를 만족하는 **비자명해(non-trivial solution)** $`\{x\} \neq \{0\}`$가 존재하면 다음과 같이 부른다.

- $`\{x\}`$를 $`[A]`$의 **고유벡터(eigenvector)** 라 한다.
- $`\lambda`$를 $`\{x\}`$에 대응하는 **고유값(eigenvalue)** 이라 한다.

> **[선형대수학]** "eigen"은 독일어로 "고유한" 또는 "특성적인"을 의미한다. 고유벡터는 행렬이 단순히 스칼라 배로 작용하는 특별한 방향이다. $`n \times n`$ 행렬은 특성 다항식의 근으로서 중복도를 포함해 정확히 $`n`$개의 (복소) 고유값을 가지며, 서로 다른 고유값은 최대 $`n`$개이다.

### 1.2 기하학적 해석

고유벡터 $`\{x\}`$는 $`[A]\{x\}`$와 **평행**하다. 즉, $`\{x\}`$에 행렬 $`[A]`$를 곱해도 원점을 지나는 같은 직선 위에 머물며 $`\lambda`$배로 스케일될 뿐이다. $`\lambda > 0`$이면 같은 방향, $`\lambda < 0`$이면 반대 방향을 가리키고, $`\lambda = 0`$이면 영벡터가 된다.

```math
\lambda \{x\} \quad \longrightarrow \quad \text{(} [A]\{x\}\text{와 같은 방향)}
```

### 1.3 예제: 투영 행렬(Projection Matrix)

$`P`$가 평면 위로의 투영 행렬이라고 하자.

**경우 1:** **평면 위에** 놓인 벡터 $`\{x\}`$에 대해서는 다음이 성립한다.

```math
[P]\{x\} = \{x\} = 1 \cdot \{x\}
```

여기서 $`\lambda = 1`$이 고유값이고, $`\{x\}`$(평면 위의 영이 아닌 임의의 벡터)가 고유벡터이다.

**경우 2:** 벡터 $`\{x\}`$가 평면에 **수직**이면 다음이 성립한다.

```math
[P]\{x\} = \{0\} = 0 \cdot \{x\}
```

여기서 $`\lambda = 0`$이 고유값이고, $`\{x\}`$(평면의 법선 방향 벡터)가 고유벡터이다.

### 1.4 예제: 반사 행렬(Reflection Matrix)

다음 반사 행렬을 생각한다.

```math
A = \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}
```

이 행렬은 직선 $`x_1 = x_2`$에 대한 반사 사상이다. 고유벡터는 $`\begin{Bmatrix} 1 \\ 1 \end{Bmatrix}`$(고유값 1)과 $`\begin{Bmatrix} 1 \\ -1 \end{Bmatrix}`$(고유값 $`-1`$)이다.

두 고유벡터가 **전체 공간을 생성(span)** 하므로 $`\mathbb{R}^2`$의 임의의 벡터를 이 고유벡터들의 선형결합으로 쓸 수 있다.

> **[선형대수학]** $`n \times n`$ 행렬의 고유벡터들이 전체 $`\mathbb{R}^n`$을 생성하면 그 행렬을 **대각화 가능(diagonalizable)** 하다고 한다. 실수 대칭 행렬은 항상 직교 행렬로 대각화할 수 있고 실수 고유값을 갖는데, 이것이 **스펙트럼 정리(spectral theorem)** 라는 중요한 결과이다.

---

<br>

## 2. 고유값과 고유벡터 계산

### 2.1 특성 방정식(Characteristic Equation)

고유값 방정식에서 출발한다.

```math
[A]\{x\} = \lambda \{x\}
```

재정렬하면 다음과 같다.

```math
[A]\{x\} - \lambda [I]\{x\} = \{0\}
```

```math
[A - \lambda I]\{x\} = \{0\}
```

**비자명해** $`\{x\} \neq \{0\}`$가 존재하려면 행렬 $`[A - \lambda I]`$가 비가역, 즉 **특이(singular)** 해야 한다. 이를 위한 조건은 다음과 같다.

```math
\det(A - \lambda I) = |A - \lambda I| = 0
```

이 식을 **특성 방정식(characteristic equation)** 이라 하고, 좌변의 행렬식을 $`\lambda`$의 다항식으로 본 것을 특성 다항식이라 한다. $`n \times n`$ 행렬에서 특성 다항식은 $`\lambda`$에 관한 $`n`$차 다항식이다.

### 2.2 2x2 경우

$`2 \times 2`$ 행렬에서는 다음과 같다.

```math
\begin{vmatrix} a_{11} - \lambda & a_{12} \\ a_{21} & a_{22} - \lambda \end{vmatrix} = 0
```

전개하면 다음과 같다.

```math
\lambda^2 - (a_{11} + a_{22})\lambda + (a_{11}a_{22} - a_{12}a_{21}) = 0
```

근의 공식을 사용하면 다음을 얻는다.

```math
\lambda = \frac{(a_{11} + a_{22}) \pm \sqrt{(a_{11} + a_{22})^2 - 4(a_{11}a_{22} - a_{12}a_{21})}}{2}
```

> **[선형대수학]** $`a_{11} + a_{22}`$는 $`A`$의 **대각합(trace)** 이고, $`a_{11}a_{22} - a_{12}a_{21}`$는 $`A`$의 **행렬식(determinant)** 이다. $`2 \times 2`$ 행렬에서 고유값의 합은 대각합과 같고 고유값의 곱은 행렬식과 같으며, 이 관계는 일반적인 $`n \times n`$ 행렬에서도 성립한다(12장 실습 1.1절).

### 2.3 예제: 고유값과 고유벡터 구하기

다음 행렬을 생각한다.

```math
A = \begin{bmatrix} 3 & 1 \\ 1 & 3 \end{bmatrix}
```

**단계 1: 특성 방정식**

```math
|A - \lambda I| = \begin{vmatrix} 3 - \lambda & 1 \\ 1 & 3 - \lambda \end{vmatrix} = (3 - \lambda)^2 - 1 = \lambda^2 - 6\lambda + 8 = (\lambda - 4)(\lambda - 2) = 0
```

따라서 $`\lambda_1 = 4`$, $`\lambda_2 = 2`$이다.

**단계 2: $`\lambda_1 = 4`$의 고유벡터**

```math
\begin{bmatrix} 3 - 4 & 1 \\ 1 & 3 - 4 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}
```

```math
-x_1 + x_2 = 0 \quad \text{와} \quad x_1 - x_2 = 0 \quad \Rightarrow \quad \text{같은 방정식}
```

고유벡터는 상수배를 제외하고 정해지므로 $`x_1 = 1`$로 놓으면 $`x_2 = 1`$이고, 고유벡터는 다음과 같다.

```math
\{x\}_1 = \begin{Bmatrix} 1 \\ 1 \end{Bmatrix}
```

**단계 3: $`\lambda_2 = 2`$의 고유벡터**

```math
\begin{bmatrix} 3 - 2 & 1 \\ 1 & 3 - 2 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}
```

```math
x_1 + x_2 = 0
```

$`x_1 = 1`$로 놓으면 $`x_2 = -1`$이고, 고유벡터는 다음과 같다.

```math
\{x\}_2 = \begin{Bmatrix} 1 \\ -1 \end{Bmatrix}
```

---

<br>

## 3. 응용: 미분방정식

### 3.1 1차 ODE: 스칼라 경우

다음 스칼라 ODE를 생각한다.

```math
\frac{dx}{dt} = ax, \quad a \in \mathbb{R}
```

해는 다음과 같다.

```math
x(t) = e^{at} C, \quad C \text{는 상수}
```

$`a`$의 값에 따른 거동은 다음 표와 같다. 마지막 행은 $`a`$를 복소수로 확장한 경우이다.

| $`a`$의 값 | 해의 거동 |
|---|---|
| $`a > 0`$ (실수, 양수) | 지수적 성장 (불안정) |
| $`a < 0`$ (실수, 음수) | 지수적 감쇠 (안정) |
| $`a = 0`$ | 상수: $`x(t) = C`$ |
| $`a = i\beta`$ (허수) | 진동: $`x(t) = e^{i\beta t}C = (\cos \beta t + i\sin \beta t)C`$ |

*표 1. 계수 a의 값에 따른 스칼라 ODE 해의 거동*

### 3.2 2차를 1차 시스템으로 변환

상수 계수를 가진 다음 2차 ODE를 생각한다.

```math
\frac{d^2 y}{dt^2} + a\frac{dy}{dt} + by = 0
```

이 식은 다음 치환으로 **1차** ODE 시스템으로 변환할 수 있다.

```math
x_1 = y, \quad x_2 = \frac{dx_1}{dt}
```

이 치환으로 다음을 얻는다.

```math
\frac{dx_1}{dt} = x_2
```

```math
\frac{dx_2}{dt} = -bx_1 - ax_2
```

행렬 형태로 쓰면 다음과 같다.

```math
\frac{d}{dt}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} 0 & 1 \\ -b & -a \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix}
```

하나의 2차 ODE를 **두 개의 1차** ODE로 변환하였다.

### 3.3 행렬 지수함수와 안정성

다음 벡터 ODE를 생각한다.

```math
\frac{d\mathbf{x}}{dt} = A\mathbf{x}, \quad \mathbf{x} \in \mathbb{R}^n, \quad A \in \mathbb{R}^{n \times n}
```

해는 다음과 같다.

```math
\mathbf{x}(t) = e^{At}\mathbf{C}, \quad \mathbf{C} \in \mathbb{R}^n
```

여기서 $`e^{At}`$는 **행렬 지수함수(matrix exponential)** 이다. $`A`$가 대칭이면 고유값의 대각 행렬 $`D`$와 정규직교 고유벡터를 열로 갖는 행렬 $`U`$로 $`A = U D U^T`$와 같이 대각화할 수 있고, 이때 다음이 성립한다(일반적인 대각화 가능 행렬에서는 $`U^T`$ 대신 $`U^{-1}`$을 쓴다).

```math
e^{At} = U \, e^{Dt} \, U^T
```

이 식은 다음을 의미한다.

- **음의 실수 고유값**은 안정적인(감쇠하는) 해를 만든다.
- **양의 실수 고유값**은 불안정한(성장하는) 해를 만든다.
- **허수 고유값**은 진동하는 해를 만든다.

### 3.4 고유값 위치와 시스템 거동

복소평면은 안정성을 직관적으로 보여 주는 지도이다.

| 영역 | 고유값 유형 | 시스템 거동 |
|---|---|---|
| 좌반면 ($`\text{Re}(\lambda) < 0`$) | 음의 실수부 | **안정** (해가 감쇠) |
| 우반면 ($`\text{Re}(\lambda) > 0`$) | 양의 실수부 | **불안정** (해가 성장) |
| 허수축 ($`\text{Re}(\lambda) = 0`$) | 순허수 | **진동** (유계 주기 운동) |
| 왼쪽 멀리 | 큰 음의 실수 | 더 빠른 수렴 |
| 오른쪽 멀리 | 큰 양의 실수 | 더 빠른 발산 |

*표 2. 복소평면에서 고유값의 위치와 시스템 거동*

---

<br>

## 4. 예제 13.2: 1차 시스템의 안정성

다음 시스템을 생각한다.

```math
\frac{dx_1}{dt} = -3x_1 + x_2
```

```math
\frac{dx_2}{dt} = x_1 - 3x_2
```

행렬 형태로 쓰면 다음과 같다.

```math
\frac{d}{dt}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} -3 & 1 \\ 1 & -3 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix}
```

**특성 방정식:**

```math
\begin{vmatrix} -3 - \lambda & 1 \\ 1 & -3 - \lambda \end{vmatrix} = (3 + \lambda)^2 - 1 = \lambda^2 + 6\lambda + 8 = (\lambda + 4)(\lambda + 2) = 0
```

```math
\therefore \lambda = -2 \text{ 또는 } -4
```

**두 고유값이 모두 음의 실수**이므로 해는 **안정적**이며, 허수부가 없으므로 진동 없이 0으로 소멸한다.

---

<br>

## 5. 고유값과 상미분방정식

### 5.1 순수 진동: 단일 2차 ODE

비감쇠 진동자를 생각한다($`a > 0`$).

```math
\ddot{y} = -ay \quad \Longleftrightarrow \quad \frac{d^2 y}{dt^2} = -ay
```

$`x_1 = y`$, $`x_2 = \frac{dx_1}{dt}`$로 놓아 1차 ODE로 바꾸면 다음과 같다.

```math
\frac{d}{dt}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} 0 & 1 \\ -a & 0 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix}
```

**$`A`$의 고유값 구하기:**

```math
|A - \lambda I| = \begin{vmatrix} -\lambda & 1 \\ -a & -\lambda \end{vmatrix} = \lambda^2 + a = 0
```

```math
\therefore \lambda = \pm \sqrt{-a} = \pm \sqrt{a}\, i
```

고유값이 **순허수**(허수축에 위치)이므로 해는 진동수 $`\sqrt{a}/2\pi`$로 **진동**한다.

오일러 공식과 주기, 진동수의 관계는 다음과 같다.

```math
e^{i\sqrt{a}\,t} = \cos(\sqrt{a}\,t) + i\sin(\sqrt{a}\,t)
```

```math
\text{주기} = \frac{2\pi}{\sqrt{a}}, \qquad \text{진동수} = \frac{1}{\text{주기}} = \frac{\sqrt{a}}{2\pi}
```

> **참고:** **고유값과 진동수.** 고유값은 **진동의 진동수**와 직접 관련된다. 이것이 진동 시스템을 분석할 때마다 물리학과 공학 전반에서 고유값 문제가 나타나는 이유이다.

### 5.2 진동 시스템의 고유벡터

$`[A - \lambda I]\{x\} = \{0\}`$를 푼다.

```math
\begin{bmatrix} -\lambda & 1 \\ -a & -\lambda \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}
```

**경우 i)** $`\lambda = \sqrt{a}\, i`$:

```math
-\sqrt{a}\,i \cdot x_1 + x_2 = 0
```

$`x_1 = 1`$로 놓으면 $`x_2 = \sqrt{a}\,i`$이다.

```math
\{x\} = \begin{Bmatrix} 1 \\ \sqrt{a}\,i \end{Bmatrix}
```

**경우 ii)** $`\lambda = -\sqrt{a}\, i`$:

```math
\sqrt{a}\,i \cdot x_1 + x_2 = 0
```

$`x_1 = 1`$로 놓으면 $`x_2 = -\sqrt{a}\,i`$이다.

```math
\{x\} = \begin{Bmatrix} 1 \\ -\sqrt{a}\,i \end{Bmatrix}
```

---

<br>

## 6. 결합 2차 ODE 시스템

### 6.1 일반 형태

결합된 두 2차 ODE의 시스템은 다음과 같다.

```math
\frac{d^2 y_1}{dt^2} = -a_{11}y_1 + a_{12}y_2
```

```math
\frac{d^2 y_2}{dt^2} = a_{21}y_1 - a_{22}y_2
```

행렬 형태로 쓰면 다음과 같다.

```math
\frac{d^2}{dt^2}\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} = \begin{bmatrix} -a_{11} & a_{12} \\ a_{21} & -a_{22} \end{bmatrix} \begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} \quad \cdots \quad (*)
```

### 6.2 주기적 해의 제안

주기적인 해의 형태를 바탕으로 다음 해를 가정한다.

```math
\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} = \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix} e^{i\omega t}
```

여기서 $`\omega`$는 각진동수이다. $`(*)`$에 대입하면 다음과 같다.

```math
-\omega^2 \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix} e^{i\omega t} = \begin{bmatrix} -a_{11} & a_{12} \\ a_{21} & -a_{22} \end{bmatrix} \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix} e^{i\omega t}
```

$`e^{i\omega t}`$를 소거하면 다음과 같다.

```math
\underbrace{-\omega^2}_{\lambda} \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix} = \underbrace{\begin{bmatrix} -a_{11} & a_{12} \\ a_{21} & -a_{22} \end{bmatrix}}_{A} \begin{Bmatrix} X_1 \\ X_2 \end{Bmatrix}
```

**이것은 고유값 문제이다.** 고유값은 $`\lambda = -\omega^2`$이므로 $`\omega^2 = -\lambda`$이다.

### 6.3 2차 시스템의 특성 방정식

$`|A - \lambda I| = 0`$에서 다음을 얻는다.

```math
\begin{vmatrix} -a_{11} - \lambda & a_{12} \\ a_{21} & -a_{22} - \lambda \end{vmatrix} = 0
```

```math
(a_{11} + \lambda)(a_{22} + \lambda) - a_{12}a_{21} = 0
```

```math
\lambda^2 + (a_{11} + a_{22})\lambda + a_{11}a_{22} - a_{12}a_{21} = 0
```

```math
\therefore \lambda = \frac{-(a_{11} + a_{22}) \pm \sqrt{(a_{11} + a_{22})^2 - 4(a_{11}a_{22} - a_{12}a_{21})}}{2}
```

$`\omega^2 = -\lambda`$로 표현하면 다음과 같다.

```math
|A + \omega^2 I| = \begin{vmatrix} -a_{11} + \omega^2 & a_{12} \\ a_{21} & -a_{22} + \omega^2 \end{vmatrix} = 0
```

```math
(\omega^2 - a_{11})(\omega^2 - a_{22}) - a_{12}a_{21} = 0
```

```math
\omega^4 - (a_{11} + a_{22})\omega^2 + a_{11}a_{22} - a_{12}a_{21} = 0
```

```math
\therefore \omega^2 = \frac{(a_{11} + a_{22}) \pm \sqrt{(a_{11} + a_{22})^2 - 4(a_{11}a_{22} - a_{12}a_{21})}}{2}
```

---

<br>

## 7. 예제 13.3: 결합 시스템의 거동

### 7.1 고유값 계산

주어진 시스템은 다음과 같다.

```math
\frac{d^2 y_1}{dt^2} = -5y_1 + 2y_2, \qquad \frac{d^2 y_2}{dt^2} = 2y_1 - 2y_2
```

행렬 형태는 다음과 같다.

```math
\frac{d^2}{dt^2}\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} = \underbrace{\begin{bmatrix} -5 & 2 \\ 2 & -2 \end{bmatrix}}_{A}\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix}
```

특성 방정식은 다음과 같다.

```math
\begin{vmatrix} -5 - \lambda & 2 \\ 2 & -2 - \lambda \end{vmatrix} = (5 + \lambda)(2 + \lambda) - 4 = \lambda^2 + 7\lambda + 6 = (\lambda + 6)(\lambda + 1) = 0
```

```math
\therefore \lambda_1 = -1, \quad \lambda_2 = -6
```

### 7.2 고유벡터 계산

**$`\lambda_1 = -1`$의 경우:**

```math
\begin{bmatrix} -5 + 1 & 2 \\ 2 & -2 + 1 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix} \quad \Rightarrow \quad -4x_1 + 2x_2 = 0 \quad \Rightarrow \quad 2x_1 - x_2 = 0
```

$`x_1 = 1`$, $`x_2 = 2`$로 놓으면 고유벡터는 다음과 같다.

```math
\{x\}_1 = \begin{Bmatrix} 1 \\ 2 \end{Bmatrix}
```

**$`\lambda_2 = -6`$의 경우:**

```math
\begin{bmatrix} -5 + 6 & 2 \\ 2 & -2 + 6 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix} \quad \Rightarrow \quad x_1 + 2x_2 = 0
```

$`x_2 = -1`$, $`x_1 = 2`$로 놓으면 고유벡터는 다음과 같다.

```math
\{x\}_2 = \begin{Bmatrix} 2 \\ -1 \end{Bmatrix}
```

### 7.3 진동수와 일반해

고유값에서 $`\omega^2 = -\lambda`$를 구한다.

```math
\omega_1^2 = 1 \quad \Rightarrow \quad \omega_1 = \pm 1
```

```math
\omega_2^2 = 6 \quad \Rightarrow \quad \omega_2 = \pm \sqrt{6}
```

해의 거동은 다음과 같다.

```math
\{x\}_1 e^{\pm it} = \{x\}_1 (\cos t \pm i\sin t)
```

```math
\{x\}_2 e^{\pm i\sqrt{6}t} = \{x\}_2 (\cos\sqrt{6}t \pm i\sin\sqrt{6}t)
```

일반해는 다음과 같다.

```math
\mathbf{y} = C_1 \{x\}_1 e^{+it} + C_2 \{x\}_2 e^{+i\sqrt{6}t} + C_3 \{x\}_1 e^{-it} + C_4 \{x\}_2 e^{-i\sqrt{6}t}
```

시간 도함수는 다음과 같다.

```math
\dot{\mathbf{y}} = iC_1 \{x\}_1 e^{+it} + i\sqrt{6}C_2 \{x\}_2 e^{+i\sqrt{6}t} - iC_3 \{x\}_1 e^{-it} - i\sqrt{6}C_4 \{x\}_2 e^{-i\sqrt{6}t}
```

### 7.4 초기 조건 적용

주어진 초기 조건은 다음과 같다.

```math
\mathbf{y}(t=0) = \begin{Bmatrix} 1 \\ -1 \end{Bmatrix}, \quad \dot{\mathbf{y}}(t=0) = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}
```

$`t = 0`$에서 속도는 $`i(C_1 - C_3)\{x\}_1 + i\sqrt{6}(C_2 - C_4)\{x\}_2 = \{0\}`$이고 두 고유벡터는 일차독립이므로, $`C_1 = C_3`$이고 $`C_2 = C_4`$이어야 한다. $`\tilde{C}_1 = C_1 = C_3`$, $`\tilde{C}_2 = C_2 = C_4`$로 놓는다.

$`t = 0`$을 위치 방정식에 대입하면 $`\mathbf{y}(0) = 2\tilde{C}_1\{x\}_1 + 2\tilde{C}_2\{x\}_2`$이다.

$`y_1`$에서: $`1 = 2\tilde{C}_1 + 4\tilde{C}_2`$

$`y_2`$에서: $`-1 = 4\tilde{C}_1 - 2\tilde{C}_2`$

두 번째 방정식에 2를 곱해 첫 번째 방정식과 나란히 쓰면 다음과 같다.

```math
-2 = 8\tilde{C}_1 - 4\tilde{C}_2
```

```math
1 = 2\tilde{C}_1 + 4\tilde{C}_2
```

두 식을 더하면 $`-1 = 10\tilde{C}_1`$이므로 $`\tilde{C}_1 = -\frac{1}{10}`$이다.

이를 $`1 = 2(-\frac{1}{10}) + 4\tilde{C}_2`$에 대입하면 $`\tilde{C}_2 = \frac{1}{4}(1 + \frac{1}{5}) = \frac{3}{10}`$이다.

**최종 해:** $`e^{it} + e^{-it} = 2\cos t`$를 이용하면 다음과 같다.

```math
\begin{Bmatrix} y_1 \\ y_2 \end{Bmatrix} = -\frac{1}{5}\begin{Bmatrix} 1 \\ 2 \end{Bmatrix}\cos(t) + \frac{3}{5}\begin{Bmatrix} 2 \\ -1 \end{Bmatrix}\cos(\sqrt{6}\,t)
```

---

<br>

## 8. 물리적 설정: 질량-스프링 시스템

### 8.1 두 질량, 세 스프링 시스템

마찰 없는 롤러 위에서 두 벽 사이에 놓여 진동하는 두 질량, 세 스프링 시스템을 생각한다.

```text
벽 ─[k]─ m₁ ─[k]─ m₂ ─[k]─ 벽
```

*그림 1. 두 벽 사이의 두 질량, 세 스프링 시스템*

두 질량 $`m_1`$과 $`m_2`$가 강성 $`k`$의 스프링으로 연결되어 있다. 각 질량에 뉴턴의 제2법칙을 적용하면 다음과 같다.

```math
m_1 \frac{d^2 x_1}{dt^2} = -kx_1 + k(x_2 - x_1)
```

```math
m_2 \frac{d^2 x_2}{dt^2} = -k(x_2 - x_1) - kx_2
```

행렬 형태는 다음과 같다.

```math
\frac{d^2}{dt^2}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} -2k/m_1 & k/m_1 \\ k/m_2 & -2k/m_2 \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix}
```

### 8.2 질량-스프링 시스템의 고유값 공식

$`|A - \lambda I| = 0`$에서 다음을 얻는다.

```math
\begin{vmatrix} -2k/m_1 - \lambda & k/m_1 \\ k/m_2 & -2k/m_2 - \lambda \end{vmatrix} = 0
```

```math
(\lambda + 2k/m_1)(\lambda + 2k/m_2) - k^2/(m_1 m_2) = 0
```

```math
\lambda^2 + \lambda(2k/m_1 + 2k/m_2) + 4k^2/(m_1 m_2) - k^2/(m_1 m_2) = 0
```

근의 공식으로 풀면 다음과 같다.

```math
-\omega^2 = \lambda = -\frac{k}{m_1 m_2}(m_1 + m_2) \pm \frac{k}{m_1 m_2}\sqrt{m_2^2 + m_1^2 - m_1 m_2}
```

---

<br>

## 9. 예제 13.4: 질량-스프링 시스템

주어진 조건은 $`m_1 = m_2 = 40`$ kg, $`k = 200`$ N/m이다.

**고유값 계산:**

```math
\lambda = -\frac{k}{m_1 m_2}(m_1 + m_2) \pm \frac{k}{m_1 m_2}\sqrt{m_2^2 + m_1^2 - m_1 m_2}
```

```math
= -\frac{200}{1600}(80) \pm \frac{1}{8}\sqrt{1600}
```

```math
= -10 \pm 5
```

```math
\therefore \lambda_1 = -5, \quad \lambda_2 = -15
```

**진동수:**

```math
\omega^2 = -\lambda \quad \Rightarrow \quad \omega_1 = \sqrt{5}, \quad f_1 = \frac{\sqrt{5}}{2\pi}
```

```math
\omega_2 = \sqrt{15}, \quad f_2 = \frac{\sqrt{15}}{2\pi}
```

여기서 $`\omega = 2\pi f`$이다.

**고유벡터:**

$`\lambda_1 = -5`$의 경우, $`A = \begin{bmatrix} -10 & 5 \\ 5 & -10 \end{bmatrix}`$에서 다음을 얻는다.

```math
\begin{bmatrix} -10 - (-5) & 5 \\ 5 & -10 - (-5) \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} -5 & 5 \\ 5 & -5 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}
```

```math
-x_1 + x_2 = 0 \quad \Rightarrow \quad \{x\}_1 = \begin{Bmatrix} 1 \\ 1 \end{Bmatrix}
```

두 질량이 **같은 방향**으로 움직인다(동위상 모드).

$`\lambda_2 = -15`$의 경우는 다음과 같다.

```math
\begin{bmatrix} -10 + 15 & 5 \\ 5 & -10 + 15 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{bmatrix} 5 & 5 \\ 5 & 5 \end{bmatrix}\begin{Bmatrix} x_1 \\ x_2 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \end{Bmatrix}
```

```math
x_1 + x_2 = 0 \quad \Rightarrow \quad \{x\}_2 = \begin{Bmatrix} 1 \\ -1 \end{Bmatrix}
```

두 질량이 **반대 방향**으로 움직인다(역위상 모드).

**일반해:**

```math
\mathbf{y} = C_1 \{x\}_1 \cos\sqrt{5}\,t + C_2 \{x\}_1 \sin\sqrt{5}\,t + C_3 \{x\}_2 \cos\sqrt{15}\,t + C_4 \{x\}_2 \sin\sqrt{15}\,t
```

계수 $`C_1, C_2, C_3, C_4`$는 초기 조건으로 정해진다.

각 성분으로 전개하면 다음과 같다.

```math
y_1 = C_1 \cdot 1 \cdot \cos\sqrt{5}\,t + C_2 \cdot 1 \cdot \sin\sqrt{5}\,t + C_3 \cdot 1 \cdot \cos\sqrt{15}\,t + C_4 \cdot 1 \cdot \sin\sqrt{15}\,t
```

```math
y_2 = C_1 \cdot 1 \cdot \cos\sqrt{5}\,t + C_2 \cdot 1 \cdot \sin\sqrt{5}\,t + C_3 \cdot (-1) \cdot \cos\sqrt{15}\,t + C_4 \cdot (-1) \cdot \sin\sqrt{15}\,t
```

---

<br>

## 10. 거듭제곱법(Power Method)

거듭제곱법은 절댓값이 **가장 큰** 지배 고유값과 이에 대응하는 고유벡터를 **반복적으로** 구한다.

### 10.1 세 질량, 네 스프링 시스템 (예제 13.5)

두 벽 사이의 세 질량, 네 스프링 시스템을 생각한다. $`m = m_1 = m_2 = m_3 = 1`$ kg, $`k = 20`$ N/m이다.

```text
벽 ─[k]─ m₁ ─[k]─ m₂ ─[k]─ m₃ ─[k]─ 벽
```

*그림 2. 두 벽 사이의 세 질량, 네 스프링 시스템*

각 질량에 뉴턴의 제2법칙을 적용하면 다음과 같다.

```math
m_1 \frac{d^2 x_1}{dt^2} = -kx_1 + k(x_2 - x_1)
```

```math
m_2 \frac{d^2 x_2}{dt^2} = -k(x_2 - x_1) + k(x_3 - x_2)
```

```math
m_3 \frac{d^2 x_3}{dt^2} = -k(x_3 - x_2) - kx_3
```

주기적인 해 $`\{x_i\} = \{X_i\}e^{i\omega t}`$를 가정하고 대입하면 다음과 같다.

```math
\begin{bmatrix} 40 - \omega^2 & -20 & 0 \\ -20 & 40 - \omega^2 & -20 \\ 0 & -20 & 40 - \omega^2 \end{bmatrix} \begin{Bmatrix} X_1 \\ X_2 \\ X_3 \end{Bmatrix} = \begin{Bmatrix} 0 \\ 0 \\ 0 \end{Bmatrix}
```

이것은 다음 $`A`$에 대한 고유값 문제 $`A\{X\} = \lambda\{X\}`$와 동등하다. 이 절에서는 부호를 바꾸어 $`\lambda = \omega^2`$로 둔다.

```math
A = \begin{bmatrix} 40 & -20 & 0 \\ -20 & 40 & -20 \\ 0 & -20 & 40 \end{bmatrix}, \quad \lambda = \omega^2
```

### 10.2 알고리즘

1. 초기 추정 벡터 $`\{x\}^{(0)}`$로 **시작**한다.
2. **곱하기:** $`\{y\}^{(k)} = [A]\{x\}^{(k-1)}`$를 계산한다.
3. **정규화:** 절댓값이 가장 큰 성분을 고유값 추정치 $`\lambda^{(k)}`$로 잡고, $`\{y\}^{(k)}`$를 이 성분으로 나누어 새 정규화 벡터 $`\{x\}^{(k)}`$를 얻는다.
4. 고유값 추정치가 수렴할 때까지 **반복**한다.

각 단계의 정규화 인자가 **지배 고유값의 추정치**이다.

### 10.3 반복 과정

**반복 1:** $`X_1 = X_2 = X_3 = 1`$을 초기 추정으로 쓴다.

```math
\begin{bmatrix} 40 & -20 & 0 \\ -20 & 40 & -20 \\ 0 & -20 & 40 \end{bmatrix} \begin{Bmatrix} 1 \\ 1 \\ 1 \end{Bmatrix} = \begin{Bmatrix} 20 \\ 0 \\ 20 \end{Bmatrix} = 20 \begin{Bmatrix} 1 \\ 0 \\ 1 \end{Bmatrix}
```

고유값 추정치는 $`\lambda \approx 20`$이다.

**반복 2:** $`\{1, 0, 1\}^T`$를 쓴다.

```math
\begin{bmatrix} 40 & -20 & 0 \\ -20 & 40 & -20 \\ 0 & -20 & 40 \end{bmatrix} \begin{Bmatrix} 1 \\ 0 \\ 1 \end{Bmatrix} = \begin{Bmatrix} 40 \\ -40 \\ 40 \end{Bmatrix} = 40 \begin{Bmatrix} 1 \\ -1 \\ 1 \end{Bmatrix}
```

고유값 추정치는 $`\lambda \approx 40`$이고, 근사 오차는 $`\varepsilon_a = \left|\frac{40 - 20}{40}\right| \times 100\% = 50\%`$이다.

**반복 3:** $`\{1, -1, 1\}^T`$를 쓴다.

```math
\begin{bmatrix} 40 & -20 & 0 \\ -20 & 40 & -20 \\ 0 & -20 & 40 \end{bmatrix} \begin{Bmatrix} 1 \\ -1 \\ 1 \end{Bmatrix} = \begin{Bmatrix} 60 \\ -80 \\ 60 \end{Bmatrix} = -80 \begin{Bmatrix} -3/4 \\ 1 \\ -3/4 \end{Bmatrix}
```

고유값 추정치는 $`\lambda \approx -80`$이고, 오차는 $`\varepsilon_a = \left|\frac{-80 - 40}{-80}\right| \times 100\% = 150\%`$이다.

**반복 4:** $`\{-3/4, 1, -3/4\}^T`$를 쓴다.

```math
A\begin{Bmatrix} -3/4 \\ 1 \\ -3/4 \end{Bmatrix} = \begin{Bmatrix} -50 \\ 70 \\ -50 \end{Bmatrix} = 70\begin{Bmatrix} -5/7 \\ 1 \\ -5/7 \end{Bmatrix}
```

고유값 추정치는 $`\lambda \approx 70`$이고, 오차는 $`\varepsilon_a = \left|\frac{70 - (-80)}{70}\right| \times 100\% \approx 214\%`$이다.

**반복 5:** $`\{-5/7, 1, -5/7\}^T`$를 쓴다.

```math
A\begin{Bmatrix} -5/7 \\ 1 \\ -5/7 \end{Bmatrix} = \begin{Bmatrix} -48.57 \\ 68.57 \\ -48.57 \end{Bmatrix} = 68.57\begin{Bmatrix} -0.71 \\ 1 \\ -0.71 \end{Bmatrix}
```

고유값 추정치는 $`\lambda \approx 68.57`$이고, 오차는 $`\varepsilon_a = \left|\frac{68.57 - 70}{68.57}\right| \times 100\% \approx 2.08\%`$이다.

### 10.4 수렴과 최소 고유값 찾기

고유값이 수렴하고 있다. 충분히 반복하면 다음을 얻는다.

```math
\lambda = 68.28427..., \quad \{x\} = \begin{Bmatrix} -0.70711 \\ 1 \\ -0.70711 \end{Bmatrix}
```

이 값은 $`40 + 20\sqrt{2}`$이며, 가장 높은 모드의 각진동수는 $`\omega = \sqrt{68.284} \approx 8.263`$ rad/s이다.

> **참고:** **거듭제곱법의 수렴 속도.** 거듭제곱법은 **절댓값이 가장 큰** 고유값으로 수렴한다. 수렴 속도는 비율 $`|\lambda_2/\lambda_1|`$에 따라 달라지며, $`\lambda_1`$은 지배 고유값, $`\lambda_2`$는 절댓값이 두 번째로 큰 고유값이다. 이 비율이 1에 가까울수록 수렴이 느려지며, 이 예에서는 $`40/68.284 \approx 0.586`$이다. 시작 벡터에 $`\lambda_1`$의 고유벡터 성분이 있어야 한다는 조건도 필요하다.

**최소 고유값을 찾으려면** $`A`$의 **역행렬**에 거듭제곱법을 적용할 수 있다.

```math
[A]^{-1}\{x\} = \frac{1}{\lambda}\{x\}
```

$`A^{-1}`$의 고유값 중 절댓값이 가장 큰 것이 $`1/\lambda_{\min}`$이므로, $`A^{-1}`$에 적용한 거듭제곱법은 $`1/\lambda_{\min}`$으로 수렴한다. 여기서 $`\lambda_{\min}`$은 절댓값이 가장 작은 고유값이다. 이것을 **역거듭제곱법(inverse power method)** 이라 한다.

---

<br>

## 11. Python 구현

### 11.1 NumPy를 이용한 고유값 계산

2.3절의 행렬에 `np.linalg.eig`를 적용하면 고유값 4와 2, 그리고 열마다 단위 길이로 정규화된 고유벡터를 얻는다.

```python
import numpy as np

# Example: 2x2 matrix from Section 2.3
A = np.array([[3, 1],
              [1, 3]])

eigenvalues, eigenvectors = np.linalg.eig(A)
print("Eigenvalues:", eigenvalues)       # [4. 2.]
print("Eigenvectors:\n", eigenvectors)   # columns are eigenvectors
```

### 11.2 거듭제곱법 구현

10.2절의 알고리즘을 그대로 구현하며, 절댓값이 가장 큰 성분으로 정규화한다. 예제 13.5의 행렬에 적용하면 11번 반복 후 지배 고유값 68.28428과 고유벡터 $`[-0.70710681, 1, -0.70710681]`$을 얻는다.

```python
import numpy as np

def power_method(A, x0, tol=1e-6, max_iter=100):
    """
    Power method to find the dominant eigenvalue and eigenvector.

    Parameters
    ----------
    A : ndarray
        Square matrix (n x n).
    x0 : ndarray
        Initial guess vector (n x 1).
    tol : float
        Convergence tolerance.
    max_iter : int
        Maximum number of iterations.

    Returns
    -------
    lam : float
        Dominant eigenvalue estimate.
    x : ndarray
        Corresponding eigenvector (normalized).
    """
    x = x0.copy().astype(float)
    lam_old = 0.0

    for k in range(max_iter):
        y = A @ x
        # Find the element with the largest absolute value
        idx = np.argmax(np.abs(y))
        lam = y[idx]
        x = y / lam

        # Check convergence
        if abs(lam - lam_old) / abs(lam) < tol:
            print(f"Converged in {k + 1} iterations")
            return lam, x
        lam_old = lam

    print("Did not converge within max_iter iterations")
    return lam, x


# Example 13.5: Three-mass four-spring system
A = np.array([[ 40, -20,   0],
              [-20,  40, -20],
              [  0, -20,  40]])

x0 = np.array([1, 1, 1])
lam, x = power_method(A, x0)
print(f"Dominant eigenvalue: {lam:.5f}")
print(f"Eigenvector: {x}")
```

### 11.3 질량-스프링 시스템 (예제 13.4)

9절의 계수 행렬로 고유값 $`-5`$, $`-15`$를 구하고, $`\omega = \sqrt{-\lambda}`$에서 각진동수 2.23607, 3.87298 rad/s와 진동수 0.35588, 0.61640 Hz를 얻는다.

```python
import numpy as np

# Parameters
m1, m2 = 40.0, 40.0  # kg
k = 200.0             # N/m

# Coefficient matrix
A = np.array([[-2*k/m1,    k/m1],
              [   k/m2, -2*k/m2]])

eigenvalues, eigenvectors = np.linalg.eig(A)
print("Eigenvalues (lambda):", eigenvalues)  # [-5. -15.]

# Natural frequencies
omega = np.sqrt(-eigenvalues)
freq = omega / (2 * np.pi)
print("Angular frequencies (omega):", omega)
print("Frequencies (Hz):", freq)
print("Eigenvectors:\n", eigenvectors)
```

### 11.4 최소 고유값을 위한 역거듭제곱법

역행렬을 만들지 않고 반복마다 `np.linalg.solve`로 $`A\mathbf{y} = \mathbf{x}`$를 푼다. 예제 13.5의 행렬에 적용하면 9번 반복 후 최소 고유값 11.71573($`= 40 - 20\sqrt{2}`$)과 고유벡터 $`[0.70710681, 1, 0.70710681]`$을 얻는다. `np.linalg.solve`는 호출할 때마다 행렬을 다시 분해하므로, 큰 행렬에서는 12장 실습처럼 LU 분해를 한 번만 해 두는 편이 효율적이다.

```python
import numpy as np

def inverse_power_method(A, x0, tol=1e-6, max_iter=100):
    """
    Inverse power method to find the smallest eigenvalue.
    """
    x = x0.copy().astype(float)
    lam_old = 0.0

    for k in range(max_iter):
        # Solve A @ y = x instead of computing A_inv @ x
        y = np.linalg.solve(A, x)
        idx = np.argmax(np.abs(y))
        lam_inv = y[idx]
        x = y / lam_inv

        lam = 1.0 / lam_inv
        if abs(lam - lam_old) / abs(lam) < tol:
            print(f"Converged in {k + 1} iterations")
            return lam, x
        lam_old = lam

    print("Did not converge within max_iter iterations")
    return lam, x


# Example: find smallest eigenvalue of the 3x3 spring system
A = np.array([[ 40, -20,   0],
              [-20,  40, -20],
              [  0, -20,  40]])

x0 = np.array([1, 1, 1])
lam_min, x_min = inverse_power_method(A, x0)
print(f"Smallest eigenvalue: {lam_min:.5f}")
print(f"Eigenvector: {x_min}")
```

---

<br>

## 요약

| 개념 | 핵심 요약 |
|---|---|
| 고유값 방정식 | $`[A]\{x\} = \lambda\{x\}`$를 만족하는 영이 아닌 $`\{x\}`$가 고유벡터, $`\lambda`$가 고유값이다. |
| 특성 방정식 | 고유값은 $`\det(A - \lambda I) = 0`$의 근이다. |
| 2x2 고유값 | $`2 \times 2`$ 행렬의 고유값은 $`\lambda = \frac{(a_{11}+a_{22}) \pm \sqrt{(a_{11}+a_{22})^2 - 4(a_{11}a_{22} - a_{12}a_{21})}}{2}`$이다. |
| 비자명해 조건 | 비자명해가 있으려면 $`[A - \lambda I]`$가 특이해야 한다. |
| 고유벡터 성질 | 고유벡터 $`\{x\}`$는 $`[A]\{x\}`$와 평행하다. |
| 스칼라 ODE $`dx/dt = ax`$ | 해는 $`x(t) = Ce^{at}`$이다. |
| 행렬 ODE $`d\mathbf{x}/dt = A\mathbf{x}`$ | 해는 $`\mathbf{x}(t) = e^{At}\mathbf{C}`$이며, $`A`$가 대칭이면 $`e^{At} = Ue^{Dt}U^T`$이다. |
| 안정성 | $`\text{Re}(\lambda) < 0`$이면 안정, $`\text{Re}(\lambda) > 0`$이면 불안정, $`\text{Re}(\lambda) = 0`$이면 진동한다. |
| 2차 ODE에서 고유값 | $`\ddot{y} = Ay`$에서 $`\lambda = -\omega^2`$이므로 $`\omega = \sqrt{-\lambda}`$이다. |
| 진동수 관계 | 주기는 $`2\pi/\omega`$이고 진동수는 $`f = \omega/(2\pi)`$이다. |
| 질량-스프링 시스템 | $`m\ddot{x} = -kx + k(x_{\text{neighbor}} - x)`$ 꼴의 운동 방정식은 고유값 문제가 된다. |
| 거듭제곱법 | 행렬을 반복해서 곱하고 정규화하여 절댓값이 가장 큰 지배 고유값을 계산한다. |
| 역거듭제곱법 | $`A^{-1}`$에 거듭제곱법을 적용하여 절댓값이 가장 작은 고유값을 찾는다. |
| 수렴 속도 | 거듭제곱법의 수렴은 $`\vert\lambda_2/\lambda_1\vert`$에 따라 달라지며, 이 비가 1에 가까울수록 느리다. |

---

<br>

## 점검 문제

1. **고유값과 고유벡터:** $`A = \begin{bmatrix} 4 & 1 \\ 2 & 3 \end{bmatrix}`$의 고유값과 고유벡터를 구하고, 대각합과 행렬식으로 확인하라.

   <details>
   <summary>정답</summary>

   특성 방정식은 $`(4-\lambda)(3-\lambda) - 2 = \lambda^2 - 7\lambda + 10 = (\lambda - 5)(\lambda - 2) = 0`$이므로 $`\lambda_1 = 5`$, $`\lambda_2 = 2`$이다. $`\lambda_1 = 5`$이면 $`-x_1 + x_2 = 0`$에서 $`\{x\}_1 = \{1, 1\}^T`$이고, $`\lambda_2 = 2`$이면 $`2x_1 + x_2 = 0`$에서 $`\{x\}_2 = \{1, -2\}^T`$이다. 고유값의 합 7은 대각합 $`4 + 3`$과 같고, 곱 10은 행렬식 $`4 \cdot 3 - 1 \cdot 2`$와 같다.

   </details>

2. **안정성 판정:** $`\frac{d\mathbf{x}}{dt} = \begin{bmatrix} -1 & 2 \\ -2 & -1 \end{bmatrix}\mathbf{x}`$의 고유값을 구하고, 해의 거동을 설명하라.

   <details>
   <summary>정답</summary>

   특성 방정식은 $`(-1-\lambda)^2 + 4 = 0`$이므로 $`\lambda = -1 \pm 2i`$이다. 실수부가 $`-1 < 0`$이므로 해는 안정하여 $`e^{-t}`$에 비례해 감쇠하고, 허수부 때문에 각진동수 2(주기 $`\pi`$)로 진동한다. 따라서 해는 진동하면서 0으로 수렴한다.

   </details>

3. **질량-스프링 시스템:** 8절의 두 질량, 세 스프링 시스템에서 $`m_1 = m_2 = 10`$ kg, $`k = 100`$ N/m일 때 고유값, 각진동수, 모드 형상을 구하라.

   <details>
   <summary>정답</summary>

   $`k/(m_1 m_2) = 100/100 = 1`$이므로 8.2절의 식에서 $`\lambda = -1 \cdot 20 \pm 1 \cdot \sqrt{100 + 100 - 100} = -20 \pm 10`$, 즉 $`\lambda_1 = -10`$, $`\lambda_2 = -30`$이다. 각진동수는 $`\omega_1 = \sqrt{10} \approx 3.162`$ rad/s, $`\omega_2 = \sqrt{30} \approx 5.477`$ rad/s이다. 계수 행렬이 $`\begin{bmatrix} -20 & 10 \\ 10 & -20 \end{bmatrix}`$이므로 예제 13.4와 같이 $`\lambda_1`$의 모드는 $`\{1, 1\}^T`$(동위상), $`\lambda_2`$의 모드는 $`\{1, -1\}^T`$(역위상)이다.

   </details>

4. **역거듭제곱법과 최저 모드:** 예제 13.5의 행렬에 역거듭제곱법을 적용하면 어떤 고유값으로 수렴하는가? 이 값에 해당하는 각진동수를 구하고, 거듭제곱법으로 구한 값과 비교하라.

   <details>
   <summary>정답</summary>

   고유값은 $`40 - 20\sqrt{2} \approx 11.716`$, $`40`$, $`40 + 20\sqrt{2} \approx 68.284`$이므로 역거듭제곱법은 절댓값이 가장 작은 11.716으로 수렴한다. 이 예제에서는 $`\lambda = \omega^2`$이므로 최저 모드의 각진동수는 $`\omega = \sqrt{11.716} \approx 3.423`$ rad/s이고, 거듭제곱법이 구한 최고 모드의 각진동수 $`\sqrt{68.284} \approx 8.263`$ rad/s보다 작다.

   </details>
