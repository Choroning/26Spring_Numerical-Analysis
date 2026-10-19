# 12장 — 선형계의 반복법

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **최종 수정일:** 2026-10-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 12

> **학습 목표:**
> 1. 가우스-자이델 방법과 야코비 방법으로 선형계를 반복해서 풀 수 있다
> 2. 대각 우위를 판정하고 수렴 조건으로서의 의미를 설명할 수 있다
> 3. 완화법이 반복법의 수렴을 어떻게 개선하는지 설명하고 적용할 수 있다
> 4. 전처리기와 반복행렬의 스펙트럼 반경으로 반복법의 수렴을 분석할 수 있다

---

<br>

## 목차

- [1. 반복법 소개](#1-반복법-소개)
- [2. 가우스-자이델 방법(Gauss-Seidel Method)](#2-가우스-자이델-방법gauss-seidel-method)
  - [2.1 문제 설정](#21-문제-설정)
  - [2.2 3x3 시스템에 대한 유도](#22-3x3-시스템에-대한-유도)
  - [2.3 위첨자를 사용한 반복 스킴](#23-위첨자를-사용한-반복-스킴)
  - [2.4 수렴 판정 기준](#24-수렴-판정-기준)
  - [2.5 예제 12.1: 가우스-자이델 반복](#25-예제-121-가우스-자이델-반복)
- [3. 수렴과 대각 우위(Diagonal Dominance)](#3-수렴과-대각-우위diagonal-dominance)
  - [3.1 대각 우위 조건](#31-대각-우위-조건)
  - [3.2 예제: 3x3 대각 우위 검사](#32-예제-3x3-대각-우위-검사)
  - [3.3 충분조건 vs. 필요조건](#33-충분조건-vs-필요조건)
- [4. 가우스-자이델 방법의 행렬 형식](#4-가우스-자이델-방법의-행렬-형식)
  - [4.1 행렬 분할: L, D, U 분해](#41-행렬-분할-l-d-u-분해)
  - [4.2 반복의 행렬 형식](#42-반복의-행렬-형식)
- [5. 야코비 반복법(Jacobi Iteration)](#5-야코비-반복법jacobi-iteration)
  - [5.1 가우스-자이델과의 비교](#51-가우스-자이델과의-비교)
  - [5.2 성분 형식](#52-성분-형식)
  - [5.3 행렬 형식](#53-행렬-형식)
- [6. 완화법(Relaxation Method)](#6-완화법relaxation-method)
  - [6.1 완화 공식](#61-완화-공식)
  - [6.2 완화의 종류](#62-완화의-종류)
  - [6.3 3x3 시스템에 적용한 완화법](#63-3x3-시스템에-적용한-완화법)
  - [6.4 완화법의 행렬 형식 (SOR)](#64-완화법의-행렬-형식-sor)
  - [6.5 예제 12.2: lambda = 1.2인 완화법](#65-예제-122-lambda--12인-완화법)
- [7. 전처리(Preconditioning)](#7-전처리preconditioning)
  - [7.1 일반 프레임워크](#71-일반-프레임워크)
  - [7.2 각 방법의 전처리기](#72-각-방법의-전처리기)
  - [7.3 오차 분석과 수렴](#73-오차-분석과-수렴)
  - [7.4 스펙트럼 반경(Spectral Radius)](#74-스펙트럼-반경spectral-radius)
  - [7.5 I - A의 고유값 관계](#75-i---a의-고유값-관계)
- [8. Python 구현](#8-python-구현)
  - [8.1 가우스-자이델 방법](#81-가우스-자이델-방법)
  - [8.2 야코비 방법](#82-야코비-방법)
  - [8.3 완화법 (SOR)](#83-완화법-sor)
  - [8.4 대각 우위 검사](#84-대각-우위-검사)
- [요약](#요약)
- [점검 문제](#점검-문제)

---

<br>

## 1. 반복법 소개

반복법은 다음과 같이 근사해의 수열을 생성하여 선형계 $`[A]\{x\} = \{b\}`$를 푼다.

```math
\mathbf{x}^{k+1} := \Psi(\mathbf{x}^k), \quad k \geq 0
```

초기 추측 $`\mathbf{x}^0`$에서 시작하여 참해에 수렴할 때까지 매핑 $`\Psi`$를 반복해서 적용한다. 직접법(가우스 소거법, LU 분해)은 반올림 오차를 무시하면 유한 번의 연산으로 정확한 답을 내지만, 반복법은 해에 점진적으로 접근하며 근사가 충분히 정확해지면 중단할 수 있다.

> **참고:** **직접법과 반복법의 비용.** 직접법은 $`O(n^3)`$의 비용이 들며 소규모에서 중규모 시스템에 잘 작동한다. 공학과 과학 계산에서 흔한 매우 크고 희소한 시스템에서는 반복법이 행렬-벡터 곱 수준의 연산만 필요하고 희소성을 활용할 수 있으므로 훨씬 효율적인 경우가 많다.

---

<br>

## 2. 가우스-자이델 방법(Gauss-Seidel Method)

### 2.1 문제 설정

다음과 같은 표준 형식의 방정식 $`n`$개가 주어졌다고 가정한다.

```math
[A]_{n \times n} \{x\}_{n \times 1} = \{b\}_{n \times 1}
```

여기서 $`[A]`$는 계수 행렬, $`\{x\}`$는 미지 벡터, $`\{b\}`$는 우변 벡터이다.

### 2.2 3x3 시스템에 대한 유도

간단하게 $`3 \times 3`$ 시스템을 생각한다.

```math
\begin{bmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \\ x_3 \end{Bmatrix} = \begin{Bmatrix} b_1 \\ b_2 \\ b_3 \end{Bmatrix}
```

각 행을 전개하면 다음과 같다.

```math
a_{11} x_1 = b_1 - a_{12} x_2 - a_{13} x_3
```

```math
a_{22} x_2 = b_2 - a_{21} x_1 - a_{23} x_3
```

```math
a_{33} x_3 = b_3 - a_{31} x_1 - a_{32} x_2
```

각 방정식을 대각 위치의 미지수에 대해 풀면 다음과 같다.

```math
x_1 = \frac{b_1 - a_{12} x_2 - a_{13} x_3}{a_{11}}
```

```math
x_2 = \frac{b_2 - a_{21} x_1 - a_{23} x_3}{a_{22}}
```

```math
x_3 = \frac{b_3 - a_{31} x_1 - a_{32} x_2}{a_{33}}
```

### 2.3 위첨자를 사용한 반복 스킴

가우스-자이델 방법은 반복 인덱스를 도입한다. 반복 $`j`$에서는 다음을 계산한다.

```math
x_1^{j} = \frac{b_1 - a_{12} x_2^{j-1} - a_{13} x_3^{j-1}}{a_{11}}
```

```math
x_2^{j} = \frac{b_2 - a_{21} x_1^{j} - a_{23} x_3^{j-1}}{a_{22}}
```

```math
x_3^{j} = \frac{b_3 - a_{31} x_1^{j} - a_{32} x_2^{j}}{a_{33}}
```

가우스-자이델의 핵심 특징은 **가장 최근에 계산된 값을 사용 가능한 즉시 사용**한다는 것이다. $`x_2^{j}`$를 계산할 때는 $`x_1^{j-1}`$이 아니라 현재 반복에서 방금 계산한 $`x_1^{j}`$를 사용한다. 마찬가지로 $`x_3^{j}`$를 계산할 때는 $`x_1^{j}`$와 $`x_2^{j}`$를 모두 사용한다.

> **참고:** **순차 갱신.** 새로 계산한 값을 같은 반복 안에서 바로 쓰는 이 성질이 가우스-자이델을 야코비 방법과 구분한다. 새로운 정보가 즉시 반영되므로 일반적으로 더 빠르게 수렴하지만, 계산이 본질적으로 순차적이어서 $`x_2^j`$는 $`x_1^j`$가 알려질 때까지 계산할 수 없다.

### 2.4 수렴 판정 기준

참해로 수렴할 때까지 절차를 반복한다. 각 변수의 **근사 상대 오차**는 다음과 같다.

```math
\varepsilon_{a,i} = \left| \frac{x_i^{j} - x_i^{j-1}}{x_i^{j}} \right| \times 100\%
```

**모든** 변수에 대해 $`\varepsilon_{a,i} \leq \varepsilon_s`$이면 반복을 멈춘다. 여기서 $`\varepsilon_s`$는 사용자가 지정한 허용오차이다.

### 2.5 예제 12.1: 가우스-자이델 반복

**주어진 시스템:**

```math
3x_1 - 0.1x_2 - 0.2x_3 = 7.85
```

```math
0.1x_1 + 7x_2 - 0.3x_3 = -19.3
```

```math
0.3x_1 - 0.2x_2 + 10x_3 = 71.4
```

**단계 1 (재배열):** 각 방정식을 대각 위치의 미지수에 대해 정리한다.

```math
x_1 = \frac{7.85 + 0.1x_2 + 0.2x_3}{3}
```

```math
x_2 = \frac{-19.3 - 0.1x_1 + 0.3x_3}{7}
```

```math
x_3 = \frac{71.4 - 0.3x_1 + 0.2x_2}{10}
```

**단계 2 (첫 번째 반복):** 초기 추측 $`\{x\}^0 = \{0, 0, 0\}`$에서 $`j = 1`$의 값을 계산한다.

```math
x_1^1 = \frac{7.85 + 0.1(0) + 0.2(0)}{3} = 2.616667
```

```math
x_2^1 = \frac{-19.3 - 0.1(2.616667) + 0.3(0)}{7} = -2.794524
```

```math
x_3^1 = \frac{71.4 - 0.3(2.616667) + 0.2(-2.794524)}{10} = 7.005610
```

$`x_2^1`$을 계산할 때 $`x_1^1`$을 바로 사용했고, $`x_3^1`$을 계산할 때 $`x_1^1`$과 $`x_2^1`$을 모두 사용했다는 점에 주목할 필요가 있다.

**단계 3 (두 번째 반복):** $`j = 2`$의 값을 계산한다.

```math
x_1^2 = \frac{7.85 + 0.1(-2.794524) + 0.2(7.005610)}{3} = 2.990557
```

```math
x_2^2 = \frac{-19.3 - 0.1(2.990557) + 0.3(7.005610)}{7} = -2.499625
```

```math
x_3^2 = \frac{71.4 - 0.3(2.990557) + 0.2(-2.499625)}{10} = 7.000291
```

**단계 4 (상대 오차 계산):** 반복 2의 상대 오차를 계산한다.

```math
\varepsilon_{a,1}^2 = \left| \frac{2.990557 - 2.616667}{2.990557} \right| \times 100\% = 12.5\%
```

$`x_2`$와 $`x_3`$의 상대 오차는 각각 $`11.8\%`$와 $`0.0760\%`$이다.

**모든** 상대 오차가 사용자 지정 허용오차 $`\varepsilon_s`$보다 작아질 때까지 이 절차를 반복한다. 세 번째 반복의 값은 $`x_1 = 3.000032`$, $`x_2 = -2.499988`$, $`x_3 = 6.999999`$이며, 참해는 $`x_1 = 3`$, $`x_2 = -2.5`$, $`x_3 = 7`$이다.

---

<br>

## 3. 수렴과 대각 우위(Diagonal Dominance)

### 3.1 대각 우위 조건

가우스-자이델 방법은 일부 시스템에서 **발산할 수 있다**. 그러나 시스템이 **대각 우위(diagonally dominant)** 이면 **반드시 수렴한다**.

행렬이 대각 우위라는 것은 모든 행 $`i`$에 대해 다음이 성립한다는 뜻이다.

```math
|a_{ii}| > \sum_{\substack{j=1 \\ j \neq i}}^{n} |a_{ij}|
```

즉, 각 대각 원소의 절댓값이 해당 행의 다른 모든 원소의 절댓값 합보다 **엄격하게 크다**. 이 조건을 엄밀한 대각 우위(strict diagonal dominance)라고도 하며, 이 장에서 대각 우위는 이 엄밀한 조건을 뜻한다.

### 3.2 예제: 3x3 대각 우위 검사

$`3 \times 3`$ 시스템을 생각한다.

```math
\begin{bmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \\ x_3 \end{Bmatrix} = \begin{Bmatrix} b_1 \\ b_2 \\ b_3 \end{Bmatrix}
```

검증할 세 가지 조건은 다음과 같다.

```math
|a_{11}| > |a_{12}| + |a_{13}|
```

```math
|a_{22}| > |a_{21}| + |a_{23}|
```

```math
|a_{33}| > |a_{31}| + |a_{32}|
```

예제 12.1에서는 $`|3| > |-0.1| + |-0.2|`$, $`|7| > |0.1| + |-0.3|`$, $`|10| > |0.3| + |-0.2|`$이다. 모든 조건이 성립하므로 수렴이 보장된다.

### 3.3 충분조건 vs. 필요조건

대각 우위는 수렴의 **충분조건**이지만 **필요조건은 아니다**.

- 행렬이 대각 우위이면 수렴이 보장된다.
- 행렬이 대각 우위가 **아니어도** 방법이 **수렴할 수 있다** (다만 보장되지는 않는다).
- 가우스-자이델 방법은 **대칭 양정치(symmetric positive definite)** 행렬에서도 수렴한다.

> **참고:** **행 교환으로 대각 우위 만들기.** 원래 형태에서 대각 우위가 아닌 시스템도 행을 교환(피벗)하면 대각 우위가 되는 경우가 있다. 6.5절의 예제 12.2가 이 전략을 사용한다.

---

<br>

## 4. 가우스-자이델 방법의 행렬 형식

### 4.1 행렬 분할: L, D, U 분해

계수 행렬 $`[A]`$는 다음과 같이 분해할 수 있다.

```math
[A] = \tilde{L} + \tilde{D} + \tilde{U}
```

여기서 각 행렬은 다음과 같다.

- $`\tilde{L}`$은 $`[A]`$의 **엄밀 하삼각** 부분이다 (대각선 아래 원소만 남기고 대각선과 그 위는 0).
- $`\tilde{D}`$는 $`[A]`$의 **대각** 행렬이다 (대각 원소만).
- $`\tilde{U}`$는 $`[A]`$의 **엄밀 상삼각** 부분이다 (대각선 위 원소만 남기고 대각선과 그 아래는 0).

$`3 \times 3`$ 경우는 다음과 같다.

```math
\tilde{L} = \begin{bmatrix} 0 & 0 & 0 \\ a_{21} & 0 & 0 \\ a_{31} & a_{32} & 0 \end{bmatrix}, \quad \tilde{D} = \begin{bmatrix} a_{11} & 0 & 0 \\ 0 & a_{22} & 0 \\ 0 & 0 & a_{33} \end{bmatrix}, \quad \tilde{U} = \begin{bmatrix} 0 & a_{12} & a_{13} \\ 0 & 0 & a_{23} \\ 0 & 0 & 0 \end{bmatrix}
```

> **참고:** **LU 분해와의 차이.** 이 $`\tilde{L} + \tilde{D} + \tilde{U}`$ 분할은 10장의 LU 분해와 다르다. 여기서 $`\tilde{L}`$, $`\tilde{D}`$, $`\tilde{U}`$는 원래 행렬 $`A`$의 부분을 그대로 떼어 낸 것이며 소거를 수행하지 않는다. LU 인수와 구별하기 위해 틸데 표기를 사용한다.

### 4.2 반복의 행렬 형식

가우스-자이델 방정식을 행렬 형식으로 재배열하면 다음과 같다.

```math
[\tilde{L} + \tilde{D}]\{x\}^{j} = \{b\} - [\tilde{U}]\{x\}^{j-1}
```

$`\tilde{L}_* = \tilde{L} + \tilde{D}`$($`A`$의 대각선을 포함한 하삼각 부분)로 놓으면 행렬-벡터 형식은 다음과 같다.

```math
\{x\}^{j} = (\tilde{L}_*)^{-1}\left(\{b\} - [\tilde{U}]\{x\}^{j-1}\right)
```

좌변은 **하삼각** 시스템이므로 **전방 대입**으로 효율적으로 풀 수 있으며, 역행렬을 명시적으로 구할 필요가 없다.

---

<br>

## 5. 야코비 반복법(Jacobi Iteration)

### 5.1 가우스-자이델과의 비교

**야코비 반복법**은 가우스-자이델 방법과 비슷하지만, $`j`$번째 반복의 모든 변수를 **오직** $`(j-1)`$번째 반복 값만으로 갱신한다. 같은 반복 안에서 새로 계산한 값은 사용하지 않는다.

| 특성 | 가우스-자이델 | 야코비 |
|---|---|---|
| 갱신 전략 | 최신 값을 즉시 사용 | 이전 반복 값만 사용 |
| 병렬성 | 순차적 (갱신 순서에 의존) | **자연적으로 병렬** (모든 갱신이 독립) |
| 수렴 속도 | 일반적으로 빠름 | 일반적으로 느림 |
| 저장 공간 | 해 벡터 1개 | 2개 필요 (이전 값과 새 값) |

*표 1. 가우스-자이델 방법과 야코비 방법의 비교*

야코비 방법도 $`[A]`$가 대각 우위이면 참해로 수렴한다.

### 5.2 성분 형식

$`3 \times 3`$ 시스템에서 야코비 반복은 다음을 계산한다.

```math
x_1^{j} = \frac{b_1 - a_{12} x_2^{j-1} - a_{13} x_3^{j-1}}{a_{11}}
```

```math
x_2^{j} = \frac{b_2 - a_{21} x_1^{j-1} - a_{23} x_3^{j-1}}{a_{22}}
```

```math
x_3^{j} = \frac{b_3 - a_{31} x_1^{j-1} - a_{32} x_2^{j-1}}{a_{33}}
```

우변의 **모든** 값이 위첨자 $`j-1`$을 사용한다는 점에 주목할 필요가 있다. 따라서 각 $`x_i^j`$를 **독립적으로** 계산할 수 있으며, 이것이 야코비 방법이 병렬 계산에 적합한 이유이다.

### 5.3 행렬 형식

```math
[\tilde{D}]\{x\}^{j} = \{b\} - [\tilde{L} + \tilde{U}]\{x\}^{j-1}
```

```math
\{x\}^{j} = [\tilde{D}]^{-1}\left(\{b\} - [\tilde{L} + \tilde{U}]\{x\}^{j-1}\right)
```

$`\tilde{D}`$가 대각행렬이므로 역행렬은 대각 원소의 역수 $`[\tilde{D}]^{-1}_{ii} = 1/a_{ii}`$로 바로 얻는다.

---

<br>

## 6. 완화법(Relaxation Method)

### 6.1 완화 공식

완화법은 **수렴을 향상**시키기 위해 가우스-자이델 방법을 약간 수정한 것이다. 이전 값과 새로 계산한 가우스-자이델 값을 가중 인수로 결합하는 것이 핵심 아이디어이다.

```math
x_i^{\text{new}} = \lambda \, x_i^{\text{new(GS)}} + (1 - \lambda) \, x_i^{\text{old}}
```

여기서 $`\lambda`$는 $`0`$과 $`2`$ 사이의 값을 갖는 가중 인수이고, $`x_i^{\text{new(GS)}}`$는 표준 가우스-자이델 방법이 계산할 값이다.

### 6.2 완화의 종류

| $`\lambda`$ 범위 | 명칭 | 효과 |
|---|---|---|
| $`\lambda = 1`$ | 수정 없음 | 표준 가우스-자이델 방법 |
| $`0 < \lambda < 1`$ | **부족 완화(Under-relaxation)** | 진동을 감쇠; 비수렴 시스템에 도움이 될 수 있음 |
| $`1 < \lambda < 2`$ | **과완화(Over-relaxation)** (SOR) | 가우스-자이델이 느릴 때 수렴을 가속 |

*표 2. 완화 계수의 범위에 따른 완화의 종류*

> **참고:** **최적 완화 계수.** $`\lambda`$의 최적값은 문제에 따라 다르며 일반적으로 미리 알 수 없다. SOR(연속 과완화) 방법에는 특수한 행렬 구조(예: 삼대각 시스템)에 대한 이론적 결과가 있지만, 실제로는 $`\lambda`$를 실험적으로 정하는 경우가 많다. $`\lambda`$가 $`0 < \lambda < 2`$를 벗어나면 SOR 반복은 수렴하지 않는다.

### 6.3 3x3 시스템에 적용한 완화법

$`3 \times 3`$ 시스템의 완화 갱신은 다음과 같다.

```math
x_1^{\text{new}} = \frac{b_1 - a_{12} x_2^{\text{old}} - a_{13} x_3^{\text{old}}}{a_{11}} \cdot \lambda + (1 - \lambda) \, x_1^{\text{old}}
```

```math
x_2^{\text{new}} = \frac{b_2 - a_{21} x_1^{\text{new}} - a_{23} x_3^{\text{old}}}{a_{22}} \cdot \lambda + (1 - \lambda) \, x_2^{\text{old}}
```

```math
x_3^{\text{new}} = \frac{b_3 - a_{31} x_1^{\text{new}} - a_{32} x_2^{\text{new}}}{a_{33}} \cdot \lambda + (1 - \lambda) \, x_3^{\text{old}}
```

가우스-자이델 패턴이 유지된다는 점에 주목할 필요가 있다. $`x_1^{\text{new}}`$은 $`x_2^{\text{new}}`$ 계산에 쓰이고, 둘 다 $`x_3^{\text{new}}`$ 계산에 쓰인다.

### 6.4 완화법의 행렬 형식 (SOR)

성분 형식을 전개하고 항을 정리하면 완화법의 행렬 형식은 다음과 같다.

```math
[\tilde{D} + \lambda \tilde{L}]\{x\}^{\text{new}} = \lambda\{b\} - \left(\lambda[\tilde{U}] + (\lambda - 1)[\tilde{D}]\right)\{x\}^{\text{old}}
```

**유도 개요:**

성분 형식에서 출발하여 양변에 $`a_{ii}`$를 곱하면 다음과 같다.

```math
a_{ii} x_i^{\text{new}} = \left(b_i - \sum_{k<i} a_{ik} x_k^{\text{new}} - \sum_{k>i} a_{ik} x_k^{\text{old}}\right) \lambda + (1-\lambda) \, a_{ii} \, x_i^{\text{old}}
```

전개하고 재배열하면 다음과 같다.

```math
a_{ii} x_i^{\text{new}} = \lambda b_i - \lambda \sum_{k>i} a_{ik} x_k^{\text{old}} + (1-\lambda) a_{ii} x_i^{\text{old}} - \lambda \sum_{k<i} a_{ik} x_k^{\text{new}}
```

위첨자가 new인 항을 좌변으로, old인 항을 우변으로 옮기고 행렬 표기로 모으면 최종 행렬 형식이 된다.

### 6.5 예제 12.2: lambda = 1.2인 완화법

**주어진 조건:** $`\lambda = 1.2`$, $`\varepsilon_s = 10\%`$

**원래 시스템:**

```math
-3x_1 + 12x_2 = 9
```

```math
10x_1 - 2x_2 = 8
```

**단계 1 (행 교환):** 대각 우위를 얻기 위해 두 행을 바꾼다.

```math
10x_1 - 2x_2 = 8
```

```math
-3x_1 + 12x_2 = 9
```

교환 후 $`|10| > |-2|`$이고 $`|12| > |-3|`$이므로 시스템이 대각 우위가 된다.

각 식을 대각 위치의 미지수에 대해 정리하면 다음과 같다.

```math
x_1^{j} = \frac{8 + 2x_2^{j-1}}{10}, \quad x_2^{j} = \frac{9 + 3x_1^{j}}{12}
```

**반복 1** (초기 추측 $`x_1^0 = 0`$, $`x_2^0 = 0`$):

가우스-자이델 단계:

```math
x_1^1 = \frac{8 + 2(0)}{10} = 0.8
```

완화 적용:

```math
x_1^{\text{new}} = \lambda \cdot x_1^1 + (1-\lambda) \cdot x_1^0 = 1.2(0.8) + (-0.2)(0) = 0.96
```

가우스-자이델 단계 (완화된 $`x_1^{\text{new}} = 0.96`$ 사용):

```math
x_2^1 = \frac{9 + 3(0.96)}{12} = 0.99
```

완화 적용:

```math
x_2^{\text{new}} = 1.2(0.99) + (-0.2)(0) = 1.188
```

**반복 2** (반복 1에서 완화한 값 $`x_1 = 0.96`$, $`x_2 = 1.188`$ 사용):

가우스-자이델 단계:

```math
x_1^2 = \frac{8 + 2(1.188)}{10} = 1.0376
```

완화 적용:

```math
x_1^{\text{new}} = 1.2(1.0376) + (-0.2)(0.96) = 1.05312
```

가우스-자이델 단계:

```math
x_2^2 = \frac{9 + 3(1.05312)}{12} = 1.01328
```

완화 적용:

```math
x_2^{\text{new}} = 1.2(1.01328) + (-0.2)(1.188) = 0.978336
```

상대 오차:

```math
\varepsilon_{a,1} = \left|\frac{1.05312 - 0.96}{1.05312}\right| \times 100\% = 8.84\%
```

```math
\varepsilon_{a,2} = \left|\frac{0.978336 - 1.188}{0.978336}\right| \times 100\% = 21.43\%
```

$`\varepsilon_{a,2} > \varepsilon_s = 10\%`$이므로 반복을 계속한다.

**반복 3** (완화 후 값):

```math
x_1^3 = 0.984177, \quad \varepsilon_{a,1} = 7.01\%
```

```math
x_2^3 = 0.999586, \quad \varepsilon_{a,2} = 2.13\%
```

두 오차가 모두 $`10\%`$ 미만이므로 반복을 멈춘다. 참해는 $`x_1 = 1`$, $`x_2 = 1`$이다.

---

<br>

## 7. 전처리(Preconditioning)

### 7.1 일반 프레임워크

가우스-자이델, 야코비, 완화법의 세 반복법은 모두 선형계 $`A\mathbf{x} = \mathbf{b}`$를 풀기 위한 **전통적 반복법**이며, 이 선형계를 식 (1)이라 한다.

$`A`$의 근사인 **전처리기(preconditioner)** $`P`$를 도입한다. $`A\mathbf{x} = \mathbf{b}`$의 양변에 $`P\mathbf{x}`$를 더하고 정리하면 다음과 같다.

```math
P\mathbf{x} + A\mathbf{x} = P\mathbf{x} + \mathbf{b}
```

```math
P\mathbf{x} = (P - A)\mathbf{x} + \mathbf{b} \quad \cdots (2)
```

해 $`\mathbf{x}^{k+1}`$을 반복해서 갱신하면 다음과 같다.

```math
P\mathbf{x}^{k+1} = (P - A)\mathbf{x}^k + \mathbf{b} \quad \cdots (3)
```

목표는 식 (3)이 식 (1)보다 **풀기 쉬운** 것이다. $`P`$를 적은 비용으로 풀 수 있는 행렬(예: 삼각행렬 또는 대각행렬)로 고르기 때문이다.

### 7.2 각 방법의 전처리기

각 반복법은 전처리기 $`P`$를 특정하게 고른 경우에 해당한다.

| 방법 | 전처리기 $`P`$ | $`P - A`$ |
|---|---|---|
| **가우스-자이델** | $`P = \tilde{L} + \tilde{D}`$ | $`P - A = -\tilde{U}`$ |
| **야코비** | $`P = \tilde{D}`$ | $`P - A = -\tilde{L} - \tilde{U}`$ |
| **완화법 (SOR)** | $`P = \frac{1}{\lambda}(\tilde{D} + \lambda\tilde{L})`$ | $`P - A = \frac{1-\lambda}{\lambda}\tilde{D} - \tilde{U}`$ |

*표 3. 각 반복법에 해당하는 전처리기와 P - A*

완화법의 $`P`$에 붙은 $`1/\lambda`$는 반복식 $`P\mathbf{x}^{k+1} = (P - A)\mathbf{x}^k + \mathbf{b}`$의 양변에 $`\lambda`$를 곱했을 때 6.4절의 행렬 형식이 그대로 나오게 한다.

결과 반복 방정식은 다음과 같다.

- **가우스-자이델:** $`[\tilde{L} + \tilde{D}]\{x\}^j = \{b\} - [\tilde{U}]\{x\}^{j-1}`$
- **야코비:** $`[\tilde{D}]\{x\}^j = \{b\} - [\tilde{L} + \tilde{U}]\{x\}^{j-1}`$
- **완화법:** $`[\tilde{D} + \lambda\tilde{L}]\{x\}^{\text{new}} = \lambda\{b\} - (\lambda[\tilde{U}] + (\lambda-1)[\tilde{D}])\{x\}^{\text{old}}`$

### 7.3 오차 분석과 수렴

반복 $`k`$에서의 오차를 다음과 같이 정의한다.

```math
\mathbf{e}^k = \mathbf{x} - \mathbf{x}^k
```

여기서 $`\mathbf{x}`$는 참해이다. 식 (2)에서 식 (3)을 빼면 다음과 같다.

```math
P(\mathbf{x} - \mathbf{x}^{k+1}) = (P - A)(\mathbf{x} - \mathbf{x}^k)
```

```math
P\mathbf{e}^{k+1} = (P - A)\mathbf{e}^k
```

```math
\mathbf{e}^{k+1} = (I - P^{-1}A)\mathbf{e}^k
```

**반복행렬(iteration matrix)** 을 다음과 같이 정의한다.

```math
M = I - P^{-1}A
```

그러면 각 단계는 오차 벡터에 $`M`$을 곱한다.

```math
\mathbf{e}^{k+1} = M\mathbf{e}^k
```

수렴은 **$`M`$의 고유값**이 결정한다. 임의의 초기 추측에서 수렴하려면 $`M`$의 모든 고유값이 $`|\lambda(M)| < 1`$을 만족해야 한다.

### 7.4 스펙트럼 반경(Spectral Radius)

$`M`$의 고유값 절댓값 중 가장 큰 값을 **스펙트럼 반경**이라 한다.

```math
\rho(M) = \max(|\lambda(M)|)
```

**수렴을 위한 필요충분조건:**

```math
\rho(M) < 1
```

오차 벡터가 $`M`$의 고유벡터이면 다음 단계의 오차는 다음과 같다.

```math
\mathbf{e}^{k+1} = M\mathbf{e}^k = \lambda \mathbf{e}^k
```

```math
\Rightarrow \mathbf{e}^k = \lambda^k \mathbf{e}^0
```

따라서 오차는 $`\lambda^k`$에 비례하여 기하급수적으로 감소한다. $`\rho(M)`$이 작을수록 수렴이 빠르다.

> **참고:** **대각 우위와 스펙트럼 반경.** 대각 우위가 수렴을 보장하는 근본적인 이유가 여기에 있다. 대각 우위 행렬에서는 가우스-자이델 방법과 야코비 방법의 반복행렬 $`M`$의 스펙트럼 반경이 1보다 작다. SOR에서 $`\lambda`$를 고르는 목적은 $`\rho(M)`$을 최소화하는 것이다.

### 7.5 I - A의 고유값 관계

전처리기를 사용하지 않으면($`P = I`$) $`M = I - A`$이다.

**정리:** $`[A]`$의 고유값이 $`\lambda_1, \lambda_2, \ldots, \lambda_n`$이면, $`[I - A]`$의 고유값은 $`1 - \lambda_1, 1 - \lambda_2, \ldots, 1 - \lambda_n`$이다.

**증명:**

$`A\mathbf{v} = \lambda\mathbf{v}`$($`\mathbf{v}`$는 고유벡터)이면 다음이 성립한다.

```math
(I - A)\mathbf{v} = \mathbf{v} - A\mathbf{v} = \mathbf{v} - \lambda\mathbf{v} = (1 - \lambda)\mathbf{v}
```

따라서 $`(1 - \lambda)`$는 같은 고유벡터 $`\mathbf{v}`$를 갖는 $`[I - A]`$의 고유값이다.

**결론:** 전처리를 하지 않으면, 반복이 수렴하기 위해 $`A`$의 모든 고유값이 복소평면에서 1을 중심으로 하는 단위원 안에 있어야 한다.

---

<br>

## 8. Python 구현

### 8.1 가우스-자이델 방법

각 반복에서 `x[i]`를 제자리에서 갱신하므로 같은 반복에서 새로 계산한 값을 바로 쓴다. 허용오차 `tol`은 퍼센트 단위이며, 예제 12.1에 `tol=0.1`을 적용하면 4번 반복 후 해 $`[3.00000035, -2.50000004, 6.99999999]`$를 얻는다.

```python
import numpy as np

def gauss_seidel(A, b, x0=None, tol=1e-6, max_iter=100):
    """
    Gauss-Seidel iterative method for solving Ax = b.

    Parameters
    ----------
    A : array_like, shape (n, n)
        Coefficient matrix.
    b : array_like, shape (n,)
        Right-hand side vector.
    x0 : array_like, shape (n,), optional
        Initial guess (default: zeros).
    tol : float
        Stopping tolerance for approximate relative error (%).
    max_iter : int
        Maximum number of iterations.

    Returns
    -------
    x : ndarray
        Solution vector.
    iterations : int
        Number of iterations performed.
    """
    n = len(b)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)

    for iteration in range(1, max_iter + 1):
        x_old = x.copy()

        for i in range(n):
            sigma = sum(A[i][j] * x[j] for j in range(n) if j != i)
            x[i] = (b[i] - sigma) / A[i][i]

        # Check convergence (max approximate relative error)
        errors = [abs((x[i] - x_old[i]) / x[i]) * 100
                  for i in range(n) if x[i] != 0]
        if max(errors) < tol:
            return x, iteration

    return x, max_iter


# Example 12.1
A = np.array([[3, -0.1, -0.2],
              [0.1, 7, -0.3],
              [0.3, -0.2, 10]], dtype=float)
b = np.array([7.85, -19.3, 71.4], dtype=float)

x, iters = gauss_seidel(A, b, tol=0.1)
print(f"Solution: {x}")
print(f"Iterations: {iters}")
```

### 8.2 야코비 방법

새 값을 별도 배열 `x_new`에 모았다가 반복이 끝난 뒤 한꺼번에 바꾸므로, 모든 갱신이 이전 반복의 값만 사용한다.

```python
def jacobi(A, b, x0=None, tol=1e-6, max_iter=100):
    """
    Jacobi iterative method for solving Ax = b.

    Parameters
    ----------
    A : array_like, shape (n, n)
        Coefficient matrix.
    b : array_like, shape (n,)
        Right-hand side vector.
    x0 : array_like, shape (n,), optional
        Initial guess (default: zeros).
    tol : float
        Stopping tolerance for approximate relative error (%).
    max_iter : int
        Maximum number of iterations.

    Returns
    -------
    x : ndarray
        Solution vector.
    iterations : int
        Number of iterations performed.
    """
    n = len(b)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)

    for iteration in range(1, max_iter + 1):
        x_new = np.zeros(n)

        for i in range(n):
            sigma = sum(A[i][j] * x[j] for j in range(n) if j != i)
            x_new[i] = (b[i] - sigma) / A[i][i]

        # Check convergence
        errors = [abs((x_new[i] - x[i]) / x_new[i]) * 100
                  for i in range(n) if x_new[i] != 0]
        if max(errors) < tol:
            return x_new, iteration

        x = x_new

    return x, max_iter
```

### 8.3 완화법 (SOR)

각 원소의 가우스-자이델 값 `x_gs`를 구한 뒤 이전 값 `x_old[i]`와 가중 평균한다. 예제 12.2에 `lam=1.2`, `tol=10`을 적용하면 3번 반복 후 $`[0.98417664, 0.99958579]`$를 얻으며, 6.5절의 손 계산과 일치한다.

```python
def sor(A, b, lam=1.2, x0=None, tol=1e-6, max_iter=100):
    """
    Successive Over-Relaxation (SOR) method for solving Ax = b.

    Parameters
    ----------
    A : array_like, shape (n, n)
        Coefficient matrix.
    b : array_like, shape (n,)
        Right-hand side vector.
    lam : float
        Relaxation factor (0 < lam < 2).
        lam = 1 => Gauss-Seidel,
        lam < 1 => under-relaxation,
        lam > 1 => over-relaxation.
    x0 : array_like, shape (n,), optional
        Initial guess (default: zeros).
    tol : float
        Stopping tolerance for approximate relative error (%).
    max_iter : int
        Maximum number of iterations.

    Returns
    -------
    x : ndarray
        Solution vector.
    iterations : int
        Number of iterations performed.
    """
    n = len(b)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)

    for iteration in range(1, max_iter + 1):
        x_old = x.copy()

        for i in range(n):
            sigma = sum(A[i][j] * x[j] for j in range(n) if j != i)
            x_gs = (b[i] - sigma) / A[i][i]          # Gauss-Seidel value
            x[i] = lam * x_gs + (1 - lam) * x_old[i]  # Relaxation

        # Check convergence
        errors = [abs((x[i] - x_old[i]) / x[i]) * 100
                  for i in range(n) if x[i] != 0]
        if max(errors) < tol:
            return x, iteration

    return x, max_iter


# Example 12.2
A2 = np.array([[10, -2],
               [-3, 12]], dtype=float)
b2 = np.array([8, 9], dtype=float)

x2, iters2 = sor(A2, b2, lam=1.2, tol=10)
print(f"Solution: {x2}")
print(f"Iterations: {iters2}")
```

### 8.4 대각 우위 검사

각 행에서 대각 원소의 절댓값이 나머지 원소의 절댓값 합보다 엄격하게 큰지 검사하며, 예제 12.1의 행렬에 대해 `True`를 출력한다.

```python
def is_diagonally_dominant(A):
    """
    Check if a matrix is strictly diagonally dominant.

    Parameters
    ----------
    A : array_like, shape (n, n)
        Matrix to check.

    Returns
    -------
    bool
        True if the matrix is strictly diagonally dominant.
    """
    n = len(A)
    for i in range(n):
        diag = abs(A[i][i])
        off_diag_sum = sum(abs(A[i][j]) for j in range(n) if j != i)
        if diag <= off_diag_sum:
            return False
    return True


# Check Example 12.1
print(is_diagonally_dominant(A))   # True
```

---

<br>

## 요약

| 개념 | 핵심 요약 |
|---|---|
| 일반 반복 | 초기 추측에서 출발해 $`\mathbf{x}^{k+1} := \Psi(\mathbf{x}^k)`$ ($`k \geq 0`$)를 반복하며 근사해의 수열을 만든다. |
| 가우스-자이델 | $`x_i^j = \frac{1}{a_{ii}}\left(b_i - \sum_{k<i} a_{ik} x_k^j - \sum_{k>i} a_{ik} x_k^{j-1}\right)`$로 갱신하며, 같은 반복에서 새로 계산한 값을 바로 사용한다. |
| 가우스-자이델 (행렬) | 행렬 형식은 $`[\tilde{L} + \tilde{D}]\{x\}^j = \{b\} - [\tilde{U}]\{x\}^{j-1}`$이며 전방 대입으로 푼다. |
| 야코비 | $`x_i^j = \frac{1}{a_{ii}}\left(b_i - \sum_{k \neq i} a_{ik} x_k^{j-1}\right)`$로 갱신하며, 모든 값을 이전 반복에서 가져오므로 병렬 계산이 가능하다. |
| 야코비 (행렬) | 행렬 형식은 $`[\tilde{D}]\{x\}^j = \{b\} - [\tilde{L} + \tilde{U}]\{x\}^{j-1}`$이다. |
| 완화법 (SOR) | $`x_i^{\text{new}} = \lambda \, x_i^{\text{GS}} + (1-\lambda) \, x_i^{\text{old}}`$로 가우스-자이델 값과 이전 값을 가중 평균하며, $`0 < \lambda < 1`$은 부족 완화, $`1 < \lambda < 2`$는 과완화이다. |
| 완화법 (행렬) | 행렬 형식은 $`[\tilde{D} + \lambda\tilde{L}]\{x\}^{\text{new}} = \lambda\{b\} - (\lambda[\tilde{U}] + (\lambda-1)[\tilde{D}])\{x\}^{\text{old}}`$이다. |
| 대각 우위 | 모든 행에서 $`\vert a_{ii}\vert > \sum_{j \neq i} \vert a_{ij}\vert`$이면 가우스-자이델 방법과 야코비 방법이 수렴하며, 이 조건은 충분조건이지만 필요조건은 아니다. |
| 수렴 판정 기준 | 모든 $`i`$에 대해 $`\varepsilon_{a,i} = \left\vert\frac{x_i^j - x_i^{j-1}}{x_i^j}\right\vert \times 100\% \leq \varepsilon_s`$이면 반복을 멈춘다. |
| 전처리 | 반복식 $`P\mathbf{x}^{k+1} = (P-A)\mathbf{x}^k + \mathbf{b}`$에서 오차는 단계마다 반복행렬 $`M = I - P^{-1}A`$가 곱해진다. |
| 스펙트럼 반경 | $`\rho(M) = \max\vert\lambda(M)\vert < 1`$은 수렴의 필요충분조건이며, $`\rho(M)`$이 작을수록 빨리 수렴한다. |
| I - A의 고유값 | $`A`$의 고유값이 $`\lambda_i`$이면 $`I-A`$의 고유값은 $`1-\lambda_i`$이다. |

---

<br>

## 점검 문제

1. **대각 우위와 행 교환:** 계수 행렬이 $`\begin{bmatrix} 2 & 5 \\ 6 & 1 \end{bmatrix}`$인 연립방정식은 대각 우위인가? 아니라면 가우스-자이델 방법의 수렴을 보장하려면 어떻게 해야 하는가?

   <details>
   <summary>정답</summary>

   대각 우위가 아니다. 첫째 행은 $`|2| < |5|`$, 둘째 행은 $`|1| < |6|`$이다. 두 행을 바꾸면 계수 행렬이 $`\begin{bmatrix} 6 & 1 \\ 2 & 5 \end{bmatrix}`$가 되고 $`|6| > |1|`$, $`|5| > |2|`$이므로 대각 우위가 되어 수렴이 보장된다.

   </details>

2. **야코비 반복:** 예제 12.1의 시스템에 초기 추측 $`\{x\}^0 = \{0, 0, 0\}`$으로 야코비 방법을 한 번 적용하여 $`\{x\}^1`$을 구하고, 가우스-자이델 방법의 $`x_2^1`$과 값이 다른 이유를 설명하라.

   <details>
   <summary>정답</summary>

   야코비 방법은 우변에 $`\{x\}^0`$만 쓰므로 $`x_1^1 = 7.85/3 = 2.616667`$, $`x_2^1 = -19.3/7 = -2.757143`$, $`x_3^1 = 71.4/10 = 7.140000`$이다. 가우스-자이델 방법은 $`x_2^1`$을 계산할 때 방금 구한 $`x_1^1 = 2.616667`$을 바로 쓰므로 $`x_2^1 = -2.794524`$가 된다.

   </details>

3. **반복행렬과 스펙트럼 반경:** $`A = \begin{bmatrix} 4 & 1 \\ 2 & 5 \end{bmatrix}`$에 야코비 방법을 적용할 때 반복행렬 $`M = I - \tilde{D}^{-1}A`$와 스펙트럼 반경을 구하고, 수렴 여부를 판정하라.

   <details>
   <summary>정답</summary>

   $`\tilde{D}^{-1}A = \begin{bmatrix} 1 & 1/4 \\ 2/5 & 1 \end{bmatrix}`$이므로 $`M = \begin{bmatrix} 0 & -1/4 \\ -2/5 & 0 \end{bmatrix}`$이다. 특성 방정식 $`\mu^2 - 1/10 = 0`$에서 고유값은 $`\mu = \pm 1/\sqrt{10} \approx \pm 0.316`$이므로 $`\rho(M) \approx 0.316 < 1`$이고 반복은 수렴한다. 오차는 대략 반복마다 0.316배로 줄어든다. 이 행렬은 대각 우위($`4 > 1`$, $`5 > 2`$)이므로 3.1절의 결과와도 일치한다.

   </details>

4. **전처리 없는 반복:** 전처리기 없이($`P = I`$) 반복할 때 $`A`$의 고유값이 0.5와 1.8이면 반복은 수렴하는가? 고유값 1.8이 2.5로 바뀌면 어떻게 되는가?

   <details>
   <summary>정답</summary>

   $`M = I - A`$의 고유값은 $`1 - 0.5 = 0.5`$와 $`1 - 1.8 = -0.8`$이므로 $`\rho(M) = 0.8 < 1`$이고 반복은 수렴한다. 1.8이 2.5로 바뀌면 $`M`$의 고유값 $`1 - 2.5 = -1.5`$가 생겨 $`\rho(M) = 1.5 > 1`$이므로 반복은 발산한다. 2.5는 1을 중심으로 하는 단위원 밖에 있다.

   </details>

5. **완화 계수:** 예제 12.2의 행을 교환한 시스템에서 초기 추측 0으로 $`\lambda = 1`$일 때 첫 반복의 $`x_1`$, $`x_2`$를 구하고, $`\lambda = 1.2`$일 때(6.5절)와 비교하라. $`0 < \lambda < 1`$은 어떤 경우에 쓰는가?

   <details>
   <summary>정답</summary>

   $`\lambda = 1`$이면 표준 가우스-자이델 방법이므로 $`x_1 = (8 + 2(0))/10 = 0.8`$, $`x_2 = (9 + 3(0.8))/12 = 0.95`$이다. $`\lambda = 1.2`$이면 각 갱신을 더 크게 움직이는 과완화가 일어나 $`x_1 = 0.96`$, $`x_2 = 1.188`$이 된다. $`0 < \lambda < 1`$인 부족 완화는 진동을 감쇠하여 수렴하지 않는 시스템을 수렴시키거나 수렴을 돕는 데 쓴다.

   </details>
