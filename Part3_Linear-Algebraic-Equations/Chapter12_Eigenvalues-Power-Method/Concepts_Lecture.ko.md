# 제12장 강의 — 선형계의 반복법

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 12

> **선수 지식**: [선형대수학] 행렬 연산 (제8-11장).
>
> **학습 목표**:
> 1. 지배 고유값 계산을 위해 거듭제곱법을 적용할 수 있다
> 2. 시프팅을 포함한 역거듭제곱법을 구현할 수 있다
> 3. 고유값 알고리즘의 수렴을 분석할 수 있다

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
- [요약](#요약)

---

<br>

## 1. 반복법 소개

반복법은 선형계 $[A]\{x\} = \{b\}$를 근사해의 수열을 생성하여 풉니다:

$$\mathbf{x}^{k+1} := \Psi(\mathbf{x}^k), \quad k \geq 0$$

초기 추측 $\mathbf{x}^0$에서 시작하여, 참해에 수렴할 때까지 매핑 $\Psi$를 반복적으로 적용합니다. 유한 번의 연산으로 정확한 답을 생성하는 직접법(가우스 소거법, LU 분해)과 달리, 반복법은 해에 점진적으로 접근하며 근사가 충분히 정확해지면 중단할 수 있습니다.

### 학습 목표

1. **가우스-자이델(Gauss-Seidel)** 방법과 **야코비(Jacobi)** 방법의 이해
2. **대각 우위(Diagonal Dominance)** 평가 방법 숙지
3. **완화법(Relaxation)**이 반복법의 수렴을 어떻게 개선하는지 이해

> **[선형대수]** 직접법은 $O(n^3)$의 비용이 들며 소규모에서 중규모 시스템에 잘 작동합니다. 매우 크고 희소한 시스템(공학 및 과학 계산에서 흔함)의 경우, 반복법은 행렬-벡터 곱만 필요하고 희소성을 활용할 수 있으므로 훨씬 효율적인 경우가 많습니다.

---

<br>

## 2. 가우스-자이델 방법(Gauss-Seidel Method)

### 2.1 문제 설정

다음과 같은 표준 형식의 $n$개 방정식이 주어졌다고 합시다:

$$[A]_{n \times n} \{x\}_{n \times 1} = \{b\}_{n \times 1}$$

여기서 $[A]$는 계수 행렬, $\{x\}$는 미지 벡터, $\{b\}$는 우변 벡터입니다.

### 2.2 3x3 시스템에 대한 유도

간단하게 $3 \times 3$ 시스템을 고려합니다:

$$\begin{bmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \\ x_3 \end{Bmatrix} = \begin{Bmatrix} b_1 \\ b_2 \\ b_3 \end{Bmatrix}$$

행을 전개하면:

$$a_{11} x_1 = b_1 - a_{12} x_2 - a_{13} x_3$$

$$a_{22} x_2 = b_2 - a_{21} x_1 - a_{23} x_3$$

$$a_{33} x_3 = b_3 - a_{31} x_1 - a_{32} x_2$$

각 방정식을 대각 미지수에 대해 풀면:

$$x_1 = \frac{b_1 - a_{12} x_2 - a_{13} x_3}{a_{11}}$$

$$x_2 = \frac{b_2 - a_{21} x_1 - a_{23} x_3}{a_{22}}$$

$$x_3 = \frac{b_3 - a_{31} x_1 - a_{32} x_2}{a_{33}}$$

### 2.3 위첨자를 사용한 반복 스킴

가우스-자이델 방법은 반복 인덱스를 도입합니다. 반복 $j$에서 다음을 계산합니다:

$$x_1^{j} = \frac{b_1 - a_{12} x_2^{j-1} - a_{13} x_3^{j-1}}{a_{11}}$$

$$x_2^{j} = \frac{b_2 - a_{21} x_1^{j} - a_{23} x_3^{j-1}}{a_{22}}$$

$$x_3^{j} = \frac{b_3 - a_{31} x_1^{j} - a_{32} x_2^{j}}{a_{33}}$$

가우스-자이델의 핵심 특징은 **가장 최근에 계산된 값을 사용 가능한 즉시 사용**한다는 것입니다. $x_2^{j}$를 계산할 때 $x_1^{j-1}$이 아닌 $x_1^{j}$(현재 반복에서 방금 계산된 값)를 이미 사용합니다. 마찬가지로 $x_3^{j}$를 계산할 때 $x_1^{j}$와 $x_2^{j}$ 모두를 사용합니다.

> **[선형대수]** 이 "순차 갱신" 성질이 가우스-자이델을 야코비 방법과 구분합니다. 새로운 정보가 즉시 반영되므로 일반적으로 더 빠른 수렴을 보이지만, 계산이 본질적으로 순차적입니다 -- $x_2^j$는 $x_1^j$가 알려질 때까지 계산할 수 없습니다.

### 2.4 수렴 판정 기준

참해로 수렴할 때까지 절차를 반복합니다. 각 변수의 **근사 상대 오차**는 다음과 같습니다:

$$\varepsilon_{a,i} = \left| \frac{x_i^{j} - x_i^{j-1}}{x_i^{j}} \right| \times 100\%$$

**모든** 변수에 대해 $\varepsilon_{a,i} \leq \varepsilon_s$일 때 반복을 중단합니다. 여기서 $\varepsilon_s$는 사용자가 지정한 허용오차입니다.

### 2.5 예제 12.1: 가우스-자이델 반복

**주어진 시스템:**

$$3x_1 - 0.1x_2 - 0.2x_3 = 7.85$$

$$0.1x_1 + 7x_2 - 0.3x_3 = -19.3$$

$$0.3x_1 - 0.2x_2 + 10x_3 = 71.4$$

**단계 1 -- 재배열**: 각 방정식을 대각 미지수에 대해 풉니다:

$$x_1 = \frac{7.85 + 0.1x_2 + 0.2x_3}{3}$$

$$x_2 = \frac{-19.3 - 0.1x_1 + 0.3x_3}{7}$$

$$x_3 = \frac{71.4 - 0.3x_1 + 0.2x_2}{10}$$

**단계 2 -- 첫 번째 반복** ($j = 1$), 초기 추측 $\{x\}^0 = \{0, 0, 0\}$:

$$x_1^1 = \frac{7.85 + 0.1(0) + 0.2(0)}{3} = 2.616667$$

$$x_2^1 = \frac{-19.3 - 0.1(2.616667) + 0.3(0)}{7} = -2.794524$$

$$x_3^1 = \frac{71.4 - 0.3(2.616667) + 0.2(-2.794524)}{10} = 7.005610$$

$x_2^1$ 계산 시 $x_1^1$이 즉시 사용되었고, $x_3^1$ 계산 시 $x_1^1$과 $x_2^1$ 모두 사용되었음에 주목하십시오.

**단계 3 -- 두 번째 반복** ($j = 2$):

$$x_1^2 = \frac{7.85 + 0.1(-2.794524) + 0.2(7.005610)}{3} = 2.990557$$

$$x_2^2 = \frac{-19.3 - 0.1(2.990557) + 0.3(7.005610)}{7} = -2.499625$$

$$x_3^2 = \frac{71.4 - 0.3(2.990557) + 0.2(-2.499625)}{10} = 7.000291$$

**단계 4 -- 상대 오차 계산** (반복 2에 대해):

$$\varepsilon_{a,1}^2 = \left| \frac{2.990557 - 2.616667}{2.990557} \right| \times 100\% = 12.5\%$$

$x_2$와 $x_3$의 상대 오차는 각각 $11.8\%$와 $0.0766\%$입니다.

**모든** 상대 오차가 사용자 지정 허용오차 $\varepsilon_s$보다 작아질 때까지 이 절차를 반복합니다.

---

<br>

## 3. 수렴과 대각 우위(Diagonal Dominance)

### 3.1 대각 우위 조건

가우스-자이델 방법은 일부 시스템에서 **발산할 수 있습니다**. 그러나 시스템이 **대각 우위(Diagonally Dominant)**이면 **반드시 수렴합니다**.

행렬이 대각 우위라 함은 모든 행 $i$에 대해:

$$|a_{ii}| > \sum_{\substack{j=1 \\ j \neq i}}^{n} |a_{ij}|$$

즉, 각 대각 원소의 절대값이 해당 행의 다른 모든 원소의 절대값 합보다 **엄격하게 큰** 것입니다.

### 3.2 예제: 3x3 대각 우위 검사

$3 \times 3$ 시스템의 경우:

$$\begin{bmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{bmatrix} \begin{Bmatrix} x_1 \\ x_2 \\ x_3 \end{Bmatrix} = \begin{Bmatrix} b_1 \\ b_2 \\ b_3 \end{Bmatrix}$$

검증할 세 가지 조건은:

$$|a_{11}| > |a_{12}| + |a_{13}|$$

$$|a_{22}| > |a_{21}| + |a_{23}|$$

$$|a_{33}| > |a_{31}| + |a_{32}|$$

예제 12.1의 경우: $|3| > |-0.1| + |-0.2|$, $|7| > |0.1| + |-0.3|$, $|10| > |0.3| + |-0.2|$. 모든 조건이 성립하므로 수렴이 보장됩니다.

### 3.3 충분조건 vs. 필요조건

대각 우위는 수렴의 **충분조건**이지만 **필요조건은 아닙니다**.

- 행렬이 대각 우위이면 수렴이 보장됨
- 행렬이 대각 우위가 **아닌** 경우에도 방법이 **수렴할 수 있음** (다만 보장되지 않을 뿐)
- 가우스-자이델 방법은 **대칭 양정치(Symmetric Positive Definite)** 행렬에도 작동함

> **[선형대수]** 시스템이 원래 형태에서 대각 우위가 아닌 경우, **행 교환**(피벗)을 통해 대각 우위를 달성할 수 있습니다. 이것이 예제 12.2에서 사용된 전략입니다.

---

<br>

## 4. 가우스-자이델 방법의 행렬 형식

### 4.1 행렬 분할: L, D, U 분해

계수 행렬 $[A]$는 다음과 같이 분해할 수 있습니다:

$$[A] = \tilde{L} + \tilde{D} + \tilde{U}$$

여기서:

- $\tilde{L}$은 $[A]$의 **엄밀 하삼각** 부분 (대각선 아래, 대각선 위와 대각선에는 0)
- $\tilde{D}$는 $[A]$의 **대각** 행렬 (대각 원소만)
- $\tilde{U}$는 $[A]$의 **엄밀 상삼각** 부분 (대각선 위, 대각선 아래와 대각선에는 0)

$3 \times 3$ 경우:

$$\tilde{L} = \begin{bmatrix} 0 & 0 & 0 \\ a_{21} & 0 & 0 \\ a_{31} & a_{32} & 0 \end{bmatrix}, \quad \tilde{D} = \begin{bmatrix} a_{11} & 0 & 0 \\ 0 & a_{22} & 0 \\ 0 & 0 & a_{33} \end{bmatrix}, \quad \tilde{U} = \begin{bmatrix} 0 & a_{12} & a_{13} \\ 0 & 0 & a_{23} \\ 0 & 0 & 0 \end{bmatrix}$$

> **[선형대수]** 이 $\tilde{L} + \tilde{D} + \tilde{U}$ 분할은 제10장의 LU 분해와 다릅니다. 여기서 $\tilde{L}$, $\tilde{D}$, $\tilde{U}$는 단순히 원래 행렬 $A$의 부분입니다 -- 소거가 수행되지 않습니다. LU 인수와 구별하기 위해 틸다 표기를 사용합니다.

### 4.2 반복의 행렬 형식

가우스-자이델 방정식을 행렬 형식으로 재배열하면:

$$[\tilde{L} + \tilde{D}]\{x\}^{j} = \{b\} - [\tilde{U}]\{x\}^{j-1}$$

$\tilde{L}_* = \tilde{L} + \tilde{D}$ ($A$의 대각선을 포함한 하삼각 부분)로 놓으면, 행렬-벡터 형식은:

$$\{x\}^{j} = (\tilde{L}_*)^{-1}\left(\{b\} - [\tilde{U}]\{x\}^{j-1}\right)$$

좌변은 **하삼각** 시스템을 포함하므로 **전방 대입**으로 효율적으로 풀 수 있습니다 -- 명시적 행렬 역전은 필요하지 않습니다.

---

<br>

## 5. 야코비 반복법(Jacobi Iteration)

### 5.1 가우스-자이델과의 비교

**야코비 반복법**은 가우스-자이델 방법과 유사하지만, $j$번째 반복의 모든 변수가 **오직** $(j-1)$번째 반복 값만을 사용하여 갱신됩니다. 같은 반복 내에서 새로 계산된 값은 사용되지 않습니다.

| 특성 | 가우스-자이델 | 야코비 |
|---|---|---|
| 갱신 전략 | 최신 값을 즉시 사용 | 이전 반복 값만 사용 |
| 병렬성 | 순차적 (갱신 순서에 의존) | **자연적으로 병렬** (모든 갱신이 독립) |
| 수렴 속도 | 일반적으로 빠름 | 일반적으로 느림 |
| 저장 공간 | 해 벡터 1개 | 2개 필요 (이전 값과 새 값) |

야코비 방법은 $[A]$가 대각 우위이면 참해로 수렴합니다.

### 5.2 성분 형식

$3 \times 3$ 시스템의 경우, 야코비 반복은 다음을 계산합니다:

$$x_1^{j} = \frac{b_1 - a_{12} x_2^{j-1} - a_{13} x_3^{j-1}}{a_{11}}$$

$$x_2^{j} = \frac{b_2 - a_{21} x_1^{j-1} - a_{23} x_3^{j-1}}{a_{22}}$$

$$x_3^{j} = \frac{b_3 - a_{31} x_1^{j-1} - a_{32} x_2^{j-1}}{a_{33}}$$

**모든** 우변 값이 위첨자 $j-1$을 사용함에 주목하십시오. 이는 각 $x_i^j$를 **독립적으로** 계산할 수 있음을 의미하며, 야코비 방법이 병렬 계산에 이상적인 이유입니다.

### 5.3 행렬 형식

$$[\tilde{D}]\{x\}^{j} = \{b\} - [\tilde{L} + \tilde{U}]\{x\}^{j-1}$$

$$\{x\}^{j} = [\tilde{D}]^{-1}\left(\{b\} - [\tilde{L} + \tilde{U}]\{x\}^{j-1}\right)$$

$\tilde{D}$가 대각행렬이므로 역행렬은 자명합니다: $[\tilde{D}]^{-1}_{ii} = 1/a_{ii}$.

---

<br>

## 6. 완화법(Relaxation Method)

### 6.1 완화 공식

완화법은 **수렴을 향상**시키기 위해 가우스-자이델 방법을 약간 수정한 것입니다. 아이디어는 이전 값과 새로 계산된 가우스-자이델 값을 가중 인수를 사용하여 결합하는 것입니다:

$$x_i^{\text{new}} = \lambda \, x_i^{\text{new(GS)}} + (1 - \lambda) \, x_i^{\text{old}}$$

여기서 $\lambda$는 $0$과 $2$ 사이의 값을 갖는 가중 인수이고, $x_i^{\text{new(GS)}}$는 표준 가우스-자이델 방법이 생성할 값입니다.

### 6.2 완화의 종류

| $\lambda$ 범위 | 명칭 | 효과 |
|---|---|---|
| $\lambda = 1$ | 수정 없음 | 표준 가우스-자이델 방법 |
| $0 < \lambda < 1$ | **부족 완화(Under-relaxation)** | 진동을 감쇠; 비수렴 시스템에 도움이 될 수 있음 |
| $1 < \lambda \leq 2$ | **과완화(Over-relaxation)** (SOR) | 가우스-자이델이 느릴 때 수렴을 가속 |

> **[선형대수]** $\lambda$의 최적값은 특정 문제에 따라 다르며 일반적으로 사전에 알 수 없습니다. SOR(연속 과완화) 방법의 경우, 특수한 행렬 구조(예: 삼대각 시스템)에 대한 이론적 결과가 있지만, 실제로는 $\lambda$를 실험적으로 결정하는 경우가 많습니다.

### 6.3 3x3 시스템에 적용한 완화법

$3 \times 3$ 시스템의 경우 완화 갱신은 다음과 같습니다:

$$x_1^{\text{new}} = \frac{b_1 - a_{12} x_2^{\text{old}} - a_{13} x_3^{\text{old}}}{a_{11}} \cdot \lambda + (1 - \lambda) \, x_1^{\text{old}}$$

$$x_2^{\text{new}} = \frac{b_2 - a_{21} x_1^{\text{new}} - a_{23} x_3^{\text{old}}}{a_{22}} \cdot \lambda + (1 - \lambda) \, x_2^{\text{old}}$$

$$x_3^{\text{new}} = \frac{b_3 - a_{31} x_1^{\text{new}} - a_{32} x_2^{\text{new}}}{a_{33}} \cdot \lambda + (1 - \lambda) \, x_3^{\text{old}}$$

가우스-자이델 패턴이 유지됨에 주목하십시오: $x_1^{\text{new}}$이 $x_2^{\text{new}}$ 계산에 사용되고, 둘 다 $x_3^{\text{new}}$에 사용됩니다.

### 6.4 완화법의 행렬 형식 (SOR)

전개하고 항을 정리하면, 완화법의 행렬 형식은:

$$[\tilde{D} + \lambda \tilde{L}]\{x\}^{\text{new}} = \lambda\{b\} - \left(\lambda[\tilde{U}] + (\lambda - 1)[\tilde{D}]\right)\{x\}^{\text{old}}$$

**유도 개요:**

성분 형식에서 출발하여 $a_{ii}$를 곱하면:

$$a_{ii} x_i^{\text{new}} = \left(b_i - \sum_{k<i} a_{ik} x_k^{\text{new}} - \sum_{k>i} a_{ik} x_k^{\text{old}}\right) \lambda + (1-\lambda) \, a_{ii} \, x_i^{\text{old}}$$

전개하고 행렬 표기법으로 재배열하면:

$$a_{ii} x_i^{\text{new}} = \lambda b_i - \lambda \sum_{k>i} a_{ik} x_k^{\text{old}} + (1-\lambda) a_{ii} x_i^{\text{old}} - \lambda \sum_{k<i} a_{ik} x_k^{\text{new}}$$

"새" 항을 좌변으로, "이전" 항을 우변으로 이동하면 최종 행렬 형식이 됩니다.

### 6.5 예제 12.2: lambda = 1.2인 완화법

**주어진 조건:** $\lambda = 1.2$, $\varepsilon_s = 10\%$

**원래 시스템:**

$$-3x_1 + 12x_2 = 9$$

$$10x_1 - 2x_2 = 8$$

**단계 1 -- 행 교환** (대각 우위를 위해):

$$10x_1 - 2x_2 = 8$$

$$-3x_1 + 12x_2 = 9$$

교환 후: $|10| > |-2|$이고 $|12| > |-3|$. 시스템이 이제 대각 우위입니다.

재배열:

$$x_1^{j} = \frac{8 + 2x_2^{j-1}}{10}, \quad x_2^{j} = \frac{9 + 3x_1^{j}}{12}$$

**반복 1** (초기 추측 $x_1^0 = 0$, $x_2^0 = 0$):

가우스-자이델 단계:

$$x_1^1 = \frac{8 + 2(0)}{10} = 0.8$$

완화 적용:

$$x_1^{\text{new}} = \lambda \cdot x_1^1 + (1-\lambda) \cdot x_1^0 = 1.2(0.8) + (-0.2)(0) = 0.96$$

가우스-자이델 단계 (완화된 $x_1^{\text{new}} = 0.96$ 사용):

$$x_2^1 = \frac{9 + 3(0.96)}{12} = 0.99$$

완화 적용:

$$x_2^{\text{new}} = 1.2(0.99) + (-0.2)(0) = 1.188$$

**반복 2** ($x_1^1 = 0.96$, $x_2^1 = 1.188$ 사용):

가우스-자이델 단계:

$$x_1^2 = \frac{8 + 2(1.188)}{10} = 1.0376$$

완화 적용:

$$x_1^{\text{new}} = 1.2(1.0376) + (-0.2)(0.96) = 1.05312$$

가우스-자이델 단계:

$$x_2^2 = \frac{9 + 3(1.05312)}{12} = 1.01328$$

완화 적용:

$$x_2^{\text{new}} = 1.2(1.01328) + (-0.2)(1.188) = 0.978336$$

상대 오차:

$$\varepsilon_{a,1} = \left|\frac{1.05312 - 0.96}{1.05312}\right| \times 100\% = 8.84\%$$

$$\varepsilon_{a,2} = \left|\frac{0.978336 - 1.188}{0.978336}\right| \times 100\% = 21.43\%$$

$\varepsilon_{a,2} > \varepsilon_s = 10\%$이므로 반복을 계속합니다.

**반복 3:**

$$x_1^3 = 0.984177, \quad \varepsilon_{a,1} = 7.01\%$$

$$x_2^3 = 0.999586, \quad \varepsilon_{a,2} = 2.13\%$$

두 오차 모두 $10\%$ 미만이므로 반복이 수렴합니다. 참해는 $x_1 = 1$, $x_2 = 1$입니다.

---

<br>

## 7. 전처리(Preconditioning)

### 7.1 일반 프레임워크

세 가지 반복법(가우스-자이델, 야코비, 완화법)은 모두 선형계 $A\mathbf{x} = \mathbf{b}$를 풀기 위한 **전통적 반복법**이며, 이를 식 (1)이라 합니다.

$A$에 대한 근사치인 **전처리기(Preconditioner)** $P$를 도입합니다. $A\mathbf{x} = \mathbf{b}$의 양변에 $P\mathbf{x}$를 더하면:

$$P\mathbf{x} + A\mathbf{x} = P\mathbf{x} + \mathbf{b}$$

$$P\mathbf{x} = (P - A)\mathbf{x} + \mathbf{b} \quad \cdots (2)$$

해 $\mathbf{x}^{k+1}$을 반복적으로 갱신하면:

$$P\mathbf{x}^{k+1} = (P - A)\mathbf{x}^k + \mathbf{b} \quad \cdots (3)$$

식 (3)이 식 (1)보다 **풀기 쉽기를** 기대합니다. $P$는 시스템을 저렴하게 풀 수 있는 행렬(예: 삼각행렬 또는 대각행렬)로 선택되기 때문입니다.

### 7.2 각 방법의 전처리기

각 반복법은 특정 전처리기 $P$의 선택에 대응합니다:

| 방법 | 전처리기 $P$ | $P - A$ |
|---|---|---|
| **가우스-자이델** | $P = \tilde{L} + \tilde{D}$ | $P - A = -\tilde{U}$ |
| **야코비** | $P = \tilde{D}$ | $P - A = -\tilde{L} - \tilde{U}$ |
| **완화법 (SOR)** | $P = \tilde{D} + \lambda\tilde{L}$ | $P - A = (\lambda - 1)\tilde{L} - \tilde{U}$ |

결과 반복 방정식:

- **가우스-자이델:** $[\tilde{L} + \tilde{D}]\{x\}^j = \{b\} - [\tilde{U}]\{x\}^{j-1}$
- **야코비:** $[\tilde{D}]\{x\}^j = \{b\} - [\tilde{L} + \tilde{U}]\{x\}^{j-1}$
- **완화법:** $[\tilde{D} + \lambda\tilde{L}]\{x\}^{\text{new}} = \lambda\{b\} - (\lambda[\tilde{U}] + (\lambda-1)[\tilde{D}])\{x\}^{\text{old}}$

### 7.3 오차 분석과 수렴

반복 $k$에서의 오차를 다음과 같이 정의합니다:

$$\mathbf{e}^k = \mathbf{x} - \mathbf{x}^k$$

여기서 $\mathbf{x}$는 참해입니다. 식 (2)에서 식 (3)을 빼면:

$$P(\mathbf{x} - \mathbf{x}^{k+1}) = (P - A)(\mathbf{x} - \mathbf{x}^k)$$

$$P\mathbf{e}^{k+1} = (P - A)\mathbf{e}^k$$

$$\mathbf{e}^{k+1} = (I - P^{-1}A)\mathbf{e}^k$$

**반복행렬(Iteration Matrix)**을 정의합니다:

$$M = I - P^{-1}A$$

그러면 각 단계는 오차 벡터에 $M$을 곱합니다:

$$\mathbf{e}^{k+1} = M\mathbf{e}^k$$

수렴은 **$M$의 고유값**에 의해 지배됩니다. 수렴을 위해서는 $M$의 모든 고유값이 $|\lambda(M)| < 1$을 만족해야 합니다.

### 7.4 스펙트럼 반경(Spectral Radius)

$M$의 (절대값에서) 가장 큰 고유값을 **스펙트럼 반경**이라 합니다:

$$\rho(M) = \max(|\lambda(M)|)$$

**수렴을 위한 필요충분조건:**

$$\rho(M) < 1$$

오차 벡터가 $M$의 고유벡터인 경우, 다음 단계의 오차는:

$$\mathbf{e}^{k+1} = M\mathbf{e}^k = \lambda \mathbf{e}^k$$

$$\Rightarrow \mathbf{e}^k = \lambda^k \mathbf{e}^0$$

따라서 오차는 $\lambda^k$로 기하급수적으로 감소합니다. $\rho(M)$이 작을수록 수렴이 빠릅니다.

> **[선형대수]** 이것이 대각 우위가 수렴을 보장하는 근본적인 이유입니다. 대각 우위 행렬의 경우, 가우스-자이델과 야코비 방법의 반복행렬 $M$의 스펙트럼 반경이 1보다 작습니다. SOR에서 $\lambda$의 선택은 $\rho(M)$을 최소화하는 것을 목표로 합니다.

### 7.5 I - A의 고유값 관계

전처리기를 사용하지 않으면 ($P = I$), $M = I - A$입니다.

**정리:** $[A]$의 고유값이 $\lambda_1, \lambda_2, \ldots, \lambda_n$이면, $[I - A]$의 고유값은 $1 - \lambda_1, 1 - \lambda_2, \ldots, 1 - \lambda_n$입니다.

**증명:**

$A\mathbf{v} = \lambda\mathbf{v}$ (여기서 $\mathbf{v}$는 고유벡터)가 주어지면:

$$(I - A)\mathbf{v} = \mathbf{v} - A\mathbf{v} = \mathbf{v} - \lambda\mathbf{v} = (1 - \lambda)\mathbf{v}$$

따라서 $(1 - \lambda)$은 같은 고유벡터 $\mathbf{v}$를 가진 $[I - A]$의 고유값입니다.

**결론:** 전처리 없이는, 반복이 수렴하려면 $A$의 모든 고유값이 복소평면에서 1을 중심으로 한 단위원 안에 있어야 합니다.

---

<br>

## 8. Python 구현

### 가우스-자이델 방법

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

### 야코비 방법

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

### 완화법 (SOR)

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

### 대각 우위 검사

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

| 주제 | 핵심 공식 / 개념 |
|---|---|
| **일반 반복** | $\mathbf{x}^{k+1} := \Psi(\mathbf{x}^k)$, $k \geq 0$ |
| **가우스-자이델** | $x_i^j = \frac{1}{a_{ii}}\left(b_i - \sum_{k<i} a_{ik} x_k^j - \sum_{k>i} a_{ik} x_k^{j-1}\right)$; 최신 값 사용 |
| **가우스-자이델 (행렬)** | $[\tilde{L} + \tilde{D}]\{x\}^j = \{b\} - [\tilde{U}]\{x\}^{j-1}$ |
| **야코비** | $x_i^j = \frac{1}{a_{ii}}\left(b_i - \sum_{k \neq i} a_{ik} x_k^{j-1}\right)$; 모두 이전 반복에서 |
| **야코비 (행렬)** | $[\tilde{D}]\{x\}^j = \{b\} - [\tilde{L} + \tilde{U}]\{x\}^{j-1}$ |
| **완화법 (SOR)** | $x_i^{\text{new}} = \lambda \, x_i^{\text{GS}} + (1-\lambda) \, x_i^{\text{old}}$ |
| **완화법 (행렬)** | $[\tilde{D} + \lambda\tilde{L}]\{x\}^{\text{new}} = \lambda\{b\} - (\lambda[\tilde{U}] + (\lambda-1)[\tilde{D}])\{x\}^{\text{old}}$ |
| **대각 우위** | $\|a_{ii}\| > \sum_{j \neq i} \|a_{ij}\|$ (충분조건, 필요조건 아님) |
| **수렴 판정 기준** | $\varepsilon_{a,i} = \left\|\frac{x_i^j - x_i^{j-1}}{x_i^j}\right\| \times 100\% \leq \varepsilon_s$ (모든 $i$에 대해) |
| **전처리** | $P\mathbf{x}^{k+1} = (P-A)\mathbf{x}^k + \mathbf{b}$; 반복행렬 $M = I - P^{-1}A$ |
| **스펙트럼 반경** | $\rho(M) = \max(\|\lambda(M)\|) < 1$ (수렴에 필요) |
| **고유값 이동** | $A$의 고유값이 $\lambda_i$이면, $I-A$의 고유값은 $1-\lambda_i$ |
