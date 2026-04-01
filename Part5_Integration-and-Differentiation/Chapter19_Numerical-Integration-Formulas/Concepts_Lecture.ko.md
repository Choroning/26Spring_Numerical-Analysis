# 제19장 강의 — 수치 적분 공식(Numerical Integration Formulas)

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 19

> **선수 지식**: [미적분학] 적분 (제1-18장).
>
> **학습 목표**:
> 1. 사다리꼴 공식과 심프슨 공식을 적용할 수 있다
> 2. 수치 적분의 오차 추정을 유도할 수 있다
> 3. 롬버그 적분과 가우스 구적법을 구현할 수 있다

---

<br>

## 목차

- [1. 수치 적분 개요](#1-수치-적분-개요)
  - [1.1 동기 — 왜 수치 적분이 필요한가?](#11-동기--왜-수치-적분이-필요한가)
  - [1.2 기본 아이디어: 다항식 근사](#12-기본-아이디어-다항식-근사)
  - [1.3 라그랑주 다항식 복습](#13-라그랑주-다항식-복습)
  - [1.4 구적법과 뉴턴-코츠 공식](#14-구적법과-뉴턴-코츠-공식)
- [2. 사다리꼴 법칙 (n = 1)](#2-사다리꼴-법칙-n--1)
  - [2.1 선형 보간으로부터의 유도](#21-선형-보간으로부터의-유도)
  - [2.2 최종 공식](#22-최종-공식)
  - [2.3 절단 오차](#23-절단-오차)
- [3. 심프슨의 1/3 법칙 (n = 2)](#3-심프슨의-13-법칙-n--2)
  - [3.1 이차 보간으로부터의 유도](#31-이차-보간으로부터의-유도)
  - [3.2 변수 치환을 통한 기저 적분 계산](#32-변수-치환을-통한-기저-적분-계산)
  - [3.3 최종 공식](#33-최종-공식)
  - [3.4 절단 오차](#34-절단-오차)
- [4. 심프슨의 3/8 법칙 (n = 3)](#4-심프슨의-38-법칙-n--3)
  - [4.1 삼차 보간으로부터의 유도](#41-삼차-보간으로부터의-유도)
  - [4.2 최종 공식](#42-최종-공식)
  - [4.3 절단 오차](#43-절단-오차)
- [5. 복합 사다리꼴 법칙](#5-복합-사다리꼴-법칙)
  - [5.1 아이디어: 구간 분할](#51-아이디어-구간-분할)
  - [5.2 유도](#52-유도)
  - [5.3 최종 공식](#53-최종-공식)
  - [5.4 오차 분석](#54-오차-분석)
- [6. 복합 심프슨의 1/3 법칙](#6-복합-심프슨의-13-법칙)
  - [6.1 유도](#61-유도)
  - [6.2 최종 공식](#62-최종-공식)
  - [6.3 오차 분석](#63-오차-분석)
- [7. 복합 심프슨의 3/8 법칙](#7-복합-심프슨의-38-법칙)
  - [7.1 오차 분석](#71-오차-분석)
- [8. 법칙의 결합 — 실용적 응용](#8-법칙의-결합--실용적-응용)
- [9. Python 구현](#9-python-구현)
- [요약](#요약)

---

<br>

## 1. 수치 적분 개요

### 1.1 동기 — 왜 수치 적분이 필요한가?

수치 적분은 구간 위에서 함수의 정적분을 근사하는 문제를 다룬다:

$$I = \int_a^b f(x)\,dx \approx I_n$$

많은 실제 적분은 **해석적으로 계산할 수 없다**. 강의에서 제시된 대표적인 동기 부여 예제는 반도체 웨이퍼 제조업체와 관련된 것이다.

**예제 — 웨이퍼 품질 관리:**

한 제조업체가 2인치 직경의 웨이퍼를 생산한다. 실제 직경은 평균 $\mu = 2$ 인치, 표준편차 $\sigma = 0.01$ 인치의 정규분포를 따른다. 규격은 직경이 1.985에서 2.02 인치 사이일 것을 요구한다. 웨이퍼가 합격할 확률은 얼마인가?

정규 확률밀도함수는 다음과 같다:

$$p(x) = \frac{1}{\sigma\sqrt{2\pi}} \exp\!\left(\frac{-(x - \mu)^2}{2\sigma^2}\right)$$

합격 웨이퍼의 확률은:

$$P(1.985 \leq x \leq 2.02) = \int_{1.985}^{2.02} p(x)\,dx = \frac{1}{\sigma\sqrt{2\pi}} \int_{1.985}^{2.02} \exp\!\left(\frac{-(x - \mu)^2}{2\sigma^2}\right) dx$$

이 적분은 **닫힌 형태의 역도함수가 존재하지 않으므로**, 수치 적분 방법이 필요하다.

### 1.2 기본 아이디어: 다항식 근사

전략은 **피적분함수를 적분하기 쉬운 다항식으로 대체**하는 것이다:

$$p(x) \approx p_n(x) = a_0 + a_1 x + \cdots + a_n x^n$$

라그랑주 보간(Lagrange Interpolation)을 이용하면, 이 다항식은 다음과 같이 쓸 수 있다:

$$p_n(x) = \sum_{j=0}^{n} p_n(x_j)\,\ell_j(x)$$

여기서 $\ell_j(x)$는 라그랑주 기저 다항식이고, $x_0, x_1, \ldots, x_n$은 보간 절점이다.

### 1.3 라그랑주 다항식 복습

두 점 $x_1, x_2$를 사용하는 1차(first-order)의 경우, 라그랑주 기저 함수는 다음과 같다:

$$N_1(x) = \frac{x - x_2}{x_1 - x_2}, \qquad N_2(x) = \frac{x_1 - x}{x_1 - x_2}$$

일반적으로 $n$차 라그랑주 기저 다항식은:

$$\ell_j(x) = \prod_{\substack{i=0 \\ j \neq i}}^{n} \frac{x - x_i}{x_j - x_i}, \qquad j = 0, 1, \ldots, n$$

### 1.4 구적법과 뉴턴-코츠 공식

$p(x) = p_n(x) + e_n(x)$이고 $e_n(x)$가 보간 오차라면:

$$\int_a^b p(x)\,dx = \int_a^b p_n(x)\,dx + \int_a^b e_n(x)\,dx$$

다항식 부분을 전개하면:

$$\int_a^b p_n(x)\,dx = \int_a^b \sum_{j=0}^{n} p_n(x_j)\,\ell_j(x)\,dx = \sum_{j=0}^{n} p_n(x_j) \underbrace{\int_a^b \ell_j(x)\,dx}_{\text{하나의 수 (가중치)}}$$

> **[미적분]** 핵심적인 통찰은 $\int_a^b \ell_j(x)\,dx$가 $p_n(x)$에 **의존하지 않는다**는 것이다 — 이것은 절점의 위치와 구간 $[a,b]$에 의해서만 결정되는 고정된 수치적 가중치이다. 이것이 이 공식들을 **구적법(Quadrature Formulas)**이라고 부르는 이유이다.

$I = \int_a^b f(x)\,dx$를 $I_n = \int_a^b f_n(x)\,dx$로 근사하는 접근법(여기서 $f_n(x)$는 $n$차 다항식)이 **뉴턴-코츠 공식(Newton-Cotes Formulas)**의 기초이다. $n$의 선택에 따라 구체적인 법칙이 결정된다:

| $n$ | 법칙 이름 | 다항식 차수 |
|-----|-----------|------------|
| 1 | 사다리꼴 법칙(Trapezoidal Rule) | 선형(Linear) |
| 2 | 심프슨의 1/3 법칙(Simpson's 1/3 Rule) | 이차(Quadratic) |
| 3 | 심프슨의 3/8 법칙(Simpson's 3/8 Rule) | 삼차(Cubic) |

---

<br>

## 2. 사다리꼴 법칙(Trapezoidal Rule) (n = 1)

### 2.1 선형 보간으로부터의 유도

피적분함수 $f(x)$를 양 끝점 $(a, f(a))$와 $(b, f(b))$를 지나는 1차 다항식 $p_1(x)$로 근사한다:

$$p_1(x) = f(a)\,N_1(x) + f(b)\,N_2(x) = f(a)\!\left(\frac{x - b}{a - b}\right) + f(b)\!\left(\frac{x - a}{b - a}\right)$$

정리하면:

$$p_1(x) = \frac{f(b) - f(a)}{b - a}\,x + \frac{b\,f(a) - a\,f(b)}{b - a}$$

이제 적분한다:

$$\int_a^b f(x)\,dx \approx \int_a^b p_1(x)\,dx = \frac{f(b) - f(a)}{b - a} \int_a^b x\,dx + \frac{b\,f(a) - a\,f(b)}{b - a} \int_a^b 1\,dx$$

기본 적분을 계산하면:

$$\int_a^b x\,dx = \frac{1}{2}(b^2 - a^2), \qquad \int_a^b 1\,dx = b - a$$

대수적 정리를 거치면:

$$= \frac{f(b) - f(a)}{2}(a + b) + b\,f(a) - a\,f(b)$$

$$= \frac{1}{2}\bigl[a\,f(b) - a\,f(a) + b\,f(b) - b\,f(a)\bigr]$$

$$= \frac{1}{2}(b - a)\,f(b) + \frac{1}{2}(b - a)\,f(a)$$

### 2.2 최종 공식

$$\boxed{I \approx I_1 = \frac{b - a}{2}\bigl[f(a) + f(b)\bigr] = (b - a) \cdot \frac{f(a) + f(b)}{2}}$$

> **[미적분]** 사다리꼴 법칙은 양 끝점의 함수값으로 형성되는 **사다리꼴**의 넓이를 계산한다. $(b - a)$는 **너비**이고, $\frac{f(a) + f(b)}{2}$는 **평균 높이**이다.

### 2.3 절단 오차(Truncation Error)

사다리꼴 법칙의 오차는 다음과 같다:

$$\boxed{E_1 = I - I_1 = -\frac{1}{12}\,f''(\xi)\,(b - a)^3, \qquad a \leq \xi \leq b}$$

**중요한 결과:** $f(x) = mx + n$ (선형 함수)이면 $f''(x) = 0$이므로 $E_1 = 0$이다. **사다리꼴 법칙은 선형 함수**(1차 이하의 다항식)에 대해 정확하다.

> **[미적분]** 오차는 $(b-a)^3$에 비례하며 이계도함수 $f''$에 의존한다. 이는 (1) 더 좁은 구간이 더 작은 오차를 주고, (2) 곡률이 작은 함수가 이 법칙으로 잘 근사됨을 알려준다.

---

<br>

## 3. 심프슨의 1/3 법칙(Simpson's 1/3 Rule) (n = 2)

### 3.1 이차 보간으로부터의 유도

$f(x)$를 등간격의 세 점 $x_1 = a$, $x_2 = \frac{a+b}{2}$, $x_3 = b$을 지나는 2차 다항식 $p_2(x)$로 근사한다. 간격 크기는 $h = \frac{b - a}{2}$이다.

다음과 같이 정의한다:

$$f_1 = f(a), \qquad f_2 = f\!\left(\frac{a+b}{2}\right), \qquad f_3 = f(b)$$

보간 다항식은:

$$p_2(x) = f_1\,\ell_1(x) + f_2\,\ell_2(x) + f_3\,\ell_3(x)$$

여기서 라그랑주 기저 다항식은 다음과 같다:

$$\ell_1(x) = \frac{(x - x_2)(x - x_3)}{(x_1 - x_2)(x_1 - x_3)}, \quad \ell_2(x) = \frac{(x - x_1)(x - x_3)}{(x_2 - x_1)(x_2 - x_3)}, \quad \ell_3(x) = \frac{(x - x_1)(x - x_2)}{(x_3 - x_1)(x_3 - x_2)}$$

따라서 적분은 다음이 된다:

$$\int_a^b p_2(x)\,dx = f_1 \int_a^b \ell_1(x)\,dx + f_2 \int_a^b \ell_2(x)\,dx + f_3 \int_a^b \ell_3(x)\,dx$$

### 3.2 변수 치환을 통한 기저 적분 계산

$\int_a^b \ell_j(x)\,dx$를 계산하기 위해 다음 변수 치환을 적용한다:

$$x = x(\zeta) = \frac{b - a}{2}\,\zeta + \frac{b + a}{2}, \qquad dx = \frac{b - a}{2}\,d\zeta$$

이는 $[a, b] \to [-1, 1]$로 사상하고, 세 절점 $x_1 = a,\; x_2 = \frac{a+b}{2},\; x_3 = b$는 각각 $\zeta = -1, 0, 1$로 사상된다.

**$\zeta$-공간에서의 기저 함수:**

$$\ell_1(\zeta) = \frac{(\zeta - 0)(\zeta - 1)}{(-1 - 0)(-1 - 1)} = \frac{\zeta(\zeta - 1)}{2}$$

$$\ell_2(\zeta) = \frac{(\zeta - (-1))(\zeta - 1)}{(0 - (-1))(0 - 1)} = (1 - \zeta)(1 + \zeta) = 1 - \zeta^2$$

$$\ell_3(\zeta) = \frac{(\zeta - (-1))(\zeta - 0)}{(1 - (-1))(1 - 0)} = \frac{\zeta(1 + \zeta)}{2}$$

**적분 계산:**

$\ell_1$의 경우:

$$\int_a^b \ell_1(x)\,dx = \int_{-1}^{1} \frac{\zeta(\zeta - 1)}{2} \cdot \frac{b - a}{2}\,d\zeta = \frac{b - a}{4} \int_{-1}^{1} (\zeta^2 - \zeta)\,d\zeta$$

$\zeta$는 $[-1, 1]$ 위의 **홀수 함수**이므로 그 적분은 소멸한다:

$$= \frac{b - a}{4} \cdot \frac{2}{3} = \frac{b - a}{6}$$

대칭성에 의해 $\int_a^b \ell_3(x)\,dx = \frac{b - a}{6}$도 마찬가지이다.

$\ell_2$의 경우:

$$\int_a^b \ell_2(x)\,dx = \int_{-1}^{1} (1 - \zeta^2) \cdot \frac{b - a}{2}\,d\zeta = \frac{b - a}{2}\left(\zeta - \frac{\zeta^3}{3}\right)\Bigg|_{-1}^{1} = \frac{b - a}{2} \cdot \frac{4}{3} = \frac{4(b - a)}{6}$$

### 3.3 최종 공식

가중치를 결합하면:

$$\int_a^b f(x)\,dx \approx \frac{b - a}{6}\bigl(f_1 + 4f_2 + f_3\bigr)$$

$h = \frac{b - a}{2}$로 놓으면:

$$\boxed{I \approx I_2 = \frac{h}{3}(f_1 + 4f_2 + f_3) = \frac{b - a}{6}\bigl[f(a) + 4f\!\left(\tfrac{a+b}{2}\right) + f(b)\bigr]}$$

> **[미적분]** "1/3 법칙"이라는 이름은 공식의 $\frac{h}{3}$ 인수에서 유래한다. 중간점은 끝점보다 **4배**의 가중치를 받는데, 이는 구간 중심 근처에서 포물선의 적합도가 더 높음을 반영한다.

### 3.4 절단 오차(Truncation Error)

심프슨의 1/3 법칙의 오차는:

$$E_2 = I - I_2 = -\frac{1}{90}\,h^5\,f^{(4)}(\xi), \qquad a \leq \xi \leq b$$

$h = \frac{b - a}{2}$이므로:

$$\boxed{E_2 = -\frac{(b - a)^5}{2880}\,f^{(4)}(\xi)}$$

**핵심 통찰:** 심프슨의 1/3 법칙은 **2차** 다항식을 사용하지만, 오차는 **4계도함수** $f^{(4)}$에 의존한다. 이는 심프슨의 1/3 법칙이 **3차 이하의 다항식에 대해 정확**함을 의미한다(삼차함수의 경우 $f^{(4)} = 0$이므로).

> **[미적분]** 이러한 "보너스 정확도"는 적분 절점이 중점에 대해 대칭이기 때문에 홀수 거듭제곱 오차항이 상쇄되어 발생한다.

---

<br>

## 4. 심프슨의 3/8 법칙(Simpson's 3/8 Rule) (n = 3)

### 4.1 삼차 보간으로부터의 유도

$f(x)$를 간격 크기 $h = \frac{b - a}{3}$인 등간격 네 점을 지나는 3차 다항식 $p_3(x)$로 근사한다:

$$x_0 = a, \quad x_1 = a + h, \quad x_2 = a + 2h, \quad x_3 = b$$

보간 다항식은:

$$p_3(x) = f_1\,\ell_1(x) + f_2\,\ell_2(x) + f_3\,\ell_3(x) + f_4\,\ell_4(x)$$

적분하면:

$$\int_a^b p_3(x)\,dx = f_1 \int_a^b \ell_1(x)\,dx + f_2 \int_a^b \ell_2(x)\,dx + f_3 \int_a^b \ell_3(x)\,dx + f_4 \int_a^b \ell_4(x)\,dx$$

### 4.2 최종 공식

$$\boxed{I \approx I_3 = \frac{3h}{8}(f_0 + 3f_1 + 3f_2 + f_3), \qquad h = \frac{b - a}{3}}$$

> **[미적분]** "3/8 법칙"이라는 이름은 $h$에 곱해지는 계수 $\frac{3}{8}$에서 유래한다. 가중치는 $1 : 3 : 3 : 1$ 패턴을 따르며, 대칭적인 삼차 보간을 반영한다.

### 4.3 절단 오차(Truncation Error)

$$E_3 = I - I_3 = -\frac{3}{80}\,h^5\,f^{(4)}(\xi), \qquad a \leq \xi \leq b$$

$h = \frac{b - a}{3}$이므로:

$$\boxed{E_3 = -\frac{(b - a)^5}{6480}\,f^{(4)}(\xi)}$$

심프슨의 1/3 법칙과 마찬가지로, 3/8 법칙도 오차가 $f^{(4)}$에 의존하므로 **3차 이하의 다항식에 대해 정확**하다.

> **[미적분]** 단일 적용 오차를 비교하면: $|E_3| = \frac{(b-a)^5}{6480}|f^{(4)}(\xi)|$ 대 $|E_2| = \frac{(b-a)^5}{2880}|f^{(4)}(\xi)|$. 3/8 법칙은 단일 적용당 **더 정확**하지만($6480 > 2880$), 하나의 함수 평가가 추가로 필요하다.

---

<br>

## 5. 복합 사다리꼴 법칙(Composite Trapezoidal Rule)

### 5.1 아이디어: 구간 분할

**Q: 사다리꼴 법칙의 오차를 어떻게 줄일 수 있는가?**

**A:** 구간 $[a, b]$를 $n$개의 소구간(세그먼트)으로 나누고, 각 세그먼트에 사다리꼴 법칙을 적용한다. 소구간이 작을수록 각 구간에서의 선형 근사가 더 정확해진다.

$n$개의 세그먼트인 경우:

$$h = \frac{b - a}{n}, \qquad x_i = a + ih \;\;(i = 0, 1, \ldots, n)$$

여기서 $x_0 = a$이고 $x_n = b$이며, $f_i = f(x_i)$이다.

### 5.2 유도

적분을 $n$개의 부분 적분으로 나누고 각각에 사다리꼴 법칙을 적용한다:

$$I = \int_a^b f(x)\,dx = \int_{x_0}^{x_1} f(x)\,dx + \int_{x_1}^{x_2} f(x)\,dx + \cdots + \int_{x_{n-1}}^{x_n} f(x)\,dx$$

$$\approx (x_1 - x_0)\frac{f_0 + f_1}{2} + (x_2 - x_1)\frac{f_1 + f_2}{2} + \cdots + (x_n - x_{n-1})\frac{f_{n-1} + f_n}{2}$$

모든 소구간의 너비가 $h$로 동일하므로:

$$= \frac{h}{2}\bigl(f_0 + f_1 + f_1 + f_2 + f_2 + \cdots + f_{n-1} + f_{n-1} + f_n\bigr)$$

항을 정리하면(내부 점들은 두 번 나타남):

$$= \frac{h}{2}\left(f_0 + 2f_1 + 2f_2 + \cdots + 2f_{n-1} + f_n\right)$$

### 5.3 최종 공식

$$\boxed{I \approx \frac{h}{2}\left(f_0 + 2\sum_{i=1}^{n-1} f_i + f_n\right) = (b - a)\,\frac{f_0 + 2\displaystyle\sum_{i=1}^{n-1} f_i + f_n}{2n}}$$

이것이 **복합 사다리꼴 법칙(Composite Trapezoidal Rule)**이다. 동일한 해석이 적용된다: **너비** $\times$ **평균 높이**, 여기서 평균 높이는 끝점은 한 번, 내부 점은 두 번 세는 가중 평균이다.

### 5.4 오차 분석

각 소구간 $[x_{i-1}, x_i]$에서의 국소 오차는:

$$-\frac{1}{12}\,f''(\xi_i)\,h^3, \qquad x_{i-1} \leq \xi_i \leq x_i$$

$n$개의 국소 오차를 모두 합하면:

$$E_{1t} = -\frac{1}{12}\,h^3 \sum_{i=1}^{n} f''(\xi_i) = -\frac{1}{12}\,h^3 \cdot n \cdot \frac{\displaystyle\sum_{i=1}^{n} f''(\xi_i)}{n}$$

$f''(\xi_i)$ 값들의 평균은 $\overline{f''}$로 근사하고, $n = \frac{b - a}{h}$이므로:

$$\boxed{E_{1a} \approx -\frac{(b - a)^3}{12\,n^2}\,\overline{f''} = -\frac{(b - a)}{12}\,\overline{f''}\,h^2}$$

오차는 **$h^2$에 비례**한다. 따라서:

> **$n$을 두 배로 늘리면, 절단 오차는 1/4로 줄어든다!** ($h^2 \to (h/2)^2 = h^2/4$이므로.)

---

<br>

## 6. 복합 심프슨의 1/3 법칙(Composite Simpson's 1/3 Rule)

### 6.1 유도

$[a, b]$를 $n$개의 소구간으로 나눈다(**$n$은 짝수여야 한다**). $h = \frac{b - a}{n}$. 연속하는 한 쌍의 소구간에 심프슨의 1/3 법칙을 적용한다:

$$I = \int_{x_0}^{x_2} f(x)\,dx + \int_{x_2}^{x_4} f(x)\,dx + \cdots + \int_{x_{n-2}}^{x_n} f(x)\,dx$$

$$\approx \frac{h}{3}(f_0 + 4f_1 + f_2) + \frac{h}{3}(f_2 + 4f_3 + f_4) + \frac{h}{3}(f_4 + 4f_5 + f_6) + \cdots + \frac{h}{3}(f_{n-2} + 4f_{n-1} + f_n)$$

항을 정리하면:

$$= \frac{h}{3}\,f_0 + \frac{4h}{3}(f_1 + f_3 + f_5 + \cdots + f_{n-1}) + \frac{2h}{3}(f_2 + f_4 + f_6 + \cdots + f_{n-2}) + \frac{h}{3}\,f_n$$

### 6.2 최종 공식

$$\boxed{I \approx \frac{h}{3}\left(f_0 + 4\!\!\sum_{\substack{i=1,3,5,\ldots}}^{n-1}\!\! f_i + 2\!\!\sum_{\substack{i=2,4,6,\ldots}}^{n-2}\!\! f_i + f_n\right)}$$

동치인 형태:

$$= (b - a)\,\frac{f_0 + 4\displaystyle\sum_{\substack{i\,\text{odd}}} f_i + 2\displaystyle\sum_{\substack{i\,\text{even}}} f_i + f_n}{3n}$$

> **[미적분]** 가중치의 패턴은 $1, 4, 2, 4, 2, 4, \ldots, 2, 4, 1$이다. 홀수 인덱스 내부 점은 가중치 4, 짝수 인덱스 내부 점은 가중치 2, 끝점은 가중치 1을 받는다. 이는 $n$이 **짝수**일 것을 요구한다.

### 6.3 오차 분석

각 소구간 쌍에 대한 국소 오차는:

$$-\frac{h^5}{90}\,f^{(4)}(\xi_k), \qquad x_{2k-2} \leq \xi_k \leq x_{2k}$$

모든 $n/2$번의 적용에 걸쳐 합산하면:

$$E_{2t} = -\frac{h^5}{90} \sum_{k=1}^{n/2} f^{(4)}(\xi_k) = -\frac{h^5}{90} \cdot \frac{n}{2} \cdot \overline{f^{(4)}}$$

$h = \frac{b-a}{n}$이므로:

$$\boxed{E_{2a} \approx -\frac{(b - a)^5}{180\,n^4}\,\overline{f^{(4)}} = -\frac{(b - a)}{180}\,\overline{f^{(4)}}\,h^4}$$

오차는 **$h^4$에 비례**한다. 따라서:

> **$n$을 두 배로 늘리면, 절단 오차는 16배 감소한다!** ($h^4 \to (h/2)^4 = h^4/16$이므로.)

---

<br>

## 7. 복합 심프슨의 3/8 법칙(Composite Simpson's 3/8 Rule)

심프슨의 3/8 법칙을 3개의 소구간 그룹에 반복 적용한다(**$n$은 3의 배수여야 한다**). 유도는 복합 1/3 법칙과 동일한 패턴을 따른다.

### 7.1 오차 분석

3개의 소구간 그룹 각각에 대한 국소 오차는:

$$-\frac{3}{80}\,h^5\,f^{(4)}(\xi_k)$$

모든 $n/3$번의 적용에 걸쳐 합산하면:

$$E_{3t} = -\frac{3}{80}\,h^5 \sum_{k=1}^{n/3} f^{(4)}(\xi_k) \approx -\frac{1}{80}\,h^5 \cdot n \cdot \overline{f^{(4)}}$$

$$\boxed{E_{3a} \approx -\frac{(b - a)}{80}\,\overline{f^{(4)}}\,h^4}$$

복합 1/3 법칙과 마찬가지로 오차는 $h^4$에 비례하지만, 상수가 약간 다르다.

> **[미적분]** 3/8 법칙이 **단일 적용**에서는 1/3 법칙보다 더 정확하지만, 복합 형태에서는 3/8 법칙이 $h^4$ 단위당 더 큰 오차 상수를 가진다($\frac{1}{80}$ 대 $\frac{1}{180}$). 이는 3/8 법칙이 적용당 더 많은 소구간을 사용하기 때문이다.

---

<br>

## 8. 법칙의 결합 — 실용적 응용

실제로 세그먼트 수 $n$이 **짝수가 아니고**(복합 심프슨의 1/3에 필요) **3의 배수도 아닌**(복합 심프슨의 3/8에 필요) 경우, 두 법칙을 **결합**할 수 있다:

- 소구간의 일부에 심프슨의 **1/3 법칙**을 적용 (세그먼트 쌍 사용)
- 나머지 소구간에 심프슨의 **3/8 법칙**을 적용 (세그먼트 3개 사용)

예를 들어, 5개의 세그먼트인 경우: 처음 2개 세그먼트에 1/3 법칙을, 나머지 3개 세그먼트에 3/8 법칙을 사용한다(또는 그 반대).

---

<br>

## 9. Python 구현

### 사다리꼴 법칙 (단일 적용)

```python
def trapezoidal(f, a, b):
    """Single-application trapezoidal rule."""
    return (b - a) / 2 * (f(a) + f(b))
```

### 복합 사다리꼴 법칙

```python
def composite_trapezoidal(f, a, b, n):
    """Composite trapezoidal rule with n subintervals."""
    h = (b - a) / n
    x = [a + i * h for i in range(n + 1)]
    result = f(x[0]) + f(x[n])
    for i in range(1, n):
        result += 2 * f(x[i])
    return h / 2 * result
```

### 심프슨의 1/3 법칙 (단일 적용)

```python
def simpson_13(f, a, b):
    """Single-application Simpson's 1/3 rule."""
    h = (b - a) / 2
    return h / 3 * (f(a) + 4 * f(a + h) + f(b))
```

### 복합 심프슨의 1/3 법칙

```python
def composite_simpson_13(f, a, b, n):
    """Composite Simpson's 1/3 rule with n subintervals (n must be even)."""
    if n % 2 != 0:
        raise ValueError("n must be even for Simpson's 1/3 rule")
    h = (b - a) / n
    x = [a + i * h for i in range(n + 1)]
    result = f(x[0]) + f(x[n])
    for i in range(1, n, 2):      # odd indices: weight 4
        result += 4 * f(x[i])
    for i in range(2, n - 1, 2):  # even indices: weight 2
        result += 2 * f(x[i])
    return h / 3 * result
```

### 심프슨의 3/8 법칙 (단일 적용)

```python
def simpson_38(f, a, b):
    """Single-application Simpson's 3/8 rule."""
    h = (b - a) / 3
    return 3 * h / 8 * (f(a) + 3 * f(a + h) + 3 * f(a + 2 * h) + f(b))
```

### 전체 예제 — 모든 방법 비교

```python
import math

def f(x):
    """Example: integrate e^x from 0 to 1 (exact answer = e - 1)."""
    return math.exp(x)

a, b = 0, 1
exact = math.e - 1

print(f"Exact value: {exact:.10f}")
print(f"Trapezoidal (n=1):       {trapezoidal(f, a, b):.10f}")
print(f"Simpson's 1/3 (n=2):     {simpson_13(f, a, b):.10f}")
print(f"Simpson's 3/8 (n=3):     {simpson_38(f, a, b):.10f}")
print(f"Composite Trap (n=10):   {composite_trapezoidal(f, a, b, 10):.10f}")
print(f"Composite S-1/3 (n=10):  {composite_simpson_13(f, a, b, 10):.10f}")
```

---

<br>

## 요약

| 법칙 | 공식 | 점 수 | 오차 (단일) | 오차 (복합) | 정확한 대상 |
|------|------|--------|------------|------------|------------|
| **사다리꼴(Trapezoidal)** | $(b-a)\,\dfrac{f(a)+f(b)}{2}$ | 2 | $O\bigl((b-a)^3\bigr)$ | $O(h^2)$ | 차수 $\leq 1$ |
| **심프슨의 1/3(Simpson's 1/3)** | $\dfrac{h}{3}(f_0 + 4f_1 + f_2)$ | 3 | $O\bigl((b-a)^5\bigr)$ | $O(h^4)$ | 차수 $\leq 3$ |
| **심프슨의 3/8(Simpson's 3/8)** | $\dfrac{3h}{8}(f_0 + 3f_1 + 3f_2 + f_3)$ | 4 | $O\bigl((b-a)^5\bigr)$ | $O(h^4)$ | 차수 $\leq 3$ |

**강의의 핵심 내용:**

1. 모든 뉴턴-코츠 공식은 **피적분함수를 다항식으로 대체**하고 그 다항식을 정확히 적분하는 방식으로 작동한다
2. **사다리꼴 법칙**은 1차 다항식을 사용하며, 복합 오차는 $O(h^2)$로 감소한다 --- $n$을 두 배로 하면 오차가 $4\times$ 줄어든다
3. **심프슨의 1/3 법칙**은 2차 다항식을 사용하지만 대칭성으로 인해 3차까지 정확하다; 복합 오차는 $O(h^4)$로 감소한다 --- $n$을 두 배로 하면 오차가 $16\times$ 줄어든다
4. **심프슨의 3/8 법칙**은 3차 다항식을 사용하며, 역시 3차까지 정확하고 $O(h^4)$ 복합 오차를 가진다; 단일 적용에서는 1/3 법칙보다 약간 더 정확하다
5. **복합 버전**은 구간을 분할하여 오차를 줄이며, $n$이 단일 법칙의 나누어짐 요건에 맞지 않을 때 서로 다른 법칙을 **결합**할 수 있다
