# 제22장 강의 — 상미분방정식의 초기값 문제

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 22

> **선수 지식**: [미적분학] 상미분방정식, 미분 (제1-21장).
>
> **학습 목표**:
> 1. 초기값 문제에 오일러 방법을 적용할 수 있다
> 2. 룽게-쿠타 방법(RK2, RK4)을 구현할 수 있다
> 3. ODE 풀이의 안정성과 정확도를 분석할 수 있다

---

<br>

## 목차

- [1. ODE와 초기값 문제 개요](#1-ode와-초기값-문제-개요)
  - [1.1 상미분방정식(Ordinary Differential Equations)](#11-상미분방정식ordinary-differential-equations)
  - [1.2 초기값 문제(IVP)](#12-초기값-문제ivp)
  - [1.3 연립 ODE](#13-연립-ode)
  - [1.4 적분 형태](#14-적분-형태)
- [2. 오일러 방법(Euler Method)](#2-오일러-방법euler-method)
  - [2.1 유도](#21-유도)
  - [2.2 기하학적 해석](#22-기하학적-해석)
- [3. 국소 및 전역 절단 오차](#3-국소-및-전역-절단-오차)
  - [3.1 전역 오차(Global Error)](#31-전역-오차global-error)
  - [3.2 오차에 대한 Big-O 표기법](#32-오차에-대한-big-o-표기법)
  - [3.3 오일러 방법의 테일러 급수 분석](#33-오일러-방법의-테일러-급수-분석)
  - [3.4 근사 국소 절단 오차](#34-근사-국소-절단-오차)
  - [3.5 국소 오차에서 전역 오차로](#35-국소-오차에서-전역-오차로)
- [4. 룽게-쿠타(RK) 방법](#4-룽게-쿠타rk-방법)
  - [4.1 핵심 아이디어](#41-핵심-아이디어)
  - [4.2 증분 함수(Increment Function)](#42-증분-함수increment-function)
  - [4.3 양적 RK 방법 — 일반 프레임워크](#43-양적-rk-방법--일반-프레임워크)
  - [4.4 부처 테이블(Butcher Tableau)](#44-부처-테이블butcher-tableau)
- [5. 2차 룽게-쿠타 방법](#5-2차-룽게-쿠타-방법)
  - [5.1 일반적인 2단계 RK 형태](#51-일반적인-2단계-rk-형태)
  - [5.2 차수 조건 유도](#52-차수-조건-유도)
  - [5.3 2차 정확도를 위한 차수 조건](#53-2차-정확도를-위한-차수-조건)
  - [5.4 호인 방법(Heun's Method) (b2 = 1/2)](#54-호인-방법heuns-method-b2--12)
  - [5.5 수정 오일러 / 중점 방법(Midpoint Method) (b2 = 1)](#55-수정-오일러--중점-방법midpoint-method-b2--1)
- [6. 고전적 4차 룽게-쿠타 방법(RK4)](#6-고전적-4차-룽게-쿠타-방법rk4)
  - [6.1 공식](#61-공식)
  - [6.2 RK4의 부처 테이블](#62-rk4의-부처-테이블)
  - [6.3 Python 구현](#63-python-구현)
- [7. 연립 ODE](#7-연립-ode)
  - [7.1 벡터 형태](#71-벡터-형태)
  - [7.2 고차 ODE의 연립 변환](#72-고차-ode의-연립-변환)
  - [7.3 연립 ODE를 위한 RK4 — Python 구현](#73-연립-ode를-위한-rk4--python-구현)
- [8. 수치적 안정성(Numerical Stability)](#8-수치적-안정성numerical-stability)
  - [8.1 모델 문제(Model Problem)](#81-모델-문제model-problem)
  - [8.2 정확해의 안정성](#82-정확해의-안정성)
  - [8.3 오일러 방법(전진 오일러)의 안정성](#83-오일러-방법전진-오일러의-안정성)
  - [8.4 전진 오일러의 안정 영역](#84-전진-오일러의-안정-영역)
  - [8.5 후진 오일러 방법의 안정성](#85-후진-오일러-방법의-안정성)
  - [8.6 사다리꼴 방법의 안정성](#86-사다리꼴-방법의-안정성)
  - [8.7 A-안정성과 L-안정성](#87-a-안정성과-l-안정성)
- [요약](#요약)

---

<br>

## 1. ODE와 초기값 문제 개요

### 1.1 상미분방정식(Ordinary Differential Equations)

$y = y(t)$를 독립변수 $t$의 미지함수라 하자. **상미분방정식(ODE)**은 $y$와 그 도함수 사이의 관계를 나타낸다:

$$\frac{dy}{dt} = f(t, y)$$

함수 $f(t, y)$는 알려져 있으며, $t$-$y$ 평면의 임의의 점 $(t, y)$에서 해의 **기울기**(경향)를 제공한다. 해 $y(t)$는 이 기울기장을 따라 **해 곡선**을 그린다.

> **[미적분]** 대수 방정식에서 미지수가 숫자인 것과 달리, 미분방정식에서 미지수는 *함수*이다. ODE $\frac{dy}{dt} = f(t,y)$는 임의의 시각 $t$에서 $y$의 변화율이 $t$와 $y$의 현재 값에 의해 결정된다는 것을 의미한다.

### 1.2 초기값 문제(IVP)

초기 조건 $y(t_0) = y_0$이 주어지면, $t > t_0$에 대해 $y(t)$를 구한다. 이것이 **초기값 문제(Initial Value Problem)**이다:

$$\begin{cases} \dfrac{dy}{dt} = f(t, y), & t > t_0 \\[6pt] y(t_0) = g \end{cases}$$

초기 조건은 ODE의 모든 해 곡선 중에서 유일한 해 곡선을 선택한다.

### 1.3 연립 ODE

연립된 두 ODE를 고려하자:

$$\frac{du}{dt} = p(t, u, v), \quad u(t_0) = u_0$$

$$\frac{dv}{dt} = q(t, u, v), \quad v(t_0) = v_0$$

**벡터 형태**로, $\mathbf{y} := (u, v)^T$ 및 $\mathbf{f} := (p, q)^T$로 정의하면:

$$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}), \quad \mathbf{y}(t_0) = (u_0, v_0)^T$$

이것이 **연립 ODE**이다. 단일 방정식에 대한 모든 방법(오일러, 룽게-쿠타 등)은 스칼라 연산을 벡터 연산으로 대체하면 연립계로 자연스럽게 확장된다.

### 1.4 적분 형태

ODE를 $t_0$에서 $t_1$까지 적분하면:

$$\int_{t_0}^{t_1} \frac{dy}{dt}\,dt = y(t_1) - y(t_0) = \int_{t_0}^{t_1} f(t, y(t))\,dt$$

따라서:

$$y(t_1) = y(t_0) + \int_{t_0}^{t_1} f(t, y(t))\,dt$$

핵심 질문은: $y(t)$를 닫힌 형태로 모를 때 **이 적분을 어떻게 수치적으로 근사할 수 있는가?**이다.

---

<br>

## 2. 오일러 방법(Euler Method)

### 2.1 유도

시간 구간 $[t_0, T]$를 균일한 스텝 크기를 갖는 $n$개의 소구간으로 분할한다:

$$h = \frac{T - t_0}{n}, \quad t_{n+1} = t_n + h$$

한 스텝에 대한 적분 형태로부터:

$$y(t_{n+1}) = y(t_n) + \int_{t_n}^{t_{n+1}} f(t, y(t))\,dt$$

가장 간단한 근사는 $f$를 $[t_n, t_{n+1}]$에서 상수로 취급하고 왼쪽 끝점에서 평가하는 것이다:

$$\int_{t_n}^{t_{n+1}} f(t, y(t))\,dt \approx f(t_n, y_n) \cdot \int_{t_n}^{t_{n+1}} dt = f(t_n, y_n) \cdot h$$

이로부터 **오일러 방법**(전진 오일러, Forward Euler)을 얻는다:

$$\boxed{y_{n+1} = y_n + h\,f_n, \quad f_n := f(t_n, y_n)}$$

### 2.2 기하학적 해석

각 시간 스텝에서 오일러 방법은 현재 점 $(t_n, y_n)$에서의 접선을 따라 해를 전진시킨다. 이 접선의 기울기는 $f(t_n, y_n)$이므로, 수평 거리 $h$만큼 이동하면 수직 변화량 $h \cdot f(t_n, y_n)$을 생성한다.

> **[미적분]** 오일러 방법은 수치 ODE 풀이법의 방대한 계열 중 가장 간단한 구성원이다. 이는 적분 $\int f\,dt$를 왼쪽 끝점 직사각형 법칙(left-endpoint rectangle rule)으로 근사하는 것에 해당한다. 더 좋은 구적법은 더 좋은 ODE 방법으로 이어진다(이것이 룽게-쿠타의 아이디어이다).

---

<br>

## 3. 국소 및 전역 절단 오차

### 3.1 전역 오차(Global Error)

스텝 $n$에서의 **전역 오차**는 참 해와 수치 근사 사이의 차이이다:

$$e_n = y(t_n) - y_n$$

초기 조건에 의해 $e_0 = y(t_0) - y_0 = 0$이다. 전역 오차는 $t_0$에서 $t_n$까지 모든 스텝에 걸친 국소 오차의 **누적**이다.

### 3.2 오차에 대한 Big-O 표기법

양 $z$가 $O(h^p)$라 함은 상수 $h_0, C > 0$가 존재하여 다음이 성립하는 것이다:

$$|z| \le C h^p \quad \forall\; 0 < h < h_0$$

이는 $z$가 $h \to 0$일 때 0으로 수렴하며, **수렴 속도**가 $p$임을 의미한다.

### 3.3 오일러 방법의 테일러 급수 분석

$h = t_{n+1} - t_n$으로 $(t_n, y_n)$ 주위에서 $y$의 테일러 전개를 상기하자:

$$y(t_n + h) = y(t_n) + y'(t_n)\,h + y''(t_n)\frac{h^2}{2} + \cdots + y^{(n)}(t_n)\frac{h^n}{n!} + R_n$$

여기서 $R_n = y^{(n+1)}(\xi)\frac{h^{n+1}}{(n+1)!}$이며, 어떤 $\xi \in [t_n, t_n + h]$이다.

$y'(t_n) = f(t_n, y_n)$이므로, $f$의 도함수를 사용하여 테일러 전개를 다시 쓸 수 있다:

$$y(t_n + h) = \underbrace{y_n + f(t_n, y_n)\,h}_{\text{오일러 근사}} + \underbrace{f'(t_n, y_n)\frac{h^2}{2} + \cdots + f^{(n)}(t_n, y_n)\frac{h^n}{n!} + O(h^{n+1})}_{\text{절단된 항 (국소 절단 오차)}}$$

오일러 방법의 **참 국소 절단 오차**는 $O(h^2)$이다 — 버려지는 항의 선행 항(leading term).

### 3.4 근사 국소 절단 오차

오일러 방법의 국소 절단 오차의 선행 항은:

$$E_a = f'(t_n, y_n)\frac{h^2}{2} = O(h^2)$$

이것이 각 스텝에서의 **근사 국소 절단 오차**이다.

### 3.5 국소 오차에서 전역 오차로

$[0, T]$에서 균일 스텝 크기 $h$를 사용하면, $n = \frac{T}{h}$ 스텝이 된다. 국소 오차를 합산하면:

$$\left|\sum_{i=1}^{n} E_{a,i}\right| = \left|(f'_0 + f'_1 + \cdots + f'_{n-1})\frac{h^2}{2}\right| \le \tilde{f}' \cdot \frac{h^2}{2} \cdot n$$

여기서 $\tilde{f}' = \max_{0 \le i \le n-1} |f'_i|$이다.

$n = \frac{T}{h}$를 대입하면:

$$\left|\sum E_{a,i}\right| \le \tilde{f}' \cdot \frac{h^2}{2} \cdot \frac{T}{h} = \frac{\tilde{f}'\,T}{2}\,h = O(h)$$

**핵심 결과: 오일러 방법은 국소 절단 오차 $O(h^2)$, 전역 절단 오차 $O(h)$를 갖는다.** 오일러 방법은 **1차 방법**이다.

> **[미적분]** 국소 오차가 전역 오차로 누적될 때 $h$의 차수가 하나 줄어드는데, 이는 스텝 수 $n = T/h$가 $h$가 감소함에 따라 증가하기 때문이다. 이것은 보편적인 패턴이다: 국소 오차가 $O(h^{p+1})$인 방법의 전역 오차는 $O(h^p)$이다.

---

<br>

## 4. 룽게-쿠타(RK) 방법

### 4.1 핵심 아이디어

룽게-쿠타 방법은 $f$의 고차 도함수를 명시적으로 계산하지 않으면서, 구간 $[t_n, t_{n+1}]$ 내의 여러 지점에서 **다수의 기울기 추정치를 조합하여 정확도를 향상**시킨다.

적분 형태에서 출발하여 $t = (t_{n+1} - t_n)z + t_n$ ($z \in [0, 1]$)으로 변수 치환을 수행하면:

$$y_{n+1} = y_n + \int_{t_n}^{t_{n+1}} f(t, y)\,dt = y_n + h \int_0^1 f(z, y(z))\,dz$$

### 4.2 증분 함수(Increment Function)

적분 $\int_0^1 f(z, y(z))\,dz$를 구적법으로 근사한다:

$$\int_0^1 f(z, y(z))\,dz \approx \phi$$

여기서 $\phi$는 **증분 함수**(평균 기울기)이다:

$$y_{n+1} = y_n + \phi \cdot h$$

증분 함수는 기울기 평가값의 가중합이다:

$$\phi = \sum_{i=1}^{s} b_i k_i = \sum_{i=1}^{s} b_i\,f(t_n + c_i h,\; y(t_n + c_i h))$$

이는 본질적으로 스텝에 걸친 $f$의 **구적법(quadrature)**이며, $b_i$는 가중치, $c_i$는 절점(node)이다.

### 4.3 양적 RK 방법 — 일반 프레임워크

일반적으로 $y(t_n + c_i h)$를 모르므로, **내부 단계(internal stage)** 값 $Y_i$로 근사한다:

**내부 단계** (순차적으로 계산):

$$Y_1 = y_n \quad \Rightarrow \quad k_1 = f(t_n, Y_1)$$

$$Y_2 = y_n + h\,a_{21}\,k_1 \quad \Rightarrow \quad k_2 = f(t_n + c_2 h, Y_2)$$

$$Y_3 = y_n + h(a_{31}\,k_1 + a_{32}\,k_2) \quad \Rightarrow \quad k_3 = f(t_n + c_3 h, Y_3)$$

$$\vdots$$

$$Y_s = y_n + h \sum_{j=1}^{s-1} a_{sj}\,k_j \quad \Rightarrow \quad k_s = f(t_n + c_s h, Y_s)$$

**스텝 완료:**

$$y_{n+1} = y_n + h \sum_{i=1}^{s} b_i\,k_i$$

여기서 $s$는 **단계 수(number of stages)**이다.

> **[미적분]** RK 방법의 핵심 통찰은, 미지의 중간 해 값에서 $f$를 평가하는 다루기 어려운 문제를 부트스트래핑(bootstrapping) 절차로 대체하는 것이다: 오일러와 유사한 예측으로 중간값을 추정하고, 그 점에서 $f$를 평가한 뒤, 그 기울기들을 결합한다.

### 4.4 부처 테이블(Butcher Tableau)

양적 RK 방법은 행렬 $A$(순 하삼각 행렬), 가중치 벡터 $\mathbf{b} = (b_1, b_2, \ldots, b_s)^T$, 절점 벡터 $\mathbf{c} = (c_1, c_2, \ldots, c_s)^T$에 의해 완전히 지정된다.

내부 단계에 대한 행렬 형태:

$$\begin{pmatrix} Y_1 \\ Y_2 \\ Y_3 \\ \vdots \\ Y_s \end{pmatrix} = \begin{pmatrix} y_n \\ y_n \\ y_n \\ \vdots \\ y_n \end{pmatrix} + h \begin{pmatrix} 0 & 0 & 0 & \cdots & 0 \\ a_{21} & 0 & 0 & \cdots & 0 \\ a_{31} & a_{32} & 0 & \cdots & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ a_{s1} & a_{s2} & a_{s3} & \cdots & 0 \end{pmatrix} \begin{pmatrix} k_1 \\ k_2 \\ k_3 \\ \vdots \\ k_s \end{pmatrix}$$

행렬 $A$가 순 하삼각 행렬인 것이 방법을 **양적(explicit)**으로 만드는 요소이다 — 각 단계를 이전에 계산된 단계들로부터 계산할 수 있다.

**부처 테이블(Butcher tableau)**은 표준적인 간결한 표기법이다:

$$\begin{array}{c|c} \mathbf{c} & A \\ \hline & \mathbf{b}^T \end{array}$$

---

<br>

## 5. 2차 룽게-쿠타 방법

### 5.1 일반적인 2단계 RK 형태

2단계($s = 2$) 양적 RK 방법은 다음의 형태를 갖는다:

$$y_{n+1} = y_n + h(b_1 k_1 + b_2 k_2)$$

여기서:

$$k_1 = f(t_n, y_n)$$

$$k_2 = f(t_n + ch,\; y_n + h\,a\,k_1)$$

결정해야 할 매개변수는 $b_1, b_2, c, a$ (4개의 미지수)이다.

### 5.2 차수 조건 유도

**참 해의 테일러 전개** ($t$에 대한 전미분 사용):

$$y(t_n + h) = y_n + h\,f + \frac{h^2}{2}(f_t + f_y f) + O(h^3)$$

여기서 모든 $f, f_t, f_y$는 $(t_n, y_n)$에서 평가된다.

> **[미적분]** 2차 도함수 $y'' = \frac{d}{dt}f(t, y) = f_t + f_y \cdot \frac{dy}{dt} = f_t + f_y f$는 전미분에 대한 **연쇄 법칙(chain rule)**을 사용한다. 이는 $f$가 $t$에 직접적으로도, $y(t)$를 통해 간접적으로도 의존하기 때문이다.

**RK2 공식의 테일러 전개:**

$k_2 = f(t_n + ch, y_n + a\,f(t_n,y_n)\,h)$를 $\Delta t = ch$, $\Delta y = a\,f(t_n,y_n)\,h$로 다변수 테일러 전개하면:

$$f(t_n + \Delta t, y_n + \Delta y) = f + f_t \cdot ch + f_y \cdot a\,f\,h + O(h^2)$$

RK2 공식에 다시 대입하면:

$$y_{n+1} = y_n + h(b_1 + b_2)f + \frac{h^2}{2} \cdot 2b_2(c\,f_t + a\,f_y f) + O(h^3)$$

### 5.3 2차 정확도를 위한 차수 조건

$O(h^2)$ 항까지 테일러 전개와 매칭하면 **4개의 미지수에 대한 3개의 방정식** (자유도 1개)을 얻는다:

$$b_1 + b_2 = 1$$

$$b_2 c = \frac{1}{2}$$

$$b_2 a = \frac{1}{2}$$

참고: 마지막 두 조건은 $c = a$를 의미한다.

### 5.4 호인 방법(Heun's Method) (b2 = 1/2)

$b_2 = 1/2$를 선택하면:

$$b_1 = \frac{1}{2}, \quad a = 1, \quad c = 1$$

$$Y_1 = y_n \quad \Rightarrow \quad k_1 = f(t_n, y_n)$$

$$Y_2 = y_n + h\,k_1 \quad \Rightarrow \quad k_2 = f(t_n + h, Y_2)$$

$$\boxed{y_{n+1} = y_n + h\left(\frac{k_1 + k_2}{2}\right)}$$

이것이 **호인 방법(Heun's method)**(개선된 오일러 또는 사다리꼴 예측-교정법이라고도 함)이다. 부처 테이블은:

$$\begin{array}{c|cc} 0 & & \\ 1 & 1 & \\ \hline & 1/2 & 1/2 \end{array}$$

> **[미적분]** 호인 방법은 먼저 전체 오일러 스텝으로 $t_{n+1}$에서 $Y_2$를 예측하고, 그 점에서의 기울기를 평가한 뒤, 양 끝점에서의 기울기를 평균한다. 이는 적분에 대한 사다리꼴 공식과 유사하다.

### 5.5 수정 오일러 / 중점 방법(Midpoint Method) (b2 = 1)

$b_2 = 1$을 선택하면:

$$b_1 = 0, \quad a = \frac{1}{2}, \quad c = \frac{1}{2}$$

$$k_1 = f(t_n, y_n)$$

$$k_2 = f\!\left(t_n + \frac{h}{2},\; y_n + \frac{h}{2}\,k_1\right)$$

$$\boxed{y_{n+1} = y_n + h\,k_2}$$

이것이 **수정 오일러(중점) 방법**이다. 부처 테이블은:

$$\begin{array}{c|cc} 0 & & \\ 1/2 & 1/2 & \\ \hline & 0 & 1 \end{array}$$

> **[미적분]** 중점 방법은 오일러 반스텝으로 구간의 중점에서의 기울기를 추정한 뒤, 그 중점 기울기를 사용하여 전체 스텝을 전진시킨다. 이는 적분에 대한 중점 법칙(midpoint rule)과 유사하다.

---

<br>

## 6. 고전적 4차 룽게-쿠타 방법(RK4)

### 6.1 공식

고전적 RK4는 국소 절단 오차 $O(h^5)$, 전역 오차 $O(h^4)$인 4단계 방법이다:

$$k_1 = f(t_n,\; y_n)$$

$$k_2 = f\!\left(t_n + \frac{h}{2},\; y_n + \frac{h}{2}\,k_1\right)$$

$$k_3 = f\!\left(t_n + \frac{h}{2},\; y_n + \frac{h}{2}\,k_2\right)$$

$$k_4 = f(t_n + h,\; y_n + h\,k_3)$$

$$\boxed{y_{n+1} = y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)}$$

가중치 $\frac{1}{6}, \frac{2}{6}, \frac{2}{6}, \frac{1}{6}$는 네 기울기 추정치에 적용된 **심프슨의 1/3 법칙(Simpson's 1/3 rule)**에 해당한다.

> **[미적분]** RK4는 ODE 풀이의 핵심 방법이다. $f$의 어떤 도함수도 필요 없이 스텝당 4번의 함수 평가만으로 4차 정확도를 달성한다. 많은 실용적 문제에서 정확도와 계산 비용 사이의 훌륭한 균형을 제공한다.

### 6.2 RK4의 부처 테이블

$$\begin{array}{c|cccc} 0 & & & & \\ 1/2 & 1/2 & & & \\ 1/2 & 0 & 1/2 & & \\ 1 & 0 & 0 & 1 & \\ \hline & 1/6 & 1/3 & 1/3 & 1/6 \end{array}$$

### 6.3 Python 구현

```python
import numpy as np

def rk4(f, t0, y0, T, h):
    """
    Classical 4th-order Runge-Kutta method.

    Parameters
    ----------
    f  : callable, f(t, y) -> dy/dt
    t0 : float, initial time
    y0 : float, initial value
    T  : float, final time
    h  : float, step size

    Returns
    -------
    t_arr : ndarray, time values
    y_arr : ndarray, solution values
    """
    n = int((T - t0) / h)
    t_arr = np.linspace(t0, T, n + 1)
    y_arr = np.zeros(n + 1)
    y_arr[0] = y0

    for i in range(n):
        t = t_arr[i]
        y = y_arr[i]
        k1 = f(t, y)
        k2 = f(t + h/2, y + h/2 * k1)
        k3 = f(t + h/2, y + h/2 * k2)
        k4 = f(t + h, y + h * k3)
        y_arr[i+1] = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)

    return t_arr, y_arr
```

---

<br>

## 7. 연립 ODE

### 7.1 벡터 형태

$m$개의 1차 ODE 연립계는 벡터 형태로 다음과 같이 쓴다:

$$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}), \quad \mathbf{y}(t_0) = \mathbf{y}_0$$

여기서 $\mathbf{y} = (y_1, y_2, \ldots, y_m)^T$이고 $\mathbf{f} = (f_1, f_2, \ldots, f_m)^T$이다.

모든 오일러 및 룽게-쿠타 공식은 **성분별로(component-by-component)** 적용된다 — 단순히 스칼라 $y$, $f$, $k$를 벡터 $\mathbf{y}$, $\mathbf{f}$, $\mathbf{k}$로 대체하면 된다.

### 7.2 고차 ODE의 연립 변환

$m$차 ODE는 $m$개의 1차 ODE 연립계로 변환할 수 있다. 예를 들어, 2차 ODE:

$$y'' = g(t, y, y')$$

$u = y$, $v = y'$를 도입하여 다시 쓰면:

$$\frac{du}{dt} = v, \quad \frac{dv}{dt} = g(t, u, v)$$

초기 조건은 $u(t_0) = y(t_0)$, $v(t_0) = y'(t_0)$이다.

> **[미적분]** 이 축소 기법은 기본적이다: 1차 연립계에 대한 *모든* 수치적 방법으로 먼저 변환한 뒤 고차 ODE를 풀 수 있다. 예를 들어, 뉴턴의 제2법칙 $m\,x'' = F(t, x, x')$는 위치 $x$와 속도 $v = x'$에 대한 1차 연립계가 된다.

### 7.3 연립 ODE를 위한 RK4 — Python 구현

```python
import numpy as np

def rk4_system(f, t0, y0, T, h):
    """
    RK4 for a system of ODEs.

    Parameters
    ----------
    f  : callable, f(t, y) -> ndarray of shape (m,)
    t0 : float, initial time
    y0 : ndarray of shape (m,), initial conditions
    T  : float, final time
    h  : float, step size

    Returns
    -------
    t_arr : ndarray of shape (n+1,)
    y_arr : ndarray of shape (n+1, m)
    """
    n = int((T - t0) / h)
    m = len(y0)
    t_arr = np.linspace(t0, T, n + 1)
    y_arr = np.zeros((n + 1, m))
    y_arr[0] = y0

    for i in range(n):
        t = t_arr[i]
        y = y_arr[i]
        k1 = f(t, y)
        k2 = f(t + h/2, y + h/2 * k1)
        k3 = f(t + h/2, y + h/2 * k2)
        k4 = f(t + h, y + h * k3)
        y_arr[i+1] = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)

    return t_arr, y_arr


# Example: Solve y'' + y = 0  (simple harmonic oscillator)
# Convert to system: u' = v, v' = -u
# with u(0) = 1, v(0) = 0  => exact solution: u(t) = cos(t)

def harmonic(t, y):
    u, v = y
    return np.array([v, -u])

t, sol = rk4_system(harmonic, 0, np.array([1.0, 0.0]), 10, 0.01)
# sol[:, 0] is u(t) ≈ cos(t)
# sol[:, 1] is v(t) ≈ -sin(t)
```

---

<br>

## 8. 수치적 안정성(Numerical Stability)

### 8.1 모델 문제(Model Problem)

수치 해가 오차가 지수적으로 증가하면 **불안정(unstable)**하다고 한다. 안정성을 분석하기 위해, **모델(시험) 문제**를 고려한다:

$$\frac{dy}{dt} = \lambda y, \quad y(0) = y_0$$

정확해는:

$$y(t) = y_0 e^{\lambda t}$$

### 8.2 정확해의 안정성

정확해의 거동은 $\lambda$에 따라 달라진다:

| $\lambda$의 값 | $y(t) = y_0 e^{\lambda t}$의 거동 |
|:---|:---|
| $\lambda > 0$ | 지수적 증가 (불안정) |
| $\lambda = 0$ | 상수 (중립 안정) |
| $\lambda < 0$ | 지수적 감소 (안정) |

$\lambda < 0$일 때 해는 점근적으로 0에 접근한다. **안정적인** 수치 방법은 $\lambda < 0$에서 이 감소 거동을 재현해야 한다.

### 8.3 오일러 방법(전진 오일러)의 안정성

오일러 방법을 $y' = \lambda y$에 적용하면:

$$y_{n+1} = y_n + h \lambda y_n = y_n(1 + \lambda h)$$

인자 $(1 + \lambda h)$는 **증폭 인자(amplification factor)** (또는 **안정성 함수** $\phi(\lambda h)$)이다. 귀납법에 의해:

$$y_1 = y_0(1 + \lambda h)$$

$$y_2 = y_0(1 + \lambda h)^2$$

$$\vdots$$

$$y_n = y_0(1 + \lambda h)^n$$

안정성(유계 해)을 위해 다음이 필요하다:

$$|1 + \lambda h| \le 1$$

### 8.4 전진 오일러의 안정 영역

$z = \lambda h \in \mathbb{C}$로 놓자. 안정 조건은:

$$|1 + z| \le 1$$

이것은 복소 $z$-평면에서 **$(-1, 0)$을 중심으로 하는 반지름 1의 원판**이다.

**실수축** 위에서, 조건은 다음으로 축소된다:

$$-1 \le 1 + z \le 1 \quad \Rightarrow \quad -2 \le z \le 0$$

$z = \lambda h$이고 $h > 0$이므로, 다음이 요구된다:

$$\lambda \le 0 \quad \text{그리고} \quad h \le \frac{-2}{\lambda} = \frac{2}{|\lambda|}$$

**전진 오일러 방법은 조건부 안정(conditionally stable)**이다 — 안정성을 보장하기 위해 스텝 크기 $h$를 충분히 작게 선택해야 한다. 강하게 음인 $\lambda$ (강성 문제)에 대해, $h$는 매우 작아야 한다.

### 8.5 후진 오일러 방법의 안정성

**후진 오일러 방법(Backward Euler method)**(음적 방법)은 *다음* 시간 스텝에서의 기울기를 평가한다:

$$x_{n+1} = x_n + h\,f(x_{n+1})$$

$\frac{dx}{dt} = \lambda x$에 적용하면:

$$x_{n+1} = x_n + \lambda h\,x_{n+1}$$

$$(1 - \lambda h)\,x_{n+1} = x_n$$

$$x_{n+1} = \frac{1}{1 - \lambda h}\,x_n$$

귀납법에 의해:

$$x_n = (1 - \lambda h)^{-n}\,x_0$$

안정성 함수는 $\phi(z) = \frac{1}{1 - z}$이다. 안정성을 위해:

$$\frac{1}{|1 - z|} < 1 \quad \Rightarrow \quad |1 - z| > 1$$

안정 영역은 복소 $z$-평면에서 $(1, 0)$을 중심으로 하는 반지름 1 원판의 **외부**이다. 이 영역은 **전체 왼쪽 반평면** $\text{Re}(z) < 0$을 포함한다.

임의의 $\lambda < 0$에 대해, 스텝 크기 $h$를 **임의로 크게** 선택해도 방법은 안정적이다. 이는 후진 오일러 방법이 안정성을 위한 스텝 크기 제한을 요구하지 않음을 의미한다 — 오직 정확도 고려만이 $h$를 제한한다.

> **[미적분]** 이 무조건적 안정성의 대가는 후진 오일러가 **음적(implicit)**이라는 것이다: 각 스텝에서 $x_{n+1}$에 대한 방정식(비선형일 수 있음)을 풀어야 한다. 선형 시험 문제에서는 단순한 대수적 풀이로 귀결되지만, 일반적인 비선형 $f$에 대해서는 뉴턴 방법이나 유사한 기법이 필요하다.

### 8.6 사다리꼴 방법의 안정성

**사다리꼴 방법(Trapezoidal method)**은 양 끝점에서의 평균 기울기를 사용한다:

$$x_{n+1} = x_n + \frac{h}{2}(f_n + f_{n+1})$$

$\frac{dx}{dt} = \lambda x$에 적용하면:

$$x_{n+1} = x_n + \frac{\lambda h}{2}(x_n + x_{n+1})$$

$$\left(1 - \frac{\lambda h}{2}\right)x_{n+1} = \left(1 + \frac{\lambda h}{2}\right)x_n$$

귀납법에 의해:

$$x_n = \left(\frac{1 + \frac{\lambda h}{2}}{1 - \frac{\lambda h}{2}}\right)^n x_0$$

안정성 함수는 $\phi(z) = \frac{1 + z/2}{1 - z/2}$이다. 안정성을 위해:

$$\left|\frac{1 + z/2}{1 - z/2}\right| < 1 \quad \Rightarrow \quad |1 + z/2| < |1 - z/2| \quad \Rightarrow \quad |2 + z| < |2 - z|$$

이 조건은 **전체 왼쪽 반평면** $\text{Re}(z) < 0$에서 만족되므로, 사다리꼴 방법도 스텝 크기에 관계없이 임의의 $\lambda < 0$에 대해 안정적이다.

### 8.7 A-안정성과 L-안정성

**A-안정성(A-stability):** 수치 방법의 안정 영역이 **전체 왼쪽 반평면** $\{z \in \mathbb{C} : \text{Re}(z) \le 0\}$을 포함하면 **A-안정**이라 한다.

- 전진 오일러: A-안정 **아님** (유한한 안정 영역)
- 후진 오일러: **A-안정**
- 사다리꼴 방법: **A-안정**

**L-안정성(L-stability):** 수치 방법이 **L-안정**이려면:
1. 방법이 A-안정이고, **그리고**
2. 안정성 함수 $\phi(z) \to 0$ ($z \to \infty$일 때)

L-안정성은 A-안정성보다 **더 엄격한** 조건이다. 이는 **점근적 거동**에 관한 것이다 — L-안정 방법은 매우 강성인 성분의 빠른 감소를 올바르게 포착한다.

- 후진 오일러: $\phi(z) = \frac{1}{1-z} \to 0$ ($z \to -\infty$일 때) — **L-안정**
- 사다리꼴 방법: $\phi(z) = \frac{1 + z/2}{1 - z/2} \to -1$ ($z \to -\infty$일 때) — A-안정이지만 L-안정은 **아님**

> **[미적분]** 강성 문제($|\lambda|$가 매우 큰 경우)에서는 L-안정성이 바람직한데, 이는 빠르게 감소하는 성분의 수치 해도 빠르게 감소하도록 보장하기 때문이다. 사다리꼴 방법은 A-안정이지만, 매우 강성인 문제에서 안정성 함수가 큰 $|z|$에 대해 $0$이 아닌 $-1$에 접근하기 때문에 가짜 진동(spurious oscillation)을 보일 수 있다.

---

<br>

## 요약

| 주제 | 핵심 공식 / 결과 |
|:------|:---------------------|
| **초기값 문제(IVP)** | $\frac{dy}{dt} = f(t,y)$, $y(t_0) = y_0$ |
| **오일러 방법** | $y_{n+1} = y_n + h\,f(t_n, y_n)$ |
| **오일러 국소 오차** | $O(h^2)$ |
| **오일러 전역 오차** | $O(h)$ — 1차 방법 |
| **RK 일반 형태** | $y_{n+1} = y_n + h \sum b_i k_i$ |
| **RK2 차수 조건** | $b_1 + b_2 = 1$, $b_2 c = 1/2$, $b_2 a = 1/2$ |
| **호인 방법** | $y_{n+1} = y_n + \frac{h}{2}(k_1 + k_2)$, $k_2 = f(t_n+h, y_n+hk_1)$ |
| **중점 방법** | $y_{n+1} = y_n + h\,k_2$, $k_2 = f(t_n + h/2, y_n + \frac{h}{2}k_1)$ |
| **RK4** | $y_{n+1} = y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$ — 4차 |
| **연립 ODE** | 스칼라를 벡터로 대체; 동일한 공식 적용 |
| **시험 방정식** | $y' = \lambda y$이고 $y(t) = y_0 e^{\lambda t}$ |
| **전진 오일러 안정성** | $\|1 + \lambda h\| \le 1$ — **조건부 안정** |
| **후진 오일러 안정성** | $\|1 - \lambda h\|^{-1} < 1$ — **무조건 안정** (A-안정, L-안정) |
| **사다리꼴 안정성** | $\left\|\frac{1+z/2}{1-z/2}\right\| < 1$ — **A-안정이나 L-안정은 아님** |
| **A-안정** | 안정 영역이 전체 왼쪽 반평면을 포함 |
| **L-안정** | A-안정이고 $\phi(z) \to 0$ ($z \to \infty$일 때) |
