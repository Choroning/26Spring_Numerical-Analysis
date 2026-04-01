# 제16장 강의 — 푸리에 해석(Fourier Analysis)

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 16

> **선수 지식**: [미적분학] 삼각함수, 적분 (제1-15장).
>
> **학습 목표**:
> 1. 신호 분석을 위해 이산 푸리에 변환(DFT)을 적용할 수 있다
> 2. 고속 푸리에 변환(FFT) 알고리즘을 구현할 수 있다
> 3. 주파수 스펙트럼과 파워 스펙트럼 밀도를 해석할 수 있다

---

<br>

## 목차

- [1. 개요 및 목표](#1-개요-및-목표)
- [2. 우함수와 기함수(Even and Odd Functions)](#2-우함수와-기함수even-and-odd-functions)
  - [2.1 우함수(Even Functions)](#21-우함수even-functions)
  - [2.2 기함수(Odd Functions)](#22-기함수odd-functions)
  - [2.3 임의 함수의 분해](#23-임의-함수의-분해)
- [3. 주기 함수(Periodic Functions)](#3-주기-함수periodic-functions)
  - [3.1 정의](#31-정의)
  - [3.2 예제: 구형파(Square Wave)](#32-예제-구형파square-wave)
- [4. 푸리에 급수(주기 $2\pi$)](#4-푸리에-급수주기-2pi)
  - [4.1 푸리에의 주장](#41-푸리에의-주장)
  - [4.2 푸리에 계수 — 유도](#42-푸리에-계수--유도)
  - [4.3 직교성 관계(Orthogonality Relations)](#43-직교성-관계orthogonality-relations)
  - [4.4 계수 공식](#44-계수-공식)
- [5. 구형파의 푸리에 급수](#5-구형파의-푸리에-급수)
- [6. $2\pi$ 이외의 주기](#6-2pi-이외의-주기)
- [7. 푸리에 적분(비주기 함수)](#7-푸리에-적분비주기-함수)
- [8. 푸리에 변환(Fourier Transform)](#8-푸리에-변환fourier-transform)
  - [8.1 정의](#81-정의)
  - [8.2 푸리에 변환의 성질](#82-푸리에-변환의-성질)
- [9. 정현파 함수를 이용한 곡선 맞춤](#9-정현파-함수를-이용한-곡선-맞춤)
  - [9.1 정현파 함수의 매개변수](#91-정현파-함수의-매개변수)
  - [9.2 위상각(Phase Angle)](#92-위상각phase-angle)
  - [9.3 대안적 형태](#93-대안적-형태)
- [10. 정현파의 최소제곱 적합](#10-정현파의-최소제곱-적합)
  - [10.1 문제 설정](#101-문제-설정)
  - [10.2 정규 방정식](#102-정규-방정식)
  - [10.3 등간격 점에서의 간소화](#103-등간격-점에서의-간소화)
  - [10.4 닫힌 형태의 해](#104-닫힌-형태의-해)
- [11. 연속 푸리에 급수(일반 모델)](#11-연속-푸리에-급수일반-모델)
  - [11.1 일반 푸리에 급수](#111-일반-푸리에-급수)
  - [11.2 오일러 공식과 복소 형태](#112-오일러-공식과-복소-형태)
- [12. 주파수 영역과 시간 영역](#12-주파수-영역과-시간-영역)
- [13. 푸리에 적분과 푸리에 변환(재방문)](#13-푸리에-적분과-푸리에-변환재방문)
- [14. 이산 푸리에 변환(DFT)](#14-이산-푸리에-변환dft)
- [15. 나이퀴스트 주파수(Nyquist Frequency)](#15-나이퀴스트-주파수nyquist-frequency)
- [16. 고속 푸리에 변환(FFT)](#16-고속-푸리에-변환fft)
- [17. 파워 스펙트럼(Power Spectrum)](#17-파워-스펙트럼power-spectrum)
- [18. 요약 표](#18-요약-표)

---

<br>

## 1. 개요 및 목표

**N17 (푸리에 급수)에서:**

- 푸리에 급수와 푸리에 적분 이해

**N18 (푸리에 해석)에서:**

- 정현파(사인파 및 코사인파) 이해
- 최소제곱 회귀를 사용하여 데이터에 정현파 적합
- 오일러 공식(Euler's Formula) 이해
- 주파수 영역에서의 신호 분석
- 이산 푸리에 변환(Discrete Fourier Transform)
- 나이퀴스트 주파수(Nyquist Frequency)
- 고속 푸리에 변환(Fast Fourier Transform)
- 파워 스펙트럼(Power Spectrum)

> **[맥락]** 푸리에 해석은 수치 해석에서 가장 강력한 도구 중 하나이다. 임의의 주기(또는 비주기) 함수를 단순한 정현파 성분의 합으로 분해할 수 있게 해준다. 신호 처리, 이미지 압축, 미분 방정식 풀이 등에 응용된다.

---

<br>

## 2. 우함수와 기함수(Even and Odd Functions)

### 2.1 우함수(Even Functions)

함수 $f(x)$가 **우함수(Even Function)** 인 조건:

$$f(-x) = f(x)$$

$f(x)$의 그래프는 $x = 0$(y축)에 대해 **대칭** 이다.

**예:** $\cos(x)$, $x^2$, $|x|$

### 2.2 기함수(Odd Functions)

함수 $f(x)$가 **기함수(Odd Function)** 인 조건:

$$f(-x) = -f(x)$$

$f(x)$의 그래프는 $x = 0$(원점)에 대해 **반대칭** 이다.

**예:** $\sin(x)$, $x^3$, $x$

### 2.3 임의 함수의 분해

임의의 함수는 우함수와 기함수의 합으로 분해할 수 있다:

$$f(x) = \frac{1}{2}f(x) + \frac{1}{2}f(x)$$

$$= \underbrace{\frac{1}{2}\bigl[f(x) + f(-x)\bigr]}_{\text{우함수 부분}} + \underbrace{\frac{1}{2}\bigl[f(x) - f(-x)\bigr]}_{\text{기함수 부분}}$$

> **[핵심 아이디어]** 이 분해는 푸리에 급수를 이해하는 데 근본적이다. 코사인 항은 함수의 우함수 부분을, 사인 항은 기함수 부분을 포착한다.

---

<br>

## 3. 주기 함수(Periodic Functions)

### 3.1 정의

함수 $f(x)$가 주기 $T$를 갖는 **주기 함수(Periodic Function)** 인 조건:

$$f(x) = f(x + T) \quad \forall x$$

항상 **가능한 가장 작은** 양의 $T$를 선택한다.

**예:** $\sin(x)$는 주기 $T = 2\pi$를 가지며, 다음이 성립한다:

$$\sin(x + 2\pi) = \sin(x)\cos(2\pi) + \cos(x)\sin(2\pi) = \sin(x)$$

$\sin(x) = \sin(x + 2\pi) = \sin(x + 4\pi) = \cdots$이지만, 기본 주기로 $T = 2\pi$를 선택한다.

### 3.2 예제: 구형파(Square Wave)

주기 $2\pi$인 구형파를 고려하자:

$$f(x) = \begin{cases} 0, & -\pi < x \le 0 \\ 1, & 0 < x \le \pi \end{cases}$$

이 함수는 주기 $2\pi$로 반복되며, 푸리에 급수를 시연하는 데 사용되는 고전적인 예이다.

---

<br>

## 4. 푸리에 급수(주기 $2\pi$)

### 4.1 푸리에의 주장

**푸리에(Fourier, 1768--1830)** 는 $[-\pi, \pi]$에서 정의된 임의의 함수가 다음 형태로 표현될 수 있다고 주장했다:

$$f(x) = \sum_{n=0}^{\infty} (a_n \cos nx + b_n \sin nx)$$

관례적인 $\frac{a_0}{2}$ 표기법으로 다시 쓰면:

$$f(x) = \frac{a_0}{2} + \sum_{n=1}^{\infty} (a_n \cos nx + b_n \sin nx)$$

> **[왜 $a_0/2$인가?]** $n = 0$일 때, $a_0 \cos(0) + b_0 \sin(0) = a_0 \cdot 1 + 0 = a_0$이다. $a_0$의 공식이 $n \ge 1$인 $a_n$과 같은 형태를 갖도록 편의상 $\frac{1}{2}$ 인자를 도입한다.

구조에 주목하자:

$$f(x) = \underbrace{\left(\frac{a_0}{2} + \sum_{n=1}^{\infty} a_n \cos nx \right)}_{\text{우함수 부분}} + \underbrace{\left(\sum_{n=1}^{\infty} b_n \sin nx \right)}_{\text{기함수 부분}}$$

### 4.2 푸리에 계수 — 유도

**$a_0$ 구하기:** 양변을 $[-\pi, \pi]$에서 적분한다:

$$\int_{-\pi}^{\pi} f(x)\,dx = \int_{-\pi}^{\pi} \frac{a_0}{2}\,dx + \sum_{n=1}^{\infty} a_n \underbrace{\int_{-\pi}^{\pi} \cos nx\,dx}_{= 0} + \sum_{n=1}^{\infty} b_n \underbrace{\int_{-\pi}^{\pi} \sin nx\,dx}_{= 0}$$

$$= \frac{a_0}{2} \cdot 2\pi = a_0 \pi$$

$$\therefore \quad a_0 = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x)\,dx$$

**$a_n$ 구하기:** 양변에 $\cos mx$를 곱하고 적분한다:

$$\int_{-\pi}^{\pi} f(x)\cos mx\,dx = \underbrace{\int_{-\pi}^{\pi} \frac{a_0}{2}\cos mx\,dx}_{= 0} + \sum_{n=1}^{\infty} a_n \int_{-\pi}^{\pi} \cos nx \cos mx\,dx + \sum_{n=1}^{\infty} b_n \underbrace{\int_{-\pi}^{\pi} \sin nx \cos mx\,dx}_{= 0}$$

**$b_n$ 구하기:** 양변에 $\sin mx$를 곱하고 적분한다:

$$\int_{-\pi}^{\pi} f(x)\sin mx\,dx = \underbrace{\int_{-\pi}^{\pi} \frac{a_0}{2}\sin mx\,dx}_{= 0} + \sum_{n=1}^{\infty} a_n \underbrace{\int_{-\pi}^{\pi} \cos nx \sin mx\,dx}_{= 0} + \sum_{n=1}^{\infty} b_n \int_{-\pi}^{\pi} \sin nx \sin mx\,dx$$

### 4.3 직교성 관계(Orthogonality Relations)

계수 유도의 핵심은 삼각 함수의 **직교성(Orthogonality)** 에 의존한다.

**곱을 합으로 바꾸는 항등식(Product-to-Sum Identities):**

$$\cos a \cos b = \frac{1}{2}[\cos(a+b) + \cos(a-b)]$$

$$\sin a \cos b = \frac{1}{2}[\sin(a+b) + \sin(a-b)]$$

$$\sin a \sin b = -\frac{1}{2}[\cos(a+b) - \cos(a-b)]$$

**$[-\pi, \pi]$에서의 직교성 결과:**

$$\int_{-\pi}^{\pi} \cos nx \cos mx\,dx = \begin{cases} 0 & \text{if } n \neq m \\ 2\pi & \text{if } n = m = 0 \\ \pi & \text{if } n = m \ge 1 \end{cases}$$

$$\int_{-\pi}^{\pi} \sin nx \cos mx\,dx = 0 \quad \text{for all } n, m$$

$$\int_{-\pi}^{\pi} \sin nx \sin mx\,dx = \begin{cases} 0 & \text{if } n \neq m \\ \pi & \text{if } n = m \ge 1 \end{cases}$$

### 4.4 계수 공식

직교성 관계를 사용하면:

$$\boxed{a_n = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x) \cos nx\,dx, \quad n = 0, 1, 2, \ldots}$$

$$\boxed{b_n = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x) \sin nx\,dx, \quad n = 1, 2, 3, \ldots}$$

---

<br>

## 5. 구형파의 푸리에 급수

구형파 $f(x) = \begin{cases} 0, & -\pi < x \le 0 \\ 1, & 0 < x \le \pi \end{cases}$에 대해:

**$a_0$ 계산:**

$$a_0 = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x)\,dx = \frac{1}{\pi} \cdot \pi = 1$$

**$a_n$ 계산:**

$$a_n = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x)\cos nx\,dx = \frac{1}{\pi} \int_{0}^{\pi} \cos nx\,dx = \frac{1}{\pi} \left[\frac{1}{n}\sin nx\right]_0^{\pi} = 0$$

**$b_n$ 계산:**

$$b_n = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x)\sin nx\,dx = \frac{1}{\pi} \int_0^{\pi} \sin nx\,dx = \frac{1}{\pi}\left[-\frac{1}{n}\cos nx\right]_0^{\pi}$$

$$= -\frac{1}{n\pi}(\cos n\pi - 1) = \begin{cases} \frac{2}{n\pi} & \text{for odd } n \\ 0 & \text{for even } n \end{cases}$$

**결과:**

$$f(x) = \frac{1}{2} + \sum_{\substack{n=1,3,5,\ldots}} \frac{2}{n\pi} \sin nx = \frac{1}{2} + \frac{2}{\pi} \sum_{n=1}^{\infty} \frac{\sin(2n-1)x}{2n-1}$$

$n \to \infty$일 때 $\frac{1}{n} \to 0$이므로, 푸리에 급수는 **수렴** 한다.

---

<br>

## 6. $2\pi$ 이외의 주기

주기 $2\ell$($2\pi$ 대신)을 갖는 함수에 대해, 다음 변수 치환을 사용한다:

$$x' = \pi \cdot \frac{x}{\ell}$$

이것은 구간 $[-\ell, \ell]$을 $[-\pi, \pi]$로 매핑한다. 푸리에 급수는 다음과 같이 된다:

$$\boxed{f(x) = \frac{a_0}{2} + \sum_{n=1}^{\infty} a_n \cos\frac{n\pi}{\ell}x + \sum_{n=1}^{\infty} b_n \sin\frac{n\pi}{\ell}x}$$

계수는:

$$a_n = \frac{1}{\ell} \int_{-\ell}^{\ell} f(x) \cos\frac{n\pi}{\ell}x\,dx, \quad n = 0, 1, 2, \ldots$$

$$b_n = \frac{1}{\ell} \int_{-\ell}^{\ell} f(x) \sin\frac{n\pi}{\ell}x\,dx, \quad n = 1, 2, 3, \ldots$$

---

<br>

## 7. 푸리에 적분(비주기 함수)

**질문:** $f(x) = e^{-x^2}$를 푸리에 급수로 전개할 수 있는가?

**답:** 아니다 — $f(x)$는 주기 함수가 아니다(NOT periodic).

**해결:** $\ell \to \infty$로 취한다. $\ell$이 증가하면, 주파수 스펙트럼 $\frac{n\pi}{\ell}$은 점점 더 조밀해져서 **연속** 스펙트럼에 수렴한다.

| $\ell$ | 주파수 $\frac{n\pi}{\ell}$ |
|:---:|:---|
| $\pi$ | $1, 2, 3, 4, \ldots$ |
| $2\pi$ | $\frac{1}{2}, 1, \frac{3}{2}, 2, \ldots$ |
| $10\pi$ | $0.1, 0.2, 0.3, 0.4, \ldots$ |

극한에서, 이산 합은 연속 적분이 된다:

$$f(x) = \frac{a_0}{2} + \sum_{n=1}^{\infty} (a_n \cos nx + b_n \sin nx)$$

$$\Downarrow$$

$$\boxed{f(x) = \int_0^{\infty} \bigl[a(\omega)\cos\omega x + b(\omega)\sin\omega x\bigr]\,d\omega}$$

이것이 **푸리에 적분(Fourier Integral)** 표현이다.

---

<br>

## 8. 푸리에 변환(Fourier Transform)

### 8.1 정의

복소 지수를 사용하여 푸리에 적분을 재표현하면:

$$f(x) = \frac{1}{2\pi} \int_{-\infty}^{\infty} e^{i\omega x} \underbrace{\int_{-\infty}^{\infty} f(\xi) e^{-i\omega\xi}\,d\xi}_{\hat{f}(\omega)}\,d\omega = \frac{1}{2\pi} \int_{-\infty}^{\infty} \hat{f}(\omega)\,e^{i\omega x}\,d\omega$$

**푸리에 변환** 과 그 **역변환**:

$$\boxed{\hat{f}(\omega) = \mathcal{F}\{f(x)\} = \int_{-\infty}^{\infty} f(\xi)\,e^{-i\omega\xi}\,d\xi}$$

$$\boxed{f(x) = \mathcal{F}^{-1}\{\hat{f}(\omega)\} = \frac{1}{2\pi} \int_{-\infty}^{\infty} \hat{f}(\omega)\,e^{i\omega x}\,d\omega}$$

### 8.2 푸리에 변환의 성질

**1) 선형성(Linearity):**

$$\mathcal{F}\{\alpha f(x) + \beta g(x)\} = \alpha\,\mathcal{F}\{f\} + \beta\,\mathcal{F}\{g\}$$

$$\mathcal{F}^{-1}\{\alpha\hat{f} + \beta\hat{g}\} = \alpha\,\mathcal{F}^{-1}\{\hat{f}\} + \beta\,\mathcal{F}^{-1}\{\hat{g}\}$$

**2) 미분 성질(Derivative Property):**

$$\mathcal{F}\left\{\frac{d}{dx}f(x)\right\} = i\omega\,\mathcal{F}\{f\}$$

단, $f(x) \to 0$ as $x \to \pm\infty$.

*증명 (부분 적분에 의해):*

$$\int_{-\infty}^{\infty} \frac{d}{d\xi}f(\xi)\,e^{-i\omega\xi}\,d\xi = \underbrace{\left[f(\xi)e^{-i\omega\xi}\right]_{-\infty}^{\infty}}_{= 0} - (-i\omega)\int_{-\infty}^{\infty} f(\xi)\,e^{-i\omega\xi}\,d\xi = i\omega\,\mathcal{F}\{f\}$$

**일반 $n$차 도함수:**

$$\mathcal{F}\left\{\frac{d^n}{dx^n}f(x)\right\} = (i\omega)^n\,\mathcal{F}\{f\}$$

> **[중요성]** 미분 성질은 공간/시간 영역에서의 미분을 주파수 영역에서의 곱셈으로 변환한다. 이것이 푸리에 변환이 미분 방정식을 풀 때 매우 유용한 이유이다.

**3) 합성곱(Convolution):**

$$f * g := \int_{-\infty}^{\infty} f(x - \xi)\,g(\xi)\,d\xi = \mathcal{F}^{-1}\{\hat{g}(\omega)\,\hat{f}(\omega)\}$$

즉, **공간 영역에서의 합성곱** 은 **주파수 영역에서의 곱셈** 에 대응한다.

**4) $x$-이동(공간 이동, Spatial Shift):**

$$\mathcal{F}^{-1}\{e^{-ia\omega}\,\hat{f}(\omega)\} = f(x - a)$$

주파수 영역에서 $e^{-ia\omega}$를 곱하면 공간 영역에서 함수를 $a$만큼 이동시킨다.

**5) $\omega$-이동(주파수 이동, Frequency Shift):**

$$\mathcal{F}^{-1}\{\hat{f}(\omega - a)\} = e^{iax}\,f(x)$$

주파수 영역에서의 이동은 공간 영역에서 $e^{iax}$에 의한 변조(Modulation)에 대응한다.

---

<br>

## 9. 정현파 함수를 이용한 곡선 맞춤

### 9.1 정현파 함수의 매개변수

주기 함수 $f(t) = f(t + T)$를 고려하자. 여기서 $T$는 상수(주기)이다.

**정현파 함수(Sinusoidal Function)** 는 사인 또는 코사인으로 기술할 수 있는 임의의 파형이다:

$$f(t) = A_0 + C_1 \cos(\omega_0 t + \theta)$$

여기서:

| 매개변수 | 의미 |
|:---:|:---|
| $A_0$ | **평균값**(수직 오프셋) |
| $C_1$ | **진폭(Amplitude)**(진동의 높이) |
| $\omega_0$ | **각주파수(Angular Frequency)**(라디안/시간) |
| $\theta$ | **위상각(Phase Angle)**(라디안) |

주파수 양들 간의 관계:

$$\omega_0 = 2\pi f, \qquad f = \frac{1}{T}$$

여기서 $f$는 사이클/시간 $[\text{Hz}]$ 단위의 **보통 주파수(Ordinary Frequency)** 이다.

### 9.2 위상각(Phase Angle)

위상각 $\theta$는 $t = 0$에서 코사인 함수가 새 주기를 시작하는 지점까지의 거리(라디안)를 나타낸다.

$\theta > 0$일 때:

- $\cos(\omega_0 t - \theta)$: **지연(Lagging)** 위상각 (오른쪽 이동)
- $\cos(\omega_0 t + \theta)$: **선행(Leading)** 위상각 (왼쪽 이동)

또한 사인과 코사인의 관계에 주목하자:

$$\sin(t + \pi/2) = \cos(t), \qquad \cos(t - \pi/2) = \sin(t)$$

### 9.3 대안적 형태

덧셈 공식(Angle Addition Formula)을 사용하여 코사인 함수에서 $\theta$를 분리할 수 있다:

$$\cos(\omega_0 t + \theta) = \cos(\omega_0 t)\cos(\theta) - \sin(\omega_0 t)\sin(\theta)$$

따라서:

$$f(t) = A_0 + C_1\cos(\omega_0 t + \theta)$$

$$= A_0 + C_1\cos\theta\,\cos(\omega_0 t) - C_1\sin\theta\,\sin(\omega_0 t)$$

$$\boxed{= A_0 + A_1\cos(\omega_0 t) + B_1\sin(\omega_0 t)}$$

여기서:

$$A_1 = C_1\cos\theta, \qquad B_1 = -C_1\sin\theta$$

이로부터 다음을 복원할 수 있다:

$$\theta = \arctan\left(-\frac{B_1}{A_1}\right), \qquad C_1 = \sqrt{A_1^2 + B_1^2}$$

---

<br>

## 10. 정현파의 최소제곱 적합

### 10.1 문제 설정

데이터 점 $(t_i, y_i)$가 주어졌을 때, 다음 모델을 적합한다:

$$y = A_0 + A_1\cos(\omega_0 t) + B_1\sin(\omega_0 t) = \hat{y}(t)$$

기저 함수를 이용한 선형 회귀 문제로 해석한다:

$$f_0 = 1, \quad f_1 = \cos(\omega_0 t), \quad f_2 = \sin(\omega_0 t)$$

$$\hat{y}(t) = \beta_0 f_0 + \beta_1 f_1 + \beta_2 f_2$$

### 10.2 정규 방정식

오차 제곱합을 최소화한다:

$$SSE = \sum (y_i - A_0 - A_1\cos(\omega_0 t_i) - B_1\sin(\omega_0 t_i))^2$$

편도함수를 취하고 0으로 놓으면:

$$\frac{\partial SSE}{\partial A_0} = 0, \quad \frac{\partial SSE}{\partial A_1} = 0, \quad \frac{\partial SSE}{\partial B_1} = 0$$

이로부터 정규 방정식 시스템을 얻는다:

$$\begin{pmatrix} \sum y_i \\ \sum y_i \cos(\omega_0 t) \\ \sum y_i \sin(\omega_0 t) \end{pmatrix} = \begin{pmatrix} n & \sum\cos(\omega_0 t) & \sum\sin(\omega_0 t) \\ \sum\cos(\omega_0 t) & \sum\cos^2(\omega_0 t) & \sum\cos(\omega_0 t)\sin(\omega_0 t) \\ \sum\sin(\omega_0 t) & \sum\sin(\omega_0 t)\cos(\omega_0 t) & \sum\sin^2(\omega_0 t) \end{pmatrix} \begin{pmatrix} A_0 \\ A_1 \\ B_1 \end{pmatrix}$$

### 10.3 등간격 점에서의 간소화

한 주기 $T$ 내에서 등간격으로 $n$개의 관측값을 수집한다: $T = n\,\Delta t$.

**등간격** 점의 경우, 다음 직교성 성질이 성립한다:

$$\sum_{i=1}^{n} \sin(\omega_0 t_i) = 0, \qquad \sum_{i=1}^{n} \cos(\omega_0 t_i) = 0$$

$$\sum_{i=1}^{n} \sin^2(\omega_0 t_i) = \frac{n}{2}, \qquad \sum_{i=1}^{n} \cos^2(\omega_0 t_i) = \frac{n}{2}$$

$$\sum_{i=1}^{n} \cos(\omega_0 t_i)\sin(\omega_0 t_i) = \frac{1}{2}\sum_{i=1}^{n}[\sin(2\omega_0 t_i) - \sin 0] = 0$$

정규 방정식은 **대각** 시스템으로 간소화된다:

$$\begin{pmatrix} \sum y_i \\ \sum y_i\cos(\omega_0 t) \\ \sum y_i\sin(\omega_0 t) \end{pmatrix} = \begin{pmatrix} n & 0 & 0 \\ 0 & \frac{n}{2} & 0 \\ 0 & 0 & \frac{n}{2} \end{pmatrix} \begin{pmatrix} A_0 \\ A_1 \\ B_1 \end{pmatrix}$$

### 10.4 닫힌 형태의 해

$$\boxed{A_0 = \frac{1}{n}\sum y_i = \bar{y}}$$

$$\boxed{A_1 = \frac{2}{n}\sum y_i \cos(\omega_0 t_i)}$$

$$\boxed{B_1 = \frac{2}{n}\sum y_i \sin(\omega_0 t_i)}$$

> **[핵심 통찰]** 등간격 데이터 점에서는 정규 방정식이 완전히 분리되어, 푸리에 계수에 대한 단순한 닫힌 형태의 표현식을 제공한다. 이것은 큰 계산적 이점이다.

---

<br>

## 11. 연속 푸리에 급수(일반 모델)

### 11.1 일반 푸리에 급수

$m$개의 고조파(Harmonics)($2m + 1$개의 계수)를 갖는 일반 모델로 확장:

$$f(t) = a_0 + \sum_{k=1}^{\infty} \bigl(a_k \cos(k\omega_0 t) + b_k \sin(k\omega_0 t)\bigr)$$

여기서:

- $\omega_0 = \frac{2\pi}{T}$는 **기본 주파수(Fundamental Frequency)**
- $2\omega_0, 3\omega_0, \ldots, m\omega_0, \ldots$는 **고조파(Harmonics)**($\omega_0$의 정수 배)

일반 모델의 계수:

$$A_0 = \frac{1}{n}\sum y_i = \bar{y}$$

$$A_j = \frac{2}{n}\sum y_i \cos(j\omega_0 t_i), \quad j = 1, 2, \ldots, m$$

$$B_j = \frac{2}{n}\sum y_i \sin(j\omega_0 t_i), \quad j = 1, 2, \ldots, m$$

> **[회귀 vs 배치]** $n > 2m + 1$(데이터 점이 계수보다 많은 경우), 계수를 **회귀(Regression)**(최소제곱) 의미로 계산한다. $n = 2m + 1$(데이터 점 수가 푸리에 계수 수와 같은 경우), 이를 **배치(Collocation)** 라 하며 — 푸리에 급수가 데이터를 정확히 보간한다.

### 11.2 오일러 공식과 복소 형태

**오일러 공식(Euler's Formula):**

$$e^{i\theta} = \cos\theta + i\sin\theta$$

이를 통해 푸리에 급수를 **복소 지수 형태** 로 표현할 수 있다:

$$f(t) = a_0 + \sum_{k=1}^{\infty}\bigl(a_k\cos(k\omega_0 t) + b_k\sin(k\omega_0 t)\bigr) = \sum_{k=-\infty}^{\infty} \tilde{C}_k\,e^{ik\omega_0 t}$$

여기서 복소 푸리에 계수는:

$$\tilde{C}_k = \frac{1}{T}\int_{-T/2}^{T/2} f(t)\,e^{-ik\omega_0 t}\,dt$$

---

<br>

## 12. 주파수 영역과 시간 영역

**주파수 영역(Frequency Domain)** 은 진동하는 함수의 동작을 특성화하기 위한 대안적 관점을 제공한다.

- **시간 영역(Time Domain)**: 일반적인 표현 $f(t)$ — 진폭 대 시간
- **주파수 영역(Frequency Domain)**: 푸리에 계수 $\hat{f}(\omega)$ — 진폭 대 주파수

푸리에 해석은 신호를 시간 영역에서 주파수 영역으로 변환하여, 어떤 주파수가 존재하고 그 상대적 강도가 어떠한지를 드러낸다.

---

<br>

## 13. 푸리에 적분과 푸리에 변환(재방문)

$T \to \infty$를 취하여 주기 함수에서 비주기 함수로 전환:

**푸리에 적분(역변환):**

$$f(t) = \frac{1}{2\pi} \int_{-\infty}^{\infty} F(\omega)\,e^{i\omega t}\,d\omega$$

**푸리에 변환:**

$$F(\omega) = \int_{-\infty}^{\infty} f(t)\,e^{-i\omega t}\,dt$$

---

<br>

## 14. 이산 푸리에 변환(DFT)

**유한 개수** 의 함수값을 수집한다. $0$에서 $T$까지의 구간을 $n$개의 등간격 부분 구간으로 나누며 $\Delta t = \frac{T}{n}$이다.

이산 표본: $f_j = f(t = t_j)$, 여기서 $t_j = t_0 + j\,\Delta t$.

**DFT (순방향):**

$$\hat{f}_k = \mathcal{F}\{f\}_k \approx \sum_{j=0}^{n-1} f_j\,e^{-ik\omega_0 j}, \quad k = 0, 1, 2, \ldots, n-1$$

**역DFT (Inverse DFT):**

$$f_j = \mathcal{F}^{-1}\{\hat{f}\}_j \approx \frac{1}{n}\sum_{k=0}^{n-1} \hat{f}_k\,e^{ik\omega_0 j}, \quad j = 0, 1, 2, \ldots, n-1$$

여기서 $\omega_0 = \frac{2\pi}{n}$.

```python
import numpy as np

def dft(f):
    """Compute the Discrete Fourier Transform of a 1D array f."""
    n = len(f)
    F = np.zeros(n, dtype=complex)
    for k in range(n):
        for j in range(n):
            F[k] += f[j] * np.exp(-2j * np.pi * k * j / n)
    return F

def idft(F):
    """Compute the Inverse Discrete Fourier Transform of a 1D array F."""
    n = len(F)
    f = np.zeros(n, dtype=complex)
    for j in range(n):
        for k in range(n):
            f[j] += F[k] * np.exp(2j * np.pi * k * j / n)
    return f / n
```

> **[복잡도]** 단순 DFT의 시간 복잡도는 $O(n^2)$이다. 이는 $n$개의 출력값 각각이 $n$개의 입력값에 대한 합을 필요로 하기 때문이다.

---

<br>

## 15. 나이퀴스트 주파수(Nyquist Frequency)

**나이퀴스트 주파수** 는 신호에서 측정할 수 있는 최대 주파수이다. 이는 **표본 추출 주파수의 절반** 과 같다:

$$f_{\max} = \frac{1}{2} f_s = \frac{1}{2\,\Delta t}$$

**예제:** 표본 추출 주파수 $f_s = 1000$ Hz로 100개의 표본($n = 100$)을 수집:

- 표본 추출 간격: $\Delta t = \frac{1}{f_s} = \frac{1}{1000} = 10^{-3}$ sec/sample
- 총 시간: $T = n \cdot \Delta t = 100 \times 10^{-3} = 0.1$ s
- **나이퀴스트 주파수:** $f_{\max} = \frac{1}{2}f_s = 500$ Hz
- **검출 가능한 최저 주파수:** $f_{\min} = \frac{1}{T} = \frac{1}{n\,\Delta t} = 10$ Hz

> **[앨리어싱(Aliasing)]** 신호에 나이퀴스트 주파수를 초과하는 주파수가 포함되면, 해당 주파수가 더 낮은 주파수로 잘못 표현되는 현상이 발생한다 — 이를 **앨리어싱(Aliasing)** 이라 한다. 이를 방지하려면, 표본 추출률이 신호에 존재하는 최대 주파수의 최소 2배여야 한다(**나이퀴스트-섀넌 표본화 정리, Nyquist-Shannon Sampling Theorem**).

---

<br>

## 16. 고속 푸리에 변환(FFT)

**고속 푸리에 변환(Fast Fourier Transform)** 은 다음을 활용하여 계산 비용을 $O(n^2)$ (DFT)에서 $O(n \log n)$으로 줄인다:

1. 삼각 함수의 **대칭성**
2. 복소 지수의 **주기성**

FFT의 핵심 성질:

- FFT는 입력을 **짝수 인덱스와 홀수 인덱스** 로 재귀적으로 분할한다
- FFT는 DFT 공식의 **대칭성과 주기성** 을 활용한다
- FFT는 반복 계산을 극적으로 **줄인다**

```python
import numpy as np

# Using NumPy's built-in FFT
t = np.linspace(0, 1, 256, endpoint=False)
signal = np.sin(2 * np.pi * 10 * t) + 0.5 * np.sin(2 * np.pi * 20 * t)

# Compute FFT
F = np.fft.fft(signal)

# Compute frequencies
freqs = np.fft.fftfreq(len(t), d=t[1] - t[0])

# Inverse FFT to recover original signal
recovered = np.fft.ifft(F)
```

> **[왜 $n = 2^p$인가?]** 쿨리-튜키(Cooley-Tukey) FFT 알고리즘은 $n$이 2의 거듭제곱일 때 가장 효율적으로 작동한다. 반으로의 재귀적 분할이 $\log_2 n$번 균등하게 수행될 수 있기 때문이다.

---

<br>

## 17. 파워 스펙트럼(Power Spectrum)

**파워 스펙트럼** 은 신호의 전력이 서로 다른 주파수에 어떻게 분포되어 있는지를 보여준다.

**이산의 경우** (이산 주파수 $k\omega_0$에서의 전력):

$$P_k = |\tilde{C}_k|^2$$

**연속의 경우** (연속 주파수 $\omega$에서의 전력):

$$P(\omega) = |\hat{F}(\omega)|^2$$

```python
import numpy as np
import matplotlib.pyplot as plt

# Generate a signal with two frequency components
n = 256
dt = 0.001
t = np.arange(n) * dt
signal = np.sin(2 * np.pi * 50 * t) + 0.5 * np.sin(2 * np.pi * 120 * t)

# Compute FFT and power spectrum
F = np.fft.fft(signal)
freqs = np.fft.fftfreq(n, d=dt)
power = np.abs(F) ** 2

# Plot only positive frequencies
pos_mask = freqs >= 0
plt.plot(freqs[pos_mask], power[pos_mask])
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power |F(w)|^2")
plt.title("Power Spectrum")
plt.show()
```

---

<br>

## 18. 요약 표

| 주제 | 핵심 수식 | 비고 |
|:---|:---|:---|
| **푸리에 급수** ($2\pi$ 주기) | $f(x) = \frac{a_0}{2} + \sum_{n=1}^{\infty}(a_n\cos nx + b_n\sin nx)$ | $[-\pi,\pi]$에서의 주기 함수 |
| **계수** ($2\pi$ 주기) | $a_n = \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\cos nx\,dx$, $b_n = \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\sin nx\,dx$ | 직교성으로부터 유도 |
| **일반 주기** $2\ell$ | $f(x) = \frac{a_0}{2} + \sum a_n\cos\frac{n\pi}{\ell}x + \sum b_n\sin\frac{n\pi}{\ell}x$ | 변수 치환 $x' = \frac{\pi x}{\ell}$ |
| **푸리에 적분** | $f(x) = \int_0^{\infty}[a(\omega)\cos\omega x + b(\omega)\sin\omega x]\,d\omega$ | 비주기; $\ell \to \infty$ |
| **푸리에 변환** | $\hat{f}(\omega) = \int_{-\infty}^{\infty}f(\xi)e^{-i\omega\xi}\,d\xi$ | 복소 지수 형태 |
| **역 푸리에 변환** | $f(x) = \frac{1}{2\pi}\int_{-\infty}^{\infty}\hat{f}(\omega)e^{i\omega x}\,d\omega$ | 원래 함수 복원 |
| **정현파 적합** | $f(t) = A_0 + A_1\cos(\omega_0 t) + B_1\sin(\omega_0 t)$ | $C_1 = \sqrt{A_1^2+B_1^2}$, $\theta = \arctan(-B_1/A_1)$ |
| **이산 계수** | $A_0 = \bar{y}$, $A_j = \frac{2}{n}\sum y_i\cos(j\omega_0 t_i)$, $B_j = \frac{2}{n}\sum y_i\sin(j\omega_0 t_i)$ | 등간격 데이터 |
| **오일러 공식** | $e^{i\theta} = \cos\theta + i\sin\theta$ | 삼각함수와 지수함수 연결 |
| **복소 푸리에 급수** | $f(t) = \sum_{k=-\infty}^{\infty}\tilde{C}_k e^{ik\omega_0 t}$ | $\tilde{C}_k = \frac{1}{T}\int_{-T/2}^{T/2}f(t)e^{-ik\omega_0 t}\,dt$ |
| **DFT** | $\hat{f}_k = \sum_{j=0}^{n-1}f_j e^{-ik\omega_0 j}$ | $O(n^2)$ 복잡도 |
| **FFT** | DFT와 같은 결과 | $O(n\log n)$ 복잡도 |
| **나이퀴스트 주파수** | $f_{\max} = \frac{1}{2}f_s = \frac{1}{2\Delta t}$ | 검출 가능한 최대 주파수 |
| **파워 스펙트럼** | $P_k = \lvert\tilde{C}_k\rvert^2$, $P(\omega) = \lvert\hat{F}(\omega)\rvert^2$ | 주파수별 에너지 분포 |
| **FT 선형성** | $\mathcal{F}\{\alpha f + \beta g\} = \alpha\mathcal{F}\{f\} + \beta\mathcal{F}\{g\}$ | 선형 연산자 |
| **FT 미분** | $\mathcal{F}\{f^{(n)}\} = (i\omega)^n\mathcal{F}\{f\}$ | 미분이 곱셈으로 변환 |
| **FT 합성곱** | $\mathcal{F}\{f*g\} = \hat{f}(\omega)\cdot\hat{g}(\omega)$ | 합성곱이 곱으로 변환 |
| **FT $x$-이동** | $\mathcal{F}^{-1}\{e^{-ia\omega}\hat{f}(\omega)\} = f(x-a)$ | 공간 이동 |
| **FT $\omega$-이동** | $\mathcal{F}^{-1}\{\hat{f}(\omega-a)\} = e^{iax}f(x)$ | 주파수 변조 |
