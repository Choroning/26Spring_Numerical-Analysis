# 강의 노트

[레포지토리 안내](../README.ko.md) | [English](README.md) | [한국어](README.ko.md)

> **최종 수정일:** 2026-10-27

## 1부 모델링, 컴퓨터, 오차 분석

| 단원 | 주제 | 문서 |
|:-----|:-----|:-----|
| [1장 수학적 모델링과 공학 문제 해결](Ch01_Mathematical-Modeling/) | 수학적 모델, 뉴턴의 제2법칙과 번지 점프 예제, 오일러 방법, 수치적으로 푸는 공학 문제의 유형, NumPy 배열과 행렬 연산, Matplotlib 시각화, NumPy의 뷰와 복사본 | [강의 노트](Ch01_Mathematical-Modeling/Concepts_Lecture.ko.md), [실습 노트](Ch01_Mathematical-Modeling/Concepts_Lab.ko.md) |
| [2장 Python 기초와 프로그래밍](Ch02_Python-Fundamentals-and-Programming/) | 프로그래밍 환경, Python 기본 문법과 NumPy, SciPy 함수, 제어 구조와 함수, Python으로 구현한 번지 점프 예제, 변수 스코프와 LEGB 규칙, 구조적 프로그래밍, 예외 처리, 데코레이터, 객체 지향 기초, 타입 힌트, 제너레이터 | [강의 노트](Ch02_Python-Fundamentals-and-Programming/Concepts_Lecture.ko.md), [실습 노트](Ch02_Python-Fundamentals-and-Programming/Concepts_Lab.ko.md) |
| [3장 반올림 오차](Ch03_Roundoff-Errors/) | 유효 숫자, 정확도와 정밀도, 오차의 정의와 스카보로 기준, 반올림 오차와 부동 소수점 산술의 문제점, 오차 지표, 절단 오차와 테일러 급수 | [강의 노트](Ch03_Roundoff-Errors/Concepts_Lecture.ko.md), [실습 노트](Ch03_Roundoff-Errors/Concepts_Lab.ko.md) |
| [4장 절단 오차](Ch04_Truncation-Errors/) | 테일러 급수 전개, 근사의 차수, 테일러 급수를 이용한 수치 미분, 총 수치 오차, 안정성과 조건수, 유한 차분의 Python 구현 | [강의 노트](Ch04_Truncation-Errors/Concepts_Lecture.ko.md) |

## 2부 근과 최적화

| 단원 | 주제 | 문서 |
|:-----|:-----|:-----|
| [5장 근: 구간법](Ch05_Roots-Bracketing-Methods/) | 근 문제와 그래프 방법, 구간법과 개방법의 비교, 증분 탐색, 이분법과 그 수렴, 가위치법과 수정 가위치법 | [강의 노트](Ch05_Roots-Bracketing-Methods/Concepts_Lecture.ko.md) |
| [6장 근: 개방법](Ch06_Roots-Open-Methods/) | 고정점 반복, Wegstein 방법, 뉴턴-랩슨법, 할선법, 역 이차 보간, 브렌트 방법 | [강의 노트](Ch06_Roots-Open-Methods/Concepts_Lecture.ko.md) |
| [7장 최적화](Ch07_Optimization/) | 국소 최적성과 전역 최적성, 최적해의 존재, 1차원 최적화, 황금 분할 탐색, 포물선 보간, 등고선과 기울기, 볼록 문제, 최급 하강법 | [강의 노트](Ch07_Optimization/Concepts_Lecture.ko.md) |

## 3부 선형계

| 단원 | 주제 | 문서 |
|:-----|:-----|:-----|
| [8장 연립 일차 대수방정식의 풀이](Ch08_Linear-Algebraic-Equations/) | 소규모 시스템의 그래프 방법, 행렬식과 크래머 법칙, 단순 가우스 소거법과 계산 비용, 단순 가우스 소거법의 함정, 질량-스프링 시스템 응용 | [강의 노트](Ch08_Linear-Algebraic-Equations/Concepts_Lecture.ko.md), [실습 노트](Ch08_Linear-Algebraic-Equations/Concepts_Lab.ko.md) |
| [9장 가우스 소거법: 개선된 방법](Ch09_Gauss-Elimination/) | 부분 피봇팅과 스케일링된 부분 피봇팅, 가우스-조르단 방법, 삼중대각 시스템과 토마스 알고리즘, 조건수, 열전도 문제 | [강의 노트](Ch09_Gauss-Elimination/Concepts_Lecture.ko.md), [실습 노트](Ch09_Gauss-Elimination/Concepts_Lab.ko.md) |
| [10장 LU 분해](Ch10_LU-Factorization/) | 가우스 소거법과 LU 분해의 관계, 전진 대입과 후진 대입, 피봇팅을 적용한 LU 분해, 촐레스키 분해, Doolittle 알고리즘, 다중 우변 문제와 열저항 네트워크 | [강의 노트](Ch10_LU-Factorization/Concepts_Lecture.ko.md), [실습 노트](Ch10_LU-Factorization/Concepts_Lab.ko.md) |
| [11장 역행렬과 조건수](Ch11_Matrix-Inverse-and-Condition/) | 역행렬과 LU 분해를 이용한 계산, 질량-스프링 시스템의 물리적 해석, 악조건의 판별, 벡터 노름과 행렬 노름, 행렬 조건수와 오차 한계 정리, 힐베르트 행렬 | [강의 노트](Ch11_Matrix-Inverse-and-Condition/Concepts_Lecture.ko.md), [실습 노트](Ch11_Matrix-Inverse-and-Condition/Concepts_Lab.ko.md) |
| [12장 선형계의 반복법](Ch12_Iterative-Methods/) | 가우스-자이델 방법, 대각 우위와 수렴, 야코비 반복법, 완화법, 전처리, 고유값 계산을 위한 거듭제곱법, 역거듭제곱법과 편이 거듭제곱법 | [강의 노트](Ch12_Iterative-Methods/Concepts_Lecture.ko.md), [실습 노트](Ch12_Iterative-Methods/Concepts_Lab.ko.md) |
| [13장 고유값 방법: 대칭 행렬](Ch13_Eigenvalues-Symmetric-Matrices/) | 고유값과 고유벡터의 계산, 고유값과 상미분방정식의 안정성, 결합 2차 ODE 시스템, 질량-스프링 시스템, 거듭제곱법, 대칭 행렬의 스펙트럼 분해, 진동 해석 | [강의 노트](Ch13_Eigenvalues-Symmetric-Matrices/Concepts_Lecture.ko.md), [실습 노트](Ch13_Eigenvalues-Symmetric-Matrices/Concepts_Lab.ko.md) |

## 4부 곡선 적합

| 단원 | 주제 | 문서 |
|:-----|:-----|:-----|
| [14장 선형 회귀](Ch14_Linear-Regression/) | 통계 복습과 정규 분포, 직선 최소제곱 회귀, 적합도의 정량화 | [강의 노트](Ch14_Linear-Regression/Concepts_Lecture.ko.md) |
| [15장 비선형 회귀](Ch15_Nonlinear-Regression/) | 표준 오차와 결정 계수, 비선형 관계의 선형화, 다항 회귀, 다중 선형 회귀, 일반 선형 최소제곱 회귀, 가우스-뉴턴 비선형 회귀 | [강의 노트](Ch15_Nonlinear-Regression/Concepts_Lecture.ko.md) |
| [16장 푸리에 해석](Ch16_Fourier-Analysis/) | 우함수와 기함수, 푸리에 급수와 푸리에 적분, 푸리에 변환, 정현파의 최소제곱 적합, 주파수 영역과 시간 영역, 이산 푸리에 변환, 나이퀴스트 주파수, 고속 푸리에 변환, 파워 스펙트럼 | [강의 노트](Ch16_Fourier-Analysis/Concepts_Lecture.ko.md) |
| [17장 다항식 보간](Ch17_Interpolation/) | 다항식 보간, 뉴턴 보간 다항식, 라그랑주 보간 다항식, 역보간, 외삽과 진동 | [강의 노트](Ch17_Interpolation/Concepts_Lecture.ko.md) |
| [18장 스플라인과 구간별 보간](Ch18_Splines/) | 연속성 조건, 선형 스플라인, 이차 스플라인, 삼차 스플라인, 쌍선형 보간 | [강의 노트](Ch18_Splines/Concepts_Lecture.ko.md) |

## 5부 적분과 미분

| 단원 | 주제 | 문서 |
|:-----|:-----|:-----|
| [19장 수치 적분 공식](Ch19_Numerical-Integration-Formulas/) | 사다리꼴 법칙, 심프슨의 1/3 법칙과 3/8 법칙, 복합 적분 공식, 적분 법칙의 결합 | [강의 노트](Ch19_Numerical-Integration-Formulas/Concepts_Lecture.ko.md) |
| [20장 함수의 수치 적분](Ch20_Numerical-Integration-of-Functions/) | 롬버그 적분, 가우스 구적법 | [강의 노트](Ch20_Numerical-Integration-of-Functions/Concepts_Lecture.ko.md) |
| [21장 수치 미분](Ch21_Numerical-Differentiation/) | 1계 도함수와 2계 도함수의 유한 차분 근사, 고차 전진 차분, 리처드슨 외삽법, 비등간격 데이터의 도함수, 수치 편도함수, 미정 계수법, 미분 가능 프로그래밍 | [강의 노트](Ch21_Numerical-Differentiation/Concepts_Lecture.ko.md) |

## 6부 상미분방정식

| 단원 | 주제 | 문서 |
|:-----|:-----|:-----|
| [22장 상미분방정식의 초기값 문제](Ch22_Initial-Value-Problems/) | 오일러 방법, 국소 절단 오차와 전역 절단 오차, 2차 룽게-쿠타 방법, 고전적 4차 룽게-쿠타 방법, 연립 ODE, 수치적 안정성 | [강의 노트](Ch22_Initial-Value-Problems/Concepts_Lecture.ko.md) |
| [23장 강성 ODE와 적응적 방법](Ch23_Stiff-ODEs/) | 적응적 룽게-쿠타 방법, 내장 RK 방법 RKF45와 스텝 크기 제어, 다단계 방법, 강성 시스템 | [강의 노트](Ch23_Stiff-ODEs/Concepts_Lecture.ko.md) |
| [24장 경계값 문제](Ch24_Boundary-Value-Problems/) | 경계값 문제와 해석해, 사격법, 유한 차분법 | [강의 노트](Ch24_Boundary-Value-Problems/Concepts_Lecture.ko.md) |
