# [2026학년도 봄학기] 수치해석

![Last Commit](https://img.shields.io/github/last-commit/Choroning/26Spring_Numerical-Analysis)
![Languages](https://img.shields.io/github/languages/top/Choroning/26Spring_Numerical-Analysis)

이 레포지토리는 대학 강의 및 과제를 위해 작성된 수치해석 기법 구현 코드를 체계적으로 정리하고 보관합니다.

*작성자: 박철원 (고려대학교(세종), 컴퓨터소프트웨어학과) - 2026년 기준 3학년*
<br><br>

## 📑 목차

- [레포지토리 소개](#about-this-repository)
- [강의 정보](#course-information)
- [사전 요구사항](#prerequisites)
- [주차별 일정](#weekly-schedule)
- [레포지토리 구조](#repository-structure)
- [라이선스](#license)

---


<br><a name="about-this-repository"></a>
## 📝 레포지토리 소개

이 레포지토리에는 대학 수준의 수치해석 과목을 위해 작성된 이중 언어 학습 자료와 코드가 포함되어 있습니다:

- 매 강의 및 실습 세션별 이중 언어 개념 정리 노트 (한국어 `.ko.md` + 영어 `.md`)
- 교재 기반 수치 알고리즘의 Python 구현
- Chapra & Canale 전체 커리큘럼을 다루는 Part/Chapter 디렉토리 구조

> **🤖 AI 에이전트 활용**
> 본 과목은 AI 에이전트 사용을 권장합니다.
> 수업 전반에 걸쳐 [Claude Code](https://claude.ai/download)와 [Gemini CLI](https://github.com/google-gemini/gemini-cli)를 코딩 어시스턴트로 활용하였습니다.

<br><a name="course-information"></a>
## 📚 강의 정보

- **학기:** 2026학년도 봄학기 (3월 - 6월)
- **소속:** 고려대학교(세종)

|학수번호      |강의명    |이수구분|교수자|개설학과|
|:----------:|:-------|:----:|:------:|:----------------|
|`DCSS305-00`|수치해석(영강)|전공선택|강신후 교수|컴퓨터소프트웨어학과|

- **📖 참고 자료**

| 유형 | 내용 |
|:----:|:---------|
|교재|"Numerical Methods for Engineers" by Steven C. Chapra and Raymond P. Canale (8th Edition, McGraw Hill)|
|강의자료|교수자 제공 슬라이드 및 온라인 영상|

<br><a name="prerequisites"></a>
## ✅ 사전 요구사항

- 선형대수 및 이산수학에 대한 이해
- Python 인터프리터 설치
- 과학 계산 라이브러리 숙지

- **💻 개발 환경**

| 도구 | 회사 |  운영체제  | 비고 |
|:-----|:-------:|:----:|:------|
|Python 3|Python Software Foundation|macOS|    |
|NumPy|NumFOCUS|macOS|수치 연산|
|SciPy|NumFOCUS|macOS|과학 연산|
|Matplotlib|NumFOCUS|macOS|시각화|
|JupyterLab|Project Jupyter|macOS|대화형 노트북|

<br><a name="weekly-schedule"></a>
## 📅 주차별 일정

| 주차 | 강의 자료 | Part / Chapter | 주제 |
|:---:|:---|:---|:---|
| 1 | N00-N03 | Part1: Ch01, Ch02, Ch03 | 모델링, Python, 반올림 오차 |
| 2 | N04 | Part1: Ch04 | 절단 오차 |
| 3 | N05-N06 | Part3: Ch08, Ch09 | 선형 시스템, 가우스 소거법 |
| 4 | N07-N08 | Part3: Ch10, Ch11 | LU 분해, 역행렬 및 조건수 |
| 5 | N09-N10 | Part3: Ch12, Ch13 | 고유값 |
| 6 | N11-N12 | Part2: Ch05, Ch06 | 근 탐색 |
| 7 | N13-N14 | Part2: Ch07 | 최적화 |
| 8 | 중간고사 | — | — |
| 9 | N15-N16 | Part4: Ch14, Ch15 | 회귀분석 |
| 10 | N17-N18 | Part4: Ch16 | 푸리에 해석 |
| 11 | N19-N20 | Part4: Ch17, Ch18 | 보간법, 스플라인 |
| 12 | N21-N22 | Part5: Ch19, Ch20 | 수치 적분 |
| 13 | N23-N24 | Part5: Ch21 / Part6: Ch22 | 수치 미분, 초기값 문제 |
| 14 | N25-N26 | Part6: Ch22 | 초기값 문제 (계속) |
| 15 | N27-N28 | Part6: Ch23, Ch24 | 강성 ODE, 경계값 문제 |
| 16 | 기말고사 | — | — |

<br><a name="repository-structure"></a>
## 🗂 레포지토리 구조

```plaintext
26Spring_Numerical-Analysis
├── Part1_Modeling-Computers-and-Error-Analysis
│   ├── Chapter01_Mathematical-Modeling
│   │   ├── Concepts_Lab.ko.md
│   │   ├── Concepts_Lab.md
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter02_Programming-and-Software
│   │   ├── Concepts_Lab.ko.md
│   │   ├── Concepts_Lab.md
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter03_Approximations-and-Round-Off-Errors
│   │   ├── Concepts_Lab.ko.md
│   │   ├── Concepts_Lab.md
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   └── Chapter04_Truncation-Errors
│       ├── Concepts_Lecture.ko.md
│       └── Concepts_Lecture.md
├── Part2_Roots-and-Optimization
│   ├── Chapter05_Roots-Bracketing-Methods
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter06_Roots-Open-Methods
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   └── Chapter07_Optimization
│       ├── Concepts_Lecture.ko.md
│       └── Concepts_Lecture.md
├── Part3_Linear-Algebraic-Equations
│   ├── Chapter08_Linear-Algebraic-Equations
│   │   ├── Concepts_Lab.ko.md
│   │   ├── Concepts_Lab.md
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter09_Gauss-Elimination
│   │   ├── Concepts_Lab.ko.md
│   │   ├── Concepts_Lab.md
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter10_LU-Factorization
│   │   ├── Concepts_Lab.ko.md
│   │   ├── Concepts_Lab.md
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter11_Matrix-Inverse-and-Condition
│   │   ├── Concepts_Lab.ko.md
│   │   ├── Concepts_Lab.md
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter12_Eigenvalues-Power-Method
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   └── Chapter13_Eigenvalues-Symmetric-Matrices
│       ├── Concepts_Lecture.ko.md
│       └── Concepts_Lecture.md
├── Part4_Curve-Fitting
│   ├── Chapter14_Linear-Regression
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter15_Nonlinear-Regression
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter16_Fourier-Analysis
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter17_Interpolation
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   └── Chapter18_Splines
│       ├── Concepts_Lecture.ko.md
│       └── Concepts_Lecture.md
├── Part5_Integration-and-Differentiation
│   ├── Chapter19_Numerical-Integration-Formulas
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter20_Numerical-Integration-of-Equations
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   └── Chapter21_Numerical-Differentiation
│       ├── Concepts_Lecture.ko.md
│       └── Concepts_Lecture.md
├── Part6_Ordinary-Differential-Equations
│   ├── Chapter22_Initial-Value-Problems
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   ├── Chapter23_Stiff-ODEs
│   │   ├── Concepts_Lecture.ko.md
│   │   └── Concepts_Lecture.md
│   └── Chapter24_Boundary-Value-Problems
│       ├── Concepts_Lecture.ko.md
│       └── Concepts_Lecture.md
├── images/
├── .gitignore
├── LICENSE
├── README.ko.md
└── README.md
```

<br><a name="license"></a>
## 🤝 라이선스

이 레포지토리는 [MIT License](LICENSE) 하에 배포됩니다.

---
