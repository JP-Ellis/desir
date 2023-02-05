Runge-Kutta methods are a family of both implicit and explicit iterative methods
to approximate solutions of simultaneous non-linear equations. They attempt to
solve an initial value problem (IVP) of the form

```math
\begin{aligned}
  \ddfrac{y}{t} &= f(t, y) \\
  y(t_0) &= y_0
\end{aligned}
```

where `$y(t)$` and `$f(t, y)$` will be, in general, vector-valued functions.

## Mathematical Background

The estimate `$y_{n+1}$` from `$y_n$` is given by:

```math
y_{n+1} = y_n + h \sum_{i=1}^s b_i k_i
```

where `$h$` is the timestep and the `$k_i$` are intermediate estimates of the
solution defined by

```math
k_i = f \left( t_n + c_i h, y_n + \sum_{j=1}^{s} a_{ij} k_j \right)
```

The coefficients `$a_{ij}$`, `$b_i$`, and `$c_i$` entirely define the
Runge-Kutta method, though specific implementations may use different algorithms
to what is naively shown above in order to exploit specific features of these
coefficients.

A number of Runge-Kutta methods will have a second (lower order) method embedded
with coefficients `$b^*_i$` which can be used to estimate the error of the
higher order method. This error estimate is given by

```math
e_{n+1} = h \sum_{i=1}^s d_i k_i,
```

where the `$d_i$` are the differences between the coefficients of the embedded
method and the coefficients of the higher order method:

```math
d_i = b^*_i - b_i.
```

All of the above coefficients are usually represented in a _Butcher tableau_

```math
\begin{array}{c|c}
  \vt c & A \\
  \hline
  & \vt b^\transpose \\
  & \vt {b^*}^\transpose
\end{array}
=
\begin{array}{c|cccc}
  c_1    & a_{11} & a_{12} & \cdots & a_{1s} \\
  c_2    & a_{21} & a_{22} & \cdots & a_{2s} \\
  \vdots & \vdots & \vdots & \ddots & \vdots \\
  c_s    & a_{s1} & a_{s2} & \cdots & a_{ss} \\
  \hline
         & b_1   & b_2   & \cdots & b_s \\
         & b^*_1 & b^*_2 & \cdots & b^*_s
 \end{array}
```

and the matrix `$A$` is called the _Runge-Kutta matrix_ will the coefficients
`$\vt b$` and `$\vt c$` are called the _weights_ and _nodes_ respectively.

### Explicit Methods

Explicit methods are those where the intermediate estimates `$k_i$` can be
computed entirely from `$y_n$` and `$k_1$` through `$k_{i-1}$`. As a result,
these methods will always have `$a_{ij} = 0$` for `$j \geq i$` (i.e., the matrix
`$A$` is strictly lower triangular).

Explicit methods are easier to compute than implicit methods, but suffer from a
small region of absolute stability.

### Implicit Methods

Implicit methods are those where the intermediate estimates `$k_i$` cannot be
computed entirely from `$y_n$` and `$k_1$` through `$k_{i-1}$`. As a result,
these methods will typically require solving a system of equations.

### Terminology

A Runge-Kutta method is said to be _of order `$p$`_ if the local truncation
error is `$O(h^{p+1})$`. The number of stages `$s$` is the number of
intermediate estimates `$k_i$` required to compute the solution estimate
`$y_{n+1}$`.
