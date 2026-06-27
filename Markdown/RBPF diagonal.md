## Hypothèses :

- supposons que l'on peut décomposer le vecteur d'état sous la forme : $$
x_k = \begin{bmatrix}
x_k^l \\
x_k^n
\end{bmatrix}
$$
- de telle sorte que l'on puisse décomposer $x_k = f_k(x_{k-1}, w_k)$ et $z_k = h_k(x_k, n_k)$ en : $$
\left\{\begin{matrix}
x_{k+1}^l & = & F_k^l(x_k^n) \: x_k^l & + & f_k^l(x_k^n) & + & B_k^l \: u_k & + & w_k^l \\
x_{k+1}^n & = & F_k^n(x_k^n) \: x_k^l & + & f_k^n(x_k^n) & + & B_k^n \: u_k & + & w_k^n \\
z_k & = & H_k(x_k^n) \: x_k^l & + & h_k(x_k^n) & & & + & n_k
\end{matrix}\right.
$$ avec $$
w_k = \begin{bmatrix}
w_k^l \\
w_k^n
\end{bmatrix} \sim \mathcal{N}\left(0, \begin{bmatrix}
Q_k^l & Q_k^{ln} \\
{Q_k^{ln}}^T & Q_k^n
\end{bmatrix}\right)
\hspace{1cm}
n_k \sim \mathcal{N}(0, R_k)
$$
- on suppose aussi que $x_0^l$ est gaussien, et que la loi de $x_0^n$ est connue
- afin d'alléger les notations, posons : $$
\tilde F_k^l = F_k^l(\hat{x}_{k|k}^n, u_k, 0)
\text{,} \quad
\tilde F_k^n = F_k^n(\hat{x}_{k|k}^n, u_k, 0)
\quad \text{et} \quad
\tilde H_k = H_k(\hat{x}_{k|k}^n)
$$

## Modèle diagonal :

- dans ce modèle, on suppose que : $$
f_k^l(\cdot) = 0
\text{,} \quad
\tilde F_k^n = 0
\quad \text{et} \quad
Q_k^{ln} = 0
$$
- posons également $z_k^{(1)} = z_k - h_k(x_k^n)$
- on a donc : $$
\left\{\begin{matrix}
x_{k+1}^l & = & F_k^l(x_k^n) \: x_k^l & & & + & B_k^l \: u_k & + & w_k^l \\
x_{k+1}^n & = & & & f_k^n(x_k^n) & + & B_k^n \: u_k & + & w_k^n \\
z_k^{(1)} & = & H_k(x_k^n) \: x_k^l & & & & & + & n_k
\end{matrix}\right.
$$
- en particulier : $$
\left\{\begin{matrix}
x_{k+1}^l & = & F_k^l(x_k^n) \: x_k^l & + & B_k^l \: u_k & + & w_k^l \\
z_k^{(1)} & = & H_k(x_k^n) \: x_k^l & & & + & n_k
\end{matrix}\right.
$$

## Algorithme :

- la partie non-linéaire est indépendante de la partie linéaire, et est évaluée par un filtre particulaire
- on retrouve exactement la forme des équations de Kalman pour la partie linéaire, et on peut donc appliquer le filtre comme d'habitude
- on réécrit le terme d'innovation pour faire apparaître la différence entre l'observation réelle et l'observation prédite
	- correction
		- $\tilde{y}_k = z_k - \left(\tilde H_k \: \hat{x}_{k|k-1}^l - h_k(\hat{x}_{k|k}^n)\right)$ ($\hat{x}_{k|k}^n$ est une donnée du FdK estimée par le FP)
		- $S_k = \tilde H_k \: P_{k|k-1} \: \tilde H_k^T + R_k$
		- $K_k = P_{k|k-1} \: \tilde H_k^T \: S_k^{-1}$
		- $\hat{x}_{k|k}^l = \hat{x}_{k|k-1}^l + K_k \: \tilde{y}_k$
		- $P_{k|k} = (I - K_k \: \tilde H_k) \: P_{k|k-1}$
	- prédiction
		- $\hat{x}_{k+1|k}^l = \tilde F_k^l \: \hat{x}_{k|k}^l + B_k^l \: u_k$
		- $P_{k+1|k} = \tilde F_k^l \: P_{k|k} \: {\tilde F_k^l}^T + Q_k^l$
