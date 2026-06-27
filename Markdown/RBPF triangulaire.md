## Modèle triangulaire :

- il s'agit d'une généralisation du [[RBPF diagonal]]
- ici, supposons uniquement que : $$
f_k^l(\cdot) = 0
\quad \text{et} \quad
Q_k^{ln} = 0
$$
- on a alors : $$
\left\{\begin{matrix}
x_{k+1}^l & = & F_k^l(x_k^n) \: x_k^l & & & + & B_k^l \: u_k & + & w_k^l \\
x_{k+1}^n & = & F_k^n(x_k^n) \: x_k^l & + & f_k^n(x_k^n) & + & B_k^n \: u_k & + & w_k^n \\
z_k^{(1)} & = & H_k(x_k^n) \: x_k^l & & & & & + & n_k
\end{matrix}\right.
$$
- posons $z_k^{(2)} = x_{k+1}^n - \left(f_k^n(x_k^n) + B_k^n \: u_k\right)$
- on retrouve un modèle diagonal, mais à deux mesures : $$
\left\{\begin{matrix}
x_{k+1}^l & = & F_k^l(x_k^n) \: x_k^l & + & B_k^l \: u_k & + & w_k^l \\
z_k^{(2)} & = & F_k^n(x_k^n) \: x_k^l & & & + & w_k^n \\
z_k^{(1)} & = & H_k(x_k^n) \: x_k^l & & & + & n_k
\end{matrix}\right.
$$
- on sépare les mesures en réécrivant le système comme suit : $$
\begin{matrix}
\left\{\begin{matrix}
{x_k^l}^* & = & x_k^l & & & & \\
z_k^{(1)} & = & H_k(x_k^n) \: x_k^l & & & + & n_k
\end{matrix}\right.
\\[5mm]
\left\{\begin{matrix}
x_{k+1}^l & = & F_k^l(x_k^n) \: {x_k^l}^* & + & B_k^l \: u_k & + & w_k^l \\
z_k^{(2)} & = & F_k^n(x_k^n) \: {x_k^l}^* & & & + & w_k^n
\end{matrix}\right.
\end{matrix}
$$

## Algorithme :

- puisqu'on a deux systèmes, on va appliquer un filtre de Kalman sur chacun
- premier système :
	- correction
		- $\hat{z}_k^{(1)} = z_k - h_k(\hat{x}_{k|k}^n)$
		- $\tilde{y}_k^{(1)} = \hat{z}_k^{(1)} - \tilde H_k \: \hat{x}_{k|k-1}^l$
		- $S_k = \tilde H_k \: P_{k|k-1} \: \tilde H_k^T + R_k$
		- $K_k = P_{k|k-1} \: \tilde H_k^T \: S_k^{-1}$
		- ${\hat{x}_{k|k}^l}^* = \hat{x}_{k|k-1}^l + K_k \: \tilde{y}_k^{(1)}$
		- $P_{k|k}^* = (I - K_k \: \tilde H_k) \: P_{k|k-1}$
	- prédiction
		- ${\hat{x}_{k+1|k}^l}^* = {\hat{x}_{k|k}^l}^*$
		- $P_{k+1|k}^* = P_{k|k}^*$
- second système :
	- correction
		- $\hat{z}_k^{(2)} = \hat{x}_{k+1|k}^n - \left(f_k^n(\hat{x}_{k|k}^n) + B_k^n \: u_k\right)$
		- $\tilde{y}_k^{(2)} = \hat{z}_k^{(2)} - \tilde F_k^n \: {\hat{x}_{k+1|k}^l}^*$
		- $N_k = \tilde F_k^n \: P_{k+1|k}^* \: {\tilde F_k^n}^T + Q_k^n$
		- $L_k = P_{k+1|k}^* \: {\tilde F_k^n}^T \: N_k^{-1}$
		- $\hat{x}_{k|k}^l = {\hat{x}_{k+1|k}^l}^* + L_k \: \tilde{y}_k^{(2)}$
		- $P_{k|k} = (I - L_k \: \tilde F_k^n) \: P_{k+1|k}^*$
	- prédiction
		- $\hat{x}_{k+1|k}^l = \tilde F_k^l \: \hat{x}_{k|k}^l + B_k^l \: u_k$
		- $P_{k+1|k} = \tilde F_k^l \: P_{k|k} \: {\tilde F_k^l}^T + Q_k^l$
- finalement, on peut supprimer la première étape de prédiction, puisque le premier système n'a pas de dynamique
- on peut donc résumer ces deux filtres en un filtre à deux mesures :
	- première correction
		- $\tilde{y}_k^{(1)} = z_k - \left(\tilde H_k \: \hat{x}_{k|k-1}^l + h_k(\hat{x}_{k|k}^n)\right)$
		- $S_k = \tilde H_k \: P_{k|k-1} \: \tilde H_k^T + R_k$
		- $K_k = P_{k|k-1} \: \tilde H_k^T \: S_k^{-1}$
		- ${\hat{x}_{k|k}^l}^* = \hat{x}_{k|k-1}^l + K_k \: \tilde{y}_k^{(1)}$
		- $P_{k|k}^* = (I - K_k \: \tilde H_k) \: P_{k|k-1}$
	- seconde correction
		- $\tilde{y}_k^{(2)} = \hat{x}_{k+1|k}^n - \left(\tilde F_k^n \: {\hat{x}_{k|k}^l}^* + f_k^n(\hat{x}_{k|k}^n) + B_k^n \: u_k\right)$
		- $N_k = \tilde F_k^n \: P_{k|k}^* \: {\tilde F_k^n}^T + Q_k^n$
		- $L_k = P_{k|k}^* \: {\tilde F_k^n}^T \: N_k^{-1}$
		- $\hat{x}_{k|k}^l = {\hat{x}_{k|k}^l}^* + L_k \: \tilde{y}_k^{(2)}$
		- $P_{k|k} = (I - L_k \: \tilde F_k^n) \: P_{k|k}^*$
	- prédiction
		- $\hat{x}_{k+1|k}^l = \tilde F_k^l \: \hat{x}_{k|k}^l + B_k^l \: u_k$
		- $P_{k+1|k} = \tilde F_k^l \: P_{k|k} \: {\tilde F_k^l}^T + Q_k^l$
- remarque : on aurait pu changer l'ordre dans lequel on traite chaque système, on aurait alors obtenu une solution différente car les termes d'innovations auraient été différents
- cet ordre est le plus naturel (au niveau de l'interprétation des innovations)
