## Modèle général à bruits corrélés :

- le modèle le plus général de navigation pour le RBPF, qu'on obtient en généralisant le [[RBPF triangulaire]]
- on ne fait plus d'hypothèse sur $f_k^l$, $\tilde F_k^n$ ni $Q_k^{ln}$
- en toute généralité, on a donc : $$
\left\{\begin{matrix}
x_{k+1}^l & = & F_k^l(x_k^n) \: x_k^l & + & f_k^l(x_k^n) & + & B_k^l \: u_k & + & w_k^l \\
z_k^{(2)} & = & F_k^n(x_k^n) \: x_k^l & & & & & + & w_k^n \\
z_k^{(1)} & = & H_k(x_k^n) \: x_k^l & & & & & + & n_k
\end{matrix}\right.
$$
- posons ${x_{k+1}^l}^{(1)} = x_{k+1}^l - f_k^l(x_k^n)$
- on retrouve alors les équations du modèle triangulaire : $$
\left\{\begin{matrix}
{x_{k+1}^l}^{(1)} & = & F_k^l(x_k^n) \: x_k^l & + & B_k^l \: u_k & + & w_k^l \\
z_k^{(2)} & = & F_k^n(x_k^n) \: x_k^l & & & + & w_k^n \\
z_k^{(1)} & = & H_k(x_k^n) \: x_k^l & & & + & n_k
\end{matrix}\right.
$$
- puisque $Q_k^{ln} \neq 0$ ici, il faut décorréler $w_k^l$ et $w_k^n$ avant de pouvoir appliquer les équations du modèle précédent
- définissons $\overline{w_k^l}$ le bruit linéaire décorrélé : $$
\overline{w_k^l} =
w_k^l - Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: w_k^n
$$
- sa matrice de covariance est : $$
\begin{align}
\overline{Q_k^l}
& = \mathbb{E}\left(\overline{w_k^l} \: {\overline{w_k^l}}^T\right) \\
& = \mathbb{E}\left((w_k^l - Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: w_k^n) \:
{(w_k^l - Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: w_k^n)}^T\right) \\
& = \mathbb{E}\left(w_k^l \: {w_k^l}^T\right) -
\mathbb{E}\left(w_k^l \: (Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: w_k^n)^T\right) -
\mathbb{E}\left((Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: w_k^n) \: {w_k^l}^T\right) +
\mathbb{E}\left((Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: w_k^n) \: (Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: w_k^n)^T\right) \\
& = \mathbb{E}\left(w_k^l \: {w_k^l}^T\right) -
\mathbb{E}\left(w_k^l \: {w_k^{n\phantom{l\!}}}^T\right) \: {Q_k^{n\phantom{l\!}}}^{-T} \: {Q_k^{ln}}^T -
Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: \mathbb{E}\left(w_k^n \: {w_k^l}^T\right) +
Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: \mathbb{E}\left(w_k^n \: {w_k^{n\phantom{l\!}}}^T\right) \: {Q_k^{n\phantom{l\!}}}^{-T} \: {Q_k^{ln}}^T \\
& = Q_k^l -
Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-T} \: {Q_k^{ln}}^T -
Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: {Q_k^{ln}}^T +
Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: Q_k^n \: {Q_k^{n\phantom{l\!}}}^{-T} \: {Q_k^{ln}}^T \\
& = Q_k^l -
Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: {Q_k^{ln}}^T -
Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: {Q_k^{ln}}^T +
Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: {Q_k^{ln}}^T \\
& = Q_k^l -
Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: {Q_k^{ln}}^T
\end{align}$$
- on a alors : $$
\begin{align}
\mathbb{E}\left(\overline{w_k^l} \: {w_k^{n\phantom{l\!}}}^T\right)
& = \mathbb{E}\left((w_k^l - Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: w_k^n) \: {w_k^{n\phantom{l\!}}}^T\right) \\
& = \mathbb{E}\left(w_k^l \: {w_k^{n\phantom{l\!}}}^T\right) -
Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: \mathbb{E}\left(w_k^n \: {w_k^{n\phantom{l\!}}}^T\right) \\
& = Q_k^{ln} - Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: Q_k^n \\
& = 0
\end{align}
$$
- on peut alors réécrire le modèle sous la forme : $$
\left\{\begin{matrix}
{x_{k+1}^l}^{(1)} & = & \left(F_k^l(x_k^n) - Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: F_k^n(x_k^n)\right) \: x_k^l & + & B_k^l \: u_k & + & \overline{w_k^l} & + & Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: \left(F_k^n(x_k^n) \: x_k^l + w_k^n\right) \\
z_k^{(2)} & = & F_k^n(x_k^n) \: x_k^l & & & + & w_k^n \\
z_k^{(1)} & = & H_k(x_k^n) \: x_k^l & & & + & n_k
\end{matrix}\right.
$$
- posons : $$
\overline{F_k^l(x_k^n)} = F_k^l(x_k^n) - Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: F_k^n(x_k^n)
$$
- d'où : $$
\left\{\begin{matrix}
{x_{k+1}^l}^{(1)} & = & \overline{F_k^l(x_k^n)} \: x_k^l & + & B_k^l \: u_k & + & \overline{w_k^l} & + & Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: z_k^{(2)} \\
z_k^{(2)} & = & F_k^n(x_k^n) \: x_k^l & & & + & w_k^n \\
z_k^{(1)} & = & H_k(x_k^n) \: x_k^l & & & + & n_k
\end{matrix}\right.
$$
- en séparant le système en deux comme précédemment, on obtient : $$
\begin{matrix}
\left\{\begin{matrix}
{x_k^l}^* & = & x_k^l & & & & \\
z_k^{(1)} & = & H_k(x_k^n) \: x_k^l & & & + & n_k
\end{matrix}\right.
\\[5mm]
\left\{\begin{matrix}
{x_{k+1}^l}^{(1)} & = & \overline{F_k^l(x_k^n)} \: {x_k^l}^* & + & B_k^l \: u_k & + & \overline{w_k^l} & + & Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: z_k^{(2)} \\
z_k^{(2)} & = & F_k^n(x_k^n) \: {x_k^l}^* & & & + & w_k^n
\end{matrix}\right.
\end{matrix}
$$
- finalement, posons ${x_{k+1}^l}^{(2)} = {x_{k+1}^l}^{(1)} - Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: z_k^{(2)}$
- on retombe sur le modèle triangulaire séparé, de la forme : $$
\begin{matrix}
\left\{\begin{matrix}
{x_k^l}^* & = & x_k^l & & & & \\
z_k^{(1)} & = & H_k(x_k^n) \: x_k^l & & & + & n_k
\end{matrix}\right.
\\[5mm]
\left\{\begin{matrix}
{x_{k+1}^l}^{(2)} & = & \overline{F_k^l(x_k^n)} \: {x_k^l}^* & + & B_k^l \: u_k & + & \overline{w_k^l} \\
z_k^{(2)} & = & F_k^n(x_k^n) \: {x_k^l}^* & & & + & w_k^n
\end{matrix}\right.
\end{matrix}
$$

## Algorithme :

- on va appliquer les résultats du modèle triangulaire à ${x_{k+1}^l}^{(2)}$
- il faut remplacer $\tilde F_k^l$ et $Q_k^l$ par $\tilde{\overline{F_k^l}}$ et $\overline{Q_k^l}$ respectivement afin de prendre les corrélations de bruits de modèle en compte
- afin de gagner en lisibilité, on remplace directement ${x_{k+1}^l}^{(2)}$ et $\overline{Q_k^l}$ par leur expression
- on a donc :
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
	- prédiction
		- $\hat{z}_k^{(2)} = \hat{x}_{k+1|k}^n - \left(f_k^n(\hat{x}_{k|k}^n) + B_k^n \: u_k\right) = \tilde{y}_k^{(2)} + \tilde F_k^n \: {\hat{x}_{k|k}^l}^*$
		- $\tilde{\overline{F_k^l}} = \tilde F_k^l - Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: \tilde F_k^n$
		- $\hat{x}_{k+1|k}^l = \tilde{\overline{F_k^l}} \: \hat{x}_{k|k}^l + f_k^l(\hat{x}_{k|k}^n) + B_k^l \: u_k + Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: z_k^{(2)}$
		- $P_{k+1|k} = \tilde{\overline{F_k^l}} \: P_{k|k} \: {\tilde{\overline{F_k^l}}}^T + \left(Q_k^l - Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: {Q_k^{ln}}^T\right)$
- en remplaçant $\overline{F_k^l}$ par son expression dans le calcul de $\hat{x}_{k+1|k}^l$, on obtient $$
\begin{align}
\hat{x}_{k+1|k}^l & = \tilde{\overline{F_k^l}} \: \hat{x}_{k|k}^l + f_k^l(\hat{x}_{k|k}^n) + B_k^l \: u_k + Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: z_k^{(2)} \\
& = \tilde{F_k^l} \: \hat{x}_{k|k}^l - Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: \tilde F_k^n \: \hat{x}_{k|k}^l + f_k^l(\hat{x}_{k|k}^n) + B_k^l \: u_k + Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: z_k^{(2)} \\
& = \tilde{F_k^l} \: \hat{x}_{k|k}^l + f_k^l(\hat{x}_{k|k}^n) + B_k^l \: u_k + Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: \left(- \tilde F_k^n \: \hat{x}_{k|k}^l + \tilde{y}_k^{(2)} + \tilde F_k^n \: {\hat{x}_{k|k}^l}^*\right) \\
& = \tilde{F_k^l} \: \hat{x}_{k|k}^l + f_k^l(\hat{x}_{k|k}^n) + B_k^l \: u_k + Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: \left(- \tilde F_k^n \: L_k \: \tilde{y}_k^{(2)} + \tilde{y}_k^{(2)}\right) \\
& = \tilde{F_k^l} \: \hat{x}_{k|k}^l + f_k^l(\hat{x}_{k|k}^n) + B_k^l \: u_k + Q_k^{ln} \: {Q_k^{n\phantom{l\!}}}^{-1} \: \left( I - \tilde F_k^n \: L_k\right) \: \tilde{y}_k^{(2)}
\end{align}
$$

## Remarques :

- puisque la seconde correction concerne la dynamique du système, on considère qu'elle fait partie de l'étape de prédiction
- on retrouve alors la structure usuelle du filtre en deux étapes :
	- correction (uniquement la première correction)
	- prédiction (seconde correction + prédiction à proprement parler)
- en toute généralité, il faut un filtre de Kalman par particule
- si $\tilde F_k^l$, $\tilde F_k^n$ et $\tilde H_k$ ne dépendent pas de $\hat{x}_{k|k}^n$, alors les gains $K_k$ et $L_k$, et $P_{k|k}$ non plus
- dans ce cas, on peut se contenter d'un unique filtre de Kalman pour gérer la partie linéaire
