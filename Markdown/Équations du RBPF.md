## Marginalisation :

- supposons que l'on peut décomposer le vecteur d'état sous la forme $$
x_k=\begin{bmatrix}x_k^n\\x_k^l\end{bmatrix}
$$ de telle sorte qu'on ait $$
(S):\left\{\begin{matrix}
x_k^n&=&F_k^n(x_{k-1}^n)\:x_{k-1}^l&+&f_k^n(x_{k-1}^n)&+&b_k^{ne}&+&w_k^n\\
x_k^l&=&F_k^l(x_{k-1}^n)\:x_{k-1}^l&+&f_k^l(x_{k-1}^n)&+&b_k^{le}&+&w_k^l\\
z_k&=&H_k(x_k^n)\:x_k^l&+&h_k(x_k^n)&+&b_k^o&+&n_k
\end{matrix}\right.
$$ avec $$
w_k=\begin{bmatrix}w_k^n\\w_k^l\end{bmatrix}\sim\mathcal{N}\left(0,\begin{bmatrix}Q_k^n&Q_k^{nl}\\(Q_k^{nl})^T&Q_k^l\end{bmatrix}\right),\quad n_k\sim\mathcal{N}(0,R_k)
$$
- on suppose aussi que $x_0^l$ est gaussien, et que la loi de $x_0^n$ est connue
- on souhaite utiliser un filtre de Kalman pour obtenir une estimation de $x_k^l$ connaissant $x_k^n$ et $z_k$

## Décorrélation des bruits :

- on ne peut pas appliquer directement Kalman car $w_k^l$ et $w_k^n$ sont corrélées ($Q_k^{nl}\neq0$ en toutes généralités)
- définissons $\overline{w_k^l}$ le bruit linéaire décorrélé :$$
\overline{w_k^l}:=w_k^l-Q_k^{nl}\:(Q_k^n)^{-1}\:w_k^n
$$
- sa matrice de covariance est $$
\begin{align}
\overline{Q_k^l}=&\mathbb{E}\left(\overline{w_k^l}\:\overline{w_k^l}^T\right)\\
=&\mathbb{E}\left((w_k^l-Q_k^{nl}\:(Q_k^n)^{-1}\:w_k^n)\:(w_k^l-Q_k^{nl}\:(Q_k^n)^{-1}\:w_k^n)^T\right)\\
\\
=&\mathbb{E}\left(w_k^l\:(w_k^l)^T\right)-\mathbb{E}\left(w_k^l\:(Q_k^{nl}\:(Q_k^n)^{-1}\:w_k^n)^T\right)-\\
&\phantom{\mathbb{E}((}\mathbb{E}\left((Q_k^{nl}\:(Q_k^n)^{-1}\:w_k^n)\:(w_k^l)^T\right)+\mathbb{E}\left((Q_k^{nl}\:(Q_k^n)^{-1}\:w_k^n)\:(Q_k^{nl}\:(Q_k^n)^{-1}\:w_k^n)^T\right)\\
\\
=&\mathbb{E}\left(w_k^l\:(w_k^l)^T\right)-\mathbb{E}\left(w_k^l\:(w_k^n)^T\right)\:(Q_k^n)^{-T}\:(Q_k^{nl})^T-\\
&\phantom{\mathbb{E}((}Q_k^{nl}\:(Q_k^n)^{-1}\:\mathbb{E}\left(w_k^n\:(w_k^l)^T\right)+Q_k^{nl}\:(Q_k^n)^{-1}\:\mathbb{E}\left(w_k^n\:(w_k^n)^T\right)\:(Q_k^n)^{-T}\:(Q_k^{nl})^T\\
\\
=&Q_k^l-Q_k^{nl}\:(Q_k^n)^{-T}\:(Q_k^{nl})^T-Q_k^{nl}\:(Q_k^n)^{-1}\:(Q_k^{nl})^T+Q_k^{nl}\:(Q_k^n)^{-1}\:Q_k^n\:(Q_k^n)^{-T}\:(Q_k^{nl})^T\\
=&Q_k^l-Q_k^{nl}\:(Q_k^n)^{-1}\:Q_k^{nl}-Q_k^{nl}\:(Q_k^n)^{-1}\:Q_k^{nl}+Q_k^{nl}\:(Q_k^n)^{-1}\:Q_k^{nl}\\
=&Q_k^l-Q_k^{nl}\:(Q_k^n)^{-1}\:Q_k^{nl}
\end{align}
$$
- on a alors $$
\begin{align}
\mathbb{E}\left(\overline{w_k^l}\:(w_k^n)^T\right)=&\mathbb{E}\left((w_k^l-Q_k^{nl}\:(Q_k^n)^{-1}\:w_k^n)\:(w_k^n)^T\right)\\
=&\mathbb{E}\left(w_k^l\:(w_k^n)^T\right)-Q_k^{nl}\:(Q_k^n)^{-1}\:\mathbb{E}\left(w_k^n\:(w_k^n)^T\right)\\
=&Q_k^{nl}-Q_k^{nl}\:(Q_k^n)^{-1}\:Q_k^n\\
=&0
\end{align}
$$
- on écrit $$
\begin{align}
x_k^l=&x_k^l-Q_k^{nl}\:(Q_k^n)^{-1}\:x_k^n+Q_k^{nl}\:(Q_k^n)^{-1}\:x_k^n\\
=&\underbrace{\left(F_k^l(x_{k-1}^n)-Q_k^{nl}\:(Q_k^n)^{-1}\:F_k^n(x_{k-1}^n)\right)}_{\overline{F_k^l(x_{k-1}^n)}:=}\:x_{k-1}^l+\underbrace{f_k^l(x_{k-1}^n)-Q_k^{nl}\:(Q_k^n)^{-1}f_k^n(x_{k-1}^n)}_{\overline{f_k^l(x_{k-1}^n)}:=}+\\
&\underbrace{b_k^{le}-Q_k^{nl}\:(Q_k^n)^{-1}b_k^{ne}}_{\overline{b_k^{le}}:=}+Q_k^{nl}\:(Q_k^n)^{-1}\:x_k^n+\overline{w_k^l}
\end{align}
$$
- $(S)$ se réécrit alors $$
(\overline{S}):\left\{\begin{matrix}
x_k^n&=&F_k^n(x_{k-1}^n)\:x_{k-1}^l&+&f_k^n(x_{k-1}^n)&+&b_k^{ne}&+&&&w_k^n\\
x_k^l&=&\overline{F_k^l(x_{k-1}^n)}\:x_{k-1}^l&+&\overline{f_k^l(x_{k-1}^n)}&+&\overline{b_k^{le}}&+&Q_k^{nl}\:(Q_k^n)^{-1}\:x_k^n&+&\overline{w_k^l}\\
z_k&=&H_k(x_k^n)\:x_k^l&+&h_k(x_k^n)&+&b_k^o&+&&&n_k
\end{matrix}\right.
$$

## Estimation de Kalman :

- $(\overline{S})$ se décompose en 2 sous-systèmes $$
\begin{matrix}
(S_1):\left\{\begin{matrix}
{x_{k-1}^l}^*&=&x_{k-1}^l&&&&\\
x_k^n&=&F_k^n(x_{k-1}^n)\:{x_{k-1}^l}^*&+&f_k^n(x_{k-1}^n)&+&b_k^{ne}&+&w_k^n\\
\end{matrix}\right.\\
\\
(S_2):\left\{\begin{matrix}
x_k^l&=&\overline{F_k^l(x_{k-1}^n)}\:{x_{k-1}^l}^*&+&\overline{f_k^l(x_{k-1}^n)}&+&\overline{b_k^{le}}&+&Q_k^{nl}\:(Q_k^n)^{-1}\:x_k^n&+&\overline{w_k^l}\\
z_k&=&H_k(x_k^n)\:x_k^l&+&h_k(x_k^n)&+&b_k^o&+&&&n_k
\end{matrix}\right.
\end{matrix}
$$
- le premier sous-système n'a pas de dynamique, et son observation correspond à la sortie du filtre particulaire de la partie non-linéaire
- le second correspond à la prédiction puis correction classiques du filtre de Kalman
- par souci de clarté, on omettra les indices des particules, et on notera $$
\left\{\begin{align}
\tilde{F}_k^n=&F_k^n\left({x_{k-1}^n}^{(i)}\right)\\
\tilde{F}_k^l=&\overline{F_k^l\left({x_{k-1}^n}^{(i)}\right)}\\
\tilde{H}_k=&H_k\left({x_k^n}^{(i)}\right)\\
\tilde{b}_k^{ne}=&f_k^n\left({x_{k-1}^n}^{(i)}\right)+b_k^{ne}\\
\tilde{b}_k^{le}=&\overline{f_k^l\left({x_{k-1}^n}^{(i)}\right)}+\overline{b_k^{le}}+Q_k^{nl}\:(Q_k^n)^{-1}\:{x_k^n}^{(i)}\\
\tilde{b}_k^o=&h_k\left({x_k^n}^{(i)}\right)+b_k^o\\
\end{align}\right.
$$
- pré-correction :
	- $y_k^{(1)}={x_k^n}^{(i)}-(\tilde{F}_k^n\:\hat{x}_{k|k-1}^l+\tilde{b}_k^{ne})$
	- $S_k=\tilde{F}_k^n\:P_{k-1|k-1}^l\:(\tilde{F}_k^n)^T+Q_k^n$
	- $K_k=P_{k-1|k-1}^l\:(\tilde{F}_k^n)^T\:S_k^{-1}$
	- ${\hat{x}_{k-1|k-1}^l}^*=\hat{x}_{k-1|k-1}^l+K_k\:y_k^{(1)}$
	- ${P_{k-1|k-1}^l}^*=(I-K_k\:\tilde{F}_k^n)\:P_{k-1|k-1}^l$

- prédiction :
	- $\hat{x}_{k|k-1}^l=\tilde{F}_k^l\:{\hat{x}_{k-1|k-1}^l}^*+\tilde{b}_k^{le}$
	- $P_{k|k-1}^l=\tilde{F}_k^l\:{P_{k-1|k-1}^l}^*\:(\tilde{F}_k^l)^T+\overline{Q_k^l}$
- correction :
	- $y_k^{(2)}=z_k-(\tilde{H}_k\:\hat{x}_{k|k-1}^l+\tilde{b}_k^o)$
	- $N_k=\tilde{H}_k\:P_{k|k-1}^l\:\tilde{H}_k^T+R_k$
	- $L_k=P_{k|k-1}^l\:\tilde{H}_k^T\:N_k^{-1}$
	- $\hat{x}_{k|k}^l=\hat{x}_{k|k-1}^l+L_k\:y_k^{(2)}$
	- $P_{k|k}^l=(I-L_k\:\tilde{H}_k)\:P_{k|k-1}^l$

## Cas d'unicité du filtre

- en général, il faut un filtre de Kalman par particule
- si $F_k^n$, $F_k^l$ et $H_k$ ne dépendent pas de $x_k^n$ ni de $x_{k-1}^n$, alors $K_k^{(1)}$ et $K_k^{(1)}$ sont uniques pour tout le nuage de particules
- on peut donc calculer les gains une fois pour toutes les particules, ce qui réduit grandement le coup calculatoire du RBPF
