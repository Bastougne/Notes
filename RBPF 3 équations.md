## Marginalisation :

- supposons que l'on peut décomposer le vecteur d'état sous la forme $$
x_k=\begin{bmatrix}x_k^l\\x_k^n\end{bmatrix}
$$ de telle sorte que l'on puisse décomposer $x_k=f_k(x_{k-1},w_k)$ et $z_k=h_k(x_k,n_k)$ en $$
(S):\left\{\begin{matrix}
x_k^l&=&F_k^l(x_{k-1}^n)\:x_{k-1}^l&+&f_k^l(x_{k-1}^n)&+&B_k^l\:u_k&+&w_k^l\\
x_k^n&=&F_k^n(x_{k-1}^n)\:x_{k-1}^l&+&f_k^n(x_{k-1}^n)&+&B_k^n\:u_k&+&w_k^n\\
z_k&=&H_k(x_k^n)\:x_k^l&+&h_k(x_k^n)&&&+&n_k
\end{matrix}\right.
$$ avec $$
w_k=\begin{bmatrix}w_k^l\\w_k^n\end{bmatrix}\sim\mathcal{N}\left(0,\begin{bmatrix}Q_k^l&Q_k^{ln}\\(Q_k^{ln})^T&Q_k^n\end{bmatrix}\right),\quad n_k\sim\mathcal{N}(0,R_k)
$$
- on suppose aussi que $x_0^l$ est gaussien, et que la loi de $x_0^n$ est connue
- on souhaite revenir aux équations de Kalman, on pose donc $$
\begin{align}
\widetilde{x_k^l}:=&x_k^l-f_k^l(x_{k-1}^n)\\
z_{k-1}^{(1)}:=&x_k^n-\left(f_k^n(x_{k-1}^n)+B_k^n\:u_k\right)\\
z_k^{(2)}:=&z_k-h_k(x_k^n)
\end{align}
$$
- on peut alors réécrire $(S)$ sous la forme $$
(\tilde{S}):\left\{\begin{matrix}
\widetilde{x_k^l}&=&F_k^l(x_{k-1}^n)\:x_{k-1}^l&+&B_k^l\:u_k&+&w_k^l\\
z_{k-1}^{(1)}&=&F_k^n(x_{k-1}^n)\:x_{k-1}^l&&&+&w_k^n\\
z_k^{(2)}&=&H_k(x_k^n)\:x_k^l&&&+&n_k
\end{matrix}\right.
$$

## Décorrélation des bruits

- on ne peut pas appliquer directement Kalman car $w_k^l$ et $w_k^n$ sont corrélées ($Q_k^{ln}\neq0$ en toutes généralités)
- définissons $\overline{w_k^l}$ le bruit linéaire décorrélé $$
\overline{w_k^l}:=w_k^l-Q_k^{ln}\:(Q_k^n)^{-1}\:w_k^n
$$
- sa matrice de covariance est $$
\begin{align}
\overline{Q_k^l}=&\mathbb{E}\left(\overline{w_k^l}\:\overline{w_k^l}^T\right)\\
=&\mathbb{E}\left((w_k^l-Q_k^{ln}\:(Q_k^n)^{-1}\:w_k^n)\:(w_k^l-Q_k^{ln}\:(Q_k^n)^{-1}\:w_k^n)^T\right)\\
\\
=&\mathbb{E}\left(w_k^l\:(w_k^l)^T\right)-\mathbb{E}\left(w_k^l\:(Q_k^{ln}\:(Q_k^n)^{-1}\:w_k^n)^T\right)-\\
&\phantom{\mathbb{E}((}\mathbb{E}\left((Q_k^{ln}\:(Q_k^n)^{-1}\:w_k^n)\:(w_k^l)^T\right)+\mathbb{E}\left((Q_k^{ln}\:(Q_k^n)^{-1}\:w_k^n)\:(Q_k^{ln}\:(Q_k^n)^{-1}\:w_k^n)^T\right)\\
\\
=&\mathbb{E}\left(w_k^l\:(w_k^l)^T\right)-\mathbb{E}\left(w_k^l\:(w_k^n)^T\right)\:(Q_k^n)^{-T}\:(Q_k^{ln})^T-\\
&\phantom{\mathbb{E}((}Q_k^{ln}\:(Q_k^n)^{-1}\:\mathbb{E}\left(w_k^n\:(w_k^l)^T\right)+Q_k^{ln}\:(Q_k^n)^{-1}\:\mathbb{E}\left(w_k^n\:(w_k^n)^T\right)\:(Q_k^n)^{-T}\:(Q_k^{ln})^T\\
\\
=&Q_k^l-Q_k^{ln}\:(Q_k^n)^{-T}\:(Q_k^{ln})^T-Q_k^{ln}\:(Q_k^n)^{-1}\:(Q_k^{ln})^T+Q_k^{ln}\:(Q_k^n)^{-1}\:Q_k^n\:(Q_k^n)^{-T}\:(Q_k^{ln})^T\\
=&Q_k^l-Q_k^{ln}\:(Q_k^n)^{-1}\:Q_k^{ln}-Q_k^{ln}\:(Q_k^n)^{-1}\:Q_k^{ln}+Q_k^{ln}\:(Q_k^n)^{-1}\:Q_k^{ln}\\
=&Q_k^l-Q_k^{ln}\:(Q_k^n)^{-1}\:Q_k^{ln}
\end{align}
$$
- on a alors $$
\begin{align}
\mathbb{E}\left(\overline{w_k^l}\:(w_k^n)^T\right)=&\mathbb{E}\left((w_k^l-Q_k^{ln}\:(Q_k^n)^{-1}\:w_k^n)\:(w_k^n)^T\right)\\
=&\mathbb{E}\left(w_k^l\:(w_k^n)^T\right)-Q_k^{ln}\:(Q_k^n)^{-1}\:\mathbb{E}\left(w_k^n\:(w_k^n)^T\right)\\
=&Q_k^{ln}-Q_k^{ln}\:(Q_k^n)^{-1}\:Q_k^n\\
=&0
\end{align}
$$
- en posant $$
\overline{F_k^l(x_{k-1}^n)}:=\left(F_k^l(x_{k-1}^n)-Q_k^{ln}\:(Q_k^n)^{-1}\:F_k^n(x_{k-1}^n)\right)
$$ $(\tilde{S})$ se réécrit $$
\left\{\begin{matrix}
\widetilde{x_k^l}&=&\overline{F_k^l(x_{k-1}^n)}\:x_{k-1}^l&+&B_k^l\:u_k&+&\overline{w_k^l}&+&Q_k^{ln}\:(Q_k^n)^{-1}\:z_{k-1}^{(1)}\\
z_{k-1}^{(1)}&=&F_k^n(x_{k-1}^n)\:x_{k-1}^l&&&+&w_k^n\\
z_k^{(2)}&=&H_k(x_k^n)\:x_k^l&&&+&n_k
\end{matrix}\right.
$$
- posons finalement $$
\overline{x_k^l}:=\widetilde{x_k^l}-Q_k^{ln}\:(Q_k^n)^{-1}\:z_{k-1}^{(1)}
$$ pour obtenir $$
(\overline{S}):\left\{\begin{matrix}
\overline{x_k^l}&=&\overline{F_k^l(x_{k-1}^n)}\:x_{k-1}^l&+&B_k^l\:u_k&+&\overline{w_k^l}\\
z_{k-1}^{(1)}&=&F_k^n(x_{k-1}^n)\:x_{k-1}^l&&&+&w_k^n\\
z_k^{(2)}&=&H_k(x_k^n)\:x_k^l&&&+&n_k
\end{matrix}\right.
$$

## Estimation de Kalman :

- le système se décompose en 2 sous-systèmes : $$
\begin{matrix}
\left\{\begin{matrix}
{x_{k-1}^l}^*&=&x_{k-1}^l&&&&\\
z_{k-1}^{(1)}&=&F_k^n(x_{k-1}^n)\:x_{k-1}^l&&&+&w_k^n\\
\end{matrix}\right.\\
\\
\left\{\begin{matrix}
\overline{x_k^l}&=&\overline{F_k^l(x_{k-1}^n)}\:{x_{k-1}^l}^*&+&B_k^l\:u_k&+&\overline{w_k^l}\\
z_k^{(2)}&=&H_k(x_k^n)\:x_k^l&&&+&n_k
\end{matrix}\right.
\end{matrix}
$$
- le premier sous-système n'a pas de dynamique, et son observation correspond à la sortie du filtre particulaire de la partie non-linéaire
- le second correspond à la prédiction puis correction classiques du filtre de Kalman
- afin d'alléger les notations, posons $$
\tilde{F}_k^l:=F_k^l(\hat{x}_{k-1|k-1}^n,u_k,0)\text{,}\quad\tilde{F}_k^n:=F_k^n(\hat{x}_{k-1|k-1}^n,u_k,0)\quad\text{et}\quad\tilde{H}_k=H_k(\hat{x}_{k|k}^n)
$$
- pré-correction
	- $\hat{z}_{k-1}^{(1)}=\hat{x}_{k|k-1}^n-\left(f_k^n(\hat{x}_{k-1|k-1}^n)+B_k^n\:u_k\right)$
	- $y_{k-1}^{(1)}=\hat{z}_{k-1}^{(1)}-\tilde{F}_k^n\:\hat{x}_{k-1|k-1}^l$
	- $S_{k-1}^{(1)}=\tilde{F}_k^n\:P_{k-1|k-1}^l\:(\tilde{F}_k^n)^T+Q_k^n$
	- $K_{k-1}^{(1)}=P_{k-1|k-1}^l\:(\tilde{F}_k^n)^T\:(S_{k-1}^{(1)})^{-1}$
	- ${\hat{x}_{k-1|k-1}^l}^*=\hat{x}_{k-1|k-1}^l+K_{k-1}^{(1)}\:y_{k-1}^{(1)}$
	- ${P_{k-1|k-1}^l}^*=(I-K_{k-1}^{(1)}\:\tilde{F}_k^n)\:P_{k-1|k-1}^l$
- prédiction :
	- $\overline{\tilde{F}_k^l}=\tilde{F}_k^l-Q_k^{ln}\:(Q_k^n)^{-1}\:\tilde{F}_k^n$
	- $\hat{x}_{k|k-1}^l=\overline{\tilde{F}_k^l}\:{\hat{x}_{k-1|k-1}^l}^*+f_k^l(\hat{x}_{k-1|k-1}^n)+B_k^l\:u_k+Q_k^{ln}\:(Q_k^n)^{-1}\:\hat{z}_{k-1}^{(1)}$
	- $P_{k|k-1}=\overline{\tilde{F}_k^l}\:{P_{k-1|k-1}^l}^*\:\left(\overline{\tilde{F}_k^l}\right)^T+\left(Q_k^l-Q_k^{ln}\:(Q_k^n)^{-1}\:Q_k^{ln}\right)$
- correction :
	- $y_k^{(2)}=z_k-\left(\tilde{H}_k\:\hat{x}_{k|k-1}^l+h_k()\right)$
	- $S_k^{(2)}=\tilde{H}_k\:P_{k|k-1}\:(\tilde{H}_k)^T+R_k$
	- $K_k^{(2)}=\tilde{H}_k\:P_{k|k-1}\:(\tilde{H}_k)^T+R_k$
	- $K_k = P_{k|k-1} \: \tilde H_k^T \: S_k^{-1}$
	- $\hat{x}_{k|k} = \hat{x}_{k|k-1} + K_k \: \tilde{y}_k$
	- $P_{k|k} = (I - K_k \: \tilde H_k) \: P_{k|k-1}$
