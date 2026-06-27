## Marginalisation :

- supposons que l'on peut décomposer le vecteur d'état sous la forme $$
x_k=\begin{bmatrix}x_k^l\\x_k^n\end{bmatrix}
$$ de telle sorte qu'on ait $$
(S):\left\{\begin{matrix}
x_k^l&=&F_k^l(x_{k-1}^n)\:x_{k-1}^l&+&f_k^l(x_{k-1}^n)&+&b_k^{le}&+&w_k^l\\
x_k^n&=&F_k^n(x_{k-1}^n)\:x_{k-1}^l&+&f_k^n(x_{k-1}^n)&+&b_k^{ne}&+&w_k^n\\
z_k&=&H_k(x_k^n)\:x_k^l&+&h_k(x_k^n)&+&b_k^o&+&n_k
\end{matrix}\right.
$$ avec $$
w_k=\begin{bmatrix}w_k^l\\w_k^n\end{bmatrix}\sim\mathcal{N}\left(0,\begin{bmatrix}Q_k^l&Q_k^{ln}\\(Q_k^{ln})^T&Q_k^n\end{bmatrix}\right),\quad n_k\sim\mathcal{N}(0,R_k)
$$
- on suppose aussi que $x_0^l$ est gaussien, et que la loi de $x_0^n$ est connue
- on souhaite utiliser un filtre de Kalman pour obtenir une estimation de $x_k^l$ connaissant $x_k^n$ et $z_k$

## Décorrélation des bruits :

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
- on écrit $$
\begin{align}
x_k^l=&x_k^l-Q_k^{ln}\:(Q_k^n)^{-1}\:x_k^n+Q_k^{ln}\:(Q_k^n)^{-1}\:x_k^n\\
=&\underbrace{\left(F_k^l(x_{k-1}^n)-Q_k^{ln}\:(Q_k^n)^{-1}\:F_k^n(x_{k-1}^n)\right)}_{\overline{F_k^l(x_{k-1}^n)}:=}\:x_{k-1}^l+f_k^l(x_{k-1}^n)+\underbrace{b_k^{le}-Q_k^{ln}\:(Q_k^n)^{-1}\left(f_k^n(x_{k-1}^n)+b_k^{ne}-x_k^n\right)}_{\overline{b_k^{le}(x_k^n,x_{k-1}^n)}:=}+\overline{w_k^l}
\end{align}
$$
- $(S)$ se réécrit alors $$
(\overline{S}):\left\{\begin{matrix}
x_k^l&=&\overline{F_k^l(x_{k-1}^n)}\:x_{k-1}^l&+&f_k^l(x_{k-1}^n)&+&\overline{b_k^{le}(x_k^n,x_{k-1}^n)}&+&\overline{w_k^l}\\
x_k^n&=&F_k^n(x_{k-1}^n)\:x_{k-1}^l&+&f_k^n(x_{k-1}^n)&+&b_k^{ne}&+&w_k^n\\
z_k&=&H_k(x_k^n)\:x_k^l&+&h_k(x_k^n)&+&b_k^o&+&n_k
\end{matrix}\right.
$$

## Passage en forme matricielle :

- Kalman n'est pas applicable à $(\overline{S})$ car il possède 3 équations
- posons $$
\textbf{x}_k:=\begin{bmatrix}x_k^l\\x_k^n\end{bmatrix}=\underbrace{\begin{bmatrix}\overline{F_k^l(x_{k-1}^n)}&0\\F_k^n(x_{k-1}^n)&0\end{bmatrix}}_{\textbf{F}_k(x_{k-1}^n):=}\begin{bmatrix}x_{k-1}^l\\x_{k-1}^n\end{bmatrix}+\underbrace{\begin{bmatrix}f_k^l(x_{k-1}^n)+\overline{b_k^{le}(x_k^n,x_{k-1}^n)}\\f_k^n(x_{k-1}^n)+b_k^{ne}\end{bmatrix}}_{\textbf{b}_k^e(x_k^n,x_{k-1}^n):=}+\underbrace{\begin{bmatrix}\overline{w_k^l}\\w_k^n\end{bmatrix}}_{\textbf{w}_k:=}
$$
- on remarque que $$
z_k=\underbrace{\begin{bmatrix}H_k(x_k^n)&0\end{bmatrix}}_{\textbf{H}_k(x_k^n):=}\begin{bmatrix}x_k^l\\x_k^n\end{bmatrix}+\underbrace{h_k(x_k^n)+b_k^o}_{\textbf{b}_k^o(x_k^n):=}+n_k
$$
- $(S)$ s'écrit donc sous forme matricielle comme suit $$
(\textbf{S}):\left\{\begin{matrix}
\textbf{x}_k&=&\textbf{F}_k(x_{k-1}^n)\:\textbf{x}_{k-1}&+&\textbf{b}_k^e(x_k^n,x_{k-1}^n)&+&\textbf{w}_k\\
z_k&=&\textbf{H}_k(x_k^n)\:\textbf{x}_k&+&\textbf{b}_k^o(x_k^n)&+&n_k
\end{matrix}\right.
$$
- sous cette forme, Kalman calcule également une estimation de l'état non-linéaire $x_k^n$, mais on l'ignore car c'est le rôle du filtre particulaire

## Estimation de Kalman :

- par souci de clarté, on notera $$
\left\{\begin{align}
\tilde{F}_k^l=&\overline{F_k^l(\hat{x}_{k-1|k-1}^n)}\\
\tilde{F}_k^n=&F_k^n(\hat{x}_{k-1|k-1}^n)\\
\tilde{H}_k=&H_k(\hat{x}_{k|k}^n)\\
\tilde{b}_k^{le}=&f_k^l(\hat{x}_{k-1|k-1}^n)+\overline{b_k^{le}(\hat{x}_{k|k-1}^n,\hat{x}_{k-1|k-1}^n)}\\
\tilde{b}_k^{ne}=&f_k^n(\hat{x}_{k-1|k-1}^n)+b_k^{ne}\\
\tilde{b}_k^o=&h_k(\hat{x}_{k|k}^n)+b_k^o\\
\end{align}\right.
$$
- prédiction
	- $$
	\begin{align}
	\begin{bmatrix}\hat{x}_{k|k-1}^l\\\hat{x}_{k|k-1}^n\end{bmatrix}=&\hat{\textbf{x}}_{k|k-1}=\textbf{F}_k(\hat{x}_{k-1|k-1}^n)\:\hat{\textbf{x}}_{k-1|k-1}+\textbf{b}_k^e(\hat{x}_{k|k-1}^n,\hat{x}_{k-1|k-1}^n)\\
	=&\begin{bmatrix}\tilde{F}_k^l&0\\\tilde{F}_k^l&0\end{bmatrix}\begin{bmatrix}\hat{x}_{k-1|k-1}^l\\\hat{x}_{k-1|k-1}^n\end{bmatrix}+\begin{bmatrix}\tilde{b}_k^{le}\\\tilde{b}_k^{ne}\end{bmatrix}\\
	=&\begin{bmatrix}\tilde{F}_k^l\:\hat{x}_{k-1|k-1}^l+\tilde{b}_k^{le}\\\tilde{F}_k^n\:\hat{x}_{k-1|k-1}^l+\tilde{b}_k^{ne}\end{bmatrix}
	\end{align}
	$$
	- $$
	\begin{align}
	\begin{bmatrix}P_{k|k-1}^l&P_{k|k-1}^{ln}\\(P_{k|k-1}^{ln})^T&P_{k|k-1}^n\end{bmatrix}=&\textbf{P}_{k|k-1}=\textbf{F}_k(\hat{x}_{k-1|k-1}^n)\:\textbf{P}_{k-1|k-1}\:\textbf{F}_k(\hat{x}_{k-1|k-1}^n)^T+\textbf{Q}_k\\
	=&\begin{bmatrix}\tilde{F}_k^l&0\\\tilde{F}_k^n&0\end{bmatrix}\begin{bmatrix}P_{k-1|k-1}^l&P_{k-1|k-1}^{ln}\\(P_{k-1|k-1}^{ln})^T&P_{k-1|k-1}^n\end{bmatrix}\begin{bmatrix}(\tilde{F}_k^l)^T&(\tilde{F}_k^n)^T\\0&0\end{bmatrix}+\begin{bmatrix}\overline{Q_k^l}&0\\0&Q_k^n\end{bmatrix}\\
	=&\begin{bmatrix}\tilde{F}_k^l\:P_{k-1|k-1}^l\:(\tilde{F}_k^l)^T+\overline{Q_k^l}&\tilde{F}_k^l\:P_{k-1|k-1}^l\:(\tilde{F}_k^n)^T\\\tilde{F}_k^n\:P_{k-1|k-1}^l\:(\tilde{F}_k^l)^T&\tilde{F}_k^n\:P_{k-1|k-1}^l\:(\tilde{F}_k^n)^T+Q_k^n\end{bmatrix}
	\end{align}
	$$
- correction
	- $$
	\begin{align}
	y_k=&z_k-\left(\textbf{H}_k(\hat{x}_{k|k}^n)\:\hat{\textbf{x}}_{k|k-1}+\textbf{b}_k^o(\hat{x}_{k|k}^n)\right)\\
	=&z_k-\left(\begin{bmatrix}\tilde{H}_k&0\end{bmatrix}\begin{bmatrix}\hat{x}_{k|k-1}^l\\\hat{x}_{k|k-1}^n\end{bmatrix}+\tilde{b}_k^o\right)\\
	=&z_k-\left(\tilde{H}_k\:\hat{x}_{k|k-1}^l+\tilde{b}_k^o\right)
	\end{align}
	$$
	- $$
	\begin{align}
	\begin{bmatrix}S_k^x&S_k^{xz}\\(S_k^{xz})^T&S_k^z\end{bmatrix}=&\textbf{S}_k=\textbf{H}_k(\hat{x}_{k|k}^n,\hat{x}_{k-1|k-1}^n)\:\textbf{P}_{k|k-1}\:\textbf{H}_k(\hat{x}_{k|k}^n,\hat{x}_{k-1|k-1}^n)^T+\textbf{R}_k\\
	=&\begin{bmatrix}0&\tilde{F}_k^n\\\tilde{H}_k&0\end{bmatrix}\begin{bmatrix}P_{k|k-1}&P_{k-1,k|k-1}\\P_{k-1,k|k-1}^T&P_{k-1|k-1}\end{bmatrix}\begin{bmatrix}0&\tilde{H}_k^T\\(\tilde{F}_k^n)^T&0\end{bmatrix}+\begin{bmatrix}Q_k^n&0\\0&R_k\end{bmatrix}\\
	=&\begin{bmatrix}\tilde{F}_k^n\:P_{k-1|k-1}\:(\tilde{F}_k^n)^T+Q_k^n&\tilde{F}_k^n\:P_{k-1,k|k-1}^T\:\tilde{H}_k^T\\\tilde{H}_k\:P_{k-1,k|k-1}\:(\tilde{F}_k^n)^T&\tilde{H}_k\:P_{k|k-1}\:\tilde{H}_k^T+R_k\end{bmatrix}
	\end{align}
	$$
	- $$
	\begin{align}
	\begin{bmatrix}\bar{A}&\bar{B}\\\bar{C}&\bar{D}\end{bmatrix}=&\textbf{S}_k^{-1}\\
	=&\begin{bmatrix}\tilde{F}_k^n\:P_{k-1|k-1}\:(\tilde{F}_k^n)^T+Q_k^n&\tilde{F}_k^n\:P_{k-1,k|k-1}^T\:\tilde{H}_k^T\\\tilde{H}_k\:P_{k-1,k|k-1}\:(\tilde{F}_k^n)^T&\tilde{H}_k\:P_{k|k-1}\:\tilde{H}_k^T+R_k\end{bmatrix}
	\end{align}
	$$
	- $K_k=P_{k|k-1}\:H_k^T\:S_k^{-1}$
	- $\hat{x}_{k|k}=\hat{x}_{k|k-1}+K_k\:y_k$
	- $P_{k|k}=(I-K_k\:H_k)\:P_{k|k-1}$
