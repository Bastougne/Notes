## But :

- on dispose d'une chaîne de Markov cachée $(x_k)_{k\geq0}$ à valeurs dans $\mathbb{R}^p$ et d'une suite d'observations $(z_k)_{k\geq0}$ à valeurs dans $\mathbb{R}^q$ ($q\leq p$) associées une à une aux états inconnus $x_k$
- formellement, il s'agit de remonter aux $x_k$ ou à leur image par une fonction $\varphi(x_k)$ dans le système $$
\left\{\begin{aligned}
x_k=&f_k(x_{k-1},u_k,w_k)\\
z_k=&h_k(x_k,n_k)
\end{aligned}\right.
$$ en connaissant $$
\left\{\begin{aligned}
p_0&\quad\text{la densité de probabilité de l'état initial }x_0\\
f_k&\quad\text{la fonction d'évolution}\\
u_k&\quad\text{la commande courante}\\
p_k^{(w)}&\quad\text{la densité de probabilité du bruit d'évolution }w_k\\
z_k&\quad\text{l'observation courante bruitée}\\
h_k&\quad\text{la fonction d'observation}\\
p_k^{(n)}&\quad\text{la densité de probabilité du bruit d'observation }n_k
\end{aligned}\right.
$$
- il faut donc s'aider des observations partielles, car $h_k$ est non-injective en générale, et bruitées $z_k$ pour déterminer les $x_k$
- 3 problèmes :
	- filtrage : estimer "au mieux" $\varphi(x_k)$ connaissant $z_{1:k}$
	- lissage : estimer "au mieux" $\varphi(x_{0:k})$ connaissant $z_{1:k}$
	- prédiction : estimer "au mieux" $\varphi(x_l)$ connaissant $z_{1:k}$ avec $l>k$

## Estimation :

- que signifie estimer "au mieux" une variable aléatoire $X$ connaissant la réalisation d'une autre variable aléatoire $Y$ ?
- c'est trouver un estimateur $\psi(Y)$ de $\varphi(X)$ dans un espace métrique fonctionnel $(F,d)$ qui minimise $$
\underset{\psi\in F}{\text{argmin}}(d(\psi(Y)-\varphi(X)))
$$
- il existe plusieurs [[Estimateurs intéressants]]
- la plupart de ces estimateurs se basent sur la densité de probabilité _a posteriori_
- en toute généralité, on va donc chercher à estimer la densité _a posteriori_ $p(x_k|\textbf{z}_{1:k})$ à chaque pas de temps discret $k$


## Algorithme :

- le filtre optimal permet son estimation récursive en 2 étapes :
	- prédiction avec l'équation de Chapman-Kolmogorov : $$
	\begin{align}
	p(x_k|\textbf{z}_{1:k-1})=&\int p(x_k,x_{k-1}|\textbf{z}_{1:k-1})\:dx_{k-1}\\
	=&\int p(x_k|x_{k-1},\textbf{z}_{1:k-1})\:p(x_{k-1}|\textbf{z}_{1:k-1})\:dx_{k-1}\\
	=&\int p(x_k|x_{k-1})\:p(x_{k-1}|\textbf{z}_{1:k-1})\:dx_{k-1}\\
	\end{align}
	$$
	- correction par théorème de Bayes : $$
	p(x_k|\textbf{z}_{1:k})=\frac{p(z_k|x_k)\:p(x_k|\textbf{z}_{1:k-1})}{\displaystyle\int p(z_k|x_k)\:p(x_k|\textbf{z}_{1:k-1})\:dx_k}
	$$
- en pratique, le calcul de l'équation de Chapman-Kolmogorov n'est pas réalisable de manière analytique

## Borne de Cramér-Rao

- [[Borne de Cramér-Rao]]
