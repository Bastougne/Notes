## But :

- on souhaite calculer $$
\mathbb{E}_{p(\cdot|\textbf{z}_{1:k})}(\varphi(\textbf{X}_{0:k}))=\int\varphi(\textbf{x}_{0:k})\:p(\textbf{x}_{0:k}|\textbf{z}_{1:k})\:d\textbf{x}_{0:k}
$$ avec $\textbf{x}_{0:k}$ la trajectoire de l'état et $\textbf{z}_{1:k}$ les observations successives
-  la [[Méthode de Monte-Carlo pondérée]] semble particulièrement adaptée au vu de sa flexibilité et de la faiblesse de ses hypothèses
- elle permet en fait d'obtenir une approximation de toute la densité $\textbf{x}_{0:k}\mapsto p(\textbf{x}_{0:k}|\textbf{z}_{1:k})$ sous la forme de la distribution empirique $S^N_\delta(\mu_{\textbf{z}_{1:k}},\nu_{\textbf{z}_{1:k}})$

## Hypothèses :

- approximation de Monte-Carlo : $$
p(\textbf{x}_{0:k}|\textbf{z}_{1:k})\:d\textbf{x}_{0:k}\approx\left(S^N_\delta(\mu_{\textbf{z}_{1:k}},\nu_{\textbf{z}_{1:k}})\right)(d\textbf{x}_{0:k})=\sum_{i=1}^N w_k^{(i)}\delta_{\textbf{X}_{0:k}^{(i)}}(d\textbf{x}_{0:k})
$$
- états markoviens : $p(x_k|\textbf{x}_{0:k-1},\textbf{z}_{1:k-1})=p(x_k|\textbf{x}_{0:k-1})=p(x_k|x_{k-1})$
- observations : $p(z_k|\textbf{x}_{0:k},\textbf{z}_{1:k-1})=p(z_k|\textbf{x}_{0:k})=p(z_k|x_k)$
- choix de $q$ : $q(\textbf{x}_{0:k-1}|\textbf{z}_{1:k})=q(\textbf{x}_{0:k-1}|\textbf{z}_{1:k-1})$
	- souvent : $q(x_k|\textbf{x}_{0:k-1},\textbf{z}_{1:k})=q(x_k|x_{k-1},z_k)$

## Relations de récursivité :

- formulation récursive de $p$ : $$
\begin{align}
p(\textbf{x}_{0:k}|\textbf{z}_{1:k})&=p(\textbf{x}_{0:k}|z_k,\textbf{z}_{1:k-1})\\
&=\frac{p(z_k|\textbf{x}_{0:k},\textbf{z}_{1:k-1})\:p(\textbf{x}_{0:k}|\textbf{z}_{1:k-1})}{p(z_k|\textbf{z}_{1:k-1})}\\
&=\frac{p(z_k|\textbf{x}_{0:k},\textbf{z}_{1:k-1})\:p(x_k,\textbf{x}_{0:k-1}|\textbf{z}_{1:k-1})}{p(z_k|\textbf{z}_{1:k-1})}\\
&=\frac{p(z_k|\textbf{x}_{0:k},\textbf{z}_{1:k-1})\:p(x_k|\textbf{x}_{0:k-1},\textbf{z}_{1:k-1})\:p(\textbf{x}_{0:k-1}|\textbf{z}_{1:k-1})}{p(z_k|\textbf{z}_{1:k-1})}\\
&=\frac{p(z_k|x_k)\:p(x_k|\textbf{x}_{0:k-1},\textbf{z}_{1:k-1})\:p(\textbf{x}_{0:k-1}|\textbf{z}_{1:k-1})}{p(z_k|\textbf{z}_{1:k-1})}\\
&=\frac{p(z_k|x_k)\:p(x_k|x_{k-1})\:p(\textbf{x}_{0:k-1}|\textbf{z}_{1:k-1})}{p(z_k|\textbf{z}_{1:k-1})}\\
&\propto p(z_k|x_k)\:p(x_k|x_{k-1})\:p(\textbf{x}_{0:k-1}|\textbf{z}_{1:k-1})
\end{align}
$$
- formulation récursive de $q$ : $$
\begin{align}
q(\textbf{x}_{0:k}|\textbf{z}_{1:k})&=q(x_k,\textbf{x}_{0:k-1}|\textbf{z}_{1:k})\\
&=q(x_k|\textbf{x}_{0:k-1},\textbf{z}_{1:k})\:q(\textbf{x}_{0:k-1}|\textbf{z}_{1:k})\\
&=q(x_k|\textbf{x}_{0:k-1},\textbf{z}_{1:k})\:q(\textbf{x}_{0:k-1}|\textbf{z}_{1:k-1})\\
\end{align}
$$
- remarque : on a alors $$
\begin{align}
q(\textbf{x}_{0:k}|\textbf{z}_{1:k})&=q(\textbf{x}_{0:k-1}|\textbf{z}_{1:k-1})\:q(x_k|\textbf{x}_{0:k-1},\textbf{z}_{1:k})\\
&=q(x_0)\:\prod_{l=1}^kq(x_l|\textbf{x}_{0:l-1},\textbf{z}_{1:l})
\end{align}
$$
- formulation récursive des poids $w_k^{(i)}$ : $$
\begin{align}
w_k^{(i)}&\propto\frac{\tilde{p}\left(\textbf{x}_{0:k}^{(i)}|\textbf{z}_{1:k}\right)}{\tilde{q}\left(\textbf{x}_{0:k}^{(i)}|\textbf{z}_{1:k}\right)}\\
&\propto\frac{\tilde{p}\left(z_k|x_k^{(i)}\right)\:\tilde{p}\left(x_k^{(i)}|x_{k-1}^{(i)}\right)\:\tilde{p}\left(\textbf{x}_{0:k-1}^{(i)}|\textbf{z}_{1:k-1}\right)}{\tilde{q}\left(x_k^{(i)}|\textbf{x}_{0:k-1}^{(i)},\textbf{z}_{1:k}\right)\:\tilde{q}\left(\textbf{x}_{0:k-1}^{(i)}|\textbf{z}_{1:k-1}\right)}\\
&\propto w_{k-1}^{(i)}\:\frac{\tilde{p}\left(z_k|x_k^{(i)}\right)\:\tilde{p}\left(x_k^{(i)}|x_{k-1}^{(i)}\right)}{\tilde{q}\left(x_k^{(i)}|\textbf{x}_{0:k-1}^{(i)},\textbf{z}_{1:k}\right)}&\text{en général}\\
&\propto w_{k-1}^{(i)}\:\frac{\tilde{p}\left(z_k|x_k^{(i)}\right)\:\tilde{p}\left(x_k^{(i)}|x_{k-1}^{(i)}\right)}{\tilde{q}\left(x_k^{(i)}|x_{k-1}^{(i)},z_k\right)}&\text{si }\tilde{q}\text{ simplifié}
\end{align}
$$

## SIS :

- on a trouvé une relation récursive entre les poids à l'instant $k$ et ceux à l'instant $k-1$
- connaissant la densité initiale $x_0\mapsto p(x_0)$, on peut donc approcher récursivement la densité $\textbf{x}_{0:k}\mapsto p(\textbf{x}_{0:k}|\textbf{z}_{1:k})$ en 2 étapes :
	- propagation des particules $x_k^{(i)}$
	- correction des poids $w_k^{(i)}$
- en particulier, en notant $\pi_k$ la projection d'un $(k+1)$-uplet sur son dernier élément $$
\pi_k:\left\{\begin{align}
\mathbb{R}^{k+1}&\rightarrow\mathbb{R}\\
(x_0,\dots,x_k)&\mapsto x_k
\end{align}\right.
$$ une estimation de l'espérance _a posteriori_ est $$
\hat{x}_k:=\left\langle S^N_\delta(\mu_{\textbf{z}_{1:k}},\nu_{\textbf{z}_{1:k}}),\varphi\circ\pi_k\right\rangle=\sum_{i=1}^Nw_k^{(i)}\:\varphi\left(x_k^{(i)}\right)\approx\mathbb{E}_{p(\cdot|\textbf{z}_{1:k})}(\varphi(X_k))
$$
- la matrice de covariance de $\hat x_k$ est donnée par $$
P_k:=\sum_{i=1}^Nw_k^{(i)}\left(\varphi\left(x_k^{(i)}\right)-\hat x_k\right)\left(\varphi\left(x_k^{(i)}\right)-\hat x_k\right)^T
$$

## Algorithme :

- $$\left(\left\{x_k^{(i)},w_k^{(i)}\right\}_{i=1}^N\right)=SIS\left(\left\{x_{k-1}^{(i)},w_{k-1}^{(i)}\right\}_{i=1}^N,z_k\right)$$
	- POUR $i=1:N$
		- propagation : $$x_k^{(i)}\sim q\left(\cdot|x_{k-1}^{(i)},z_k\right)$$
		- correction : $$\tilde{w}_k^{(i)}=w_{k-1}^{(i)}\:\frac{\tilde{p}\left(z_k|x_k^{(i)}\right)\:\tilde{p}\left(x_k^{(i)}|x_{k-1}^{(i)}\right)}{\tilde{q}\left(x_k^{(i)}|x_{k-1}^{(i)},z_k\right)}$$
	- FIN POUR
	- POUR $i=1:N$
		- normalisation : $$w_k^{(i)}=\frac{\tilde{w}_k^{(i)}}{\displaystyle{\sum_{j=1}^N\tilde{w}_k^{(j)}}}$$
	- FIN POUR

## Résultats :

- pour T=600 itérations du système, N_s=500 particules et N_e=1000 itérations de Monte-Carlo :
	- taux de succès=13.80 %, RMSE=0.67 km
	- dégénérescence des poids : ![[SIS MATLAB 1.png]]
	- dégénérescence + divergence : ![[SIS MATLAB 2.png]]

## Remarques :

- [[Correction numérique de la matrice de covariance]]
- pour que l'algorithme fonctionne bien, on veut un fort bruit de modèle, et un faible bruit de mesure
