## But :

- le rééchantillonnage du filtre particulaire usuel a été proposé pour réduire la dégénérescence particulaire
- mais il induit une perte de diversité particulaire (_particle impoverishment_ en anglais)
- [[Improving Regularised Particle Filters.pdf]] introduit le RPF comme piste de solution pour éviter ce problème

## Principe :

- l'approximation discrète $$
p(\textbf{x}_{0:k}|\textbf{z}_{1:k})\:d\textbf{x}_{0:k}\approx\left(S^N_\delta(\mu_{\textbf{z}_{1:k}},\nu_{\textbf{z}_{1:k}})\right)(d\textbf{x}_{0:k})=\sum_{i=1}^N w_k^{(i)}\delta_{\textbf{X}_{0:k}^{(i)}}(d\textbf{x}_{0:k})
$$ limite les capacités d'exploration du filtre, et est à l'origine de la perte de diversité
- on va approcher la densité par une mixture de noyaux régularisés $S^N_{K_h}$ telle que : $$
\frac{dS^N_{K_h}(\mu_{\textbf{z}_{1:k}},\nu_{\textbf{z}_{1:k}})}{d\textbf{x}_{0:k}}=\sum_{i=1}^N w_k^{(i)}{K_h}\left(\textbf{x}_{0:k}-\textbf{X}_{0:k}^{(i)}\right)
$$
- par propriétés de la convolution, on remarque que $$
S^N_{K_h}(\mu_{\textbf{z}_{1:k}},\nu_{\textbf{z}_{1:k}})=K_h\ast S^N_\delta(\mu_{\textbf{z}_{1:k}},\nu_{\textbf{z}_{1:k}})
$$
- l'approximation aura un support de mesure non-nulle et sera donc plus fidèle à la vraie densité

![[Régularisation par noyau.png]]

## Noyau de régularisation :

- on remplace la mesure de Dirac $\delta$ par un noyau de régularisation $K$ vérifiant :
	- $K$ est une densité de probabilité $$
	K\underset{\text{p.p.}}{\geq}0\quad\text{et}\quad\int{K(x)\:dx=1}$$
	- $K$ est homogène isotrope $$
	\exists\:\phi:\mathbb{R}_+\rightarrow\mathbb{R}\:|\:\forall x\in\mathbb{R}^d,\:K(x)=\phi(||x||_2)
	$$ donc d'espérance nulle $$
	\int{x\:K(x)\:dx}=0
	$$
	- $K$ admet un moment d'ordre 2 $$
	\int{x^2\:K(x)\:dx}=0
	$$
- $K$ est ensuite régularisé : $$
\forall x\in\mathbb{R}^d,\:K_h(\textbf{x}):=\frac{1}{h^d}K\left(\frac{x}{h}\right)
$$ avec $h$ le facteur de dilatation (_bandwidth_ en anglais)
- si $h$ est petit, $K_h$ est proche d'un Dirac
- si $h$ est grand, $K_h$ se rapproche d'une densité uniforme
- [[Choix du noyau de régularisation du RPF]]

## RPF :

- on peut effectuer la convolution de régularisation avant ou après l'étape de correction
- il y a donc 2 versions du RPF :
	- RPF post-régularisé (le plus courant car plus simple à implémenter) : $$
	\begin{matrix*}[ccl]
	S^N_\delta(\mu_{\textbf{z}_{1:k-1}},\nu_{\textbf{z}_{1:k-1}})&\underset{\text{correction}}{\longrightarrow}&S^N_\delta(\mu_{\textbf{z}_{1:k}},\nu_{\textbf{z}_{1:k}})=\frac{d\mu_{\textbf{z}_{1:k}}}{d\nu_{\textbf{z}_{1:k}}}S^N_\delta(\nu_{\textbf{z}_{1:k}})\\
	&\underset{\text{régularisation}}{\longrightarrow}&S^N_{K_h}(\mu_{\textbf{z}_{1:k}},\nu_{\textbf{z}_{1:k}})=K_h\ast S^N_\delta(\mu_{\textbf{z}_{1:k}},\nu_{\textbf{z}_{1:k}})
	\end{matrix*}
	$$

## Pre-RPF :
 (**TODO à vérifier**)
- $$\left(\left\{x_k^{(i)}, w_k^{(i)}\right\}_{i=1}^N\right) = Pre\textit{-}RPF\left(\left\{x_{k-1}^{(i)}, w_{k-1}^{(i)}\right\}_{i=1}^N,z_k\right)$$
	- $S_k = \textbf{Cov}\left(\textbf{x}_k^{(1:N)}\right)$
	- $D_k$ tel que $D_k \: D_k^T = S_k$
	- POUR $i = 1:N$ 
		- $$x_k^{(i)} \sim q\left(x_k|x_{k-1}^{(i)},z_k\right)$$
		- $$\varepsilon^{(i)} \sim K$$
		- $$x_k^{(i)} = x_k^{(i)} + h \: D_k \: \varepsilon^{(i)}$$
		- $$w_k^{(i)} = w_{k-1}^{(i)} \: \frac{\tilde{p}\left(z_k|x_k^{(i)}\right) \: \tilde{p}\left(x_k^{(i)}|x_{k-1}^{(i)}\right)}{\tilde{q}\left(x_k^{(i)}|x_{k-1}^{(i)},z_k\right)}$$
	- FIN POUR
	- POUR $i = 1:N$
		- $$w_k^{(i)} = \frac{w_k^{(i)}}{\displaystyle{\sum_{j=1}^Nw_k^{(j)}}}$$
	- FIN POUR
	- $$\widehat{N_{eff}} = \frac{1}{\displaystyle{\sum_{i=1}^N\left(w_k^{(i)}\right)^2}}$$
	- SI $\widehat{N_{eff}} < N_{seuil}$
		- $$\left(\left\{x_k^{(i)}, w_k^{(i)},\_\right\}_{i=1}^N\right) = Resampling\left(\left\{x_k^{(i)}, w_k^{(i)}\right\}_{i=1}^N\right)$$
	- FIN SI

## Post-RPF

- $$\left(\left\{x_k^{(i)}, w_k^{(i)}\right\}_{i=1}^N\right) = Post\textit{-}RPF\left(\left\{x_{k-1}^{(i)}, w_{k-1}^{(i)}\right\}_{i=1}^N,z_k\right)$$
	- $$\left(\left\{x_k^{(i)}, w_k^{(i)}\right\}_{i=1}^N\right) = SIS\left(\left\{x_{k-1}^{(i)}, w_{k-1}^{(i)}\right\}_{i=1}^N,z_k\right)$$
	- $$\widehat{N_{eff}} = \frac{1}{\displaystyle{\sum_{i=1}^N\left(w_k^{(i)}\right)^2}}$$
	- SI $\widehat{N_{eff}} < N_{seuil}$
		- $$S_k = \sum_{i=1}^N\left(w_k^{(i)}\left(x_k^{(i)} - \sum_{j=1}^Nw_k^{(j)}x_k^{(j)}\right)\left(x_k^{(i)} - \sum_{j=1}^Nw_k^{(j)}x_k^{(j)}\right)^T\right)$$
		- $$D_k \text{ tel que } D_k \: D_k^T = S_k \quad \text{(factorisation de Cholesky)}$$
		- $$\left(\left\{x_k^{(i)}, w_k^{(i)},\_\right\}_{i=1}^N\right) = Resampling\left(\left\{x_k^{(i)}, w_k^{(i)}\right\}_{i=1}^N\right)$$
		- POUR $i = 1:N$
			- $$\varepsilon^{(i)} \sim K$$
			- $$x_k^{(i)} = x_k^{(i)} + h \: D_k \: \varepsilon^{(i)}$$
		- FIN POUR
	- FIN SI

## Résultats :

- pour T = 600, N_s = 500, N_e = 1000, N_seuil = 0.5 * N_s, mu = 0.3 :
	- taux de succès = 84.00%, RMSE = 4.63km
	- convergence : ![[PostRPF MATLAB 1.png]]
	- convergence : ![[PostRPF MATLAB 2.png]]
