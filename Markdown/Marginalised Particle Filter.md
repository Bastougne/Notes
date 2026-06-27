## But :

- aussi appelé filtre Rao-Blackwellisé (_RBPF_ en anglais)
- le filtrage particulaire est sous-optimal par rapport au filtre de Kalman
- le filtre de Kalman impose des hypothèses de linéarité sur les fonctions de modèle et d'observation
- si ces fonctions sont linéaires en un sous-vecteur de l'état, on peut utiliser Kalman pour trouver une estimation optimale de ce sous-vecteur et estimer la partie non-linéaire par filtrage particulaire
- [[Marginalized Particle Filters for Mixed Linear Nonlinear State-Space Models.pdf]]
- [[A New Formulation of the Rao-Blackwellised Particle Filter.pdf]]
- [[Marginalized Particle Filter for Accurate and Reliable Terrain-Aided Navigation.pdf]]

## Principe :

- on peut alors marginaliser l'état : $$
p(x_k^l,x_k^n|z_k)=\underbrace{p(x_k^l|x_k^n, z_k)}_{\text{estimé par KF}} \: \underbrace{p(x_k^n|z_k)}_{\text{estimé par PF}}
$$
- on a donc deux filtres entrelacés :
	- le filtre particulaire traite les non-linéarités et les multimodalités dans un sous-espace de l'espace d'état
	- le filtre de Kalman traite optimalement le reste en prenant comme hypothèses les estimations du filtre particulaire
- réduire la dimension de l'espace traité par le filtre particulaire (sous-optimal) permet d'améliorer les performances globales du filtre et de réduire le nombre de particules nécessaires
- cependant, il faut un filtre de Kalman par particule en principe
- [[Équations du RBPF]]

## Algorithme :

- $$\left(\left\{\begin{bmatrix}{x_k^n}^{(i)}\\{x_k^l}^{(i)}\end{bmatrix},w_k^{(i)},{P_k^l}^{(i)}\right\}_{i=1}^N\right)=RBPF\left(\left\{\begin{bmatrix}{x_{k-1}^n}^{(i)}\\{x_{k-1}^l}^{(i)}\end{bmatrix},w_{k-1}^{(i)},{P_{k-1}^l}^{(i)}\right\}_{i=1}^N,z_k\right)$$
	- POUR $i=1:N$
		- propagation particulaire : $${x_k^n}^{(i)}\sim q\left(\cdot|{x_{k-1}^n}^{(i)},z_k\right)$$
		- pré-correction de Kalman : $$\left\{\begin{align}
		K_k^{(i)}=&{P_{k-1}^l}^{(i)}\:(\tilde{F}_k^n)^T\:\left(\tilde{F}_k^n\:{P_{k-1}^l}^{(i)}\:(\tilde{F}_k^n)^T+Q_k^n\right)^{-1}\\
		{{x_{k-1}^l}^*}^{(i)}=&{x_{k-1}^l}^{(i)}+K_k^{(i)}\:\left({x_k^n}^{(i)}-\left(\tilde{F}_k^n\:{x_{k-1}^l}^{(i)}+\tilde{b}_k^{ne}\right)\right)\\
		{{P_{k-1}^l}^*}^{(i)}=&(I-K_k^{(i)}\:\tilde{F}_k^n)\:{P_{k-1}^l}^{(i)}
		\end{align}\right.$$
		- prédiction de Kalman : $$\left\{\begin{align}
		{x_{k|k-1}^l}^{(i)}=&\tilde{F}_k^l\:{{x_{k-1}^l}^*}^{(i)}+\tilde{b}_k^{le}\\
		{P_{k|k-1}^l}^{(i)}=&\tilde{F}_k^l\:{{P_{k-1}^l}^*}^{(i)}\:(\tilde{F}_k^l)^T+\overline{Q_k^l}\end{align}\right.$$
		- correction particulaire : $$\tilde{w}_k^{(i)}=w_{k-1}^{(i)}\:\frac{\tilde{p}\left(z_k|{x_k^n}^{(i)}\right)\:\tilde{p}\left({x_k^n}^{(i)}|{x_{k-1}^n}^{(i)}\right)}{\tilde{q}\left({x_k^n}^{(i)}|{x_{k-1}^n}^{(i)},z_k\right)}$$
		- correction de Kalman : $$\left\{\begin{align}
		L_k^{(i)}=&{P_{k|k-1}^l}^{(i)}\:\tilde{H}_k^T\:\left(\tilde{H}_k\:{P_{k|k-1}^l}^{(i)}\:\tilde{H}_k^T+R_k\right)^{-1}\\
		{x_k^l}^{(i)}=&{x_{k|k-1}^l}^{(i)}+L_k^{(i)}\:\left(z_k-\left(\tilde{H}_k\:{x_{k|k-1}^l}^{(i)}+\tilde{b}_k^o\right)\right)\\
		{P_k^l}^{(i)}=&(I-L_k^{(i)}\:\tilde{H}_k)\:{P_{k|k-1}^l}^{(i)}
		\end{align}\right.$$
	- FIN POUR
	- POUR $i=1:N$
		- normalisation : $$w_k^{(i)}=\frac{\tilde{w}_k^{(i)}}{\displaystyle{\sum_{j=1}^N\tilde{w}_k^{(j)}}}$$
	- FIN POUR
- on peut naturellement incorporer des étapes de différents types de filtres particulaires ou de Kalman aux étapes correspondantes du RBPF
