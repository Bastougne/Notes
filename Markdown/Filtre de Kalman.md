## But :

- on souhaite estimer la densité _a posteriori_ $p(x_k|\textbf{z}_{1:k})$ à chaque pas de temps discret $k$
- on se place dans une approche paramétrique en supposant que la densité _a posteriori_ est gaussienne à chaque pas de temps
- on peut trouver une solution analytique récursive

## Hypothèses :

- $x_k=F_k\:x_{k-1}+b_k^e+w_k$ état courant
	- $F_k$ matrice d'évolution
	- $b_k^e$ biais d'évolution
	- $w_k\sim\mathcal{N}(0,Q_k)$ bruit de modèle ou d'évolution (_process noise_) gaussien
	- $Q_k$ matrice de covariance du bruit de modèle
- $z_k=H_k\:x_k+b_k^o+n_k$ observation
	- $H_k$ matrice d'observation
	- $b_k^o$ biais d'observation
	- $n_k\sim\mathcal{N}(0,R_k)$ bruit de mesure ou d'observation (_observation noise_) gaussien
	- $R_k$ matrice de covariance du bruit de mesure
- $x_0$ gaussien

## Algorithme :

- notation : $m|n$ signifie "l'état à l'instant $m$ sachant les observations aux instants $1$ à $n$"
![[Étapes du filtre de Kalman.png]]

**TODO : à démontrer**
- prédiction de l'état après déplacement par odométrie :
	- $\hat{x}_{k|k-1}=F_k\:\hat{x}_{k-1|k-1}+b_k^e$ estimateur de l'état courant
	- $P_{k|k-1}=F_k\:P_{k-1|k-1}\:F_k^T+Q_k$ matrice de covariance de l'erreur
- correction optimale de l'état (au sens des moindres carrés) (équations de Riccati) :
	- $y_k=z_k-\left(H_k\:\hat{x}_{k|k-1}+b_k^o\right)$ innovation
	- $S_k=H_k\:P_{k|k-1}\:H_k^T+R_k$ covariance de l'innovation
	- $K_k=P_{k|k-1}\:H_k^T\:S_k^{-1}$ gain de Kalman
	- $\hat{x}_{k|k}=\hat{x}_{k|k-1}+K_k\:y_k$ mise à jour de l'état
	- $P_{k|k}=(I-K_k\:H_k)\:P_{k|k-1}$ mise à jour de la covariance

## Forme de Joseph :

- à cause d'approximations numériques, $P_{k|k}$ n'est pas forcément SDP
- on cherche à reformuler pour garantir ces propriétés :$$
\begin{align}
P_{k|k}=&(I-K_k\:H_k)\:P_{k|k-1}\\
=&(I-K_k\:H_k)\:P_{k|k-1}-K_k\:S_k
\:K_k^T+K_k\:S_k\:K_k^T\\
=&(I-K_k\:H_k)\:P_{k|k-1}-(P_{k|k-1}\:H_k^T\:S_k^{-1})\:S_k
\:K_k^T+K_k\:(H_k\:P_{k|k-1}\:H_k^T+R_k)\:K_k^T\\
=&(I-K_k\:H_k)\:P_{k|k-1}-P_{k|k-1}\:H_k^T
\:K_k^T+K_k\:H_k\:P_{k|k-1}\:H_k^T\:K_k^T+K_k\:R_k\:K_k^T\\
=&(I-K_k\:H_k)\:P_{k|k-1}-(I-K_k\:H_k)\:P_{k|k-1}\:H_k^T\:K_k^T+K_k\:R_k\:K_k^T\\
=&(I-K_k\:H_k)\:P_{k|k-1}(I-K_k\:H_k)^T+K_k\:R_k\:K_k^T\\
\end{align}
$$
- sous cette forme plus robuste, $P_{k|k}$ est SDP par somme de 2 termes SDP

## Remarques :

- $P_{k|k-1}$, $K_k$ et $P_{k|k}$ ne dépendent pas des observations, on peut donc les précalculer pour gagner du temps si on est en temps réel
- optimal au sens des moindres carrés si évolution et commande linéaires et bruits gaussiens
- les covariances successives $P_{k|k}$ correspondent donc à la borne de Cramér-Rao dans le cas linéaire déterministe à bruit additif gaussien
- fortes hypothèses, souvent fausses en pratique
- [[Extended Kalman Filter]] permet de passer outre les non-linéarités
- [[Unscented Kalman Filter]]

## Limites :

- impose toujours l'hypothèse gaussienne même dans ses variantes
- gère mal les densités multimodales ou à queues lourdes
- perd son optimalité dans le cas de fonctions d'évolution ou d'observation non-linéaires
- a du mal avec la non-injectivité de la fonction d'observation
