## But :

- estimer la densité de probabilité _a posteriori_ $p(\textbf{x}_{0:k}|\textbf{z}_{1:k})$ à chaque pas de temps $k$ de façon non-paramétrique
- méthode de Monte-Carlo par chaînes de Markov pour le filtrage non-linéaire et non-gaussien
- [[A Tutorial on Particle Filters for Online Nonlinear or Non-Gaussian Bayesian Tracking.pdf]]
- [[Filtrage Bayésien et Approximation Particulaire.pdf]]
- [[An Introduction to Sequential Monte Carlo Methods.pdf]]

## Principe :

- le problème de dégénérescence adressé par [[Sampling Importance Resampling]] ne se présente pas à chaque itération
- on peut donc se contenter de rééchantillonner lorsqu'un [[Critères de rééchantillonnage]] est validé
- en pratique, on utilise le SIR comme un filtre de démarrage (_bootstrap filter_ en anglais) car simple et peu restrictif
	
## Algorithme :

- $$\left(\left\{x_k^{(i)}, w_k^{(i)}\right\}_{i=1}^N\right) = PF\left(\left\{x_{k-1}^{(i)}, w_{k-1}^{(i)}\right\}_{i=1}^N,z_k\right)$$
	- $$\left(\left\{x_k^{(i)}, w_k^{(i)}\right\}_{i=1}^N\right) = SIS\left(\left\{x_{k-1}^{(i)}, w_{k-1}^{(i)}\right\}_{i=1}^N,z_k\right)$$
	- SI critère de rééchantillonnage
		- $$\left(\left\{x_k^{(i)}, w_k^{(i)},\_\right\}_{i=1}^N\right) = Resampling\left(\left\{x_k^{(i)}, w_k^{(i)}\right\}_{i=1}^N\right)$$
	- FIN SI

## Résultats :

- pour T = 600, N_s = 500, N_e = 1000, N_seuil = 0.5 * N_s :
	- taux de succès = 26.90%, RMSE = 0.34 km
	- convergence : ![[PF MATLAB 1.png]]
	- appauvrissement des particules + divergence : ![[PF MATLAB 2.png]]
	- appauvrissement + divergence : ![[PF MATLAB 3.png]]

## Remarques :

- [[Choix de la fonction d'importance du SIS]]
- on a résolu le problème de dégénérescence particulaires
- mais les nouvelles particules sont crées à l'endroit des anciennes particules de poids fort
- donc perte de diversité particulaires (appauvrissement) : ![[Particle impoverishment.png]]
- pistes d'améliorations :
	- [[Regularised Particle Filter]]
	- [[Auxiliary Particle Filter]]
	- [[Marginalised Particle Filter]]
