## But :

- [[Sequential Importance Sampling]] simple à mettre en œuvre
- [[On Sequential Monte Carlo Sampling Methods for Bayesian Filtering.pdf]] montre que la variance des poids de particules échantillonnées par $q$ ne peut qu'augmenter avec le temps
- donc sujet à une dégénérescence de ses particules (_sample degeneracy_ en anglais) : le poids d'un sous-ensemble de particules tend vers 0

## Principe :

- on va donc rééchantillonner toutes les particules de sorte que :
	- une particule de poids fort soit transformée en plusieurs particules de poids plus faibles
	- une particule de poids faible soit éliminée
	- de façon générale, plus une particule est de poids important, et plus il y aura de nouvelles particules pour la remplacer
- plusieurs algorithmes de [[Rééchantillonnage]] possibles
**TODO : impact du rééchantillonnage ?**

## Simplifications :

- on simplifie l'algorithme en choisissant $q\left(x_k|x_{k-1}^{(i)},z_k\right)$ comme étant la densité de transition a priori $p\left(x_k|x_{k-1}^{(i)}\right)$
- puisqu'on rééchantillonne à chaque itération, $w_{k-1}^{(i)} = 1/N$, et donc $$
\begin{align}
w_k^{(i)}\propto&\:w_{k-1}^{(i)}\:\frac{\tilde{p}\left(z_k|x_k^{(i)}\right)\:\tilde{p}\left(x_k^{(i)}|x_{k-1}^{(i)}\right)}{\tilde{q}\left(x_k^{(i)}|x_{k-1}^{(i)},z_k\right)}\\
\propto&\:\frac{1}{N}\:\frac{\tilde{p}\left(z_k|x_k^{(i)}\right)\:\tilde{p}\left(x_k^{(i)}|x_{k-1}^{(i)}\right)}{\tilde{p}\left(x_k^{(i)}|x_{k-1}^{(i)}\right)}\\
\propto&\:\tilde{p}\left(z_k|x_k^{(i)}\right)
\end{align}
$$

## Algorithme :

- $$\left(\left\{x_k^{(i)},w_k^{(i)}\right\}_{i=1}^N\right)=SIR\left(\left\{x_{k-1}^{(i)},w_{k-1}^{(i)}\right\}_{i=1}^N,z_k\right)$$
	- $$\left(\left\{x_k^{(i)},w_k^{(i)}\right\}_{i=1}^N\right)=SIS\left(\left\{x_{k-1}^{(i)}, w_{k-1}^{(i)}\right\}_{i=1}^N,z_k\right)$$
	- $$\left(\left\{x_k^{(i)},w_k^{(i)},\_\right\}_{i=1}^N\right)=Resampling\left(\left\{x_k^{(i)}, w_k^{(i)}\right\}_{i=1}^N\right)$$

## Résultats :

- pour T = 600, N_s = 500 et N_e = 1000 :
	- taux de succès = 16.60 %, RMSE = 0.38 km
	- appauvrissement des particules : ![[exemple SIR 1.png]]
	- appauvrissement + divergence : ![[exemple SIR 2.png]]
	- appauvrissement + divergence : ![[exemple SIR 3.png]]

## Remarques :

- les hypothèses du SIR sont très faibles :
	- fonction d'état $f$ et de mesure $h$ connues
	- possibilité d'échantillonner depuis $p\left(\cdot|x_{k-1}^{(i)}\right)$
