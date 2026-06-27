- [[Gaussian Processes for Machine Learning Covariance Functions.pdf]]

## Noyau exponentiel quadratique :

- _squared exponential kernel_ en anglais : $$k(x, x') = \sigma_f^2 \: \text{exp}\left(-\frac{1}{2} \: (x'-x)^T \: M \: (x'-x)\right) \quad \text{avec} \quad M = \text{diag}(l)^{-2}$$
- $\sigma_f$ caractérise la corrélation maximale entre les valeurs prises par le GP en $x$ et en $x'$:
	- si $x$ et $x'$ sont proches, la valeur de $k(x, x')$ se rapproche de $\sigma_f^2$, signifiant que les valeurs de du GP en $x$ et en $x'$ sont fortement corrélées
	- si $x$ et $x'$ sont loins, $k(x, x') \rightarrow 0$, signifiant que les les valeurs de du GP en $x$ et en $x'$ ne sont presque pas corrélées
- le paramètre $l$ correspond à la distance caractéristique du noyau, qui donne l'ordre de grandeur de la "proximité" entre deux points

## Mixture de noyaux exponentiels quadratiques :

- on peut considérer une somme de noyaux exponentiels quadratiques de la forme $$k(x, x') = \sum_{i = 0}^n \: \sigma_{f_i}^2 \: \text{exp}\left(-\frac{1}{2} \: (x'-x)^T \: M_i \: (x'-x)\right) \quad \text{avec} \quad M_i = \text{diag}(l_i)^{-2}$$
- des $l_i$ de différentes échelles permettent de prendre en compte les fluctuations locales et la tendance globale des données

## Noyaux de Matérn

**TODO**
