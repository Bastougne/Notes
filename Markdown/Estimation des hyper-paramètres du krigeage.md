## Estimation par maximum de vraisemblance :

- [[Gaussian Processes for Machine Learning Model Selection and Adaptation of Hyperparameters.pdf]]
- posons $\theta$ les hyper-paramètres du processus gaussien
- notons $K_\theta = k_\theta(X, X)$ pour simplifier les notations
- puisque $f \sim \mathcal{GP}(\mu, k)$, on a, pour $z := f(x)$, : $$
p( Z \: | \: X, \theta ) = \frac{
\displaystyle{\text{exp}\left(-\frac{1}{2} \: (Z - \mu_\theta(X))^T \: K_\theta^{-1} \: (Z - \mu_\theta(X))\right)}
}{
\displaystyle{\sqrt{{(2 \pi)}^n \: \text{det}(K_\theta)}}}
$$
- puisqu'on n'a pas d'_a priori_ sur le modèle, on suppose que $\mu_\theta$ est identiquement nul
- la log-vraisemblance s'écrit alors : $$
\begin{align}
l(\theta, Z, X) & = \text{log}(p( Z \: | \: X, \theta)) \\
& = -\frac{1}{2} \: Z^T \: K_\theta^{-1} \: Z -
\frac{1}{2} \text{log}\left({(2 \pi)}^n \: \text{det}(K_\theta)\right) \\
& = -\frac{1}{2} \left(
Z^T \: K_\theta^{-1} \: Z + \text{log}(\text{det}(K_\theta))
\right) -
\frac{n}{2} \text{log}(2 \pi)
\end{align}
$$
- les conditions nécessaires du maximum de vraisemblance sur chaque coordonnée de $\theta$ s'écrivent : $$
\begin{align}
& & \frac{\partial}{\partial \theta_i} l(\theta, Z, X) & = 0 \\
\Leftrightarrow & & -\frac{1}{2} \: \frac{\partial}{\partial \theta_i}\left(
Z^T \: K_\theta^{-1} \: Z + \text{log}(\text{det}(K_\theta))
\right) & = 0 \\
\Leftrightarrow & & - Z^T \: \frac{\partial}{\partial \theta_i}\left(K_\theta^{-1}\right) 
\: Z - \frac{\frac{\partial}{\partial \theta_i}(\text{det}(K_\theta))}{\text{det}(K_\theta)} & = 0 \\
\Leftrightarrow & & Z^T \: K_\theta^{-1} \: \frac{\partial}{\partial \theta_i}(K_\theta) \: K_\theta^{-1} \: Z - \frac{\frac{\partial}{\partial \theta_i}(\text{det}(K_\theta))}{\text{det}(K_\theta)} & = 0 & (1) \\
\Leftrightarrow & & Z^T \: K_\theta^{-1} \: \frac{\partial}{\partial \theta_i}(K_\theta) \: K_\theta^{-1} \: Z - \text{tr}\left(
K_\theta^{-1} \: \frac{\partial}{\partial \theta_i}(K_\theta)
\right) & = 0 & (2) \\
\Leftrightarrow & & \text{tr}\left(Z^T \: K_\theta^{-T} \: \frac{\partial}{\partial \theta_i}(K_\theta) \: K_\theta^{-1} \: Z -
K_\theta^{-1} \: \frac{\partial}{\partial \theta_i}(K_\theta)
\right) & = 0 \\
\Leftrightarrow & & \text{tr}\left(\left(\left(K_\theta^{-1} \: Z\right) \: \left(K_\theta^{-1} \: Z\right)^T - K_\theta^{-1}\right)
\: \frac{\partial}{\partial \theta_i}(K_\theta) \right) & = 0 \\
\end{align}
$$
- la ligne $(1)$ se démontre par composition avec la formule de la différentielle de l'inverse matricielle pour une matrice inversible $$
D(\cdot ^{-1})(X) : H \mapsto - X^{-1} \: H \: X^{-1}
$$
- la ligne  $(2)$ se démontre avec la [formule de Jacobi](https://en.wikipedia.org/wiki/Jacobi%27s_formula) dans le cas d'une matrice inversible
