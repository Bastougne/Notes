## Densité a posteriori :

**TODO : à réécrire**
- prenons $q(x_k|\textbf{x}_{0:k-1},\textbf{z}_{1:k}) = p(x_k|x_{k-1},z_k)$ la densité a posteriori
- calcul de $w_k^{(i)}$ : $$\begin{align}
w_k^{(i)} & \propto w_{k-1}^{(i)} \:
\frac{\tilde p\left(z_k|x_k^{(i)}\right) \:
\tilde p\left(x_k^{(i)}|x_{k-1}^{(i)}\right)}
{\tilde q\left(x_k^{(i)}|\textbf{x}_{0:k-1}^{(i)},\textbf{z}_{1:k}\right)} \\
& \propto w_{k-1}^{(i)} \: \frac{\tilde p\left(z_k|x_k^{(i)}\right) \:
\tilde p\left(x_k^{(i)}|x_{k-1}^{(i)}\right)}{\tilde p(x_k|x_{k-1},z_k)} \\
& \propto w_{k-1}^{(i)} \: \frac{\tilde p\left(z_k|x_k^{(i)}\right) \:
\tilde p\left(x_k^{(i)}|x_{k-1}^{(i)}\right)}
{\displaystyle{\frac{\tilde p\left(z_k|x_k^{(i)},x_{k-1}^{(i)}\right) \:
\tilde  p\left(x_k^{(i)}|x_{k-1}^{(i)}\right)}
{\tilde p\left(z_k|x_{k-1}^{(i)}\right)}}} \\
& \propto w_{k-1}^{(i)} \: \tilde p\left(z_k|x_{k-1}^{(i)}\right) \:
\frac{\tilde p\left(z_k|x_k^{(i)}\right)}
{\tilde p\left(z_k|x_k^{(i)},x_{k-1}^{(i)}\right)} \\
& \propto w_{k-1}^{(i)} \: \tilde p\left(z_k|x_{k-1}^{(i)}\right)
& \text{(observations)} 
\end{align}$$

- par récurrence, $w_k^{(i)}$ ne dépend pas de $x_k^{(i)}$, car pour $l \in [\![0,k-1]\!]$, $w_l^{(i)}$ ne dépend pas de $x_k^{(i)}$

## Variance des poids :

- $$
\begin{align}
\text{Var}_{q\left(\cdot|x_{k-1}^{(i)},z_k\right)}\left(w_k^{(i)}\right)
& = \mathbb{E}_{q\left(\cdot|x_{k-1}^{(i)},z_k\right)}
\left({w_k^{(i)}}^2\right) -
\mathbb{E}_{q\left(\cdot|x_{k-1}^{(i)},z_k\right)}
\left(w_k^{(i)}\right)^2 \\
& \propto \int {w_k^{(i)}}^2
q\left(dx_k^{(i)}|x_{k-1}^{(i)},z_k\right) dx_k^{(i)} -
\left(\int w_k^{(i)}
q\left(dx_k^{(i)}|x_{k-1}^{(i)},z_k\right) dx_k^{(i)} \right)^2 \\
& \propto {w_k^{(i)}}^2
\int q\left(dx_k^{(i)}|x_{k-1}^{(i)},z_k\right) dx_k^{(i)} -
{w_k^{(i)}}^2
\left(\int q\left(dx_k^{(i)}|x_{k-1}^{(i)},z_k\right) dx_k^{(i)} \right)^2 \\
& \propto {w_k^{(i)}}^2 - {w_k^{(i)}}^2 \\
& = 0
\end{align}
$$

## Taille effective :

- par Monte-Carlo pondéré, on a $$
\frac{p\left(\textbf{x}_{0:k}^{(i)}|\textbf{z}_{1:k}\right)}
{q\left(\textbf{x}_{0:k}^{(i)}|\textbf{z}_{1:k}\right)} \propto w_k^{(i)}
$$
- donc $$
\text{Var}\left(\frac{p\left(\textbf{x}_{0:k}^{(i)}|\textbf{z}_{1:k}\right)}
{q\left(\textbf{x}_{0:k}^{(i)}|\textbf{z}_{1:k}\right)}\right) \propto
\text{Var}\left(w_k^{(i)}\right) = 0
$$
- finalement $$
N_{eff} = \frac{N}{1 +
\text{Var}\left(\frac{p\left(\textbf{x}_{0:k}^{(i)}|\textbf{z}_{1:k}\right)}
{q\left(\textbf{x}_{0:k}^{(i)}|\textbf{z}_{1:k}\right)}\right)} = N
$$
- donc il n'y a pas de dégénérescence particulaire pour $$
q(x_k|\textbf{x}_{0:k-1},\textbf{z}_{1:k}) = p(x_k|x_{k-1},z_k)
$$ et le filtre particulaire est équivalent au SIS
- il s'agit donc de la fonction d'importance optimale pour le filtrage particulaire

## Calcul de la loi a posteriori

- calcul de $\tilde p\left(z_k|x_{k-1}^{(i)}\right)$ : $$\begin{align}
\tilde p\left(z_k|x_{k-1}^{(i)}\right) & = \int \tilde p\left(z_k, x_k|x_{k-1}^{(i)}\right) \: dx_k \\
& = \int \tilde p\left(z_k|x_k,x_{k-1}^{(i)}\right) \: \tilde p\left(x_k|x_{k-1}^{(i)}\right) \: dx_k \\
& = \int \tilde p\left(z_k|x_k\right) \: \tilde p\left(x_k|x_{k-1}^{(i)}\right) \: dx_k
\end{align}
$$
- en général, on ne sait pas calculer $\tilde p\left(z_k|x_{k-1}^{(i)}\right)$, ce choix de densité d'importance est le plus souvent impossible à mettre en place en pratique

## Densité a priori

- dans le SIR, on choisi $q\left(x_k|x_{k-1}^{(i)},z_k\right) = p\left(x_k|x_{k-1}^{(i)}\right)$ la densité a priori
- sous-optimale mais plus simple à calculer

**TODO : autres densités ?**
