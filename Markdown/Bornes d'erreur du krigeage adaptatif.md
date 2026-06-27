## Principe :
- la sélection de certains échantillons en fonction de leur proximité avec les particules dégrade l'estimation de la fonction d'observation par krigeage
- on cherche à quantifier l'erreur commise par cette sélection en distinguant les échantillons $\tilde X$ comme suit :
	1) une part $1-\alpha$ de $\tilde X$ est proche des particules dans $X$, et forme $\tilde X_p$
	2) une part $\alpha$ de $\tilde X$ est donc loin des particules, et forme $\tilde X_l$
	3) une part $\beta$ de $\tilde X_l$ est loin de tous les échantillons de $\tilde X_p$, et interagit donc peu avec les échantillons proches des particules

## Équations :

- [[Équations des bornes d'erreur du krigeage adaptatif]]
- en posant $$
\left\{\begin{align}
&\theta:=\frac{\alpha(1-\alpha)(\beta k_l^2+(1-\beta)k_p^2)}{\lambda_{min}(A)^2/\tilde n^2}\\
\\
&\psi_{\mu^*}:x\mapsto\frac{\left((1-\alpha)k_p^2x+\alpha k_l^2\right)(x+1)}{\lambda_{min}(D-CA^{-1}B)^2/(n\tilde n)}\\
\\
&\psi_{R^*}:x\mapsto\frac{\left((1-\alpha)k_p^2x+
\alpha k_l^2\right)^2}{\lambda_{min}(D-CA^{-1}B)^2/(n\tilde n)^2}\\
\end{align}\right.
$$ on trouve $$
\begin{align}
e_{\mu^*,sélection}&\leq\sqrt{\psi_{\mu^*}(\theta)}\:||\tilde Z-\mu(\tilde X)||_2\\
e_{R^*,sélection}&\leq\sqrt{\psi_{R^*}(\theta)}
\end{align}
$$
- répartition des données d'entraînement en 3 zones :
	- zone rouge : les $\#\mathcal{R}:=(1-\alpha)\tilde n$ échantillons proches des particules
	- zone verte : les $\#\mathcal{V}:=\alpha\beta\tilde n$ échantillons loin des particules et des échantillons proches
	- zone jaune : les $\#\mathcal{J}:=\alpha(1-\beta)\tilde n$ échantillons à l'interface

![[Sélection des données de krigeage.png]]
- en utilisant les cardinaux des 3 zones, on peut réécrire $$
\left\{\begin{align}
&\theta=\frac{\#\mathcal{R}(\#\mathcal{V}\:k_l^2+\#\mathcal{J}\:k_p^2)}{\lambda_{min}(A)^2}\\
\\
&\psi_{\mu^*}(x)=\frac{\left(\#\mathcal{R}\:k_p^2x+(\#\mathcal{V}+\#\mathcal{J})k_l^2\right)(x+1)}{\lambda_{min}(D-CA^{-1}B)^2/n}\\
\\
&\psi_{R^*}(x)=\frac{\left(\#\mathcal{R}\:k_p^2x+(\#\mathcal{V}+\#\mathcal{J})k_l^2\right)^2}{\lambda_{min}(D-CA^{-1}B)^2/n^2}\\
\end{align}\right.
$$

## Limites :

- idéalement, on voudrait que $k_l$ soit petit devant $k_p$ pour pouvoir négliger les échantillons dans la zone verte, et que $\beta$ soit proche de $1$ pour minimiser la taille de la zone jaune
- en faisant tendre $k_l/k_p$ vers $0$, on trouve $$
\begin{align}
\theta\underset{k_l/k_p\rightarrow0}{\longrightarrow}&\frac{\#\mathcal{R}\:\#\mathcal{J}\:k_p^2}{\lambda_{min}(A)^2}\\
\\
\\
\psi_{\mu^*}(\theta)\underset{k_l/k_p\rightarrow0}{\longrightarrow}&\frac{\#\mathcal{R}\:k_p^2\frac{\#\mathcal{R}\:\#\mathcal{J}\:k_p^2}{\lambda_{min}(A)^2}\left(\frac{\#\mathcal{R}\:\#\mathcal{J}\:k_p^2}{\lambda_{min}(A)^2}+1\right)}{\lambda_{min}(D-CA^{-1}B)^2/n}\\
\\
\underset{k_l/k_p\rightarrow0}{\longrightarrow}&\frac{\#\mathcal{R}^2\:\#\mathcal{J}\:k_p^4(\#\mathcal{R}\:\#\mathcal{J}\:k_p^2+\lambda_{min}(A)^2)}{\lambda_{min}(A)^4\:\lambda_{min}(D-CA^{-1}B)^2/n}\\
\\
\\
\psi_{R^*}(\theta)\underset{k_l/k_p\rightarrow0}{\longrightarrow}&\frac{\left(\#\mathcal{R}\:k_p^2\frac{\#\mathcal{R}\:\#\mathcal{J}\:k_p^2}{\lambda_{min}(A)^2}\right)^2}{\lambda_{min}(D-CA^{-1}B)^2/n^2}\\
\\
\underset{k_l/k_p\rightarrow0}{\longrightarrow}&\frac{\#\mathcal{R}^4\:\#\mathcal{J}^2\:k_p^8}{\lambda_{min}(A)^4\:\lambda_{min}(D-CA^{-1}B)^2/n^2}\\
\end{align}
$$
- de même, en faisant tendre $\beta$ vers $1$, on a $$
\begin{align}
\theta\underset{\beta\rightarrow1}{\longrightarrow}&\frac{\#\mathcal{R}\:\#\mathcal{V}\:k_l^2}{\lambda_{min}(A)^2}\\
\\
\\
\psi_{\mu^*}(\theta)\underset{\beta\rightarrow1}{\longrightarrow}&\frac{\left(\#\mathcal{R}\:k_p^2\frac{\#\mathcal{R}\:\#\mathcal{V}\:k_l^2}{\lambda_{min}(A)^2}+\#\mathcal{V}\:k_l^2\right)\left(\frac{\#\mathcal{R}\:\#\mathcal{V}\:k_l^2}{\lambda_{min}(A)^2}+1\right)}{\lambda_{min}(D-CA^{-1}B)^2/n}\\
\\
\underset{\beta\rightarrow1}{\longrightarrow}&\frac{\#\mathcal{V}\:k_l^2\left(\#\mathcal{R}^2\:k_p^2+\lambda_{min}(A)^2\right)\left(\#\mathcal{R}\:\#\mathcal{V}\:k_l^2+\lambda_{min}(A)^2\right)}{\lambda_{min}(A)^4\:\lambda_{min}(D-CA^{-1}B)^2/n}\\
\\
\\
\psi_{R^*}(\theta)\underset{\beta\rightarrow1}{\longrightarrow}&\frac{\left(\#\mathcal{R}\:k_p^2\frac{\#\mathcal{R}\:\#\mathcal{V}\:k_l^2}{\lambda_{min}(A)^2}+\#\mathcal{V}\:k_l^2\right)^2}{\lambda_{min}(D-CA^{-1}B)^2/n^2}\\
\\
\underset{\beta\rightarrow1}{\longrightarrow}&\frac{\#\mathcal{V}^2\:k_l^4\left(\#\mathcal{R}^2\:k_p^2+\lambda_{min}(A)^2\right)^2}{\lambda_{min}(A)^4\:\lambda_{min}(D-CA^{-1}B)^2/n^2}\\
\end{align}
$$
- finalement, en faisant tendre simultanément $k_l/k_p$ vers $0$ et $\beta$ vers $1$, il vient $$
\begin{align}
\theta\underset{\begin{matrix}k_l/k_p\rightarrow0\\\beta\rightarrow1\end{matrix}}{\longrightarrow}&0\\
\\
\psi_{\mu^*}(\theta)\underset{\begin{matrix}k_l/k_p\rightarrow0\\\beta\rightarrow1\end{matrix}}{\longrightarrow}&0\\
\\
\psi_{R^*}(\theta)\underset{\begin{matrix}k_l/k_p\rightarrow0\\\beta\rightarrow1\end{matrix}}{\longrightarrow}&0\\
\end{align}
$$
- d'où $$
\begin{align}
e_{\mu^*,sélection}\underset{\begin{matrix}k_l/k_p\rightarrow0\\\beta\rightarrow1\end{matrix}}{\longrightarrow}&0\\
\\
e_{R^*,sélection}\underset{\begin{matrix}k_l/k_p\rightarrow0\\\beta\rightarrow1\end{matrix}}{\longrightarrow}&0\\
\end{align}
$$
