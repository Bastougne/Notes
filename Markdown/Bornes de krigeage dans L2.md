## Principe :

- [[Theoretical Analysis of Kriging Error under Training Data Subsampling.pdf]]
- connaissant les données d'entraînement $\tilde{X}=\begin{bmatrix}\tilde{X}_p&\tilde{X}_l\end{bmatrix}^T$
- on cherche à déterminer l'erreur que l'omission des $\tilde{X}_l$ causent sur le krigeage
- on veut calculer $$e(x)=||\tilde{f}_{pl}(x)-\tilde{f}_p(x)||_{L^2(\Omega)}=\sqrt{\mathbb{E}\left(\left(\tilde{f}_{pl}(x)-\tilde{f}_{pl}(x)\right)^2\right)}$$

## Factorisation de l'espérance :

- en remarquant que $$
f(x)-\tilde{f}_p(x)=\left(f(x)-\tilde{f}_{pl}(x)\right)+\left(\tilde{f}_{pl}(x)-\tilde{f}_p(x)\right)
$$
- et en décomposant par orthogonalité dans l'espace des variables aléatoires connaissant $\tilde{X}_p$, on trouve que $$\begin{align}
e(x)^2=&||\tilde{f}_{pl}(x)-\tilde{f}_p(x)||_{L^2(\Omega)}^2\\
=&||f(x)-\tilde{f}_p(x)||_{L^2(\Omega)}^2-||f(x)-\tilde{f}_{pl}(x)||_{L^2(\Omega)}^2\\
=&R^*_s(x,\tilde{X}_p)-R^*_s(x,\tilde{X})\\
=&-k(x,\tilde{X}_p)\left(k(\tilde{X}_p,\tilde{X}_p)+\tilde{R}I_{\tilde{n}_p}\right)^{-1}k(\tilde{X}_p,x)+k(x,\tilde{X})\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)^{-1}k(\tilde{X},x)\\
=&k(x,\tilde{X})\left(\begin{bmatrix}A&B\\B^T&D\end{bmatrix}^{-1}-\begin{bmatrix}A^{-1}&0\\0&0\end{bmatrix}\right)k(\tilde{X},x)\\
=&k(x,\tilde{X})\begin{bmatrix}A^{-1}BS^{-1}B^TA^{-1}&A^{-1}BS^{-1}\\-S^{-1}B^TA^{-1}&S^{-1}\end{bmatrix}k(\tilde{X},x)\\
=&\begin{bmatrix}k(x,\tilde{X}_p)&k(x,\tilde{X}_l)\end{bmatrix}\begin{bmatrix}-A^{-1}B\\I\end{bmatrix}S^{-1}\begin{bmatrix}-A^{-1}B\\I\end{bmatrix}^T\begin{bmatrix}k(\tilde{X}_p, x)\\k(\tilde{X}_l,x)\end{bmatrix}
\end{align}$$
- posons $$
\varphi:x\mapsto k(\tilde{X}_l,x)-k(\tilde{X}_l,\tilde{X}_p)\left(k(\tilde{X}_p,\tilde{X}_p)+\tilde{R}I_{\tilde{n}_p}\right)^{-1}k(\tilde{X}_p,x)
$$
- on a alors $$\begin{align}
\varphi(x)=&k(\tilde{X}_l,x)-k(\tilde{X}_l,\tilde{X}_p)\left(k(\tilde{X}_p,\tilde{X}_p)+\tilde{R}I_{\tilde{n}_p}\right)^{-1}k(\tilde{X}_p,x)\\
=&\begin{bmatrix}-k(\tilde{X}_l,\tilde{X}_p)\left(k(\tilde{X}_p,\tilde{X}_p)+\tilde{R}I_{\tilde{n}_p}\right)^{-1}\\I\end{bmatrix}^T\begin{bmatrix}k(\tilde{X}_p, x)\\k(\tilde{X}_l,x)\end{bmatrix}\\
=&\begin{bmatrix}-A^{-1}B\\I\end{bmatrix}^T\begin{bmatrix}k(\tilde{X}_p, x)\\k(\tilde{X}_l,x)\end{bmatrix}\\
\end{align}$$
- de même $$\begin{align}
\varphi(\tilde{X}_l)=&k(\tilde{X}_l,\tilde{X}_l)-k(\tilde{X}_l,\tilde{X}_p)\left(k(\tilde{X}_p,\tilde{X}_p)+\tilde{R}I_{\tilde{n}_p}\right)^{-1}k(\tilde{X}_p,\tilde{X}_l)\\
=&D-\tilde{R}I_{\tilde{n}_l}-B^TA^{-1}B\\
=&S-\tilde{R}I_{\tilde{n}_l}\\
\end{align}$$
- finalement, il vient que $$
e(x)=\sqrt{\varphi(x)^T\left(\varphi(\tilde{X}_l)+\tilde{R}I_{\tilde{n}_l}\right)^{-1}\varphi(x)}
$$

## Étude asymptotique de $\varphi$ :

- quitte à multiplier $\varphi$ par $\sigma_{max}$, on se ramène à $\text{max}_x(k(x,x))=1$
- à l'ordre 1, on a $$
\begin{align}
k(\tilde{X}_l,\tilde{X}_p)\approx&\begin{bmatrix}\alpha&\cdots&\cdots&\alpha\\\vdots&\ddots&\ddots&\vdots\\\alpha&\cdots&\cdots&\alpha\end{bmatrix}=\alpha\:J_{\tilde{n}_l,\tilde{n}_p}\\
k(\tilde{X}_i,\tilde{X}_i)+\tilde{R}I_{\tilde{n}_i}\approx&\begin{bmatrix}1+\tilde{R}&\kappa&\cdots&\kappa\\\kappa&\ddots&\ddots&\vdots\\\vdots&\ddots&\ddots&\vdots\\\kappa&\cdots&\kappa&1+\tilde{R}\end{bmatrix}=(1+\tilde{R}-\kappa)\:I_{\tilde{n}_i}+\kappa\:J_{\tilde{n}_i},\quad i\in\{p,l\}\\
k(\tilde{X}_p,x)\approx&\begin{bmatrix}\kappa\\\vdots\\\kappa\end{bmatrix}=\kappa\:\mathbb{1}_{\tilde{n}_p}\\
k(\tilde{X}_l,x)\approx&\begin{bmatrix}\varepsilon\\\vdots\\\varepsilon\end{bmatrix}=\varepsilon\:\mathbb{1}_{\tilde{n}_l}\\
\end{align}
$$
- puisque $A$ est de la forme $aI+bJ$, on sait par [Sherman-Morrison](https://fr.wikipedia.org/wiki/Formule_de_Sherman-Morrison) que son inverse s'écrit sous la forme $$
A^{-1}=\frac{1}{a}\:I-\frac{b}{a(a+nb)}\:J
$$
- d'où $$\begin{align}
k(\tilde{X}_l,\tilde{X}_p)\left(k(\tilde{X}_p,\tilde{X}_p)+\tilde{R}I_{\tilde{n}_p}\right)^{-1}\approx&\alpha\:J_{\tilde{n}_l,\tilde{n}_p}\times\left(\frac{1}{a}\:I_{\tilde{n}_p}-\frac{b}{a(a+nb)}\:J_{\tilde{n}_p}\right)\\
=&\frac{\alpha}{a}\:J_{\tilde{n}_l,\tilde{n}_p}-\frac{\alpha b\tilde{n}_p}{a(a+\tilde{n}_pb)}\:J_{\tilde{n}_l,\tilde{n}_p}\\
=&\frac{\alpha(a+b\tilde{n}_p)-\alpha b\tilde{n}_p}{a(a+b\tilde{n}_p)}\:J_{\tilde{n}_l,\tilde{n}_p}\\
=&\frac{\alpha}{a+b\tilde{n}_p}\:J_{\tilde{n}_l,\tilde{n}_p}\\
=&\frac{\alpha}{\underbrace{1+\tilde{R}+(\tilde{n}_p-1)\kappa}_{\lambda:=}}\:J_{\tilde{n}_l,\tilde{n}_p}\\
\end{align}$$
- on a donc $$\begin{align}
\varphi(x)\approx&\varepsilon\:\mathbb{1}_{\tilde{n}_l}-\frac{\alpha}{\lambda}\:J_{\tilde{n}_l,\tilde{n}_p}\times \kappa\:\mathbb{1}_{\tilde{n}_p}\\
=&\frac{\lambda\varepsilon-\alpha\tilde{n}_p\kappa}{\lambda}\:\mathbb{1}_{\tilde{n}_l}
\end{align}$$
- de même $$\begin{align}
\varphi(\tilde{X}_l)+\tilde{R}I_{\tilde{n}_l}\approx&(1+\tilde{R}-\kappa)\:I_{\tilde{n}_l}+\kappa\:J_{\tilde{n}_l}-\frac{\alpha}{\lambda}\:J_{\tilde{n}_l,\tilde{n}_p}\times\alpha\:J_{\tilde{n}_p,\tilde{n}_l}\\
=&(\lambda-\tilde{n}_p\kappa)\:I_{\tilde{n}_l}+\frac{\lambda\kappa-\alpha^2\tilde{n}_p}{\lambda}\:J_{\tilde{n}_l}\\
\end{align}$$
- la formule de Sherman-Morrison s'applique également ici pour donner $$\begin{align}
\left(\varphi(\tilde{X}_l)+\tilde{R}I_{\tilde{n}_l}\right)^{-1}\approx&\frac{1}{\lambda-\tilde{n}_p\kappa}\:I_{\tilde{n}_l}-\frac{\frac{\lambda\kappa-\alpha^2\tilde{n}_p}{\lambda}}{(\lambda-\tilde{n}_p\kappa)(\lambda-\tilde{n}_p\kappa+\tilde{n}_l\frac{\lambda\kappa-\alpha^2\tilde{n}_p}{\lambda})}\:J_{\tilde{n}_l}\\
=&\frac{1}{\lambda-\tilde{n}_p\kappa}\:I_{\tilde{n}_l}-\frac{\lambda\kappa-\alpha^2\tilde{n}_p}{(\lambda-\tilde{n}_p\kappa)(\lambda^2-\lambda\tilde{n}_p\kappa+\tilde{n}_l(\lambda\kappa-\alpha^2\tilde{n}_p))}\:J_{\tilde{n}_l}\\
=&\frac{1}{\lambda-\tilde{n}_p\kappa}\left(I_{\tilde{n}_l}-\frac{\lambda\kappa-\alpha^2\tilde{n}_p}{\lambda(\lambda-\tilde{n}_p\kappa)+\tilde{n}_l(\lambda\kappa-\alpha^2\tilde{n}_p)}\:J_{\tilde{n}_l}\right)\\
\end{align}$$
- finalement, il vient que $$\begin{align}
e(x)^2\approx&\sigma_{max}^2\frac{(\lambda\varepsilon-\alpha\tilde{n}_p\kappa)^2}{\lambda^2(\lambda-\tilde{n}_p\kappa)}\:\mathbb{1}_{\tilde{n}_l}^T\left(I_{\tilde{n}_l}-\frac{\lambda\kappa-\alpha^2\tilde{n}_p}{\lambda(\lambda-\tilde{n}_p\kappa)+\tilde{n}_l(\lambda\kappa-\alpha^2\tilde{n}_p)}\:J_{\tilde{n}_l}\right)\mathbb{1}_{\tilde{n}_l}\\
=&\sigma_{max}^2\frac{(\lambda\varepsilon-\alpha\tilde{n}_p\kappa)^2}{\lambda^2(\lambda-\tilde{n}_p\kappa)}\left(\tilde{n}_l-\tilde{n}_l^2\:\frac{\lambda\kappa-\alpha^2\tilde{n}_p}{\lambda(\lambda-\tilde{n}_p\kappa)+\tilde{n}_l(\lambda\kappa-\alpha^2\tilde{n}_p)}\right)\\
=&\sigma_{max}^2\frac{(\lambda\varepsilon-\alpha\tilde{n}_p\kappa)^2}{\lambda^2(\lambda-\tilde{n}_p\kappa)}\:\tilde{n}_l\:\frac{\lambda(\lambda-\tilde{n}_p\kappa)+\tilde{n}_l(\lambda\kappa-\alpha^2\tilde{n}_p)-\tilde{n}_l(\lambda\kappa-\alpha^2\tilde{n}_p)}{\lambda(\lambda-\tilde{n}_p\kappa)+\tilde{n}_l(\lambda\kappa-\alpha^2\tilde{n}_p)}\\
=&\sigma_{max}^2\frac{(\lambda\varepsilon-\alpha\tilde{n}_p\kappa)^2}{\lambda}\:\tilde{n}_l\:\frac{1}{\lambda(\lambda-\tilde{n}_p\kappa)+\tilde{n}_l(\lambda\kappa-\alpha^2\tilde{n}_p)}\\
=&\sigma_{max}^2\frac{\tilde{n}_l(\lambda\varepsilon-\alpha\tilde{n}_p\kappa)^2}{\lambda(\lambda^2+\lambda(\tilde{n}_l-\tilde{n}_p)\kappa-\alpha^2\tilde{n}_p\tilde{n}_l)},\quad\quad\lambda=1+\tilde{R}+(\tilde{n}_p-1)\kappa\\
\end{align}$$
- dans le cas particulier où $\kappa\rightarrow1$ et $\tilde{R}\rightarrow0$, on a $\lambda\rightarrow\tilde{n}_p$, et donc $$\begin{align}
e(x)\approx&\sigma_{max}\sqrt{\frac{\tilde{n}_l(\tilde{n}_p\varepsilon-\alpha\tilde{n}_p)^2}{\tilde{n}_p(\tilde{n}_p^2+\tilde{n}_p(\tilde{n}_l-\tilde{n}_p)-\alpha^2\tilde{n}_p\tilde{n}_l)}}\\
=&\sigma_{max}\frac{|\varepsilon-\alpha|}{\sqrt{1-\alpha^2}}\\
\end{align}$$
