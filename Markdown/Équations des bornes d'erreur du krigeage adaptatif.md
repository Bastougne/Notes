## Notations :
- notons $||\cdot||$ la [[Norme d'opérateur]] sur les matrices
- décomposons $$
\tilde X=\begin{bmatrix}
\tilde X_p\\\tilde X_l
\end{bmatrix}
\text{,}\quad
\tilde Z=\begin{bmatrix}
\tilde Z_p\\\tilde Z_l
\end{bmatrix}
\quad\text{et}\quad
\mu(\tilde X)=\begin{bmatrix}
\mu(\tilde X_p)\\\mu(\tilde X_l)
\end{bmatrix}
$$
- les matrices de covariance s'écrivent alors $$
k(X, \tilde X)=\begin{bmatrix}
k(X, \tilde X_p)&k(X, \tilde X_l)
\end{bmatrix}
$$ et $$
k(\tilde X,\tilde X)+I_{\tilde n}\otimes\tilde R=\begin{bmatrix}
k(\tilde X_p,\tilde X_p)+I_{(1-\alpha)\tilde n}\otimes\tilde R&k(\tilde X_l,\tilde X_p)\\k(\tilde X_p,\tilde X_l)&k(\tilde X_l,\tilde X_l)+I_{\alpha\tilde n}\otimes\tilde R
\end{bmatrix}
$$
- décomposons de plus la matrice de Gram et son inverse en blocs $$
G=k(\tilde X,\tilde X)+I_{\tilde n}\otimes\tilde R=\begin{bmatrix}
A&B\\ C&D
\end{bmatrix}
\quad\text{et}\quad
G^{-1}=\begin{bmatrix}
\bar A&\bar B\\\bar C&\bar D
\end{bmatrix}
$$
- les blocs de $G$ et $G^{-1}$ sont reliés par les égalités suivantes $$
\left\{\begin{align}
\bar A=&A^{-1}+A^{-1}BS^{-1}CA^{-1}\\
\bar B=&-A^{-1}BS^{-1}\\
\bar C=&-S^{-1}CA^{-1}\\
\bar D=&S^{-1}\\
S=&D-CA^{-1}B\\
B=&C^T\\
\end{align}\right.
$$ où $S$ est le complément de Schur de $G$ par rapport à $A$

## Majorations préliminaires :
- par hypothèse, $k$ est bornée et on note $$
k_p:=\sup_{x_1,x_2}\{||k(x_1,x_2)||_F\}>0
$$
- on suppose que les covariances particules/échantillons éloignés sont bornées en norme de Frobenius par une constante $k_l>0$
- de même, en décomposant $\tilde X_l$ selon $\beta$, on a $$
||B||^2\leq\alpha(1-\alpha)\tilde n^2(\beta k_l^2+(1-\beta)k_p^2)
$$ et de même pour $||C||^2$

## Erreur sur la moyenne _a posteriori_ :
- la moyenne _a posteriori_ du krigeage s'écrit $$
\begin{align}
\mu^*(X,\tilde X,\tilde Z)=&\mu(X)+k(X,\tilde X)G^{-1}\left(\tilde Z-\mu(\tilde X)\right)\\
=&\mu(X)+\begin{bmatrix}k(X, \tilde X_p)&k(X, \tilde X_l)\end{bmatrix}\begin{bmatrix}\bar A&\bar B\\\bar C&\bar D\end{bmatrix}\left(\begin{bmatrix}\tilde Z_p\\\tilde Z_l\end{bmatrix}-\begin{bmatrix}\mu(\tilde X_p)\\\mu(\tilde X_l)\end{bmatrix}\right)\\
\end{align}
$$
- on remarque que $$
\begin{align}
\mu(X)+k(X,\tilde X)\begin{bmatrix}A^{-1}&0\\0&0\end{bmatrix}\left(\tilde Z-\mu(\tilde X)\right)=&\mu(X)+\begin{bmatrix}k(X, \tilde X_p)&k(X, \tilde X_l)\end{bmatrix}\begin{bmatrix}A^{-1}&0\\0&0\end{bmatrix}\left(\begin{bmatrix}\tilde Z_p\\\tilde Z_l\end{bmatrix}-\begin{bmatrix}
\mu(\tilde X_p)\\\mu(\tilde X_l)
\end{bmatrix}\right)\\
=&\mu(X)+k(X,\tilde X_p)A^{-1}\left(\tilde Z_p-\mu(\tilde X_p)\right)\\
=&\mu^*(X,\tilde X_p,\tilde Z_p)\\
\end{align}
$$
- posons $$
e_{\mu^*,sélection}:=||\mu^*(X,\tilde X,\tilde Z)-\mu^*(X,\tilde X_p,\tilde Z_p)||_2
$$
- on a alors $$
\begin{align}
\frac{{e_{\mu^*,sélection}}^2}{||\tilde Z-\mu(\tilde X)||_2^2}=&\frac{||\mu^*(X,\tilde X,\tilde Z)-\mu^*(X,\tilde X_p,\tilde Z_p)||_2^2}{||\tilde Z-\mu(\tilde X)||_2^2}\\
\leq&\left|\left|\begin{bmatrix}k(X, \tilde X_p)&k(X, \tilde X_l)\end{bmatrix}\begin{bmatrix}\bar A-A^{-1}&\bar B\\\bar C&\bar D\end{bmatrix}\right|\right|^2\\
\leq&\left|\left|\begin{bmatrix}k(X, \tilde X_p)(\bar A-A^{-1})+k(X, \tilde X_l)\bar C\\k(X, \tilde X_p)\bar B+k(X, \tilde X_l)\bar D\end{bmatrix}^T\right|\right|^2\\
\leq&||k(X, \tilde X_p)(\bar A-A^{-1})+k(X, \tilde X_l)\bar C||^2+||k(X, \tilde X_p)\bar B+k(X, \tilde X_l)\bar D||^2\\
\leq&||k(X, \tilde X_p)||^2(||\bar A-A^{-1}||^2+||\bar B||^2)+||k(X, \tilde X_l)||^2(||\bar C||^2+||\bar D||^2)\\
\leq&n(1-\alpha)\tilde nk_p^2(||\bar A-A^{-1}||^2+||\bar B||^2)+n\alpha\tilde nk_l^2(||\bar C||^2+||\bar D||^2)\\
\\
\leq&n(1-\alpha)\tilde nk_p^2(||A^{-1}BS^{-1}CA^{-1}||^2+||-A^{-1}BS^{-1}||^2)+\\
&n\alpha\tilde nk_l^2(||-S^{-1}CA^{-1}||^2+||S^{-1}||^2)\\
\\
\leq&n(1-\alpha)\tilde nk_p^2(||A^{-1}||^4\:||B||^2\:||C||^2\:||S^{-1}||^2+||A^{-1}||^2\:||B||^2\:||S^{-1}||^2)+\\
&n\alpha\tilde nk_l^2(||A^{-1}||^2\:||C||^2\:||S^{-1}||^2+||S^{-1}||^2)\\
\\
\leq&n(1-\alpha)\tilde nk_p^2\left(\frac{(\alpha(1-\alpha)\tilde n^2(\beta k_l^2+(1-\beta)k_p^2))^2}{\lambda_{min}(A)^4\:\lambda_{min}(S)^2}+\right.\\
&\phantom{n(1-\alpha)\tilde nk_p^2(\:\:}\left.\frac{\alpha(1-\alpha)\tilde n^2(\beta k_l^2+(1-\beta)k_p^2)}{\lambda_{min}(A)^2\:\lambda_{min}(S)^2}\right)+\\
&n\alpha\tilde nk_l^2\left(\frac{\alpha(1-\alpha)\tilde n^2(\beta k_l^2+(1-\beta)k_p^2)}{\lambda_{min}(A)^2\:\lambda_{min}(S)^2}+\frac{1}{\lambda_{min}(S)^2}\right)\\
\end{align}
$$


## Erreur sur la matrice de covariance :
- la covariance du krigeage s'écrit $$
\begin{align}
R^*(X,\tilde X)=&k(X, X)+I_n\otimes R-k(X,\tilde X)G^{-1}k(X,\tilde X)^T\\
=&k(X, X)+I_n\otimes R-\begin{bmatrix}k(X, \tilde X_p)&k(X, \tilde X_l)\end{bmatrix}\begin{bmatrix}\bar A&\bar B\\\bar C&\bar D\end{bmatrix}\begin{bmatrix}k(X, \tilde X_p)^T\\k(X, \tilde X_l)^T\end{bmatrix}\\
\end{align}
$$
- de même que précédemment, on remarque que $$
k(X, X)+I_n\otimes R-k(X,\tilde X)\begin{bmatrix}A^{-1}&0\\0&0\end{bmatrix}k(X,\tilde X)^T=R^*(X,\tilde X_p)
$$
- posons $$
e_{R^*,sélection}:=||R^*(X,\tilde X)-R^*(X,\tilde X_p)||
$$
- d'où $$
\begin{align}
e_{R^*,sélection}^2=&\left|\left|\begin{bmatrix}k(X, \tilde X_p)&k(X, \tilde X_l)\end{bmatrix}\begin{bmatrix}\bar A-A^{-1}&\bar B\\\bar C&\bar D\end{bmatrix}\begin{bmatrix}k(X, \tilde X_p)^T\\k(X, \tilde X_l)^T\end{bmatrix}\right|\right|^2\\
\\
=&\left|\left|k(X, \tilde X_p)\left((\bar A-A^{-1})k(X, \tilde X_p)+\bar Bk(X, \tilde X_l)\right)+\right.\right.\\
&\phantom{||}\left.\left.k(X, \tilde X_p)\left(\bar Ck(X, \tilde X_p)+\bar Dk(X, \tilde X_l)\right)\right|\right|^2\\
\\
\leq&||k(X, \tilde X_p)||^2(||\bar A-A^{-1}||^2\:||k(X, \tilde X_p)||^2+||\bar B||^2\:||k(X, \tilde X_l)||^2)+\\
&||k(X, \tilde X_l)||^2(||\bar C||^2\:||k(X, \tilde X_p)||^2+||\bar D||^2\:||k(X, \tilde X_l)||^2)\\
\\
\leq&n(1-\alpha)\tilde nk_p^2(n(1-\alpha)\tilde nk_p^2||\bar A-A^{-1}||^2+n\alpha\tilde nk_l^2||\bar B||^2)+\\
&n\alpha\tilde nk_l^2(n(1-\alpha)\tilde nk_p^2||\bar C||^2+n\alpha\tilde nk_l^2||\bar D||^2)\\
\\
\leq&n(1-\alpha)\tilde nk_p^2\left(n(1-\alpha)\tilde nk_p^2\frac{(\alpha(1-\alpha)\tilde n^2(\beta k_l^2+(1-\beta)k_p^2))^2}{\lambda_{min}(A)^4\:\lambda_{min}(S)^2}+\right.\\
&\phantom{n(1-\alpha)\tilde nk_p^2(\:\:}\left.n\alpha\tilde nk_l^2\frac{\alpha(1-\alpha)\tilde n^2(\beta k_l^2+(1-\beta)k_p^2)}{\lambda_{min}(A)^2\:\lambda_{min}(S)^2}\right)+\\
&n\alpha\tilde nk_l^2\left(n(1-\alpha)\tilde nk_p^2\frac{\alpha(1-\alpha)\tilde n^2(\beta k_l^2+(1-\beta)k_p^2)}{\lambda_{min}(A)^2\:\lambda_{min}(S)^2}+n\alpha\tilde nk_l^2\frac{1}{\lambda_{min}(S)^2}\right)\\
\\
\leq&\frac{(n(1-\alpha)\tilde nk_p^2)^2(\alpha(1-\alpha)\tilde n^2(\beta k_l^2+(1-\beta)k_p^2))^2}{\lambda_{min}(A)^4\:\lambda_{min}(S)^2}+\\
&\frac{2(n(1-\alpha)\tilde nk_p^2)(n\alpha\tilde nk_l^2)(\alpha(1-\alpha)\tilde n^2(\beta k_l^2+(1-\beta)k_p^2))}{\lambda_{min}(A)^2\:\lambda_{min}(S)^2}+\frac{(n\alpha\tilde nk_l^2)^2}{\lambda_{min}(S)^2}\\
\\
\leq&\frac{\left(n(1-\alpha)\tilde nk_p^2\frac{\alpha(1-\alpha)\tilde n^2(\beta k_l^2+(1-\beta)k_p^2)}{\lambda_{min}(A)^2}+
n\alpha\tilde nk_l^2\right)^2}{\lambda_{min}(S)^2}\\
\end{align}
$$
