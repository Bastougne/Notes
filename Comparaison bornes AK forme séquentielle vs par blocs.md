Les résultats de §4.5 sont les mêmes dans les deux formulations : mêmes matrices, mêmes identités, mêmes bornes. Seules changent la façon de définir les objets, la route vers les identités centrales, et ce qu'on peut en lire.

Dans le tableau, la colonne A est la formulation actuelle (retenues = ancien lot, écartées = nouveau lot, via la 3.14). La colonne B est la formulation par blocs, telle qu'elle était dans le manuscrit ce matin (ex-4.2 à ex-4.7).

| Résultat | A : séquentielle | B : par blocs, avec $\mathfrak D_G$ |
|---|---|---|
| $k_{d\mid r}$ (déf. 4.3) | $k_{\mid r}(\tilde{\mathbf x}^d,\cdot)$, où $k_{\mid r}$ est le noyau conditionnel UK de la 3.14 | $k(\tilde{\mathbf x}^d,x)-\bar k(\tilde{\mathbf x}^d,\tilde{\mathbf x}^r)\,\mathcal G_r^{-1}\bar k(\tilde{\mathbf x}^r,x)$, posé directement, puis $=P^\top\bar k(\tilde{\mathbf x},x)$ (identités résiduelles). **Même fonction.** |
| $\tilde{\mathbf z}_{d\mid r}$ (déf. 4.4) | Innovation UK du lot écarté : $\tilde{\mathbf z}^d-\hat z_{\mathrm{UK}}(\tilde{\mathbf x}^d,\tilde{\mathbf x}^r,\tilde{\mathbf z}^r,g)$ | $\tilde{\mathbf z}^d-\bar k(\tilde{\mathbf x}^d,\tilde{\mathbf x}^r)\,\mathcal G_r^{-1}\bar{\mathbf z}^r=P^\top\bar{\mathbf z}$. **Même vecteur** (par la 3.11). |
| $S_{d\mid r}$ | Covariance de l'innovation, $k_{d\mid r}(\tilde{\mathbf x}^d)+I\otimes\tilde R$ ; SPD car c'est une covariance (4.2) | Complément de Schur de $\mathcal G(\tilde{\mathbf x})$ par rapport à $\mathcal G(\tilde{\mathbf x}^r)$, inversible car les deux le sont ; le lien avec $k_{d\mid r}(\tilde{\mathbf x}^d)$ passe par les identités résiduelles. **Même matrice.** |
| Processus résiduel (4.2, ex-4.6) | Identique. Il est le GP conditionné aux données retenues, à la limite de l'a priori vague (preuve de la 3.14) | Identique, même preuve. La lecture « conditionnelle » ne s'appuie que sur la 3.10, prise en un point |
| Estimateur complet = sous-échantillonné + correction (4.3, ex-4.7) | C'est la 3.14 avec ancien = $r$ et nouveau = $d$ | Obtenu par $\bar k^\top\mathfrak D_G\bar{\mathbf z}$, la factorisation de Schur et les identités résiduelles |
| $e_r$, $e_f$ en formes résiduelles (4.3, ex-4.4 puis ex-4.7) | $e_f$ est la norme de la correction ; $e_r$ est sa covariance, via $\mathbb V[\tilde{\mathbf z}_{d\mid r}]=S_{d\mid r}$ | D'abord $e_r=\bar k^\top\mathfrak D_G\bar k$ et $e_f=\mid\bar k^\top\mathfrak D_G\bar{\mathbf z}\mid$ (ex-4.4). Puis $\mathfrak D_G=P\,S_{d\mid r}^{-1}P^\top$ (Schur) et les identités résiduelles redonnent les mêmes formes |
| $e_r=R_{\mathrm{UK}}(x,\tilde{\mathbf x}^r,g)-R_{\mathrm{UK}}(x,\tilde{\mathbf x},g)$ et la borne gratuite (rem. 4.8) | Identité sur la covariance de la 3.14 | Argument d'orthogonalité du BLUE (ex-B.15) |
| Lecture de Kalman | Littérale : la 4.3 est l'étape de correction d'une mise à jour séquentielle | Une analogie de forme seulement |
| 4.4 (normes de Mahalanobis), 4.5 (lien avec le krigeage simple), rem. 4.9 | Identiques. B.16 fait le même calcul par blocs | Identiques |
| §4.5.3 (4.6, 4.7, loi du $\chi^2$) | Identiques. Le caractère gaussien de la différence se lit sur la 4.3 avec la 4.2 | Identiques. Il vient de la linéarité en $\tilde{\mathbf z}$ |
| §4.5.4 (4.8 à 4.12) | Identiques. L'emboîtement se lit comme deux mises à jour successives, $r\to e\to\tilde{\mathbf x}$ | Identiques, avec la même preuve par la remarque 4.8 |

**Ce qui diffère vraiment**
- **La route vers les identités centrales :**
  - A y arrive en une étape : la 3.14, appliquée aux retenues et aux écartées.
  - B enchaîne quatre énoncés intermédiaires : le pseudo-inverse de $\mathcal G_h$, la définition de $\mathfrak D_G$, l'ex-4.4, puis la décomposition de Schur et les identités résiduelles.
  - Cela fait environ trois pages et quatre propositions de plus.
- **Ce qui porte la preuve :**
  - B est de l'algèbre linéaire autonome (Schur, orthogonalité).
  - A repose sur la 3.14, donc sur un passage à la limite (B.9, puis B.7). Ce coût n'est pas propre à §4.5, puisque la 3.14 sert aussi au ch.3.
- **Ce qu'on peut en lire :**
  - A donne un sens à chaque objet : noyau conditionnel, innovation, covariance d'innovation. Cela justifie la notation $d|r$ et vaut pour n'importe quel découpage en lots.
  - B donne une lecture matricielle : toute l'erreur, pour tout $x$, vient d'une seule matrice indépendante de $x$, SPSD de rang $\tilde N^dd_z$.
  - Cette lecture matricielle n'a pas d'avantage de calcul : la forme résiduelle travaille en dimension $\tilde N^dd_z$, alors que $\mathfrak D_G$ vit en dimension $(\tilde N+r_{\mathrm{UK}})d_z$.
- **L'alignement avec ton préprint :** B est la formulation de l'article, A s'en écarte.