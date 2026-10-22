# Revue des notations mathématiques (hors chapitre 5)

*Rédigée le 2026-10-21. Ce document propose, tu tranches (§11), puis j'applique par étapes (§10). Dans cette note, « annexe A » et « annexes B » désignent celles du manuscrit ; les annexes de la note elle-même sont numérotées I et II.*

**Avancement (2026-10-22, étapes 1 et 2 faites, puis lots F à Y).** Matrices en gras majuscule (D1) ; bruits $\eta_k$, $\varepsilon_k$ et $\tilde\varepsilon_i$, seuils $\epsilon$ ; APF ($a$, $v$) ; mesures de Monte-Carlo $\mu_p$, $\mu_q$ (gardées) ; processus gaussien $h$ ; erreurs en $\zeta$ ; coefficients multiplicateurs en $\alpha_\bullet$, taux $p_\mathrm{RS}$ ; largeurs de bande et échelle de l'UKF en $\ell$ ; coefficients de tendance et de carte $\boldsymbol\gamma$ ; **noyaux $\kappa$**, multiplicateur $\boldsymbol\Xi$, Gram par blocs $\bar{\mathbf G}$ ; nomenclature à une lettre de base par ligne (D1 à D8, D10, D11, D15 à D21) : appliqués hors chapitre 5, voir le journal en fin de note ; puis, l'après-midi du 22-10, ses neuf retours sur la nomenclature, l'identité $\operatorname{id}$ et la paire $(\boldsymbol\Omega,\boldsymbol\Upsilon)$ de la décomposition polaire (D23 à D25), puis ses décisions du soir, $N_\mathrm{MC}$, $N_\mathrm{steps}$, $\ln$ et $\mathcal B$ gardé (D26), **et la dernière passe : chapitre 5 aligné, nomenclature complète (D12, D27)** (trois dernières sections du journal). Le corps de la note (tableaux §3, §5 à §10, annexes) date du 21-10 ; seuls les paragraphes C1 à C10 et les décisions ont été remis à jour. **Un état antérieur de cette note a été remis en place le 22-10 vers 14 h (cause inconnue : elle avait perdu mes mises à jour du soir et de la nuit, et c'est un état de la soirée du 21-10 que j'ai retrouvé) : je les ai réécrites le même jour, le journal ci-dessous en fait foi.** Tout le reste est encore à l'état de proposition.

Tes règles, reprises telles quelles : (R1) deux objets ne portent pas le même nom ; (R2) s'ils le portent, ils sont de nature similaire (`h` et `h_MS`, oui ; `β` et `β_RLS`, non) ; (R3) on reconnaît la nature d'un objet à son nom (matrices en gras majuscule, par exemple). « Nom » = police + lettre + accents + indices.

---

## 1. En bref

**Ce qui est déjà solide.** Les n-uplets sont en gras minuscule ($\mathbf x$, $\tilde{\mathbf z}$, $\mathbf z_{1:k}$, $\mathbf w_{k\vert k}$), les estimations en chapeau, les ensembles en calligraphique, les dimensions en $d_\cdot$, les effectifs en $N$ et $\tilde N$, les rangs en $r_\cdot$ ; les marqueurs d'indice sont stables (`(i)` particule, `n`/`l` partie non linéaire/linéaire, `r`/`d` retenu/écarté). La plupart des homonymes sont de même nature (§7).

**Ce qui ne l'est pas.**

1. La majuscule italique désignait à la fois des matrices, des variables aléatoires, des fonctions, des mesures, des ensembles et des entiers : devant `K`, `S`, `T`, `J`, `F`, `C`, `U`, `V`, `G`, `E`, on ne savait pas ce que c'était. C'était la cause de la moitié des collisions ; la charte du §3 la résout sans changer de lettre, et **elle est appliquée** (étape 1 : 2 613 occurrences passées en gras, D1).
2. **Neuf collisions globales de nature** (C1 à C9, §4) exigent un changement de lettre ou de forme, dont deux touchent presque tout le manuscrit : `k` (noyau / temps) et `w` (bruit / poids). **Levées le 2026-10-21 (étape 2) :** `w` (C2), `ν` (C4) et la moitié de `μ` (C3, particule et poids auxiliaires de l'APF) ; restent `k` (C1), C5 à C9 et les moyennes de C3. Un dixième cas (C10, `h`) est à arbitrer.
3. Onze collisions locales dans des énoncés, des preuves et l'annexe A (§6).
4. Une nomenclature qui oublie une vingtaine de symboles et mélange des objets de nature différente (§8).

**Ce que je propose.** Une charte typographique (§3, déjà appliquée aux matrices), neuf renommages (§4) dont les lettres d'arrivée sont libres dans le manuscrit **et** dans le chapitre 5, onze retouches locales (§6), une nomenclature réécrite par nature (§8, encart en annexe II de cette note). Aucun autre alphabet n'est nécessaire (§3.5).

**Ce qu'il me faut de toi.** Plus aucune décision ◆ ouverte (D1, D3 et D4 sont tranchées et appliquées) ; tout le reste est coché ✔ par défaut, tu peux valider en bloc. Le seul point de ce reste que je te signale est $h$ (C10, D10) : je propose de garder les deux sens (fonction d'observation et processus gaussien d'un côté, largeur de bande de l'autre).

---

## 2. Périmètre et méthode

- **Fichiers.** `2_bayesian_estimation`, `2_Monte_Carlo`, `2_particle_filtering`, `3_Gaussian_processes`, `3_kriging`, `3_reduced_rank`, `3_RKHS`, `3_KKF`, `4_kriging_navigation`, `4_adaptive_kriging`, `4_subsampling`, `9_appendix`, `9_proofs_1`, `9_proofs_2`, plus la nomenclature de `0_administratif`. Les blocs `comment` (v2 : §4.5, remerciements…) sont exclus de l'inventaire, mais **ils réutilisent des lettres que je propose de prendre** (§9) : ils sont à repasser dans la charte au moment de leur retour (les matrices y sont déjà en gras, étape 1). `5_SLAM` a été inventorié à part, uniquement pour vérifier que les noms proposés y sont libres.
- **Inventaire.** Un analyseur dédié relève chaque symbole mathématique avec fichier, ligne, police (italique, gras, calligraphique, ajouré), lettre, accents, indice, exposant : **20 125 relevés**, dont 13 053 lettres de niveau 0 (108 lettres de base) et 7 072 dans des indices ou exposants. Pour chaque lettre ambiguë, j'ai relu les contextes dans le source (≈ 150 extraits) et rangé les objets par nature (scalaire, vecteur, n-uplet, matrice, variable aléatoire, fonction, ensemble, opérateur, mesure, entier, estimateur).
- **Nature « similaire ».** Même type mathématique ET même rôle (deux largeurs de bande ; pas un vecteur de coefficients et un facteur d'oubli).
- **Limites.** La détection se fait par lettre : un objet qui change de nom d'un chapitre à l'autre sans collision n'est pas compté. Quelques lettres fantômes viennent de `\overset{a.s.}` et de `n_\mathrm{MC}` écrit sans accolades (écartées). Les comptes d'occurrences sont des ordres de grandeur. Je n'ai pas vérifié les unités ni les domaines de définition. Les renvois `BE:64`, `B1:418`… abrègent le fichier : BE `2_bayesian_estimation`, MC `2_Monte_Carlo`, PF `2_particle_filtering`, GP `3_Gaussian_processes`, KR `3_kriging`, RR `3_reduced_rank`, RK `3_RKHS`, KK `3_KKF`, KN `4_kriging_navigation`, AK `4_adaptive_kriging`, B1 et B2 `9_proofs_1` et `9_proofs_2` ; les numéros de ligne sont ceux de l'inventaire, pris avant l'étape 1, et ont pu bouger de quelques lignes depuis (trois retours à la ligne de l'étape 1, éditions en cours).

---

## 3. Charte typographique proposée

### 3.1 Forme par nature

*La dernière colonne donne l'état du manuscrit après l'étape 1 (2026-10-21).*

| Nature | Forme | Exemples | État actuel |
|---|---|---|---|
| Scalaire, vecteur, point | minuscule italique | $x_k$, $z_k$, $\theta$, $a$ | ✔ |
| Fonction (même à valeurs matricielles) | minuscule italique | $f_k$, $h$, $\kappa$, $g$, $m$, $q$ | ✔ (le noyau est encore $k$ : C1) |
| N-uplet, vecteur empilé, vecteur de poids ou de constantes | minuscule grasse | $\mathbf x$, $\tilde{\mathbf z}$, $\mathbf w_{k\vert k}$, $\mathbf 1_n$, $\mathbf e_i$ | ✔ ($e_i$ n'est pas encore en gras : L10) |
| **Matrice** (covariance, gain, Gram, jacobienne, facteur) | **majuscule grasse** (`\mathbf`, `\boldsymbol` pour le grec) | $\mathbf P,\mathbf Q,\mathbf R,\mathbf F,\mathbf H,\mathbf K,\mathbf S,\mathbf G,\mathbf I,\boldsymbol\Sigma,\boldsymbol\Lambda,\boldsymbol\Psi,\boldsymbol\Pi$ | ✔ **fait le 2026-10-21** (2 613 occurrences) ; exceptions voulues en italique, listées au journal ; restent $\mathcal G$ (C8) et $\mathfrak D_G$ (nomenclature, retour prévu en v2) |
| Variable aléatoire | majuscule italique | $X$, $Z$, $U$, $V$, $Y$ | ✔ (ne se confond plus avec une matrice) |
| Effectif, constante | majuscule italique | $N$, $\tilde N$, $C$ | ✔ |
| Fonction nommée | majuscule italique | $K$, $L$, $J$, $K_\nu$, $S_\nu$ | ✔ ($K$, $J$, $S_\nu$ ne se confondent plus avec les matrices) ; restent $E$ et $G$ du théorème du représentant (L1) |
| Ensemble, espace | calligraphique | $\mathcal U$, $\mathcal H$, $\mathcal F$, $\mathcal L^p$, $\mathcal M$, $\mathcal S$ | ✔ sauf $E,F$ (BE:64) et $E_\mathrm{in}$ (KN:148) ; $\Theta$ toléré (espace des paramètres) |
| Opérateur | calligraphique | $\mathcal K$, $\mathcal Q_k$, $\mathcal T$ | ~ ($T$ en RR:12) |
| Mesure, loi, plongement | minuscule grecque | $\mu$, $\eta$, $\delta_x$, $\mu_X$ | ✔ |
| Corps, opérateurs probabilistes | ajouré | $\mathbb R$, $\mathbb N$, $\mathbb E$, $\mathbb V$, $\mathbb P$, $\mathbb 1$ | ✔ |
| Fonction ou opérateur usuel | droit | $\operatorname{tr}$, $\det$, $\operatorname{argmin}$ | ✔ |

Les fonctions à valeurs matricielles (noyau $\kappa$, tendance $g$) restent en minuscule ; leurs évaluations sur des n-uplets sont des matrices par application. Une matrice qui dépend d'arguments ($\mathbf G(\tilde{\mathbf x})$, $\mathbf W(x)$, $\mathbf H_k(x_k^n)$) garde son symbole de matrice. Les n-uplets de variables aléatoires de l'EnKF ($\mathbf X$, $\mathbf Z$) se lisent comme des matrices d'échantillons : ils restent.

**Avant / après** (le gras est appliqué ; les changements de lettre attendent leur étape).

- Kalman : $K_k=P_{k\vert k-1}H_k^\top S_k^{-1}$ devient $\mathbf K_k=\mathbf P_{k\vert k-1}\mathbf H_k^\top\mathbf S_k^{-1}$ (**fait**).
- Modèle : $x_k=f_k(x_{k-1},w_k)$ devient $x_k=f_k(x_{k-1},\eta_k)$, avec $z_k=h_k(x_k,\varepsilon_k)$ (C2, **fait**).
- Carte : $h\sim\mathcal{GP}(m,k)$ devient $h\sim\mathcal{GP}(m,\kappa)$ (C1, à faire).
- Importance : $\langle\mu,\varphi\rangle=\langle\nu,\frac{d\mu}{d\nu}\varphi\rangle$ devient $\langle\mu_p,\varphi\rangle=\langle\mu_q,\frac{d\mu_p}{d\mu_q}\varphi\rangle$ (C4, **fait**).
- Gram : $G(\tilde{\mathbf x})=k(\tilde{\mathbf x},\tilde{\mathbf x})+I_{\tilde N}\otimes\tilde R$ devient $\mathbf G(\tilde{\mathbf x})=\kappa(\tilde{\mathbf x},\tilde{\mathbf x})+\mathbf I_{\tilde N}\otimes\tilde{\mathbf R}$ : le gras est **fait**, $k\to\kappa$ reste à faire (C1).

### 3.2 Accents

| Accent | Sens actuel | Verdict |
|---|---|---|
| chapeau $\hat\cdot$ (≈ 420) | estimation, estimateur ($\hat z,\hat x,\hat\theta,\hat m,\hat{\mathbf P},\hat\psi,\hat N_\mathrm{eff}$) | ✔ aucun écart |
| tilde $\tilde\cdot$ (≈ 1 620) | (a) données d'apprentissage ($\tilde x,\tilde z,\tilde N,\tilde{\mathbf R}$ : ≈ 1 400) ; (b) non normalisé ($\tilde p,\tilde q,\tilde w,\tilde\mu,\tilde\lambda$) ; (c) évalué au point de linéarisation ($\tilde{\mathbf F}_k,\tilde{\mathbf H}_k$ en EKF et RBPF) ; (d) blanchi ($\tilde K_h$) ; (e) rang réduit ($\tilde{\mathbf G}$, RR:128) ; (f) auto-normalisé ($\tilde S^N$, sens **inverse** de (b)) ; (g) réalisations centrées ($\tilde u,\tilde v$, annexe A) ; (h) profil radial ($\tilde r$, B1:175) | (a), (b) et (c) sont des familles claires, à décrire dans la nomenclature ; (d) à (h) sont locaux. Aucun renommage |
| barre $\bar\cdot$, $\overline\cdot$ (≈ 120) | moyenne ($\overline x_0$) ; forme par blocs ($\bar k,\bar{\mathbf z},\mathcal G$) ; blocs de l'inverse ($\bar{\mathbf A},\dots,\bar{\mathbf M}$, annexe A) ; grandeurs décorrélées du RBPF ($\overline{\mathbf F},\overline f,\overline{\mathbf Q},\overline w$) | à décrire ; $\mathcal G$ devient $\bar{\mathbf G}$ (C8) |
| brève $\breve\cdot$ (9) | centre de classe ($\breve x^{(c)}$, AK) | ✔ |
| caron $\check\cdot$ | libre | n'est plus proposé (l'APF prend $a$ et $v$) |
| point $\dot\cdot$ | vitesse ($\dot p_k^{(x)}$, `4_kriging_navigation`) | ✔ appliqué le 2026-10-21 (L6) |

### 3.3 Indices et exposants

`_k` temps ; `_{k\vert l}` conditionnement ; `^{(i)}`, `^{(j)}` particule ; `^{(c)}` classe ; `^{(e)}`, `^{(o)}` évolution/observation ; `^n`, `^l`, `^{nl}` partie non linéaire/linéaire ; `^r`, `^d`, `_{d\vert r}`, `_{\vert r}` retenus/écartés ; `^+`, `^-` rééchantillonné/précédent ; `_\mathrm{old}`, `_\mathrm{new}` ; noms de méthode en droit (`_\mathrm{AK}`). Cohérent. Écart toléré : `^{(i)}` désigne aussi la i-ème composante dans la définition de la jacobienne (BE:419–425), passage générique qui précède toute particule.

### 3.4 Noms locaux

Un symbole défini dans un énoncé ou une preuve peut reprendre une lettre globale **seulement si** l'objet local est de même nature que l'objet global ; sinon il prend une lettre libre dans le chapitre ou un nom droit (`\operatorname{loss}`).

### 3.5 Pas d'autre alphabet

Libres dans tout le texte visible (chapitre 5 inclus ; blocs `comment` : voir §9) : minuscules $\gamma$, $\ell$ (ell cursif) ; majuscules $\Xi$, $\Upsilon$, $Y$ ; calligraphiques $\mathcal A$, $\mathcal D$, $\mathcal J$, $\mathcal R$, $\mathcal T$, $\mathcal W$, $\mathcal Y$, $\mathcal Z$ ; accent $\check\cdot$. Depuis le 2026-10-21, $\eta$ (bruit d'évolution) et $\dot\cdot$ (vitesse) sont pris. Je laisse de côté $\iota,\vartheta,\varkappa,\varpi,\varrho,\varsigma$, trop proches de $\iota,\theta,\kappa,\pi,\rho,\sigma$. Le préambule charge déjà `amsmath`, `amsfonts`, `amssymb` : `\boldsymbol` et `\check` sont disponibles, ni `bm` ni `mathrsfs` ne sont nécessaires.

---

## 4. Collisions à lever (nature différente, objets globaux)

### C1 — `k` : noyau / indice de temps (≈ 700 / ≈ 2 500 occurrences)

- **Noyau** (fonction à deux arguments) : `3_Gaussian_processes:52`, `3_kriging`, `3_reduced_rank:28`, `3_RKHS:26`, `4_kriging_navigation:196` ($k_\theta$), `4_subsampling:53–60` ($k_{\vert r}$, $k_{d\vert r}$), preuves des annexes B du manuscrit.
- **Indice de temps** : tout le filtrage (ch.2, KKF, ch.4) et les phrases « at time $k$ ».
- Ils se côtoient dans la même ligne : `z_k=h(x_k)+n_k, h\sim\mathcal{GP}(m,k)` (`4_kriging_navigation:24`).
- **Proposition ✔** : noyau → **$\kappa$** ($\kappa_\mathrm{SE},\kappa_\nu,\kappa_\mathrm{sep},\kappa_{\mathbf B},\kappa_\theta,\kappa_x,\kappa_z,\kappa_X,\kappa_Z,\bar\kappa,\kappa_{\vert\mathrm{old}},\kappa_{\mathrm{UK}\vert\mathrm{old}},\kappa_{\vert r},\kappa_{d\vert r},\kappa_{r_\mathrm{KL}}$ ; $\mathcal H_\kappa$ ; $\mathcal{GP}(m,\kappa)$). Le $\kappa$ actuel (profil radial, `3_Gaussian_processes:172–177`) devient $\kappa_\mathrm{rad}$. « $k$-means » ne change pas (nom propre). **Fait le 2026-10-22 sur son « kappa pour les noyaux de covariance »** : lot W, le profil radial de la définition 3.3 passe d'abord en $\kappa_\mathrm{rad}$ (3 atomes : « on identifie $\kappa_\mathrm{rad}$ avec $\kappa$ ») ; lot X, **1 056 atomes** $k\to\kappa$ dans les chapitres 3, 4 (§4.2, §4.5) et les annexes B : appels $\kappa(x,x')$, noyaux indicés ($\kappa_\nu$, $\kappa_\mathrm{SE}$, $\kappa_\mathbf B$, $\kappa_\theta$, $\kappa_{|r}$, $\kappa_{|\mathrm{old}}$, $\kappa_{\mathrm{UK}|\mathrm{old}}$, $\kappa_{r_\mathrm{KL}}$, $\kappa_x$, $\kappa_z$, $\kappa_X$, $\kappa_Z$, $\kappa_{3/2}$ …), formes par blocs $\bar\kappa$, RKHS $\mathcal H_\kappa$, $\mathcal H_{\kappa_X}$ ; l'indice de temps $k$ (2 467 atomes) n'est pas touché, ni « $k$-means ». **Noyau ou fonction de covariance (sa consigne : « les noyaux sont des fonctions à 2 paramètres, sauf les cas radiaux par abus de langage, donc $k_{d|r}$ n'est pas un noyau mais une fonction de covariance »)** : $\kappa_{d|r}$ (et $\kappa^0_{d|r}$, $\kappa^{\max}_{d|r}$) garde la lettre $\kappa$ mais se dit « residual covariance function », comme le fait déjà le §4.5 (définition 4.1) ; $\kappa_{|r}$ est « the residual kernel » ; aucune phrase n'a été touchée, et la recherche de « kernel function » ne trouve rien.
- Pourquoi dans ce sens : renommer le noyau coûte près de 4 fois moins que renommer l'indice, et $k$ = temps est la convention de tous les algorithmes de filtrage. Chapitre 5 : $k_\nu$ (10 occurrences), $k_x$, $k_z$ à suivre.

### C2 — `w` : bruit d'évolution / poids de particules (≈ 70 / ≈ 190)

- **Bruit d'évolution** (vecteur) : $w_k$, $w_k^n$, $w_k^l$, $\overline w_k^l$, $d_w$ (`2_bayesian_estimation:15,19,316,339`, `2_particle_filtering:300–329`, preuves KF/RBPF de `9_proofs_1`).
- **Poids** (scalaires) : $w_k^{(i)}$, $\tilde w$, $w^{(i)}$ (MC, PF, AK). $\mathbf w$ (poids plongés du KKF) est de même nature : ✔.
- Deux collisions « de proximité » : dans l'EnKF, un échantillon de bruit s'écrit $w_k^{(i)}$ (`2_bayesian_estimation:507`), exactement comme un poids ; dans le RBPF, le bruit $w_k^n$ (`2_particle_filtering:300`) et le poids $w_k^{n(i)}$ (`:393`) ne diffèrent que par « (i) ».
- **Décidé et appliqué le 2026-10-21** (ta règle : un bruit est une lettre grecque, $w$ ne peut pas être à la fois un poids et un bruit) : bruit d'évolution → **$\eta_k$** ($\eta_k^n$, $\eta_k^l$, $\overline\eta_k^l$, $\eta_k^{(i)}$, $d_\eta$) ; bruit d'observation $n_k$ → **$\varepsilon_k$** (`\varepsilon`, $d_\varepsilon$, $\boldsymbol\varepsilon_{1:k}$), qui retrouve le $\varepsilon$ du krigeage du ch.3 et libère $n$ (entier, marqueur) ; seuils et correction numérique → **$\epsilon$** (`\epsilon`) ; $w$ ne désigne plus que des poids ; la vitesse de l'état s'écrit $\dot p$ (L6). La mesure instrumentale ne peut plus prendre $\eta$ : voir C4. *(Proposition initiale, remplacée : $v_k$, avec $\eta_k$ en alternative.)*
- Chapitre 5 : `w` y est un poids, rien à changer.

### C3 — `μ` : mesure / plongement / moyenne / particule auxiliaire

| Objet | Nature | Où | Occ. | Sort |
|---|---|---|---|---|
| mesure cible $\mu$, $\tilde\mu$, $\mu_k$ | mesure | MC, PF | ≈ 84 | → **$\mu_p$**, $\tilde\mu_p$, $\mu_{p,k}$ ✔ **appliqué le 2026-10-21** (voir C4) ; **gardé le 2026-10-22** quand il a repris la surcharge $\mu$ mesure / KME : la note de Le Gland (source du chapitre 2) écrit $\mu$, $\nu$, $S^N(\mu)$, $\langle\mu,\varphi\rangle$ comme nous, et $\mathbb P,\mathbb Q$ avait déjà été écarté (C4) ; ne pas reproposer |
| plongement $\mu_X,\mu_{XZ},\mu_{Z\vert X},\hat\mu_{XX},\mu_{k\vert k}$ | plongement d'une loi (représentant de $\langle\mu,\cdot\rangle$ dans $\mathcal H$) | RK, KK | ≈ 112 | ✔ reste : c'est la **même** famille que la mesure (le KME de $\mu$ en est le représentant de Riesz), R2 est respectée |
| moyenne d'un vecteur | vecteur | `3_Gaussian_processes:35–37`, `2_bayesian_estimation:473–490` (UKF), `9_appendix:11–67` ($\mu_U,\mu_V,\mu_{U\vert V}$) | ≈ 18 | **reste $\mu$** ✔ **décidé le 2026-10-22** (« gardons mu comme ça » : la paire $(\mu,\Sigma)$ remplace $(x,P)$, comme les plongements ; la nomenclature écrit $\mathcal N(\mu,\boldsymbol\Sigma)$) ; GP:35 voir L9 |
| particule auxiliaire $\mu_k^{(i)}$ | vecteur | `2_particle_filtering:236–276` | ≈ 14 | → **$a_k^{(i)}$** ✔ **appliqué le 2026-10-21** (tirage de la densité de transition ou sa moyenne : « auxiliaire ») |
| poids auxiliaires $\lambda_k^{(i)},\tilde\lambda_k^{(i)}$ | scalaires | `2_particle_filtering:247–277` | ≈ 8 | → **$v_k^{(i)}$**, $\tilde v_k^{(i)}$ ✔ **appliqué le 2026-10-21** (même convention de tilde que $w$) ; libère $\lambda$ de ce sens |

### C4 — `ν` : mesure instrumentale / régularité de Matérn (≈ 65 / ≈ 50)

- **Mesure instrumentale** : `2_Monte_Carlo:85–173`, PF ($\nu$, $\tilde\nu$, $\nu_k$, $C_\nu$, $\mathbb E_\nu$, $\mathbb V_\nu$, $\mathrm d\mu/\mathrm d\nu$).
- **Régularité** (scalaire) : $k_\nu$, $K_\nu$, $S_\nu$, $p_\nu$ (`3_Gaussian_processes:184–231`, `9_proofs_1:172–230`, chapitre 5).
- **Décidé et appliqué le 2026-10-21** (puisque $\eta$ est pris par le bruit) : dans le Monte-Carlo (MC et PF), mesure cible $\mu\to$ **$\mu_p$** (densité $p$) et mesure instrumentale $\nu\to$ **$\mu_q$** (densité $q$) ; constantes $C_\mu,C_\nu\to C_p,C_q$ ($\tilde\mu_p=C_p\mu_p$, $\tilde q=C_q q$) ; $\nu$ ne désigne plus que la régularité de Matérn, canonique et présente jusqu'au ch.5. Même lettre et rôle en indice, comme $h$ / $h_\mathrm{MS}$. *(Proposition initiale, remplacée : $\eta$.)*
- *Alternatives écartées :* (b) cible et instrumentale deviennent $(\mathbb P,\mathbb Q)$ (réécrit aussi $\mu$ et empiète sur le $\mathbb P$ des événements) ; (c) la régularité devient $s$ (peu naturel : $K_s$ pour Bessel). Variances alignées sur les espérances (lot E) : $\mathbb V_p$, $\mathbb V_q$ ; seuls restent en indice de mesure $\mathbb E_{\mu_{q,k}}$ et $\mathbb V_{\mu_{q,k}}$ (variance des poids de trajectoire, §2.4). *Rouvert le 2026-10-22 (« $\mu$ surchargé entre les mesures de proba et le KME ») et clos : il garde $\mu$ pour les trois (voir C3).*

### C5 — `α`, `β` : vecteurs / facteurs scalaires

| Objet | Nature | Où |
|---|---|---|
| $\alpha$ (coefficients KL ; coefficients de la carte au ch.5) | vecteur | `3_reduced_rank:203–210`, ch.5 |
| $\alpha_\mathrm{reg}$, $\alpha_\mathrm{cov}$ | scalaires | PF:187–201, AK:104–171, KN:349–382 |
| $\beta$ (coefficients de tendance UK) | vecteur | `3_kriging:262–271`, B1 |
| $\beta$ (facteur d'oubli du RLS, ch.5) | scalaire | `5_SLAM` ; la nomenclature a deux entrées $\beta$ |

**Décidé le 2026-10-21 et 22 (la proposition initiale, $\gamma_\bullet$ pour les facteurs, est remplacée) :** les coefficients multiplicateurs sont tous en **$\alpha_\bullet$** : $\alpha_\mathrm{reg}\to\alpha_\mathrm{RPF}$ (9 atomes), $\alpha_\mathrm{cov}\to\alpha_\mathrm{AK}$ (7), $\lambda_\mathrm{RLS}\to\alpha_\mathrm{RLS}$ (nomenclature ; ch.5 plus tard), $\lambda_x,\lambda_z\to\alpha_{\mathrm{RKHS},x},\alpha_{\mathrm{RKHS},z}$ (30) et, par analogie, le $\lambda$ nu de régularisation du RKHS et du KKF $\to\alpha_\mathrm{RKHS}$ (31) ; le facteur d'oubli du RLS est $\beta_\mathrm{RLS}$ (nomenclature ; ch.5 plus tard) ; le taux du tirage aléatoire est **$p_\mathrm{RS}$** (« pas fan de alpha pour RS », 4 atomes) et l'échelle de l'UKF **$\ell_\mathrm{UKF}$** (« un coefficient de dilatation », 5). Les **coefficients des fonctions de base** prennent un même symbole **$\boldsymbol\gamma$** (« va pour gamma » : $\beta$ de la tendance de l'UK, 60 atomes, et $\alpha$ de Karhunen–Loève, 4 ; composantes $\gamma_i$, $\boldsymbol\gamma_0$, $\boldsymbol\gamma^\top$ ; la covariance a priori reste $\mathbf B$). **Reste :** les coefficients de la carte du ch.5 ($\alpha$, $\hat\alpha_k$, $\mathbf P_k^\alpha\to\boldsymbol\gamma$, $\hat{\boldsymbol\gamma}_k$, $\mathbf P_k^\gamma$, déjà annoncés dans la nomenclature) et le scalaire $\alpha$ de l'annexe A, en collision avec la famille $\alpha_\bullet$ ($s$, L2).

### C6 — `ξ` : variables gaussiennes / multiplicateur de Lagrange

- $\xi_i\sim\mathcal N(0,1)$ (coefficients de Karhunen–Loève, `3_reduced_rank:165`, `9_proofs_2`) et $\xi\sim\mathcal N(0,\mathbf I_{d_z})$ (`4_subsampling:477`) : même nature ✔.
- $\xi(x)$ multiplicateur de Lagrange, **matrice** (`4_subsampling:92–111`, `9_proofs_1:220–228`), dans la même section que l'autre $\xi$.
- **Proposition ✔** : multiplicateur → **$\boldsymbol\Xi(x)$** (capitale grecque libre ; matrice donc grasse). **Fait le 2026-10-22 (lot Y, 28 atomes)** : `4_subsampling` 6 ($\boldsymbol\Xi(x)$ du §4.5), `9_proofs_1` 18 (preuve du krigeage ordinaire et de l'UK), `9_proofs_2` 4 ; les $\xi_i$ de Karhunen–Loève et le vecteur $\xi\sim\mathcal N(0,\mathbf I_{d_z})$ de la borne en $\chi^2$ restent. La nomenclature a maintenant $\xi_i$ et $\boldsymbol\Xi$.

### C7 — `𝒰` : ensemble d'entrée / loi uniforme

- $\mathcal U$, $\mathcal U_X$, $\mathcal U_Z$ : ensembles d'entrée (≈ 160 occurrences, ch.3 et annexes).
- $\mathcal U([a,b])$ : loi uniforme (`2_Monte_Carlo:436–454`, nomenclature, ch.5).
- **Proposition ✔** : loi uniforme → $\operatorname{Unif}([a,b])$ (3 occurrences + ch.5). **Remplacée le 2026-10-22 sur sa question « que penses-tu de $\mathfrak U$ pour les ensembles de départ des GP ? ou alors $\Omega$ ? » puis « go pour Omega »** : ce sont les ensembles qui changent, $\mathcal U$ → **$\Omega$** (lot Z, 215 atomes hors chapitre 5 : `3_RKHS` 41, `3_reduced_rank` 40, `9_proofs_2` 88, `3_Gaussian_processes` 26, `4_subsampling` 8, `3_kriging` 6, `9_proofs_1` 5, `3_KKF` 1 ; $\mathcal U_X,\mathcal U_Z\to\Omega_X,\Omega_Z$ ; la cellule de partition $\mathcal U_i^{(n)}$ de la preuve de Mercer devient $\Omega_i^{(n)}$), la loi uniforme garde $\mathcal U([a,b])$ ; raisons : $\Omega$ est la lettre de la méthode de l'espace de Hilbert qu'il suit (Solin et Särkkä, « Let $\Omega\subset\mathbb R^d$ », 22 occurrences dans son dossier), elle se lit bien dans $\int_\Omega$, $\mathcal L^2(\Omega,\mathbb R^{d_z})$ et $\Omega^n$, le gothique $\mathfrak U$ est peu lisible en indice, et le manuscrit n'écrit jamais d'espace probabilisé $(\Omega,\mathcal F,\mathbb P)$ ; seule la matrice polaire $\boldsymbol\Omega$ de `2_Monte_Carlo:292` (en gras) porte déjà la lettre. **Chapitre 5 :** ses 11 $\mathcal U$ (ensembles) suivront à sa passe ; la loi uniforme de `5_SLAM:189` ne change pas.

### C8 — `𝒢` : matrice écrite en calligraphique

- La matrice de Gram par blocs de l'UK (`3_kriging`, B1, `4_subsampling`) est calligraphique alors que le calligraphique désigne des ensembles et des opérateurs ; la nomenclature la range avec $\bar k$ et $\bar{\mathbf z}$ (« forme par blocs »).
- **Proposition ✔** : $\mathcal G\to\bar{\mathbf G}$ (≈ 27 occ.), cohérent avec $\bar\kappa$ et $\bar{\mathbf z}$. $\mathcal{GP}$ n'est pas touché. **Fait le 2026-10-22 (lot V, 66 atomes, nomenclature comprise)** quand il l'a rappelé en voyant la ligne $\mathcal G$ à côté de $\bar k$ et $\bar{\mathbf z}$ (« tu n'aimais pas le cursif pour cette matrice non ? ») : `3_kriging` 6, `4_subsampling` 34 (dont $\mathcal G_h\to\bar{\mathbf G}_h$, blocs parqués), `9_proofs_1` 18 (dont $\mathcal G_{\mathbf B}$), `9_proofs_2` 7.

### C9 — `l` : longueur de corrélation / indice / marqueur

- $l$ longueur de corrélation de Matérn (`3_Gaussian_processes:184–231,283`, `9_proofs_1`, AK:60) ; $l$ aussi indice de lissage (`2_bayesian_estimation:40`), composante (`3_Gaussian_processes:276`) et marqueur de la partie linéaire ($x_k^l$, $\mathbf Q_k^l$…).
- **Proposition ✔** : longueur de corrélation → **$\ell$** (≈ 40 occ., convention des processus gaussiens). Le marqueur $l$ ne change pas. **Fait le 2026-10-21 (lot F, 26 atomes : « tu peux prendre \ell si tu veux »)**, puis étendu le 22-10 aux largeurs de bande et à l'échelle de l'UKF (C10, C5).

### C10 ◆ — `h` : fonction d'observation / largeur de bande

Tu cites $h$ / $h_\mathrm{MS}$ comme exemple acceptable (deux largeurs de bande, je suppose). Reste que $h$ est aussi la fonction d'observation $h_k$ (≈ 70 occ.) et la carte du ch.4 ($h\sim\mathcal{GP}$). R2 voudrait qu'un des deux change. Le contexte lève presque toujours le doute (la fonction porte un indice $k$, la largeur de bande jamais) mais je n'ai aucune lettre libre qui s'impose pour la largeur de bande. **Décidé le 2026-10-22 (reco initiale, garder les deux, remplacée) :** la largeur de bande passe en **$\ell$** (lot S, 65 atomes : $K_\ell$, $S_{K_\ell}^N$, $\tilde K_\ell$, $\ell_\mathrm{Epa}$, $\ell_\mathrm{Gauss}$, le $h$ des quatre boîtes d'algorithme du chapitre 4, le mean shift $h_\mathrm{MS}\to\ell_\mathrm{MS}$), sur sa proposition (« tu peux appeler la bandwidth \ell ? il me semble que c'est standard dans nos notations ») ; $\alpha_\mathrm{RPF}$ reste le coefficient dans $[0,1]$ qui la multiplie ; $h$ n'est plus que la fonction d'observation et le GP. Son message précédent (« tu n'as pas renommé h en alpha_RPF dans la présentation du RPF ») ne pouvait pas s'appliquer tel quel : $\alpha_\mathrm{RPF}\,h$ serait devenu $\alpha_\mathrm{RPF}\,\alpha_\mathrm{RPF}$.

---

## 5. Collisions levées par la charte seule (aucun changement de lettre)

Depuis l'étape 1 (matrices en gras), ces homonymes sont des noms distincts (police différente, nature différente).

| Lettre | Matrice (gras) | Autres objets (italique / calligraphique) |
|---|---|---|
| K | gain de Kalman $\mathbf K_k$ (BE, PF, KKF) | noyau de régularisation $K,K_h$ (PF, AK) ; Bessel $K_\nu$ (GP, B1) ; opérateur $\mathcal K$ |
| S | covariance d'innovation $\mathbf S_k,\mathbf S_{d\vert r},\mathbf S_{\mathrm{new}\vert\mathrm{old}}$ ; complément de Schur (annexe A) | mesure de Monte Carlo $S^N(\mu)$ ; densité spectrale $S_\nu$ ; ensembles $\mathcal S$ |
| T | transition plongée $\mathbf T_k$ (KKF) | opérateur $T$ (RR:12 → $\mathcal T$, L5) ; itérations $T$ (AK:370 → $n_\mathrm{it}$) |
| O | observation plongée $\mathbf O$ (KKF) | Landau $\mathcal O(\cdot)$, $o(\cdot)$ |
| J | jacobienne $\mathbf J_f$ (BE) | coût $J$ (RK:169) |
| F | évolution $\mathbf F_k$ (+ RBPF) | espace normé $F$ (BE:64 → L7) ; vecteur aléatoire $F$ (GP:63) ; $\mathcal F$ |
| C | contrainte $\mathbf C=\mathbf 1\otimes\mathbf I$ (KR:212) | nombre de classes $C,C_k$ (AK) ; constantes $C_\mu,C_\eta$ ; ensembles $\mathcal C$ |
| U, V | poids UK $\mathbf U$ (B1:418) ; SVD (MC:287) ; $\mathbf V$ (KKF) | vecteurs aléatoires $U,V$ (annexe A) ; $\mathbb V$, $\mathcal V$ |
| G | Gram $\mathbf G,\mathbf G_\theta,\mathbf G^z,\mathbf G^x$ | pénalité du représentant $G$ (RK:97 → L1) |
| W | poids du krigeage $\mathbf W(x)$ | poids de l'UKF $W^{(i)}$ (BE:476 → $w^{(i)}$, L3) |
| A, B | blocs, facteurs $\mathbf A_k$ | paire d'estimateurs $A,B\in\mathcal F$ (BE:81,151,229 → L8) |
| Q | covariance $\mathbf Q_k$ | opérateur gain $\mathcal Q_k$ (KKF) |
| E | — | $\mathbb E$ ; $E$ espace mesuré (BE:64) ; $E_\mathrm{in}$ ; $E$ perte (RK:97) |

---

## 6. Collisions locales (énoncés, preuves, annexes)

Onze retouches, toutes à ✔ par défaut.

| # | Lieu | Collision | Proposition |
|---|---|---|---|
| L1 | `3_RKHS:96–143` théorème du représentant | $E:\mathbb R^{\tilde N}\to\mathbb R$ (perte) et $G:\mathbb R_+\to\mathbb R$ (pénalité) : $G$ est la Gram, $E$ l'espace mesuré | noms droits $\operatorname{loss}$ et $\operatorname{pen}$ |
| L2 | `9_appendix:11–69` et `186–195` | $\mu_U,\mu_V,\mu_{U\vert V}$ (moyennes) ; $\alpha$ (valeur scalaire d'une forme quadratique, `:44`) ; $\mathbf P,\mathbf Q\in\mathcal M_{n,p}$ pour l'identité de commutation | $m_U,m_V,m_{U\vert V}$ ; $s$ ; $\mathbf X,\mathbf Y$ |
| L3 | `2_bayesian_estimation:473–513` UKF et EnKF | $W^{(i)}$ (poids, majuscule) ; $w_k^{(i)}$ (échantillon de bruit) ; $\mathcal X^{(i)}$ (points écrits en calligraphique : toléré, notation canonique de l'UKF) | $w^{(i)}$ (à faire) ; $\eta_k^{(i)}$ (**fait**, C2) |
| L4 | `2_particle_filtering:299–362` RBPF | systèmes $(S),(\bar S),(S_1),(S_2)$ à côté des covariances d'innovation $\mathbf S_k^{n(i)},\mathbf S_k^{l(i)}$ | labels droits $(\mathrm{RB})$, $(\overline{\mathrm{RB}})$, $(\mathrm{RB}_1)$, $(\mathrm{RB}_2)$ |
| L5 | `3_reduced_rank:12`, `4_adaptive_kriging:370` | $T$ opérateur et $T$ nombre d'itérations, face à $\mathbf T$ (transition du KKF) | $\mathcal T$ ; $n_\mathrm{it}$ |
| L6 | `4_kriging_navigation:324,373` | vitesse $v_k^{(x)},v_k^{(y)}$ (face aux poids auxiliaires $v$) ; position $p_k^{(x)}$ face à la densité $p$ | $\dot p_k^{(x)},\dot p_k^{(y)}$ ; la position reste $p$ (locale au tableau) (**fait**) |
| L7 | `2_bayesian_estimation:64–71` | $E$, $F$, $m$ : espace mesuré, espace normé, mesure ($m$ est la moyenne ; $E$ et $F$ sont ailleurs une perte et une matrice) | $\mathcal L^p(\mathcal D,\mathcal A)$, la mesure de $\mathcal D$ étant notée $\mathrm dx$ |
| L8 | `2_bayesian_estimation:81,151,229` | $(A,B)\in\mathcal F^2$ : des fonctions nommées comme des matrices | $(\psi_1,\psi_2)$ |
| L9 | `3_Gaussian_processes:35–40` | $\mu\in\mathbb R$, $\Sigma\in\mathbb R_+$ : $\Sigma$ est ailleurs une matrice, $\mu$ une mesure | $Y\triangleq\sum_i a_i^\top f(x_i)\sim\mathcal N(m_Y,\sigma_Y^2)$ |
| L10 | `3_kriging:136`, `4_subsampling:28–44`, `3_RKHS:52` | $e_i$ (base canonique), $e_r,e_f$ (erreurs), $e_x$ (évaluation) | $\mathbf e_i$ ; $e_r,e_f$ restent ; $e_x\to\operatorname{ev}_x$ |
| L11 | `4_adaptive_kriging:9` | $p\in[0,1]$ probabilité de Bernoulli de l'échantillonnage aléatoire, face à la densité $p$ | $\gamma_\mathrm{RS}$ (famille des facteurs de réglage, C5) |

Deux cas laissés tels quels : $\delta\triangleq\max(0,-\lambda_{\min})+\varepsilon$ de la correction SPD (`2_Monte_Carlo:294`, scalaire local) face au Dirac $\delta_x$ (toujours indicé) ; $y$ variable d'intégration (`2_particle_filtering:143`, `3_reduced_rank:30`) face à l'innovation $y_k$ (toujours indicée).

---

## 7. Homonymes de même nature et points à surveiller

**Acceptés** (R2 respectée, aucune action) : $x$ (état, point d'entrée du GP, point de requête) ; $z$ (observation, valeur de requête, valeur d'apprentissage) ; $f$ (fonction d'évolution $f_k$ ; élément générique de $\mathcal H$) ; $h$ (fonction d'observation et, depuis le 2026-10-21, processus gaussien du ch.3 : un seul nom pour la carte du ch.2 au ch.4) ; $m$ (fonction moyenne, moyenne constante $\hat m_\mathrm{OK}$, paramètre de $\mathcal N(m,\mathbf P)$) ; $\boldsymbol\Sigma$ (covariances) ; $\mathbf G$ (Grams) ; $\mathbf R$ (covariances de bruit d'observation) ; $\sigma$ (écarts-types) ; $N,\tilde N$ (effectifs) ; $\varphi$ (fonctions test) ; $\mathbf A$ (racines carrées d'une covariance : UKF, RPF, GP) ; $\delta_x$ (Dirac) ; $\Delta_g,\Delta_t$ (pas).

**À surveiller** (même type, concepts différents ; aucune modification proposée, dis-moi si tu veux être strict, D14) :

| Lettre | Concepts | Remarque |
|---|---|---|
| $\lambda$ | valeurs propres $\lambda_i$ ; coefficients de régularisation $\lambda_x,\lambda_z,\lambda_\mathrm{RLS}$ ; échelle de l'UKF | en RK, $\lambda$ (ridge) et $\lambda_i$ (Mercer) cohabitent. Option : UKF → $\lambda_\mathrm{UT}$ |
| $\psi,\phi,\varphi$ | estimateurs ($\hat\psi_\mathrm{MMSE}$), fonctions propres $\psi_i$, éléments de $\mathcal H$, fonctions test, applications d'état | $\phi$ et $\varphi$ sont deux glyphes d'une même lettre (RR vs ch.2) |
| $r$ | rayon ($r_\mathrm{cov},r_{\min},r_\mathrm{DB}$, argument de Matérn) ; entiers ($r_\mathrm{UK},r_\mathrm{KL},r_\mathrm{RFF}$) ; marqueur « retenu » | rayon vs entier : même cas que $h$ |
| $n,m,p$ | entiers génériques ; mesure ou moyenne ; densité | le bruit $n_k$ est devenu $\varepsilon_k$ (2026-10-21) |
| $\varepsilon$, $\epsilon$ | bruits $\varepsilon$ (observation $\varepsilon_k$, entraînement $\varepsilon_i$, tirage du RPF $\varepsilon^{(j)}$) ; seuils et correction numérique $\epsilon$ (arrêt de l'IEKF BE:463, MAP BE:229, décalage SPD MC:294) | depuis le 2026-10-21 ; deux glyphes d'une même lettre que le contexte sépare (vecteur bruit / scalaire positif) |
| $\mathbf R_\bullet$ | bruit d'observation $\mathbf R$ ; covariances d'erreur du krigeage $\mathbf R_\mathrm{SK}$ | même nature (covariances) |
| $\tilde S^N$ | auto-normalisé, alors que $\tilde p,\tilde w$ = non normalisé | sens inverse du tilde |

---

## 8. Nomenclature

**Symboles utilisés et absents de la liste** (vérifié par lettre) : $\mathcal H_\kappa$ (RKHS), $\mathcal F$, $\mathcal V$, $\mathcal C$, $\mathcal Q$, $\xi_i$ et $\boldsymbol\Xi$, $\chi^2$, $\tau$, $\Delta_g$, $\Delta_t$, $\sigma_f$ (seul $\sigma_r$ y est), $\ell$, la largeur de bande $h$, $L$ (log-vraisemblance), $\Theta$, $\mathbf A_k$, $c_{d_x}$, $\mathbf W$, $\Omega$, $\mathcal X$, $W^{(i)}$.

**Entrées à corriger** : les deux lignes $\beta$ (C5) ; $\alpha$ et $\alpha_\mathrm{reg},\alpha_\mathrm{cov}$ ; $\mu$ (moyennes seules à corriger, C3) ; $\mathcal U([a,b])$ (C7) ; $\mathcal G,\bar k,\bar{\mathbf z}$ (C1, C8) ; $\mathfrak D_G$ reste (retour prévu en v2). **Étape 2 :** les symboles $\eta_k,d_\eta$, $\varepsilon_k,d_\varepsilon$, $a_k^{(i)}$, $v_k^{(i)}$, $\mu_p$, $\mu_q$, $C_p,C_q$ ont été remplacés sur place ; leur rang alphabétique, le doublon $\varepsilon,\varepsilon_i$ / $\varepsilon_k$ et l'entrée manquante de $\epsilon$ sont à reprendre avec la réécriture.

**Structure proposée** : regrouper par nature plutôt que par ordre alphabétique mêlant grec et latin, avec en tête un encart « Conventions » (annexe II de cette note) :
1. Conventions typographiques (une demi-page).
2. Opérateurs et notations générales.
3. Ensembles et espaces.
4. Entiers, dimensions, indices.
5. Scalaires et paramètres.
6. Vecteurs et n-uplets.
7. Matrices.
8. Fonctions, noyaux, opérateurs.
9. Mesures, lois, plongements.
10. Estimateurs.

---

## 9. Chapitre 5, conclusion, code

- **Chapitre 5.** Les renommages C1 ($k_\nu,k_x,k_z$), C5 ($\beta\to\gamma_\mathrm{RLS}$), C7 ($\mathcal U\to\operatorname{Unif}$) et les gras ($\mathbf P^\alpha,\mathbf K^\alpha,\mathbf S,\boldsymbol\Psi,\boldsymbol\Lambda$) s'y appliquent, mais ne sont **pas faits** : l'étape 1 s'est arrêtée avant ce chapitre ; C2 ne le touche pas. Aucun nom proposé n'y est pris : `v`, `η`, `γ`, `ℓ`, `Ξ`, `\check`, `\dot` y sont absents (inventaire du chapitre 5 : 1 533 relevés, 61 lettres).
- **Blocs `comment` de la v2.** Le gras des matrices y est déjà appliqué (429 occurrences). Les bornes calculables de §4.5 (`4_subsampling:506–662`) et leurs preuves (`9_proofs_2:435–628`) emploient comme **constantes scalaires** $\alpha$, $\kappa$, $\eta$ (bornes sur $\|\kappa(\tilde{\mathbf x}^d,\tilde{\mathbf x}^r)\|_{\max}$, $\|\kappa(\tilde{\mathbf x}^r,x)\|_{\max}$, $\|\kappa(\tilde{\mathbf x}^d,x)\|_{\max}$), ainsi que $\gamma$ (échelle du bruit, `:551`), $\rho_r$, $k_{\max}$, $g_{\max}$, $\lambda_j^\pm$, $\lambda_r^g$ : $\kappa$ et $\eta$ entreraient en collision avec C1 et C4, et $\alpha$ avec C5 (vecteur). À renommer au retour de la v2 (par exemple une famille $\mathrm c_{dr},\mathrm c_{\mathrm{in}},\mathrm c_{dx},\mathrm c_{\max}$ pour les bornes d'entrée, et $s$ pour l'échelle du bruit) ; je peux le faire dans les `comment` dès maintenant si tu préfères que tout soit déjà cohérent. Même remarque pour le bloc commenté de `5_SLAM:466–509`, où $\gamma$ est le niveau de dégradation de la carte (à distinguer de $\gamma_\mathrm{RLS}$).
- **Conclusion et introduction** : sans formules, rien à faire.
- **MATLAB** : les noms de variables sont indépendants de la notation ; seuls les libellés LaTeX des figures (légendes, axes) seraient à relire.

---

## 10. Plan de migration

*Les points 0 à 6 sont les lots du plan ; l'« étape 1 » du journal (audit des marqueurs, puis matrices en gras) correspond au point 2, exécuté en premier à ta demande.*

0. **Gel.** Tu commites l'état actuel (je ne commite jamais). Le renommage se fait dans la copie de travail, un lot à la fois. (L'étape 1 a été appliquée sans commit de ma part : voir « Retour arrière » dans le journal.)
1. **Renommages de lettre**, du moins risqué au plus risqué, chacun par script d'après l'inventaire (on cible l'atome, jamais le texte), diff relu, compilation, 0 référence non définie, comparaison du texte extrait du PDF avant/après (seuls les symboles changent) : $\operatorname{Unif}$ (C7) ; $\ell$ (C9) ; $\boldsymbol\Xi$ (C6) ; $\bar{\mathbf G}$ (C8) ; $\gamma$ (C5, L11) ; moyennes et APF (C3) ; $\eta$ (C4) ; $v$ (C2, L6) ; $\kappa$ (C1).
2. **Gras des matrices** ($\mathbf R$, $\mathbf P$, $\mathbf Q$, $\mathbf F$, $\mathbf H$, $\mathbf K$, $\mathbf S$, $\mathbf G$, $\mathbf I$, $\mathbf A$–$\mathbf D$, $\mathbf M$, $\mathbf U$, $\mathbf W$, $\boldsymbol\Sigma$, $\boldsymbol\Lambda$, $\boldsymbol\Psi$, $\boldsymbol\Pi$, $\mathbf O$, $\mathbf T$, $\mathbf J$) : **fait le 2026-10-21**, d'un seul tenant plutôt que famille par famille, avant les renommages de lettre du point 1 (qui se feront donc sur du texte déjà en gras). Contrôles réalisés : chapeau/tilde sur gras (`\hat{\mathbf P}`), largeur des algorithmes (trois boîtes trop larges, résorbées), pages types relues en image (journal).
3. **Retouches locales** L1 à L11.
4. **Nomenclature** réécrite (§8) et vérifiée contre l'inventaire (aucun symbole manquant parmi ceux qui apparaissent dans au moins deux sections).
5. **Chapitre 5**, relecture Word en-GB (0 nouveau signalement), notes de mémoire du projet mises à jour.
6. **Contrôle final.** Je relance l'inventaire et vérifie programmatiquement qu'aucune des neuf collisions n'existe plus.

Risques des lots restants : un renommage peut déplacer des coupures de ligne dans les algorithmes et les équations, et la pagination peut bouger d'une ou deux pages. Pour le gras (point 2) : trois boîtes trop larges, résorbées ; aucun défaut visible sur les chapeaux des pages relues ; 223 pages avant et après. Rien d'irréversible tant que tu commites avant chaque lot.

---

## 11. Décisions

À cocher dans la note, ou à me répondre par numéro.

- [x] **D1 ◆ Matrices en gras majuscule** (2 613 occurrences dans le texte visible, 429 de plus dans les blocs `comment`). **Fait le 2026-10-21**, voir le journal.
- [x] **D2 Noyau $k\to\kappa$**, profil radial $\kappa_\mathrm{rad}$ (C1). **Fait le 2026-10-22** (lots W et X).
- [x] **D3 ◆ Bruit d'évolution** : $\eta_k$ (avec $\varepsilon_k$ pour l'observation et $\epsilon$ pour les seuils, D15). **Fait le 2026-10-21**, voir le journal.
- [x] **D4 ◆ Mesure instrumentale** : $\mu_q$, cible $\mu_p$ dans le Monte-Carlo, $\nu$ gardé pour Matérn (C4). **Fait le 2026-10-21.**
- [x] **D5 Moyennes** $\to m$ : **abandonnée le 2026-10-22**, elles restent $\mu$ (paire $(\mu,\Sigma)$ ; « 2) oui » à $\mathcal N(\mu,\boldsymbol\Sigma)$ dans la nomenclature). Particule auxiliaire $a_k^{(i)}$ et poids auxiliaires $v_k^{(i)}$ : **faits le 2026-10-21**.
- [x] **D6 Coefficients multiplicateurs en $\alpha_\bullet$** : $\alpha_\mathrm{RPF}$, $\alpha_\mathrm{AK}$, $\alpha_\mathrm{RLS}$, $\alpha_\mathrm{RKHS}$, $\alpha_{\mathrm{RKHS},x/z}$ (C5). **Fait le 2026-10-21** (ch.5 plus tard) ; le taux du tirage aléatoire est $p_\mathrm{RS}$ et l'échelle de l'UKF $\ell_\mathrm{UKF}$ (D18).
- [x] **D7 Longueur de corrélation** $l\to\ell$ (C9). **Fait le 2026-10-21 (lot F)**, étendu aux largeurs de bande et à l'UKF le 22-10 (D18).
- [ ] **D8 ✔ Multiplicateur** $\boldsymbol\Xi$, **Gram par blocs** $\bar{\mathbf G}$, **loi uniforme** $\operatorname{Unif}$ (C6 à C8). **Faits le 2026-10-22 : $\bar{\mathbf G}$ (lot V) et $\boldsymbol\Xi$ (lot Y)** ; $\operatorname{Unif}$ abandonné au profit de $\Omega$ pour les ensembles (C7, lot Z).
- [ ] **D9 ✔ Retouches locales** L1 à L11.
- [x] **D10 $h$ et largeur de bande** (C10) : la largeur de bande passe en $\ell$ (D18) ; $h$, $h_k$ restent la fonction d'observation et le GP. **Fait le 2026-10-22.**
- [x] **D11 Nomenclature** : alphabétique (opérateurs, grec, latin), fusions, descriptions concises, une lettre de base par ligne (D20). **Fait le 2026-10-21 et 22** ; l'encart « Conventions » de l'annexe II n'est pas ajouté.
- [x] **D12 ✔ Chapitre 5 et libellés de figures** : mêmes règles, dans un second temps (§9). **Fait le 2026-10-22** (lot C5, 421 symboles ; étiquettes des figures déjà conformes ; le bloc `comment` de l'ancien scénario reste à l'ancienne notation).
- [ ] **D13 (facultatif)** Différentielle droite $\mathrm dx$ (≈ 200 occurrences) pour séparer $d$ dimension et $d$ différentielle.
- [ ] **D14 (facultatif)** Aller jusqu'au bout sur $\lambda$, $r$, $n/m/p$ (§7) ($\varepsilon$ est réglé par D15).
- [x] **D15 Bruits et seuils** : $\eta_k$, $\varepsilon_k$, $\epsilon$ (voir D3). **Fait le 2026-10-21.**
- [x] **D16 Processus gaussien du ch.3** : $f\to h$ ($\sigma_f\to\sigma_h$), pour que la carte garde un seul nom et que la famille $(h,\varepsilon)$ soit celle du modèle d'observation. **Fait le 2026-10-21.**
- [x] **D17 Erreurs en $\zeta$ et évaluation $\operatorname{ev}_x$** : $e_{\mathbf R}\to\zeta_{\mathbf R}$, $e_{\hat z}\to\zeta_{\hat z}$ (et les erreurs des preuves du §4.5), pour rejoindre $\zeta_h$, $\zeta_\varepsilon$ ; $e_x\to\operatorname{ev}_x$ (RKHS) ; $e_i$ (base canonique) reste. **Fait le 2026-10-21.**
- [x] **D18 Échelles et coefficients, le 2026-10-22** : $p_\mathrm{RS}$ ; $\ell_\mathrm{UKF}$ ; largeurs de bande $\ell$ avec $\alpha_\mathrm{RPF}$ pour le coefficient dans $[0,1]$ ; coefficients de tendance et de carte $\boldsymbol\gamma$ ; $\mu$ gardé pour mesures, plongements et moyennes. **Fait** (lots R, S, T ; ch.5 plus tard).
- [x] **D19 Bruits d'entraînement** $\varepsilon_i\to\tilde\varepsilon_i$ (tilde = entraînement, comme $\tilde{\mathbf x}$, $\tilde{\mathbf z}$, $\tilde{\mathbf R}$). **Fait le 2026-10-22** (lot U).
- [x] **D20 Nomenclature à une lettre de base par ligne**, plus les fusions demandées ($\top,{-1},\dagger$ ; les deux $\sim$ ; les deux limites). **Fait le 2026-10-22.**
- [x] **D21 Noyaux et multiplicateur, le 2026-10-22** : $k\to\kappa$ (lots W et X, $\kappa_\mathrm{rad}$ pour le profil radial), $\boldsymbol\Xi$ pour le multiplicateur de Lagrange (lot Y) ; $\kappa_{d|r}$ est une « fonction de covariance », pas un noyau. **Fait.** *(Reste $\operatorname{Unif}$ : abandonné, voir D22.)*
- [x] **D22 Ensembles de départ en $\Omega$ et nomenclature complète, le 2026-10-22** : $\mathcal U\to\Omega$ (C7, lot Z) ; **règle de sa part : « si un symbole apparaît dans le manuscrit, il doit être dans la nomenclature » (« pour les retraits »)** : $\hat{\boldsymbol\gamma}_k$ remis (ligne séparée « Estimated map coefficients »), intervalles continus et tous les symboles manquants ajoutés, voir le journal. **Fait.**
- [x] **D23 Paire de la décomposition polaire** (`2_Monte_Carlo:292–293`) : $(\boldsymbol\Omega,\boldsymbol\Upsilon)$, $\boldsymbol\Upsilon$ à la place de $\boldsymbol\Sigma$ (en collision avec les $\boldsymbol\Sigma_{k|k-1}$ et le $\Sigma$ du GP) ; une ligne de nomenclature chacun. Chemin : $(\boldsymbol\Gamma,\boldsymbol\Delta)$ refusé (« clashe avec la fonction Gamma »), $(\mathbf Y_\mathrm{orth},\mathbf Y_\mathrm{SPSD})$, « oups tu peux garder Omega et Gamma », puis « au lieu de Gamma tu vas prendre Upsilon (on en est là) » : **état final**. **Fait le 2026-10-22.**
- [x] **D24 Identité** : $\operatorname{id}$ pour l'identité d'un espace de Hilbert (`3_KKF`, lot AA, 10 atomes), $\mathbf I$ pour les matrices, une seule ligne de nomenclature. **Fait le 2026-10-22.**
- [x] **D25 Nomenclature, après-midi du 22-10** : $\boldsymbol\gamma$ sur une ligne ; intervalles employés seulement ; pas de routines d'algorithme ; ligne des plongements de moyenne ; $\mathbf C$ « Trend matrix of a constant mean » ; $\tilde{\mathbf F}_k$ et $\tilde{\mathbf H}_k$ sur les lignes de $f$ et de $h$ ; symboles nommés oubliés (inventaire strict). **Fait.**
- [x] **D27 Dernière passe, le 2026-10-22, corrigée le soir** : chapitre 5 aligné sur la notation du reste (D12) ; nomenclature complète sur le texte visible seulement (aucune grandeur commentée, aucun symbole inventé, pas de ligne de marques), 189 entrées de symboles. **Fait.** À trancher : l'indice $\mathrm{MS}$ à deux sens, $a_d,b_d,L_d$ au ch.5, et les grandeurs $n$/$l$ du filtre de Kalman du RBPF (voir le journal).
- [x] **D26 Décisions du soir du 22-10** : $n_\mathrm{MC},n_\mathrm{steps}\to N_\mathrm{MC},N_\mathrm{steps}$ (comme les autres effectifs en $N$) ; $\mathcal B(x,r)$ **gardé** comme ensemble des indices des particules de la boule, « par abus de notation » (aucun symbole de boule ailleurs : seulement le mot « ball ») ; $\log\to\ln$ ; lettres locales $\mathbf M,\mathbf U,\mathbf D,\mathbf V$ de la décomposition polaire acceptées telles quelles. **Fait** (sauf la phrase éventuelle sur l'abus de notation, qu'il écrira lui-même).
- [x] **D28 Noyaux des espaces d'état et d'observation, la nuit du 22-10** : une seule graphie, celle de la théorie, $\kappa_X,\kappa_Z,\Phi_X,\Phi_Z,\mathcal H_{\kappa_X},\mathcal H_{\kappa_Z}$, jusque dans le KKF, l'AKKF et le ch.5 (lot AC, 100 atomes). **Fait.** À trancher : $\alpha_{\mathrm{RKHS},x/z}$ et $\mathbf G^{x/z}$ (minuscules gardées) ; quels $\mu$ il visait.

---

## Journal d'application

### 2026-10-21, étape 1 : matrices en gras majuscule (D1)

**Fait.** `\mathbf{P}` pour le latin, `\boldsymbol{\Sigma}` pour le grec, atome par atome (un analyseur qui conserve la position de chaque symbole, jamais une recherche de texte). 2 613 remplacements dans le texte visible (quinze fichiers, nomenclature comprise : 47, et la table `kriging_covariance`), 429 dans les blocs `comment` de la v2 (§4.5 et ses preuves, fin du krigeage séquentiel). Chapitre 5 et `resultats_slam.tex` non touchés, comme prévu.

**Exceptions voulues** (restées en italique) : noyau de régularisation $K$, $K_h$ et Bessel $K_\nu$ ; mesures $S^N$ et densité spectrale $S_\nu$ ; systèmes $(S)$, $(\bar S)$, $(S_1)$, $(S_2)$ du RBPF ; opérateur $T$ de RR et RK, et nombre d'itérations $T$ de AK ; identité de l'espace de Hilbert $I$ du KKF ; perte $E$ et pénalité $G$ du représentant, coût $J$ de RK ; espaces $E$, $F$ de $\mathcal L^p$ et vecteur aléatoire $F$ du GP ; paire d'estimateurs $A$, $B$ de la section 2.2 ; poids $W^{(i)}$ de l'UKF ; nombre de classes $C$, constantes $C_\mu$, $C_\nu$ et constante $C$ de la preuve de Matérn ; vecteurs aléatoires $U$, $V$ de l'annexe A ; variance scalaire $\Sigma$ de la définition d'un GP ; entrées $A_{i,j}$ d'une matrice. Deux matrices ne sont ni en italique ni en gras majuscule : $\mathcal G$ (C8, à traiter) et $\mathfrak D_G$ (fraktur, dans la nomenclature et les blocs de la v2, à revoir avec leur retour). **Les indices qui désignent une matrice sont en gras eux aussi** ($k_{\mathbf B}$, $\mathbf G_{\mathbf B}$, $\|\cdot\|_{\mathbf S_{d|r}}$, $\|x\|_{\mathbf M}$) ; dis-moi si c'est trop lourd.

**Contrôles.** Compilation : 223 pages comme avant, 0 erreur, 0 référence ou citation non définie, mêmes avertissements qu'avant ; trois boîtes trop larges apparues (les gras élargissent) et résorbées par un retour à la ligne `&\qquad+` dans `2_bayesian_estimation.tex` (formule de Joseph), `9_proofs_1.tex` (preuves B.2 et B.7). Le texte extrait du PDF contient exactement les mêmes caractères qu'avant, aux points de conduite de la nomenclature et d'un sommaire près (sept points de moins en tout, les symboles étant plus larges). La mise en page n'a bougé que sur six pages du PDF (43, 44, 109, 115, 138, 163) : quatre césures déplacées et, aux pages 43–44, le premier point de la liste des trois défauts de l'EKF passé en haut de la page suivante (effet de la ligne ajoutée à la formule de Joseph), si bien que « three reasons: » termine la page 43. Les blocs de la v2 ont été compilés une fois activés, avant et après : même résultat (243 pages, 0 erreur). Pages relues en image : nomenclature, algorithme de Kalman, formule de Joseph, annexes, §4.5 v2.

**Hors manuscrit.** Figure 2.1 régénérée avec $\mathbf P$ en gras (`figures_chapitre_2.m`, seule la figure 2.1 a été relancée) ; l'en-tête de `kriging_covariance.tex` est aussi passé en gras dans son générateur (`figures_chapitre_4.m`, ligne 527) pour qu'une relance ne l'annule pas.

**Reste de l'étape 1.** $\mathcal G\to\bar{\mathbf G}$ (C8) n'est pas fait : $\bar{\mathbf G}$ est une décision de nom, pas de typographie. La nomenclature est en gras mais n'a pas encore son encart « Conventions ». Le chapitre 5 et `resultats_slam.tex` n'ont pas reçu le gras : ils sont traités à part (§9).

**Retour arrière.** Je n'ai rien commité. Les fichiers d'avant l'étape sont copiés dans le dossier temporaire de la session (`D:\bhubert\AppData\Local\Temp\claude\D--bhubert-Desktop-Notes\737225cb-2766-4eb3-9fd4-901f1616bd01\scratchpad\backup_bold`, à copier ailleurs si tu veux les garder). Git ne l'isole pas : depuis le dernier commit (19-10) les fichiers portent aussi d'autres modifications.

### 2026-10-21, audit des marqueurs $r/d$, old/new, $n/l$

| Marqueur | Place | Forme | Occurrences | Où |
|---|---|---|---|---|
| $n$, $l$ (non linéaire, linéaire) | exposant | lettre italique | ≈ 330 | RBPF (§2) et ses preuves (B.2) |
| $r$, $d$ sur des données | exposant | lettre italique | ≈ 140 | $\tilde{\mathbf x}^r$, $\tilde{\mathbf z}^d$, $\tilde N^r$ (§4.3, §4.5) |
| $d\vert r$, $\vert r$ dans un conditionnement | indice | lettres italiques | ≈ 100 | $k_{d\vert r}$, $\tilde{\mathbf z}_{d\vert r}$, $\mathbf S_{d\vert r}$, $\mathbf R_{\mathrm{SK}\vert r}$ (§4.5) |
| old, new sur des données | indice | mots droits | ≈ 130 | $\tilde{\mathbf x}_{\mathrm{old}}$, $\tilde{\mathbf x}_{\mathrm{new}}$, $\tilde N_{\mathrm{new}}$ (§3.3) |
| $\mathrm{new}\vert\mathrm{old}$, $\vert\mathrm{old}$ dans un conditionnement | indice | mots droits | ≈ 60 | $k_{\vert\mathrm{old}}$, $\mathbf S_{\mathrm{new}\vert\mathrm{old}}$, $\hat z_{\mathrm{SK}\vert\mathrm{old}}$ (§3.3) |
| $r$ en indice nu | indice | lettre italique | 47 | $e_r$ (19), $\sigma_r$ (13), $\hat z_r$ (6), $\mathbf G_r$ (4), $g_r$ (5), §4.5 |
| $n(i)$, $l(i)$ avec indice de particule | exposant | juxtaposés | ≈ 110 | RBPF |
| $r,(c)$, $r,(i)$ avec indice de classe ou de particule | exposant | séparés par une virgule | 17 | ICAK, NNAK (§4.4) |

Lecture. Sur les données, $r/d$ sont au même niveau que $n/l$ (exposant d'une lettre) et $\mathrm{old}/\mathrm{new}$ suivent ta règle du 19-10 (mots droits en indice) ; sous une barre de conditionnement, $d\vert r$ et $\mathrm{new}\vert\mathrm{old}$ sont tous deux en indice, ce qui les met au même niveau. Deux écarts : le $r$ en indice nu de §4.5 (47 occurrences, des noms et des raccourcis plutôt que des repères de sous-ensemble, où $e^r$ ou $\mathbf G^r$ se liraient comme des puissances), et la façon de coller le repère à l'indice de particule ou de classe (`n(i)` dans le RBPF, `r,(c)` dans les méthodes AK). Le second se corrige à moindre coût en passant les 17 occurrences de AK à `r(c)` et `r(i)`.

### 2026-10-21 (soir), étape 2 : bruits, APF, mesures de Monte-Carlo, processus gaussien

**Fait** (D3, D4, D5 en partie, D15, D16), hors chapitre 5, par le même analyseur (atome par atome : seul l'intérieur des formules change, jamais le texte). Dans l'ordre : lot C (APF, vitesse, vecteurs muets) 58 remplacements ; lot A (mesures de Monte-Carlo) 170 ; lot B (bruits et seuils) 124 ; lot D (processus gaussien) 140 ; plus deux atomes à la main (le $n$ de la preuve du filtre de Kalman, `9_proofs_1:64`) : **494 au total**, blocs `comment` compris, puis 9 pour le lot E (variances) : **503**.

- **Lot C.** $\mu_k^{(i)}\to a_k^{(i)}$ (11 dans la section de l'APF, 1 dans la nomenclature) ; $\lambda,\tilde\lambda\to v,\tilde v$ (7 + 1) ; vitesse $v_k^{(x)},v_k^{(y)}\to\dot p_k^{(x)},\dot p_k^{(y)}$ (4) ; vecteurs muets libérant $v$ : $u,v\to u,y$ dans $\langle u|v\rangle_{\mathbf S_{d|r}}$ (§4.5, 3) et réalisation de $V$, $v,\tilde v\to y,\tilde y$ (annexe A, 31).
- **Lot A.** Monte-Carlo et filtrage particulaire (et la nomenclature) : $\mu\to\mu_p$, $\nu\to\mu_q$, $\mu_k\to\mu_{p,k}$, $\nu_k\to\mu_{q,k}$ ; constantes $C_\mu,C_\nu\to C_p,C_q$ ; indices de variance et d'espérance $\mathbb V_\mu\to\mathbb V_{\mu_p}$, $\mathbb V_\nu\to\mathbb V_{\mu_q}$, $\mathbb V_{\nu_k},\mathbb E_{\nu_k}\to\mathbb V_{\mu_{q,k}},\mathbb E_{\mu_{q,k}}$ (les $\mathbb E_p$ et $\mathbb E_q$ du texte sont inchangés). **Lot E (même soirée, à sa demande, 9 atomes) :** $\mathbb V_{\mu_p},\mathbb V_{\mu_q}\to\mathbb V_p,\mathbb V_q$, alignés sur $\mathbb E_p,\mathbb E_q$ ; l'entrée $\mathbb V_\mu$ de la nomenclature devient $\mathbb V_p$ (et le $\mu$ de sa description $p$ : le mot « measure » reste, alors que l'entrée de $\mathbb E$ dit « law ») ; restent $\mathbb V_{\mu_{q,k}}$ et $\mathbb E_{\mu_{q,k}}$ de la variance des poids de trajectoire (§2.4), faute d'un symbole pour la densité de la trajectoire. $C_p$ et $C_q$ validés. Les entrées génériques $\langle\mu,\varphi\rangle$ et $\mathbb V_\mu$ de la nomenclature restent : mesure quelconque.
- **Lot B.** $w\to\eta$ dans le modèle du §2.1, le RBPF, les preuves B.1 et B.2, le §4.1 et la nomenclature (78 atomes, dont $d_w\to d_\eta$) ; $n_k\to\varepsilon_k$ (35 atomes, dont $d_n$, $\mathbf n_{1:k}\to\boldsymbol\varepsilon_{1:k}$ et les échantillons de l'EnKF et de l'AKKF) ; seuils $\varepsilon\to\epsilon$ (11 : MAP, IEKF, correction SPD). Le $\varepsilon^{(j)}$ du RPF et le $\varepsilon$ du krigeage restent.
- **Lot D.** $f\to h$ pour le processus du ch.3, du §4.5 et de leurs preuves (124), $\sigma_f\to\sigma_h$ (16). Restent $f_k$ (évolution), les éléments $f\in\mathcal H$ des définitions de §3.5, le vecteur $F$ des valeurs du GP et l'indice $f$ de $e_f$ (« full »).

**Contrôles.** Pour chaque lot : texte hors formules identique, même nombre de lignes, comptes d'atomes égaux aux comptes attendus ; les lots C et D rejoués à blanc ne trouvent plus rien. Compilation : 223 pages, 0 erreur, 0 non défini, **une boîte trop large de 1,7 pt** (définition 2.3, ligne 100 de `2_Monte_Carlo.tex` : les indices de $\mu_p$ et $\mu_q$ élargissent la ligne). Pages relues en image : algorithme de l'APF, §2.4.2, RBPF, nomenclature, annexe A.1.

**À voir ensemble (aucun texte touché).** (1) `4_kriging_navigation:26` : « Here, $\varepsilon_k$ plays the role of the query noise $\varepsilon$ of … » est devenue redondante. (2) Nomenclature : $\varepsilon,\varepsilon_i$ (bruits de requête et d'entraînement) et $\varepsilon_k$ (bruit d'observation) sont deux entrées ; $\eta_k$, $a_k^{(i)}$, $v_k^{(i)}$, $\mu_p$, $\mu_q$ sont restés aux rangs de $w_k$, $\mu$, $\lambda$, $\nu$ ; $\epsilon$ n'a pas d'entrée. (3) $\mathbb E_{\mu_{q,k}}$ et $\mathbb V_{\mu_{q,k}}$ de la variance des poids de trajectoire (§2.4) : seuls indices de mesure restants ; on pourrait écrire $\mathbb V_{q(\mathbf x_{0:k}|\mathbf z_{1:k})}$, comme l'introduction du filtrage particulaire écrit $\mathbb V_{q(\cdot|\mathbf z_{1:k})}$. (4) La boîte trop large de la définition 2.3. (5) Les blocs `comment` de la v2 emploient encore $\eta$ comme constante scalaire (bornes de §4.5 et leurs preuves) : collision avec le bruit, à renommer au retour de la v2. (6) Chapitre 5 : $n_k$, $\varepsilon$, $\sigma_f$, $f$, $\mu$, $\nu$, le poids $w$ et les gras restent à passer dans le second temps. (7) **En attente** : $\sigma_h\to\sigma$ (il a proposé le $\sigma$ nu « s'il n'y a pas de conflit ») ; conflit avec le $\sigma$ de $\tilde{\mathbf R}=\sigma^2\mathbf I$ (§3.4 sur le rang réduit et sa preuve, 18 atomes) : proposition $\sigma_\varepsilon$ pour le bruit et $\sigma$ nu pour le noyau, ou rien. (8) **Questions posées, rien changé :** $F$ (vecteur des valeurs du GP, 5 atomes, ne fait plus écho à $h$ : le laisser ou le passer à $H$), indice $f$ de $e_f$ (« full » ou « field » ? s'il désigne l'ancien $f$, $e_h$), éléments génériques $f\in\mathcal H$ du RKHS (les laisser ou passer à $\psi$).

**Retour arrière.** Copie des fichiers d'avant l'étape 2 : `D:\bhubert\AppData\Local\Temp\claude\D--bhubert-Desktop-Notes\737225cb-2766-4eb3-9fd4-901f1616bd01\scratchpad\backup_notation` (dossier temporaire de la session, à copier ailleurs si tu veux la garder). Rien n'est commité.

### 2026-10-21 (nuit) et 2026-10-22 : lots F à Y et nomenclature

**Fait** (D6 à D8, D10, D11, D17 à D21), hors chapitre 5, atome par atome (texte hors formules identique, même nombre de lignes, comptes d'atomes égaux aux comptes attendus ; chaque lot passé à blanc avant d'être écrit, blocs `comment` compris) : 2 157 atomes depuis le début de l'étape 2 (503 après le lot E, 858 après le lot Q, 1 004 après le lot U, 1 070 après le lot V, 2 157 après le lot Y).

- **Lots F à I (93).** $l\to\ell$, longueur de corrélation (F, 26) ; vecteur des valeurs du GP $F\to H$ (G, 5) ; $e_r\to e_{\mathbf R}$ et $e_f\to e_{\hat z}$ (H, 55) ; variances et espérances sous la loi du chemin $q(\cdot|\mathbf z_{1:k})$ (I, 7). $\sigma_h$ gardé (le $\sigma$ nu aurait été en conflit avec celui de $\tilde{\mathbf R}=\sigma^2\mathbf I$), $C_p$ et $C_q$ validés.
- **Lots J à N (152).** $\sigma_r\to\sigma_{\mathbf R}$ (J, 28) ; éléments génériques du RKHS $f\to\psi$ (K, 8) ; argument explicite $\mathbf x_{0:k}$ des variances et espérances sous la loi du chemin (L, 12 ; l'équation de la variance totale est coupée en deux lignes) ; nombre de classes $C\to N_\mathrm{cl}$ et matrice de contrainte du §4.2 en gras (M, 26) ; coefficients multiplicateurs en $\alpha_\bullet$ (N, 78).
- **Lots O à Q (110).** Toutes les erreurs en $\zeta$ (O, 99) ; $e_x\to\operatorname{ev}_x$ (P, 1) ; $\alpha_\mathrm{UKF}$ et $\alpha_\mathrm{RS}$ (Q, 10), remplacés le lendemain par $\ell_\mathrm{UKF}$ et $p_\mathrm{RS}$.
- **22-10, lots R à U (146).** R : $p_\mathrm{RS}$ (4 atomes, §4.3) et $\ell_\mathrm{UKF}$ (5, UKF), nomenclature comprise. S : largeurs de bande $h\to\ell$ (65). T : coefficients des fonctions de base $\boldsymbol\gamma$ (64 : $\beta\to\boldsymbol\gamma$ 33, $\beta_0\to\boldsymbol\gamma_0$ 21, $\beta^\top$ 2, composantes $\gamma_i$ 4, $\alpha\to\boldsymbol\gamma$ de Karhunen–Loève 4). U : $\varepsilon_i\to\tilde\varepsilon_i$ (5, `3_kriging:14–17`) ; le bruit de requête $\varepsilon$ et celui de la navigation $\varepsilon_k$ restent.
- **Lot V (66).** Matrice de Gram par blocs $\mathcal G\to\bar{\mathbf G}$ (C8), sur sa remarque.
- **Lots W, X, Y (1 087).** W : profil radial $\kappa\to\kappa_\mathrm{rad}$ (3) ; X : noyaux $k\to\kappa$ (1 056 : `3_RKHS` 213, `9_proofs_1` 216, `9_proofs_2` 169, `4_subsampling` 159, `3_kriging` 122, `3_KKF` 70, `3_Gaussian_processes` 66, `3_reduced_rank` 31, `4_kriging_navigation` 9, `4_adaptive_kriging` 1) ; Y : multiplicateur de Lagrange $\xi\to\boldsymbol\Xi$ (28). Le temps $k$ ne bouge pas (2 467 atomes restants, dont les « $k$-means »).
- **Nomenclature (21 et 22-10).** Réordonnée (opérateurs, grec, latin, chaque groupe alphabétique), descriptions raccourcies ; fusions demandées ($\top,{-1},\dagger$ : « les 3 premiers symboles » ; les deux $\sim$ ; les deux limites) ; **une lettre de base par ligne** (« sépare les lignes qui contiennent des symboles qui ne sont pas dans l'ordre alphabétique, par exemple S sur la même ligne que k et z ») : $\mathcal G/\bar k/\bar{\mathbf z}$, $k_{d|r}/\tilde{\mathbf z}_{d|r}/\mathbf S_{d|r}$, $\tilde{\mathbf z}_{\mathrm{new}|\mathrm{old}}/\mathbf S_{\mathrm{new}|\mathrm{old}}$, $S_\nu/p_\nu$, $\psi,\boldsymbol\Psi/\boldsymbol\Lambda$, $\omega/b$, $a/v$, contrôlé par script. Les noyaux passent dans le bloc grec ($\kappa$, $\bar\kappa$, $\kappa_{d|r}$ « Residual covariance function ») avec $\nu$, $\xi_i$, $\boldsymbol\Xi$ ; $\omega_j$ devient « Frequencies of the map basis functions » et $b_j$ « Phases of the RFF » (il avait demandé « ω_j ce sont les fonctions de base de quoi ? » : les sinus de la méthode de l'espace de Hilbert, `5_SLAM:63`, et les RFF, `5_SLAM:204`). **Relecture de la fin de la nomenclature (sa remarque « bonne au moins jusqu'à cdf inclus ») :** $g_l\to g_i$ (le texte écrit $g_i$) ; $x_{k|k-1},x_{k|k}\to\hat x_{k|k-1},\hat x_{k|k}$ (le texte n'écrit que le chapeau, sauf les échantillons de l'EnKF) ; $\alpha_\mathrm{AK}$ « dilation factor » (le texte dit ainsi) ; entrées ajoutées : $\Delta_g$ (résolution de la grille, 25 occurrences), $\epsilon$ (seuils et corrections numériques, 11), $\nu$ (régularité de Matérn, 49 avec ses indices), $\chi^2_d(\cdot)$ (quantile, 17), $\kappa_\theta$ dans la liste des noyaux ; non ajoutés car locaux à une section : $\tau$ (pas du gradient, §4.2.2), $\Theta$ (ensemble des hyperparamètres admissibles), $\Omega$ (facteur polaire de §2), $\phi$ (élément générique de $\mathcal L^2$ en annexes). Lui-même fusionnait des lignes pendant ce temps ($G,G_\theta,\bar G$ ; $f_k$ et $\mathbf F_k$ ; $h_k$ et $\mathbf H_k$ ; $\lambda$ et $\boldsymbol\Lambda$ ; entrée de $\mathfrak D_G$ supprimée) : mes scripts retrouvent les lignes par leur symbole, donc ses retouches sont conservées.

**Retours sur la nomenclature (même jour : « je valide jusqu'à Monte Carlo measures inclus, à l'exception de… »).** Le texte de $\boldsymbol\gamma$ prend la forme de celui de $\boldsymbol\Lambda$ et $\boldsymbol\Psi$ (« Trend and map coefficients and their vector form », symboles $\gamma_i,\boldsymbol\gamma$ ; l'estimation $\hat{\boldsymbol\gamma}_k$ n'y figure plus, comme les chapeaux de $\lambda$) ; **toutes les dimensions sur une ligne** ($d_x,d_z,d_\varphi,d_\eta,d_\varepsilon,d_\theta,d_\mathrm{pos}$ : « Dimensions of vectors and subspaces » ; $d_\theta$, utilisé au §4.2.2, n'y était pas) ; **les nombres sur deux lignes** (« dans une moindre mesure ») : $N,\tilde N,N_\mathrm{cl},N_{\mathrm{cl},k},\tilde N^r,\tilde N_{\min}$ « Numbers of particles, training points and clusters » et $\hat N_\mathrm{eff},N_\mathrm{eff},N_\mathrm{th}$ « Effective sample size and resampling threshold » ; $r_\mathrm{UK},r_\mathrm{KL},r_\mathrm{RFF}$ : « Trend, KL and RFF ranks » (le mot *rank*, sur le modèle de « Kriging and subsampling errors »). **Question de lui, sans réponse appliquée :** « que penses-tu de $\mathfrak U$ pour les ensembles de départ des GP ? ou alors $\Omega$ ? » (il a déjà mis $\mathfrak U$ dans la nomenclature, le texte est encore en $\mathcal U$ : 215 atomes hors chapitre 5, 11 au chapitre 5, dont 13 pour $\mathcal U_X,\mathcal U_Z$) ; ma réponse : $\Omega$, comme la méthode de l'espace de Hilbert qu'il suit (Solin et Särkkä, « Let $\Omega\subset\mathbb R^d$ », 22 occurrences), $\mathcal U$ restant à la loi uniforme (4 atomes, C7 clos sans $\operatorname{Unif}$) ; à appliquer sur son accord.

**Retours suivants (même jour : « pour les retraits : si un symbole apparait dans le manuscrit, il doit être dans la nomenclature […] go pour Omega »).** Règle appliquée par un inventaire (`inv.pl`, atome par atome : chaque police, lettre et accent du texte visible comparé à la nomenclature) : **$\hat{\boldsymbol\gamma}_k$ revient** sur une ligne à part (« Estimated map coefficients », pour que la ligne de $\boldsymbol\gamma$ garde la forme de celle de $\boldsymbol\Lambda$ et $\boldsymbol\Psi$) ; **intervalles continus** $[\cdot,\cdot],(\cdot,\cdot],[\cdot,\cdot),(\cdot,\cdot)$ (« Continuous intervals », à côté de l'intervalle discret) ; **opérateurs** $\partial$, $\operatorname{Card}$, $\operatorname{dim}$, $\operatorname{erf}$, $\operatorname{ev}_x$, $\lim$, $\log$, $\max,\min$, $\operatorname{span}$ et les routines des algorithmes $\operatorname{Clustering}$, $\operatorname{MLE}$, $\operatorname{Resampling}$ ; **lettres** : $k$ (temps), $e_i$, $c_d$, $L$, $\phi$, $\tau$, $\Theta$ et $\hat\theta_k$, $\tilde\mu_p,\tilde\mu_q$, $\hat\mu_X$…, $\hat\varphi$, $\boldsymbol\varepsilon_{1:k}$, $\boldsymbol\Omega$, $\mathbf A$, $\mathbf B$, $\mathcal B(x,r)$, $\mathbf C$, $\mathcal C^{(c)}$, $\mathcal F$, $\tilde{\mathbf F}_k$, $\tilde{\mathbf G}$, $\tilde{\mathbf H}_k$, $\mathcal H,\mathcal H_\kappa$, $I$, $\tilde K_\ell$, $\hat{\mathbf K}_k$, $\mathcal Q_k$, $\tilde p$, $\tilde q$, $p_\mathrm{pos}$, les composantes de position et de vitesse, $\hat{\mathbf P}_k$, $W^{(i)}$ et $\mathbf W$, $\tilde x_i$, $\breve x^{(c)}$, $\overline x_0$, $X$ et $Z$ avec leurs ensembles d'échantillons et leurs moyennes, $\mathcal X^{(i)}$, $\tilde v$, $\tilde z_i$ : de 128 à 169 entrées, quatre pages. **Ligne des dimensions :** « Dimensions of vectors and subspaces » faisait fourre-tout et $d_\mathrm{pos}$ est la dimension du vecteur de position ; la phrase naturelle (« … of the state, observation, mapping, noise, hyperparameter and position vectors ») ne tient pas sur une ligne (une description dépasse une cinquantaine de caractères avec sept symboles à gauche), j'ai donc nommé les vecteurs par leurs symboles : « Dimensions of $x$, $z$, $\varphi(x)$, $\eta$, $\varepsilon$, $\theta$ and $\boldsymbol\Pi_\mathrm{pos}x$ » (dans le style de « Value sets of $X$ and $Z$ »). Les lignes $N_\bullet$ et $r_\bullet$ sont validées. **Non ajoutés, car locaux à une preuve, un lemme ou un algorithme (à sa décision) :** annexe A ($U,V,u,\tilde u,\tilde y,\bar{\mathbf A},\bar{\mathbf B},\bar{\mathbf C},\bar{\mathbf D},\bar{\mathbf M}$), preuve de Matérn ($\tilde r$), meilleure approximation de rang réduit ($\tilde\kappa$), théorème du représentant ($J$, $G$, $\mathcal V_{\tilde{\mathbf x}}$), l'opérateur générique $T$, les systèmes $(S),(\bar S)$ du RBPF, la paire d'estimateurs $A,B$ du chapitre 2, les indices et entiers génériques ($u$, $n$, $j$, $l$) ; chapitre 5 : ses formes à l'ancienne notation s'aligneront sur les entrées existantes à sa passe, plus $t$, $\bar x_k$ et $L_d$ qui n'ont pas d'entrée.

**Hors manuscrit.** La figure 3.1 (« Noyaux de Matérn ») portait `k_ν(r)` en ordonnée : `figures_chapitre_3.m:61` passe à `\kappa_\nu(r)` et seule cette figure a été régénérée (MATLAB R2023b en lot, copie temporaire du script dans le dossier de la session, `return` après le premier `exporter`) ; les figures dessinées à la main n'ont pas de symbole de noyau (le diagramme de régularisation du RPF n'a que des mots).

**Contrôles.** Compilation propre après chaque série : 223 pages, 0 erreur, 0 non défini, 0 boîte trop large (un seul dépassement de 0,22 pt, apparu avec $\kappa$ dans une ligne déjà serrée par un `\hspace{-7pt}` de `9_proofs_2:211`, résorbé en passant à `-8pt`). Pages relues en image : nomenclature (p. xx–xxii), présentation du RPF (p. 36–37), bruits d'entraînement (p. 53), tendance de l'UK (p. 58), forme par blocs de l'UK (p. 60), définition du noyau radial (p. 49). **Incident de build, trois fois dans la journée :** des fichiers auxiliaires tronqués (`File ended while scanning use of \contentsline`, caractères nuls dans `main.aux`, `\ab` inachevé), des centaines de « undefined » ensuite, sans rapport avec les modifications ; cause probable : une compilation d'un autre outil (LaTeX Workshop) lancée en même temps que la mienne ; remède : attendre la fin de tout processus LaTeX, supprimer `*.aux`, `*.lof`, `*.lot`, `*.toc`, `*.loa`, `*.out`, `*.fdb_latexmk`, `*.run.xml` de `build/`, relancer ; le PDF de la racine n'est copié qu'après un build propre.

**À voir ensemble (aucun texte touché).** (a) Les constantes $\alpha$, $\kappa$, $\eta$, $\gamma$ des blocs `comment` de la v2 (`4_subsampling` ≈ 506–662, `9_proofs_2` ≈ 435–628) entrent en collision avec les symboles actuels : $\kappa$ est maintenant le noyau (et $k_{\max}\to\kappa_{\max}$ y est la borne du noyau, à côté de la constante $\kappa$), $\eta$ le bruit, $\alpha$ la famille des coefficients, $\gamma$ les coefficients de base ; à renommer au retour des blocs. (b) $\alpha_\mathrm{AK}$ est lui aussi un coefficient de dilatation (le rayon est $r_{\mathrm{cov},k}=\alpha_\mathrm{AK}\sqrt{\chi^2_{d_\mathrm{pos}}(99\%)}\sqrt{\lambda_{\max}}$, « a dilation factor ») : famille $\alpha_\bullet$ (décision du 21-10) ou $\ell_\bullet$ ? (c) $\kappa_X,\kappa_Z$ (`3_RKHS` et début de `3_KKF`) et $\kappa_x,\kappa_z$ (fin de `3_KKF`) nomment les mêmes noyaux avec deux graphies, comme $\Phi_X,\Phi_Z$ et $\Phi_x,\Phi_z$ : à unifier ? (d) La définition d'un GP écrit $\mathcal N(\mu,\Sigma)$ avec un $\Sigma$ **scalaire** (`3_Gaussian_processes:35–37`) et l'UKF $(\mu,\mathbf P)$, alors que la nomenclature annonce $\mathcal N(\mu,\boldsymbol\Sigma)$ : passer le scalaire en $\sigma^2$ (L1) ? (e) Chapitre 5 : $\alpha,\hat\alpha_k,\mathbf P_k^\alpha\to\boldsymbol\gamma,\hat{\boldsymbol\gamma}_k,\mathbf P_k^\gamma$, $\alpha_\mathrm{reg}h\to\alpha_\mathrm{RPF}\ell$, $\beta\to\beta_\mathrm{RLS}$, $\lambda_\mathrm{RLS}\to\alpha_\mathrm{RLS}$, $\lambda_{x,z}\to\alpha_{\mathrm{RKHS},x/z}$, $k_\nu\to\kappa_\nu$, $k_x,k_z\to\kappa_x,\kappa_z$, $w$, $n_k$, $f$, $\sigma_f$, gras (la nomenclature est déjà à jour pour ceux qui y figurent).

**Retour arrière.** Sauvegardes dans le dossier temporaire de la session : `backup_notation`, `…2`, `…3`, `…4` (avant l'étape 2 et après les lots E, I, N), `backup_notation5` (avant le lot R), `…6` (avant T et U), `…7` (avant V), `…8` (avant W, X et Y) ; états de la nomenclature `0_administratif_avant_*.tex`. Rien n'est commité.

### 2026-10-22 (après-midi) : neuf retours sur la nomenclature, identité $\operatorname{id}$, paire de la décomposition polaire

**Fait** (D23 à D25). Nomenclature par script (`nomen_v6` à `nomen_v9` : lignes retrouvées par leur symbole, donc ses retouches conservées ; il a fusionné entre-temps la ligne de l'intervalle discret et celle des intervalles continus en « Discrete, closed and left-open intervals »), texte par le même analyseur (lot AA) ou à la main. Ses neuf points, dans l'ordre :

1. **$\boldsymbol\gamma$ sur une ligne** : $\gamma_i,\boldsymbol\gamma,\hat{\boldsymbol\gamma}_k$, « Trend and map coefficients, their vector form and estimate » (la ligne « Estimated map coefficients » disparaît).
2. **Intervalles** : seules les formes employées. Le texte écrit $[a,b]$ (une douzaine : $[0,1]$, $[0,1/N]$, $[0,2\pi]$, les bornes des scénarios) et $(0,1]$ (3 : `4_subsampling:475` et `:498`, `5_SLAM:230`) ; ni $[a,b)$ (pas même $[0,+\infty)$) ni $(a,b)$.
3. **Routines d'algorithme** : $\operatorname{Clustering}$, $\operatorname{MLE}$, $\operatorname{Resampling}$ retirées de la nomenclature (règle : les noms de routines des algorithmes n'y figurent pas).
4. **Décomposition polaire** (`2_Monte_Carlo:292–293`), voir plus bas.
5. **Plongements de moyenne** : « Kernel mean embeddings and their estimates » (une ligne).
6. **$\mathcal B(x,r)$ n'est pas une boule** : le texte le définit comme l'ensemble des indices des particules à distance au plus $r$ de $x$ (`4_adaptive_kriging:349`, repris par le mean shift et DBSCAN), donc un ensemble d'indices ; la lettre $\mathcal B$ fait lire une boule de $\mathbb R^d$ (R3). La description de la nomenclature, fidèle au texte, est inchangée ; j'avais proposé $\mathcal I(x,r)$ (le $\mathcal I$ calligraphique n'étant utilisé nulle part ; $\mathcal C^{(c)}$ est lui aussi un ensemble d'indices, ceux de la classe $c$) ; **il a tranché le soir : $\mathcal B$ est gardé**, voir la section suivante.
7. **$\mathbf C$** : « Trend matrix of a constant mean » (le texte dit : « the trend matrices reduce to $g(x)=\mathbf I_{d_z}$ and $g(\tilde{\mathbf x})=\mathbf C$ », `3_kriging:309`, et la preuve de la proposition 3.9 écrit $\mathbf C\to g(\tilde{\mathbf x})$) ; l'ancien libellé « Replication matrix » venait de `3_kriging:248` (« its replication for each of the training points », $\mathbf C\hat m_\mathrm{OK}$). Si « nom » visait la lettre et non le libellé, à reprendre avec lui.
8. **Identité** : l'italique $I$ était l'identité de l'espace de Hilbert $\mathcal H_{\kappa_X}$ (ou $\mathcal H_{\kappa_z}$) dans les inverses régularisés de `3_KKF` ($(\hat\mu_{XX}+\alpha_\mathrm{RKHS}I)^{-1}$, « where $I$ denotes the identity of … », 10 atomes, 8 lignes) ; il passe à $\operatorname{id}$ (lot AA, comme $\operatorname{ev}_x$, $\operatorname{rank}$, $\operatorname{tr}$), l'identité matricielle restant $\mathbf I$, $\mathbf I_n$ : j'avais mis une seule ligne « Identity operator and matrix, subscript: of size $n$ » ; **il l'a séparée lui-même** ($\mathbf I,\mathbf I_n$ « Identity matrix, subscript: of size $n$ » et $\operatorname{id}$ « Identity operator »), état final. Les $I$ à indice du chapitre 5 sont des matrices, à passer en gras à sa passe.
9. **« Natural logarithm »** : c'est le logarithme en base $e$ (le logarithme décimal est $\log_{10}$) ; les sept $\log$ du manuscrit (six lignes) étaient tous naturels (log-vraisemblance gaussienne $-\tfrac12\log\det\mathbf G_\theta-\tfrac{\tilde Nd_z}{2}\log2\pi$, sa dérivée et la dérivée en $\log\theta_i$, $\mathcal O(N^2\log N)$ et $\mathcal O(N\log N)$ qui ne dépendent pas de la base, $\sqrt{2\log(2/\delta)}$) ; il n'y avait ni $\ln$ ni $\log_{10}$. J'avais proposé `\ln` pour lever l'ambiguïté avec la convention française ; **il a dit oui le soir** (« passe log en ln »), voir la section suivante.

**Décomposition polaire.** Le $\boldsymbol\Sigma$ du facteur symétrique entrait en collision avec les $\boldsymbol\Sigma_{k|k-1}$ de l'AKKF et le $\Sigma$ scalaire du GP, et seul le facteur orthogonal figurait dans la nomenclature (mon inventaire compare les lettres sans leur indice : $\boldsymbol\Sigma$ y était donc déjà présent). Chemin parcouru, par honnêteté : $(\boldsymbol\Gamma,\boldsymbol\Delta)$ d'abord, refusé (« clashe avec la fonction Gamma ») ; recensement de toutes les majuscules par police (`libres.pl`, chapitre 5 compris) : **seules $Y$ (latin) et $\Upsilon$ (grec) ne sont utilisées nulle part**, et se ressemblent trop pour faire une paire ; $(\mathbf Y_\mathrm{orth},\mathbf Y_\mathrm{SPSD})$ ensuite, avec un retour à la ligne disgracieux dans la puce ; puis sa réponse (« oups tu peux garder Omega et Gamma pour la DP »), aussitôt corrigée (« au lieu de Gamma tu vas prendre Upsilon (on en est là) »). **État final : $\mathbf M=\boldsymbol\Omega\boldsymbol\Upsilon$, $\boldsymbol\Omega\triangleq\mathbf U\mathbf V^\top$, $\boldsymbol\Upsilon\triangleq\mathbf V\mathbf D\mathbf V^\top$, $\mathbf M_+\triangleq(\mathbf M+\boldsymbol\Upsilon)/2$**, une ligne de nomenclature chacun (« Orthogonal factor » après $\Omega,\Omega_X,\Omega_Z$, « SPSD factor » entre $\tau$ et $\varphi$ ; une lettre de base par ligne). $\Upsilon$ est la seule majuscule grecque libre (aucune collision) ; $\boldsymbol\Omega$ garde son homonyme maigre, l'ensemble $\Omega$, qu'il a accepté en le retenant. Le texte alentour (« the SPSD term of the polar decomposition ») est inchangé. **Reste signalé, non appliqué :** $\mathbf M$, $\mathbf U$, $\mathbf D$, $\mathbf V$ de la SVD du même bloc sont des noms locaux qui entrent en collision avec les $\mathbf M$ (covariance des coefficients de tendance estimés) et $\mathbf U$ (leurs poids) de l'UK, $\mathbf D$ et $\mathbf V$ (transition plongée) du KKF, au ch.3 ; aucune majuscule n'étant libre, le plus court est de les garder locaux (comme ceux de l'annexe A) ou de les noter $\hat{\mathbf P}_k^{\mathrm{sym}}$ et autres dérivés. **Accepté le soir : « pas de problème pour les matrices locales, elles sont loin les unes des autres anyway »** ; rien à changer.

**Constante de la preuve de Matérn.** Erreur de mon lot « matrices en gras » (21-10) : la constante scalaire $C>0$ de la transformée de Fourier de $\kappa_\nu$ (`9_proofs_1:199, 202, 212`, « a normalising constant which only depends on $d_x$ and the parameters of $\kappa_\nu$ ») avait été passée en $\mathbf C$ (la règle d'exclusion de ce passage n'a pas opéré). **Remise en $C$ maigre (3 atomes)**, remarquée par lui à la relecture. Contrôle étendu : les 147 groupes (fichier, lettre) en gras du texte visible sont listés (`boldscan.pl`) et ceux qui pouvaient cacher un scalaire (une soixantaine, dont tous les gros groupes) relus en contexte : plus aucun scalaire, ensemble ni constante en gras.

**Lignes de $f$ et de $h$.** $\tilde{\mathbf F}_k$ et $\tilde{\mathbf H}_k$ passent sur la ligne de $f$ et sur celle de $h$ : « Evolution functions and matrices, linearised » et « Observation function and matrix, linearised » (la première formulation, « …, matrices and linearised matrix », débordait d'un mot sur la ligne de $f$, déjà à neuf symboles, et se retrouvait collée à la marge). **Il les a ensuite ramenées lui-même à « Evolution functions and matrices » et « Observation function and matrix »** (les versions linéarisées y sont des matrices) : état final.

**Inventaire strict (nouveau).** `inv.pl` compare la police, la lettre et les accents, pas l'indice : il laissait passer $\Delta_t$ (sous $\Delta_g$), $n_\mathrm{MC}$ (sous $n$), $\sigma_h$, etc. `inv2.pl` compare aussi l'indice nommé. Entrées ajoutées : $\Delta_t$ (avec $\Delta_g$ : « Grid resolution and sampling period »), $n_\mathrm{MC},n_\mathrm{steps}$ (« Numbers of Monte Carlo simulations and time steps »), $r_\mathrm{DB}$ (avec $r_{\mathrm{cov},k},r_{\min,k}$ : « Covariance, minimum-count and DBSCAN radii »), $w_{\min}$ (« DBSCAN weight threshold »), $\sigma_h$ (avec $\sigma_{\mathbf R}$ : « Kernel and subsampling standard deviations »), $E_\mathrm{in},E_\mathrm{gap},E_\mathrm{out}$ (« Subsets of the prediction domain »). Restent locaux, non ajoutés : $\theta_\mathrm{init}$ (valeur initiale de l'ascension), $\hat m_\mathrm{ML}$, les $\psi_\bullet$ du théorème du représentant, les covariances croisées $\hat{\mathbf P}_{\mathbf{XZ}}$ de l'EnKF, $u_j$ du rééchantillonnage, et ceux de la liste du 22-10 (annexe A, Matérn, $T$, $A,B$…).

**Contrôles.** Compilation propre après chaque série : 223 pages, 0 erreur, 0 non défini, 0 boîte trop large. Pages relues en image : nomenclature (p. xx à xxiii, sans ligne qui déborde ; la page xxiv est la page blanche avant le chapitre suivant), bloc de la décomposition polaire (ch.2, p. 27 : la puce tient sur une ligne), identité (ch.3, p. 76), preuve de Matérn (annexe B.4, p. 171). `inv.pl` : plus que les locaux de la liste ; `nomen_check.pl` (une lettre de base par ligne, ordre) : rien à signaler en dehors des symboles composés ($\mathcal{GP}$, $S^N(\mu_p)$, $\mathcal U([a,b])$…). **Nomenclature : 166 entrées de symboles** (après ses retouches de 16 h 17 à 16 h 33 : « Subsets of the example kriging domain », descriptions de $f$ et $h$, $\operatorname{id}$ séparé ; mes scripts retrouvent les lignes par leur symbole et ne les ont pas touchées).

**Retour arrière.** Dossier temporaire de la session : `backup_notation10` (avant cette série : `0_administratif`, `2_Monte_Carlo`, `3_KKF`, `9_proofs_1`), états `0_administratif_avant_v7.tex`, `…_v8.tex`, `…_v9.tex`. Rien n'est commité.

### 2026-10-22 (soir) : $N_\mathrm{MC}$ et $N_\mathrm{steps}$, $\ln$, $\mathcal B$ gardé

**Fait** (D26), sur ses quatre décisions :

- **$n_\mathrm{MC}$ et $n_\mathrm{steps}\to N_\mathrm{MC}$ et $N_\mathrm{steps}$** (« pour garder la cohérence avec les autres $N$ ») : lot AB, 16 atomes (10 dans `4_kriging_navigation`, 6 dans `5_SLAM`), mêmes contrôles que les autres lots. J'ai inclus le chapitre 5, pourtant en attente, parce que l'entrée de la nomenclature est globale : sans cela son texte aurait contredit la nomenclature. L'entrée passe dans le bloc $N$ (« Numbers of Monte Carlo simulations and time steps », après $\hat N_\mathrm{eff},N_\mathrm{eff},N_\mathrm{th}$). **Hors texte :** l'étiquette d'ordonnée de la figure « Erreur finale par maille » portait `$n_\mathrm{MC}$` (`figures_campagnes.m:591`) ; seule cette figure a été régénérée (script temporaire : préambule et campagne, section de la figure, fonctions locales ; sortie dans le dossier de la session, comparée à l'ancien PDF au pixel près : 1 276 pixels différents, tous dans la bande de l'étiquette ; PDF remplacé dans `Images/`). La ligne de commentaire que le script écrit en tête de `ak_vs_ok.tex` et `cak_vs_ak.tex` suit aussi (`N_MC`, sans effet sur le rendu).
- **$\ln$** (« passe log en ln ») : sept occurrences sur six lignes, `4_kriging_navigation` 209 (deux), 255, 283 ; `4_adaptive_kriging` 313, 380 ; `4_subsampling` 498 ; l'entrée de la nomenclature devient « $\ln$, Natural logarithm » (`log_to_ln.ps1` ne remplace que la macro `\log`).
- **$\mathcal B(x,r)$ gardé** : « à moins qu'on n'utilise la boule explicitement ailleurs, je pense qu'on peut identifier une particule à son indice et garder $\mathcal B$ pour désigner les indices des particules qui sont dans la boule […] par abus de notation ». Vérifié : le manuscrit n'a aucun symbole de boule, seulement le mot (« within a ball of radius », `4_adaptive_kriging:101` et `:107`, pour la fenêtre de l'AK, et « the unit ball » pour $c_d$). $\mathcal B$ n'est introduit qu'une fois (`4_adaptive_kriging:347–351`) et repris dans l'itération du mean shift (`:355`) et dans DBSCAN (`:375`). **Question à lui, où préciser l'abus ?** Ma réponse : dans la phrase qui suit la définition (`:351`, « contains the indices of the particles within a distance $r$ of $x$ »), une subordonnée suffit ; il l'écrit lui-même.
- **Matrices locales** $\mathbf M,\mathbf U,\mathbf D,\mathbf V$ du bloc de la décomposition polaire : acceptées (« elles sont loin les unes des autres anyway »).

**Contrôles.** Compilation propre : 223 pages, 0 erreur, 0 non défini, 0 boîte trop large. Pages relues en image : nomenclature p. xxii (lignes $N$, $\ln$), §4.2 (log-vraisemblance, p. 95), tableau 4.5 et NEES (p. 100), figure régénérée (« then over the $N_\mathrm{MC}$ runs »).

**Retour arrière.** `backup_notation11` (avant ces changements : `0_administratif`, `4_kriging_navigation`, `4_adaptive_kriging`, `4_subsampling`, `5_SLAM`), `0_administratif_avant_v11.tex`, `efm_ancien_sauvegarde.pdf` (ancienne figure), dossier temporaire de la session. Rien n'est commité.

### 2026-10-22 (fin d'après-midi) : dernière passe, chapitre 5 et nomenclature complète

Sa demande : « fais une dernière passe sur la nomenclature, vérifie que toutes les notations y soient et que le manuscrit en entier (chapitre 5 inclus) s'y conforme ».

**Chapitre 5 (lot C5, 421 symboles, texte visible).** Mêmes contrôles que les autres lots : texte hors formules identique, mêmes lignes, comptes d'atomes égaux aux comptes attendus. Seules les correspondances déjà décidées sont appliquées, aucune lettre nouvelle :
- matrices en gras : $\mathbf P$, $\mathbf Q$, $\mathbf R$, $\mathbf H$, $\mathbf F$, $\mathbf I$, $\mathbf K$, $\mathbf S$, $\mathbf A$, $\mathbf D$, $\mathbf G$, $\mathbf T$, $\boldsymbol\Sigma$, $\boldsymbol\Lambda$, $\boldsymbol\Psi$ et $\boldsymbol\Pi_\mathrm{pos}$ ;
- bruit d'observation $n_k\to\varepsilon_k$ (et $\varepsilon_k^{(i)}$) ;
- noyaux $k_\nu$, $k_{3/2}$, $k_x$, $k_z\to\kappa$ ($\mathcal H_{\kappa_\nu}$) ; $\sigma_f\to\sigma_h$ ; longueur de corrélation $l\to\ell$ (le $l$ des indices de temps du RLS reste) ;
- coefficients de la carte $\alpha\to\boldsymbol\gamma$ ($\hat{\boldsymbol\gamma}_k$, $\mathbf P_k^\gamma$, $\mathbf K_k^\gamma$) ;
- $\lambda_x,\lambda_z\to\alpha_{\mathrm{RKHS},x},\alpha_{\mathrm{RKHS},z}$, $\lambda_\mathrm{RLS}\to\alpha_\mathrm{RLS}$, $\beta\to\beta_\mathrm{RLS}$ ;
- domaine $\mathcal U\to\Omega$ (la loi uniforme $\mathcal U([0,2\pi])$ reste) ; $\alpha_\mathrm{reg}h\to\alpha_\mathrm{RPF}\ell$.

L'algorithme 5.2 reprend désormais l'algorithme 3.2 (AKKF) symbole pour symbole. **Le bloc `comment` de l'ancien scénario (l. 466–509) n'est pas touché** : il emploie $\gamma$ pour son niveau de dégradation, ce qui heurterait $\boldsymbol\gamma$ ; à reprendre s'il revient.

**Nomenclature (166 → 184 entrées de symboles).** Inventaires stricts sur tout le texte visible, chapitre 5, introduction, conclusion et tableaux inclus : `inv2.pl` pour les indices nommés, `inv_sup.pl` pour les exposants nommés, `inv_rev.pl` pour le contrôle inverse (entrées sans occurrence). Ajouts :
- **variantes de familles déjà présentes** : $\kappa_{d|r}^0$, $\mathbf S_{d|r}^0$, $\tilde{\mathbf z}_{d|r}^0$ (« superscript 0: for SK ») ; $\boldsymbol\Sigma_{k|k}^+$ ; $\sigma_{\mathbf R,0}$, $\sigma_{\mathbf R}^{\max}$ ; $\hat\theta_\mathrm{ML}$, $\theta_\mathrm{init}$ ; $\hat m_\mathrm{ML}$ ; $\bar{\mathbf G}_h$, $\mathbf G^z$, $\mathbf G_k^x$, $\mathbf G_k^z$ ; $\Phi_X,\Phi_Z,\Phi_x,\Phi_z$ ; $w_k^{(i,c)}$ ; $x_k^{+(i)}$ ; $\hat\mu_{k|k}$ ; $\hat z_{\mathrm{SK},\mathbf B}$, $\hat z_{\mathrm{SK}|\mathrm{old}}$, $\hat z_{\mathrm{SK}|r}$, et les covariances correspondantes ;
- **lignes nouvelles** : $\kappa_X,\kappa_Z,\kappa_x,\kappa_z$ ; $\kappa_\mathrm{rad}$ ; plongements joints ($\mu_{XX},\mu_{XZ},\mu_{ZX}$ et estimations) ; plongements conditionnels ($\mu_{Z|x}$ et estimations) ; opérateurs de transition et d'observation du KKF ; $\boldsymbol\Sigma,\boldsymbol\Sigma_i$ des noyaux séparable et multi-échelle ; $\psi$ (estimateur ou élément générique d'un RKHS) ; $\psi_\mathrm{RKHS}$ ; $\varepsilon^{(j)}$ du RPF ; $i^{(j)}$ ; $\mathbf K_k^\gamma$ ; $\hat{\mathbf P}_{\mathrm{pos},k|k-1}$ ; $\hat{\mathbf P}_{\mathbf{XZ}}$ ; $\mathbf R_{\hat m},\mathbf R_{\hat m_\mathrm{OK}},\mathbf R_{\hat m_\mathrm{UK}}$ (séparés des covariances de krigeage) ; $\tilde N^d,\tilde N_\mathrm{new},\tilde N_{\cdot,k}^r$ ; $\tilde{\mathbf x}^r,\tilde{\mathbf x}^d,\tilde{\mathbf x}^r_{\cdot,k}$ (« subscript: method ») et de même pour $\tilde{\mathbf z}$ ;
- **une ligne de marques** $\cdot^n,\cdot^l$ (« Non-linear and linear parts of an RBPF quantity »), qui couvre $\mathbf P^l$, $\mathbf K^l$, $\mathbf S^l$, $\hat x^l$, $y^l$ du RBPF sans les énumérer ;
- **descriptions** : $\psi_i,\boldsymbol\Psi$ devient « Map basis functions and their matrix form » (fonctions propres de Mercer au ch.3, base de Laplace et RFF au ch.5) ; l'entrée $\tilde{\mathbf x}_k^r$, qui n'apparaissait nulle part, est remplacée.

Deux lignes trop longues ont été scindées (plongements marginaux et joints ; covariances de krigeage et des moyennes estimées) : aucune entrée ne passe plus à la ligne, la nomenclature va de la p. xix à la p. xxiv. **Hors liste, car locaux** : symboles des preuves et de l'annexe A, abréviations de preuve ($\mathbf G_{\mathbf B}$, $\mathbf M_{\mathbf B}$, $\hat z_\mathrm{old}$, $\mathbf G_r$…), indices et composantes génériques, et au chapitre 5 les bornes et longueurs du domaine $a_d,b_d,L_d$, la translation $t$, la moyenne empirique $\bar x_k$ et le critère $J_k$.

**Figures et tableaux.** Les étiquettes des scripts MATLAB sont conformes ; les figures dessinées à la main n'ont pas de symbole. `fenetrage.tex` (tableau généré, non inclus dans le manuscrit) et son générateur portaient encore $\alpha_{\mathrm{cov}}$ : passés à $\alpha_{\mathrm{AK}}$.

**Contrôles.** Compilation propre après chaque étape : 223 pages, 0 erreur, 0 non défini, 0 boîte trop large. Pages relues en image : algorithmes 5.1 et 5.2, équations 5.10 à 5.13, tableau 5.1, remarque 5.3, nomenclature entière. Il n'y a plus d'ancienne graphie dans le texte visible ; il en reste dans le bloc `comment` du chapitre 5.

**À trancher (rien d'appliqué).**
1. L'indice $\mathrm{MS}$ a deux sens : multi-échelle dans $\kappa_\mathrm{MS}$ (§3.1) et mean shift dans $\ell_\mathrm{MS}$ et $K_{\ell_\mathrm{MS}}$ (§4.4).
2. Au chapitre 5, $a_d,b_d,L_d$ partagent leurs lettres avec $a_k^{(i)}$, les phases $b_j$ et la log-vraisemblance $L$ (R2), bien qu'ils soient locaux.
3. Toujours ouverts : les deux graphies $\kappa_X,\kappa_Z$ et $\kappa_x,\kappa_z$ (listées toutes deux), le $\Sigma$ scalaire de la définition d'un GP, la famille de $\alpha_\mathrm{AK}$ ($\alpha_\bullet$ ou $\ell_\bullet$).

**Retour arrière.** Dossier temporaire de cette session (`…\35cf4a95-…\scratchpad`) : `backup_c5` (avant le lot C5 : `5_SLAM`, `0_administratif`), `0_administratif_avant_v12.tex`, `…_v13.tex`, `…_v14.tex`. Rien n'est commité.

### 2026-10-22 (soir) : correction de la dernière passe

Ses retours : « non pour la ligne .^n et .^l », « il ne faut pas faire apparaître des grandeurs qui sont commentées », « globalement tu as fait n'importe quoi, revérifie ». **La section précédente du journal n'est plus exacte sur ces points** :
- la ligne $\cdot^n,\cdot^l$ est retirée (par lui) ;
- les exposants 0 et les symboles génériques $\tilde{\mathbf x}^r_{\cdot,k}$ sont abandonnés ;
- le décompte passe de 184 à 189 entrées.

**Ce qui était faux.**
1. **Grandeurs commentées.** J'avais repéré une partie des ajouts par une recherche de texte, qui ne distingue pas les blocs `comment`. $\bar{\mathbf G}_h$, $\kappa_{d|r}^0$, $\mathbf S_{d|r}^0$, $\tilde{\mathbf z}_{d|r}^0$, $\sigma_{\mathbf R,0}$ et $\sigma_{\mathbf R}^{\max}$ n'existent que dans les blocs `comment` de §4.5 (`4_subsampling` l. 217–380, 416–468, 516–672, et `9_proofs_2` l. 435–628). Il en avait retiré trois ; j'ai retiré les trois autres et rétabli « Residual covariance » et « Residual observation vector ».
2. **Symboles inventés.** $\tilde{\mathbf x}^r_{\cdot,k}$, $\tilde{\mathbf z}^r_{\cdot,k}$ et $\tilde N^r_{\cdot,k}$ (« subscript: method ») n'apparaissent pas dans le texte. Ils sont remplacés par les symboles écrits : $\tilde{\mathbf x}^r_{\mathrm{RS},k}$, $\tilde{\mathbf x}^r_{\mathrm{AK},k}$, $\tilde{\mathbf x}^r_{\mathrm{NNAK},k}$, $\tilde{\mathbf x}^r_{\mathrm{MCAK},k}$ (« Training inputs retained by each method ») ; $\tilde{\mathbf x}^{r,(i)}_{\mathrm{NNAK},k}$, $\tilde{\mathbf x}^{r,(c)}_{\mathrm{ICAK},k}$ (« Windows of a particle and of a cluster ») ; les $\tilde{\mathbf z}$ correspondants ; $\tilde N^r_{\mathrm{RS},k}$ et $\tilde N^r_{\mathrm{AK},k}$.
3. **Libellés inexacts**, corrigés d'après les définitions du texte :
   - $\varepsilon^{(j)}$ : « Draws from the regularisation kernel » (tirés selon $K$) ;
   - $\kappa_X,\kappa_Z$ sur les ensembles de valeurs de $X$, $Z$ et $\kappa_x,\kappa_z$ sur les espaces d'état et d'observation, en deux lignes ;
   - $\mu_{x_k|x_{k-1}}$, $\mu_{z_k|x_k}$ : « Transition and likelihood embedding operators » (le texte : « their respective conditional embedding operators ») ;
   - $\boldsymbol\Sigma,\boldsymbol\Sigma_i$ : « Matrices of the separable and multi-scale kernels » ;
   - $\psi_\mathrm{RKHS}$ : « estimator » (le texte : « the least-squares estimator ») ;
   - $\hat{\mathbf P}_{\mathrm{pos},k|k-1}$ : « Covariance estimate of the predicted position ».
4. **Oubli.** $\tilde{\mathbf z}_\mathrm{old},\tilde{\mathbf z}_\mathrm{new}$ (lots d'observations du krigeage séquentiel, ch.3 et B.2), signalés par l'inventaire mais laissés de côté : ajoutés en miroir de $\tilde{\mathbf x}_\mathrm{old},\tilde{\mathbf x}_\mathrm{new}$.

**Revérification.** `where.pl` dit, pour chaque symbole, s'il apparaît dans le texte visible ou seulement en commentaire. `inv_rev.pl` ne trouve plus d'entrée sans occurrence visible, hors formes génériques dont les instances sont visibles : $\chi^2_d$, $c_d$, $d_\eta$, $d_\varepsilon$ (en exposant), $\ell_i$, $p_{A|B}$, $\mathcal S_d^+$, $w_k^{(i,c)}$. Les inventaires directs ne laissent hors liste que les symboles locaux (annexe A, preuves et leurs abréviations, $x_\infty^{(i)}$ du mean shift, $a_d,b_d,L_d$, $t$, $\bar x_k$ et $J_k$ au ch.5), les indices génériques et les instances de familles. Compilation propre, 223 pages, nomenclature de la p. xix à la p. xxiv sans ligne qui déborde. **Question laissée à lui** : sans la ligne de marques, les grandeurs du filtre de Kalman des parties linéaire et non linéaire du RBPF ne sont pas listées une à une : $\hat x^{l}$, $\mathbf P^{l}$, $\mathbf K^{n}$, $\mathbf K^{l}$, $\mathbf S^{n}$, $\mathbf S^{l}$, $y^{n}$, $y^{l}$, $\tilde{\mathbf F}^{n}$, $\tilde{\mathbf F}^{l}$ (§2.5).

**Retour arrière.** Dossier temporaire de la session `…\0debc8f0-…\scratchpad` : `0_administratif_avant_v15.tex`, `…_v16.tex`.

### 2026-10-22 (nuit) : une seule graphie pour les noyaux des espaces d'état et d'observation

Sa demande : « si possible il faut uniformiser kappa_x/X et kappa_z/Z, et je pense qu'on peut faire de même pour certains mu ».

**Choix : les majuscules.** La théorie des plongements (§3.5) définit $\kappa_X,\kappa_Z$ sur $\Omega_X,\Omega_Z$, avec $\Phi_X,\Phi_Z$, $\mathcal H_{\kappa_X},\mathcal H_{\kappa_Z}$, $\mu_X$ et $p_{XZ}$. Le ch.2 (`2_bayesian_estimation` l. 59) fait de $X$ et $Z$ l'état et l'observation génériques, par opposition aux processus indicés $x_k,z_k$. Le KKF, l'AKKF et le ch.5 instancient donc ces noyaux sur $\mathbb R^{d_x}$ et $\mathbb R^{d_z}$. La minuscule aurait mis $\kappa_x$ face à $\Omega_X$ dans la théorie, et $\kappa_x$ se lit souvent comme la section $\kappa(x,\cdot)$.

**Lot AC** (`rename.pl`, texte visible ; aucune occurrence dans les blocs `comment`) : $\kappa_x\to\kappa_X$, $\kappa_z\to\kappa_Z$, $\Phi_x\to\Phi_X$, $\Phi_z\to\Phi_Z$ ; $\mathcal H_{\kappa_x},\mathcal H_{\kappa_z}$ suivent. `3_KKF` : 95 remplacements (37 $\Phi_x$, 28 $\kappa_x$, 16 $\kappa_z$, 14 $\Phi_z$). `5_SLAM` : 5 ($\kappa_x$ ×3, $\kappa_z$ ×2 : algorithme 5.2, remarque 5.3, réglages de la campagne). Contrôles : prose, comptes d'atomes et nombre de lignes identiques.

**Nomenclature** (`nomen_v17.ps1`) : la ligne $\kappa_x,\kappa_z$ disparaît, il reste « $\kappa_X,\kappa_Z$ : Kernels on the value sets of $X$ and $Z$ » ; la ligne des applications de caractéristiques devient $\Phi,\Phi_X,\Phi_Z$. 185 entrées de symboles, en comptant ses fusions de la soirée ($p,\tilde p,p_0$ ; $q,\tilde q$ ; les gains $\mathbf K$).

**Laissés en minuscules, à trancher par lui** : $\alpha_{\mathrm{RKHS},x},\alpha_{\mathrm{RKHS},z}$ (35 atomes) et $\mathbf G^z,\mathbf G_k^x,\mathbf G_k^z,\mathbf G_{k-1}^x$ (40), qui côtoient maintenant $\kappa_X$ dans les mêmes formules (proposition 3.32, algorithme 5.2). $d_x,d_z$ restent dans tous les cas (dimensions de $x_k$ et $z_k$).

**$\mu$ : rien d'appliqué.** Aucun $\mu$ n'est écrit de deux façons pour un même objet. $\mu_{Z|X}$ (opérateur de $\mathcal H_{\kappa_X}$ dans $\mathcal H_{\kappa_Z}$) et $\mu_{Z|x}=\mu_{Z|X}\Phi_X(x)$ (plongement conditionnel au point $x$) sont deux objets, et la casse les distingue comme dans $p_{Z|X}(\cdot|x)$. $\mu_{x_k|x_{k-1}}$ et $\mu_{z_k|x_k}$ suivent les densités $p(x_k|x_{k-1})$ et $p(z_k|x_k)$ du ch.2. Question posée : lesquels visait-il ?

**Contrôles.** Compilation : 223 pages, 0 erreur, 0 non défini. Une boîte trop large de 0,29 pt, nouvelle, dans la proposition 3.32 (`3_KKF.tex` l. 297–298), parce que $\Phi_X$ est un peu plus large que $\Phi_x$ ; elle est invisible. Pages 77 et 82 relues en image.

**Retour arrière.** `…\0debc8f0-…\scratchpad\backup_AC` (`3_KKF`, `5_SLAM`) et `0_administratif_avant_v17.tex`. Attention : ses fusions $p,\tilde p$ et $q,\tilde q$ sont postérieures à cette copie de la nomenclature.

## Annexe I. Inventaire des lettres à objets multiples

*Mise à jour du 2026-10-21 : les verdicts « §5 » sont levés par l'étape 1 (matrices en gras) ; restent ceux qui renvoient à C1–C10 et L1–L11.*

| Lettre | Objets (nature) | Verdict |
|---|---|---|
| `a` | coefficients (vecteur), point de Taylor, $a=\mathbf G^{-1}(\tilde{\mathbf z}-\mathbf Cm)$ ; particule auxiliaire $a_k^{(i)}$ (APF) | ✔ même nature |
| `b` | biais $b_k^{(e)},b_k^{(o)}$ (vecteurs, BE) ; phases $b_j$ (scalaires, ch.5) ; $b_i$ local (GP:268) | ~ ch.5 |
| `c` | $c_{d_x}$ (volume de la boule), indice de classe, vecteur $c$ (KK:71) | local |
| `d` | dimensions ; différentielle ; marqueur « écarté » | ✔ (D13) |
| `e` | erreurs $e_r,e_f$ ; base canonique $e_i$ ; évaluation $e_x$ | L10 |
| `f` | évolution $f_k$ ; élément générique de $\mathcal H$ (le processus et $f_r$ sont devenus $h$, $h_r$) | ✔ fonctions, fait |
| `g` | tendance $g,g_l$ ; $g_r$ (matrice) ; $\Delta_g$ | ✔ |
| `h` | observation $h_k$ ; largeur de bande $h$ | C10 |
| `k` | noyau / temps | **C1** |
| `l` | longueur de corrélation / indice / marqueur | **C9** |
| `m` | moyenne ; mesure (BE:64) ; entiers ; $m_+$ | C3, L7 |
| `n` | entiers ; marqueur ; $n_\mathrm{MC}$ (le bruit est devenu $\varepsilon_k$) | ✔ fait |
| `p` | densité ; entier ; probabilité de Bernoulli ; position | L6, L11 |
| `r` | rayon ; entiers ; marqueur | §7 |
| `u`, `v` | $v$ : poids auxiliaires de l'APF ; $u$ et $y$ : vecteurs muets et réalisations (la vitesse est devenue $\dot p$) | ✔ fait |
| `w` | poids (le bruit est devenu $\eta$) | ✔ fait |
| `y` | innovation ; variable d'intégration | §6 |
| `α` | vecteur / facteurs | **C5** |
| `β` | vecteur / facteur d'oubli (ch.5) | **C5** |
| `δ` | Dirac ; Kronecker ; confiance ; décalage | §6 |
| `ε` | bruits $\varepsilon$ ; seuils $\epsilon$ | ✔ fait |
| `λ` | valeurs propres ; régularisation ; échelle (les poids auxiliaires sont devenus $v$) | §7 |
| `μ` | mesures $\mu_p$, $\mu_q$ ; plongement ; moyenne (→ $m$ à faire) ; l'auxiliaire est devenue $a$ | C3 (moyennes) |
| `ν` | Matérn (l'instrumentale est devenue $\mu_q$) | ✔ fait |
| `ξ` | gaussiennes / multiplicateur | **C6** |
| `ψ`, `φ`, `ϕ` | fonctions de rôles divers | §7 |
| `A`, `B` | matrices ; paire d'estimateurs | L8 |
| `C` | contrainte ; classes ; constantes ; bloc | §5 |
| `D` | SVD ; résidu KKF ; bloc | ✔ matrices |
| `E` | espace mesuré ; domaines ; perte | L1, L7 |
| `F` | évolution ; espace normé ; vecteur aléatoire | §5 |
| `G` | Gram ; pénalité ; $\mathcal G$ | §5, **C8** |
| `J` | jacobienne ; coût | §5 |
| `K` | gain ; noyau de régularisation ; Bessel ; $\mathcal K$ | §5 |
| `M` | covariance de tendance ; $\mathbf M$ générique ; SVD | ✔ matrices |
| `P`, `Q` | covariances ; génériques (annexe A) | L2 |
| `S` | innovation ; Schur ; mesure ; densité spectrale ; système | §5, L4 |
| `T` | transition ; opérateur ; itérations | L5 |
| `U`, `V` | matrices / vecteurs aléatoires | §5 |
| `W` | poids du krigeage / poids UKF | §5, L3 |
| `𝒰` | ensemble / loi | **C7** |
| `𝒢` | matrice en calligraphique | **C8** |

Lettres à objet unique (aucun problème) : $\theta,\Theta,\rho$ (distance), $\sigma$, $\tau$, $\chi$, $\Gamma$, $\Delta$, $\Phi$, $\boldsymbol\Lambda$, $\boldsymbol\Psi$, $\boldsymbol\Pi$, $\pi$, $\omega$ (fréquences), $q$, $s$, $t$, $x$, $z$, $N$, $L$, $\mathbf H$, $\mathbf I$, $\mathbf O$, $\mathbf R$, $X$, $Z$, $\mathcal H$, $\mathcal L$, $\mathcal M$, $\mathcal N$, $\mathcal O$.

---

## Annexe II. Encart « Conventions » proposé pour la nomenclature (anglais, tel qu'il irait dans le manuscrit)

> **Notation conventions.** Lowercase italic letters denote scalars, vectors and functions ($x_k$, $\theta$, $h$, $\kappa$); lowercase bold letters denote tuples and stacked vectors ($\mathbf x=(x_1,\dots,x_n)$, $\mathbf 1_n$); uppercase bold letters denote matrices ($\mathbf P$, $\mathbf G$, $\boldsymbol\Sigma$); uppercase italic letters denote random variables ($X$, $Z$), sizes ($N$, $\tilde N$) and named functions ($K$, $L$); calligraphic letters denote sets, spaces and operators ($\mathcal U$, $\mathcal H$, $\mathcal K$); lowercase Greek letters denote measures and their embeddings ($\mu_p$, $\mu_q$, $\mu_X$), noises ($\eta_k$, $\varepsilon_k$) and parameters. A hat denotes an estimate or an estimator; a tilde denotes training data ($\tilde{\mathbf x}$, $\tilde{\mathbf z}$, $\tilde N$), an unnormalised quantity ($\tilde p$, $\tilde w$), or a quantity evaluated at the current linearisation point ($\tilde{\mathbf F}_k$, $\tilde{\mathbf H}_k$); a bar denotes a mean, a block form ($\bar\kappa$, $\bar{\mathbf G}$) or a decorrelated quantity of the RBPF. Subscripts $k$ and $k\vert l$ index time, superscripts $(i)$ and $(c)$ index particles and clusters, superscripts $n$ and $l$ denote the non-linear and linear parts of the state, and superscripts $r$ and $d$ denote retained and discarded training data.
