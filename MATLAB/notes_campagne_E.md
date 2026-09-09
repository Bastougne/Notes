# Campagne E — les deux contrôles négatifs de la sélection

Note de passation, écrite le 2026-09-08. Elle explique **pourquoi** ces deux modèles existent, ce
qu'on attend d'eux, ce qui a été ajouté au code, et comment lancer. Rien n'a encore été mesuré
sérieusement : seul un test à blanc à 2 essais a tourné.

## Pourquoi

Le chapitre 4 propose deux contributions successives, le krigeage adaptatif (AK, FUSION 2024) et sa
version par clusters (CA-OK, EUSIPCO 2025). Chacune repose sur un choix de conception qui, tel quel,
ressemble à une astuce :

- l'AK borne par le bas le nombre d'échantillons retenus, $n_{\min}$ ;
- le CA-OK partitionne le nuage avant de sélectionner, au lieu d'une seule fenêtre.

Ces deux modèles sont les **contrôles négatifs** qui rendent ces choix nécessaires plutôt
qu'arbitraires. Chacun échoue précisément par le défaut que l'un des deux choix corrige. C'est aussi
la réponse à une inquiétude de Bastien sur la maigreur des deux sections : un contrôle qui échoue
pour la raison exacte qui motive un choix vaut plus qu'un paragraphe d'exposition.

### `ak_alea` — la place des points contre leur nombre

Tire un sous-ensemble uniforme du relevé, **une fois pour toutes**, et krige dessus toute la mission
avec des hyperparamètres ajustés une seule fois. C'est l'analogue honnête d'un relevé éclairci.

Il sert **deux expériences distinctes**, que Bastien avait d'abord confondues :

1. **À budget égal à la fenêtre** (`n_alea = n_min = 25`) : est-ce que l'AK gagne parce qu'il prend
   peu de points, ou parce qu'il les prend *au bon endroit* ? Un tirage de 25 points sur un corridor
   de $158\times187$ km donne un espacement typique de $34$ km, très au-delà de la longueur de
   corrélation ajustée ($\ell \approx 5{,}5$ km) : le krigeage rend la moyenne partout.
2. **À effectif égal à un réseau plus grossier** (`n_alea` = l'effectif d'une maille du balayage,
   par exemple $676$ pour $7{,}5$ km) : un tirage uniforme est-il équivalent à agrandir la maille ?
   **Non**, et c'est le point à démontrer. À densité moyenne égale, un tirage uniforme laisse des
   trous — le plus grand cercle vide d'un processus de Poisson croît en $p\sqrt{\log n}$ — là où un
   réseau garantit une distance maximale de $p/\sqrt{2}$ au plus proche échantillon. Le côté régulier
   du balayage de maille est **déjà mesuré** dans la campagne A, donc cette moitié est gratuite.

### `ak_boules` — la sélection « évidente »

La fenêtre est la réunion des boules de rayon $\rho$ centrées sur **chaque particule**, sans plancher
et sans clustering. C'est ce qu'un lecteur proposerait spontanément. On attend qu'elle échoue par les
deux bouts, et ce sont exactement les deux défauts que $n_{\min}$ et le clustering corrigent :

1. **En régime multimodal**, un seul ajustement d'hyperparamètres sur une réunion qui enjambe les
   modes revient à réapprendre $\theta$ sur des régions séparées de dizaines de kilomètres — le
   défaut de moyennage global auquel le fenêtrage devait échapper. C'est la motivation du CA-OK.
   (Variante à mentionner si un rapporteur la soulève : ajuster par composante connexe de la réunion
   serait un clustering géométrique, donc un CA-OK du pauvre. Le mean-shift lui est préféré pour le
   coût et pour le contrôle explicite de l'échelle par la bande passante.)
2. **Quand le nuage se resserre**, la réunion se réduit à une boule de rayon $\approx\rho$ et ne
   contient plus qu'une poignée d'échantillons, voire aucun. Chiffres du scénario de référence :
   RMSE finale $\approx 380$ m, donc un nuage d'extension de l'ordre du kilomètre, contre une maille
   de $3{,}5$ km — le nuage entier tient dans une maille. À $\rho = 3{,}5$ km la réunion tient
   $\approx 4$ échantillons, à $\rho = 1{,}75$ km environ $1$. C'est la démonstration empirique que
   $n_{\min}$ est nécessaire, et non une précaution de confort.

## Ce que le test à blanc a déjà montré

Deux essais, 40 pas, un seul relevé — **à ne citer sous aucune forme**, c'est un test de
non-régression, pas une mesure. Les deux mécanismes y sont néanmoins visibles :

| modèle | RMSE $k=10$ | $k=30$ | $k=40$ | conv. seuil | conv. Mahalanobis | fenêtre |
|---|---|---|---|---|---|---|
| OK adaptatif | $945$ | $587$ | $608$ | $100\%$ | $100\%$ | $130$ |
| OK aléatoire | $15722$ | $15193$ | $15148$ | $0\%$ | $100\%$ | $25$ |
| OK par boules | $765$ | $654$ | $1317$ | $100\%$ | $100\%$ | $39$ |

- L'aléatoire ne converge jamais et reste à son incertitude initiale, comme prévu.
- Les boules partent bien puis **se dégradent** quand le nuage se resserre — l'effondrement annoncé.
- Bonus inattendu et utile au chapitre : l'aléatoire est à $0\%$ sous le critère de seuil mais
  $100\%$ sous Mahalanobis. Il sait parfaitement qu'il ne sait rien. C'est **l'illustration la plus
  nette possible de la mise en garde déjà écrite** sur le critère de Mahalanobis — un filtre très
  faux le satisfait s'il annonce une incertitude assez grande. À rapprocher de
  `4_kriging_navigation.tex:331`.

## Le code ajouté

Tout est dans `krigeage.m` et `campagne_chapitre_4.m`, dans l'idiome existant.

### `krigeage.m`

Deux cas dans `build_model`, insérés avant `case 'cak'` :

- **`'ak_alea'`** → nom de modèle `OK aléatoire`. Tire `idx_a` par `randperm` sur un `RandStream`
  `threefry` de graine `par.krig.survey_seed + 7`, donc **reproductible et indépendant** des autres
  flux. Ajuste les hyperparamètres une fois sur ce sous-ensemble, à la construction.
- **`'ak_boules'`** → nom de modèle `OK par boules`. Construit un `KDTreeSearcher` sur le relevé une
  seule fois.

Deux fonctions de requête, insérées avant `query_cak` :

- **`query_alea`** : résolution exacte sur le sous-ensemble figé. Coût `[n; n^3; 0; 0]` — zéro
  ajustement par pas, puisqu'ils sont figés. La factorisation s'amortirait sur la mission ; on la
  compte par pas pour rester comparable aux autres lignes, et la remarque `th:tabulation` vaut ici.
- **`query_boules`** : `rangesearch` sur l'arbre, réunion des indices, puis réajustement et
  résolution comme l'AK. Coût `[n; n^3; 1; n_dist]` avec `n_dist = n_part * log2(n_relevé)`, la
  recherche de voisinage logée dans la **quatrième case**, celle du mean-shift, pour rester
  comparable au CA-OK. **Cas dégénéré traité et non levé** : si la réunion tient moins de deux
  points, la requête rend la loi a priori du relevé (moyenne du relevé, variance
  $\sigma_f^2+\sigma_{\mathrm{map}}^2+\sigma_{\mathrm{obs}}^2$) et enregistre `n = 0`. L'échec doit
  être *mesurable*, pas fatal — sans quoi la campagne s'arrête au lieu de documenter l'effondrement.

### `campagne_chapitre_4.m`

Deux paramètres par défaut, après `par.krig.bandwidth` :

```matlab
par.krig.n_alea      = 25;      % budget du tirage uniforme
par.krig.rayon_boule = 3500;    % m, rayon des boules
```

Et une **campagne `'E'`**, six modèles, sans balayage :
`{'carte', 'ok', 'ak', 'cak', 'ak_alea', 'ak_boules'}`.

Le balayage du rayon est délibérément **sorti** de la campagne : seul `ak_boules` en dépend, et
rejouer la tabulation de `ok` à chaque valeur coûterait bien plus cher que la mesure elle-même.

## Comment lancer

La comparaison principale :

```bash
matlab -batch "cd('D:/bhubert/Desktop/Notes/MATLAB'); over.campagne='E'; over.mc.n_runs=1000; over.krig.n_surveys=100; campagne_chapitre_4"
```

Le balayage du rayon des boules, seul :

```bash
matlab -batch "cd('D:/bhubert/Desktop/Notes/MATLAB'); over.compare.models={'ak_boules'}; over.balayage=struct('champ','krig.rayon_boule','valeurs',[1750 3500 7000 14000]); over.mc.n_runs=1000; campagne_chapitre_4"
```

Le second volet de `ak_alea`, tirage contre réseau à effectif égal — une exécution par effectif, à
comparer aux mailles déjà mesurées de la campagne A :

```bash
matlab -batch "cd('D:/bhubert/Desktop/Notes/MATLAB'); over.compare.models={'ak_alea'}; over.balayage=struct('champ','krig.n_alea','valeurs',[676 900 1296 1936 3136 6241]); over.mc.n_runs=1000; campagne_chapitre_4"
```

Les six valeurs sont les effectifs $\tilde n$ des six mailles du balayage, de $7{,}5$ à $2{,}5$ km.

## Points de vigilance

- **Le temps de calcul n'est pas une contrainte** (décision du 2026-09-03) : lancer à 1000 essais
  d'emblée, et se méfier des projections de durée, qui se sont trompées de $-15$ à $+73\%$.
- `ak_boules` fait une `rangesearch` de $5000$ particules par pas. C'est un coût réel du modèle, à
  compter dans l'analyse et non à optimiser : il fait partie de ce qui le disqualifie.
- Les hyperparamètres de `ak_alea` sont ajustés sur $25$ points quand `n_alea = 25`. L'ajustement est
  alors mauvais, et **c'est voulu** : c'est une conséquence du budget, pas un défaut d'implémentation.
- Le seuil de convergence du manuscrit est à $2000$ m dans les figures, contre `par.conv.seuil` à
  $3500$ dans la campagne. `figures_campagnes.m` le redéfinit ; ne pas s'étonner de l'écart.
- Modèles appariés par sous-chaîne : `'ak_alea'` s'appelle `OK aléatoire` et `'ak_boules'`
  `OK par boules`. Les deux contiennent `OK`, comme `OK carte totale` et `OK adaptatif`. **Ne jamais
  apparier sur `'OK'` ni sur `'carte'`** — un piège déjà tombé une fois, qui avait silencieusement
  faussé deux figures.
