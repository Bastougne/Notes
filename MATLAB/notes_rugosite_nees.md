# Rugosité locale du champ, et ce qu'elle explique de la NEES

Mesures du 2026-09-21, tirées des campagnes déjà enregistrées et de la carte. **Aucune
simulation nouvelle.** Elles répondent à une question de §4.3.3 : pourquoi la NEES du krigeage
statique monte vers le pas 70 alors que celle de l'AK reste plate.

## La trajectoire

Le vol de référence fait 150 pas : ligne droite jusqu'au 32, **premier demi-tour des pas 33 à
64**, ligne droite jusqu'au 96, **second demi-tour des pas 97 à 128**, puis ligne droite. Chaque
demi-tour dure 32 pas et tourne de 180°.

## La rugosité du champ le long du vol

Écart-type du champ vrai dans un carré de $\pm10$ km autour de la position vraie, échantillonné
tous les 2 km, lu avec le `lire` de `krigeage.m` :

| $k$ | 20 | 33 | **45** | 64 | **70** | **85** | 97 | 110 | 125 | 140 |
|---|---|---|---|---|---|---|---|---|---|---|
| écart-type (nT) | 113 | 59 | **268** | 98 | **177** | **71** | 121 | 120 | 120 | 133 |
| gradient (nT/km) | 39 | 23 | **99** | 37 | **71** | **34** | 54 | 49 | 37 | 47 |

- Médiane sur le vol : **125 nT** et 41,5 nT/km.
- **Maximum de tout le vol au pas 47** (281 nT).
- Par segment, les cinq tronçons se valent : 131, 112, 122, 117 et 132 nT. **Le découpage
  droite/virage ne porte donc aucune information de rugosité.**

## La NEES décomposée

Médianes sur les essais convergés (seuil 2000 m), campagne A pour le statique
(`20260904_084456`), lot des quatre modèles pour l'AK (`20260912_030842`). L'« incertitude
déclarée » est `err / √d²`, en mètres : un indicateur, pas un écart-type, puisque $d^2$ porte sur
les quatre composantes de l'état et l'erreur sur les deux positions.

**Maille 4,5 km**

| $k$ | 20 | 33 | 45 | 64 | **70** | 85 | 97 | 110 | 125 | 140 |
|---|---|---|---|---|---|---|---|---|---|---|
| statique, erreur (m) | 513 | 622 | 666 | 819 | **1048** | 792 | 628 | 559 | 459 | 799 |
| statique, incertitude (m) | 401 | 452 | 379 | 395 | **293** | 435 | 365 | 372 | 507 | 463 |
| statique, NEES | 0,46 | 0,47 | 0,74 | 1,01 | **2,29** | 0,80 | 0,72 | 0,61 | 0,22 | 0,69 |
| AK, erreur (m) | 533 | 608 | 624 | 724 | **721** | 687 | 618 | 567 | 460 | 616 |
| AK, incertitude (m) | 418 | 464 | 537 | 536 | **572** | 459 | 437 | 409 | 509 | 489 |
| AK, NEES | 0,42 | 0,40 | 0,32 | 0,41 | **0,38** | 0,52 | 0,42 | 0,49 | 0,20 | 0,43 |

À 7,5 km, au pas 70 : statique 1670 m d'erreur pour 503 m déclarés, AK 1342 pour 688.

**Au pas 70, le statique atteint sa plus grande erreur et sa plus petite incertitude déclarée de
tout le vol** : il devient faux et confiant à la fois. L'AK, au même endroit, élargit la sienne.

## Corrélations avec la rugosité

Spearman, phase de suivi ($k\geq35$), sur les courbes médianes :

| | statique 4,5 km | AK 4,5 km | statique 7,5 km | AK 7,5 km |
|---|---|---|---|---|
| NEES | $+0{,}29$ | $-0{,}03$ | $+0{,}31$ | $+0{,}33$ |
| erreur | $+0{,}07$ | $+0{,}01$ | $+0{,}15$ | $+0{,}00$ |
| **incertitude déclarée** | $\mathbf{-0{,}51}$ | $\mathbf{-0{,}01}$ | $\mathbf{-0{,}48}$ | $-0{,}21$ |

Et le rapport des deux incertitudes déclarées, AK sur statique, corrèle $+0{,}59$ à 4,5 km et
$+0{,}39$ à 7,5 km.

## Ce que cela établit

**Le virage n'explique rien.** Deux arguments indépendants. Les deux filtres volent la même
trajectoire avec les mêmes bruits, donc un effet de dynamique se verrait dans les deux : or ils
divergent au pas 70, et c'est bien le cas des traits communs, comme le creux du pas 125, présent
chez les deux. Et l'incertitude du statique se resserre pendant le premier virage (452 → 379 m)
mais s'élargit pendant le second (365 → 507 m) : un virage ne produit pas deux fois le même
effet. Objection de Bastien qui a tué cette piste : il y a **quatre** changements de dynamique,
et la NEES n'en marque que deux.

**La sélection locale seule n'explique rien non plus.** À 4,5 km la fenêtre est au plancher de 25
points sur toute la phase de suivi : même effectif et même géométrie au pas 70 qu'au pas 85,
alors que le rapport des incertitudes passe de 1,95 à 1,06. Surtout, **à noyau figé la variance de
krigeage ne dépend que de la position des points retenus, jamais des valeurs du champ** : une
sélection ne peut pas produire un rapport corrélé à $+0{,}59$ avec la rugosité. Seul le
réajustement des hyperparamètres le peut.

**Reste le champ, et la mesure le soutient.** L'incertitude déclarée du statique est anticorrélée
à $-0{,}51$ avec la rugosité locale : là où le champ varie le plus, l'observation discrimine
mieux et le nuage se resserre, mais l'erreur de reconstruction y est aussi la plus grande, et le
noyau global, aveugle à l'amplitude locale, ne l'élargit pas. Celle de l'AK est plate ($-0{,}01$).
La divergence commence au pas 45, et le maximum de rugosité du vol est au pas 47.

## Réserves

- Les corrélations portent sur des séries temporelles autocorrélées de 116 pas : ce sont des
  statistiques descriptives, pas des tests.
- À 7,5 km le réajustement ne compense plus entièrement : la NEES de l'AK y corrèle aussi
  $+0{,}33$ avec la rugosité et son incertitude $-0{,}21$. Vingt-cinq points sur une maille lâche
  ne suffisent pas.
- L'essai `refit = 'global'` — fenêtre conservée, hyperparamètres figés — confirmerait
  directement ce que l'élimination établit ici. Il n'est plus nécessaire à l'argument.

## Où cela va dans le manuscrit

Le TODO est posé en `4_kriging_navigation.tex:372`, dans un bloc `comment` après la figure
d'ARMSE ; **Bastien le déplacera autour de `:288`**, avec la figure de la carte, la rugosité étant
une propriété du scénario. §4.3.3 pourra alors s'y appuyer pour justifier son explication de la
NEES par la non-stationnarité, plutôt que de la laisser au conditionnel.

## Reproduire

Rugosité le long du vol :

```matlab
A = load(fullfile('resultats','20260904_084456_reference_A_1000runs.mat'), 'scen', 'par');
cm = load(A.par.map.file);
map = struct('h', double(cm.mesures).', 'step', double(cm.pas(:))' .* [1 1]);
krig = krigeage();  X = A.scen.x_true(1:2, :);
[dx, dy] = meshgrid(-10e3:2e3:10e3);  off = [dx(:)'; dy(:)'];
rug = arrayfun(@(k) std(krig.lire(map.h, map.step, X(:, k) + off)), 1:150);
```

Décomposition, pour une entrée `r` de `res` : `ok = median(r.err(end-9:end, :), 1) <= 2000`, puis
`median(r.err(:, ok), 2)`, `median(r.d2(:, ok), 2, 'omitnan') / d_x` et
`median(r.err(:, ok) ./ sqrt(r.d2(:, ok)), 2, 'omitnan')`.
